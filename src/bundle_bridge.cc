/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle_bridge.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"
#include "bridger.h"

bundle_bridge::bundle_bridge(bundle_base &b)
	: bb(b)
{
	//compute_strand();
}

bundle_bridge::~bundle_bridge()
{}

int bundle_bridge::build()
{
	build_junctions();
	extend_junctions();
	build_regions();

	align_hits_transcripts();
	index_references();

	remove_tiny_boundary();
	build_fragments();
	//group_fragments();

	bridger bdg(this);
	bdg.bridge();
	return 0;
}

int bundle_bridge::build_junctions()
{
	int min_max_boundary_quality = min_mapping_quality;
	map< int64_t, vector<int> > m;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		vector<int64_t> v = bb.hits[i].spos;
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		int32_t p1 = high32(it->first);
		int32_t p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bb.hits[v[k]];
			nm += h.nm;
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", bb.chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(it->first, v.size());
		jc.nm = nm;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);

	}
	return 0;
}

int bundle_bridge::extend_junctions()
{
	map< int64_t, vector<int> > m;
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		vector<PI32> v = ref_trsts[i].get_intron_chain();
		for(int k = 0; k < v.size(); k++)
		{
			assert(v[k].first < v[k].second);
			// TODO, TODO
			if(v[k].first <= bb.lpos) continue;
			if(v[k].second >= bb.rpos) continue;
			int64_t p = pack(v[k].first, v[k].second);

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		for(int k = 0; k < v.size(); k++)
		{
			char c = ref_trsts[v[k]].strand;
			if(c == '.') s0++;
			if(c == '+') s1++;
			if(c == '-') s2++;
		}

		junction jc(it->first, 0 - v.size());
		jc.nm = 0;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);
	}
	return 0;
}

int bundle_bridge::build_regions()
{
	MPI s;
	s.insert(PI(bb.lpos, START_BOUNDARY));
	s.insert(PI(bb.rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	for(int k = 0; k < v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second; 
		int rtype = v[k + 1].second; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		region rr(l, r, ltype, rtype);
		evaluate_rectangle(bb.mmap, l, r, rr.ave, rr.dev, rr.max);
		regions.push_back(rr);
	}

	return 0;
}

int bundle_bridge::align_hits_transcripts()
{
	map<int32_t, int> m;
	for(int k = 0; k < regions.size(); k++)
	{
		if(k >= 1) assert(regions[k - 1].rpos == regions[k].lpos);
		m.insert(pair<int32_t, int>(regions[k].lpos, k));
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		align_hit(m, bb.hits[i], bb.hits[i].vlist);
		bb.hits[i].vlist = encode_vlist(bb.hits[i].vlist);
	}

	ref_phase.resize(ref_trsts.size());
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		align_transcript(m, ref_trsts[i], ref_phase[i]);
	}

	return 0;
}

int bundle_bridge::align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &vv)
{
	vv.clear();
	vector<int64_t> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return 0;

	vector<PI> sp;
	sp.resize(v.size());

	int32_t p1 = high32(v.front());
	int32_t p2 = low32(v.back());

	sp[0].first = locate_region(p1);
	for(int k = 1; k < v.size(); k++)
	{
		p1 = high32(v[k]);

		map<int32_t, int>::const_iterator it = m.find(p1);

		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[sp.size() - 1].second = locate_region(p2 - 1);
	for(int k = 0; k < v.size() - 1; k++)
	{
		p2 = low32(v[k]);
		map<int32_t, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = 0; k < sp.size(); k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > 0) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle_bridge::align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &vv)
{
	vv.clear();
	int k1 = -1;
	int k2 = -1;
	for(int k = 0; k < t.exons.size(); k++)
	{
		if(t.exons[k].second > bb.lpos)
		{
			k1 = k;
			break;
		}
	}
	for(int k = t.exons.size() - 1; k >= 0; k--)
	{
		if(t.exons[k].first < bb.rpos)
		{
			k2 = k;
			break;
		}
	}

	if(k1 > k2) return 0;
	if(k1 == -1 || k2 == -1) return 0;

	vector<PI> sp;
	sp.resize(k2 + 1);

	int32_t p1 = t.exons[k1].first > bb.lpos ? t.exons[k1].first : bb.lpos;
	int32_t p2 = t.exons[k2].second < bb.rpos ? t.exons[k2].second : bb.rpos;

	sp[k1].first = locate_region(p1);
	for(int k = k1 + 1; k <= k2; k++)
	{
		p1 = t.exons[k].first;
		map<int32_t, int>::const_iterator it = m.find(p1);
		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[k2].second = locate_region(p2 - 1);
	for(int k = k1; k < k2; k++)
	{
		p2 = t.exons[k].second;
		map<int32_t, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = k1; k <= k2; k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > k1) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle_bridge::index_references()
{
	ref_index.clear();
	ref_index.resize(regions.size());
	for(int k = 0; k < ref_phase.size(); k++)
	{
		vector<int> &v = ref_phase[k];
		for(int j = 0; j < v.size(); j++)
		{
			int x = v[j];
			ref_index[x].push_back(PI(k, j));
		}
	}
	return 0;
}

int bundle_bridge::locate_region(int32_t x)
{
	if(regions.size() == 0) return -1;

	int k1 = 0;
	int k2 = regions.size();
	while(k1 < k2)
	{
		int m = (k1 + k2) / 2;
		region &r = regions[m];
		if(x >= r.lpos && x < r.rpos) return m;
		else if(x < r.lpos) k2 = m;
		else k1 = m;
	}
	return -1;
}

int bundle_bridge::remove_tiny_boundary()
{
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		vector<int> vv = decode_vlist(h.vlist);
		int n = vv.size();
		if(n >= 2 && vv[n - 2] + 1 == vv[n - 1])
		{
			int k = vv[n - 1];
			int32_t total = regions[k].rpos - regions[k].lpos;
			int32_t flank = h.rpos - regions[k].lpos;

			if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
			{
				vector<int> v(vv.begin(), vv.begin() + n - 1);
				assert(v.size() + 1 == vv.size());
				h.vlist = encode_vlist(v);
				h.rpos = regions[k].lpos;
			}
		}

		vv = decode_vlist(h.vlist);
		n = vv.size();
		if(n >= 2 && vv[0] + 1 == vv[1])
		{
			int k = vv[0];
			int32_t total = regions[k].rpos - regions[k].lpos;
			int32_t flank = regions[k].rpos - h.pos;

			if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
			{
				vector<int> v(vv.begin() + 1, vv.end());
				assert(v.size() + 1 == vv.size());
				h.vlist = encode_vlist(v);
				h.pos = regions[k].rpos;
			}
		}
	}
	sort(bb.hits.begin(), bb.hits.end(), compare_hit_pos);
	return 0;
}

int bundle_bridge::build_fragments()
{

	int ctp = 0;// count fragments number from paired-end reads
	int ctu = 0;// count fragments number from UMI
	int ctb = 0;// count fragments number for both

	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	fragments.clear();
	if(bb.hits.size() == 0) return 0;

	int max_index = bb.hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index);

	// first build index
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		if(h.isize >= 0) continue;
		if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.qhash % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i);
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		if(h.paired == true) continue;
		if(h.isize <= 0) continue;
		if(h.vlist.size() == 0) continue;

		int k = (h.qhash % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			//if(z.hi != h.hi) continue;
			if(z.paired == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = vv[k][j];
			break;
		}

		/*
		SI si(h.qname, h.hi);
		MSI::iterator it = vv[k].find(si);
		if(it == vv[k].end()) continue;
		int x = it->second;
		*/

		//printf("HIT: i = %d, x = %d, bb.hits[i].vlist = %lu | ", i, x, bb.hits[i].vlist.size(), bb.hits[i].qname.c_str()); bb.hits[i].print();

		if(x == -1) continue;
		if(bb.hits[x].vlist.size() == 0) continue;

		fragment fr(&bb.hits[i], &bb.hits[x]);

		// ===============================
		// TODO: dit for UMI
		bb.hits[i].pi = x;
		bb.hits[x].pi = i;
		bb.hits[i].fidx = fragments.size();
		bb.hits[x].fidx = fragments.size();
		ctp += 1;
		fr.type = 0; 
		// ================================
		fr.lpos = h.pos;
		fr.rpos = bb.hits[x].rpos;

		vector<int> v1 = decode_vlist(bb.hits[i].vlist);
		vector<int> v2 = decode_vlist(bb.hits[x].vlist);
		fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
		fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
		fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
		fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
		}

		fragments.push_back(fr);

		bb.hits[i].paired = true;
		bb.hits[x].paired = true;

	}

	//printf("total bb.hits = %lu, total fragments = %lu\n", bb.hits.size(), fragments.size());
	

	// TODO
	// ===============================================
	// build fragments based on UMI
	// ===============================================
	
	vector<string> ub;
	vector< vector<int> > hlist;

	int sp = 0;

	// create ub and corresponding hlist
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		assert(h.pos >= sp);
		sp = h.pos;

		if((h.flag & 0x4) >= 1 || h.umi == "") continue;

		bool new_umi = true;
		int ubidx = -1;
		for(int j = 0; j < ub.size(); j++)
		{
			if(h.umi == ub[j])
			{
				new_umi = false;
				ubidx = j;
				break;
			}
		}

		if(new_umi)
		{
			ub.push_back(h.umi);
			vector<int> head;
			head.push_back(i);
			hlist.push_back(head);
		}
		else hlist[ubidx].push_back(i);
	}


	/*	
	// print ub and hlist
        assert(ub.size() == hlist.size());
        for(int k = 0; k < ub.size(); k++)
        {
                printf("ub: %s, hlist # %d: (  ", ub[k].c_str(), k);
                for(int kk = 0; kk < hlist[k].size(); kk++)
                {
                        printf("hit %d: ", hlist[k][kk]);
			vector<int> v1 = decode_vlist(bb.hits[(hlist[k][kk])].vlist);
			for(int kkk = 0; kkk < v1.size(); kkk++)
			{
				printf("%d ", v1[kkk]);
			}
			printf(";");

                }
                printf(")\n");
        }
	*/

	

	/*
	//print ub
	for(int k = 0; k < ub.size(); k++)
        {
                printf("check repeat UMI %s\n", ub[k].c_str());
        }
	*/

	

	// build fragments based on hlist
	// every two consecutive hits in hlist stored as one fragment
	int hidx1 = -1;
	int hidx2 = -1;

	umiLink.clear();

	for(int i = 0; i < hlist.size(); i++)
	{
		assert(hlist[i].size() > 0);
		if(hlist[i].size() == 1) continue;

		vector<int> flist;
		flist.clear();

		for(int j = 0; j < hlist[i].size() - 1; j++)
		{
			hidx1 = hlist[i][j];
			hidx2 = hlist[i][j+1];

			// check whether has been paired
			if ( bb.hits[hidx1].pi == hidx2 && bb.hits[hidx2].pi == hidx1 && bb.hits[hidx1].paired == true && bb.hits[hidx2].paired == true)
			{
				//printf("Already exist paired-end fr: (%d, %d), fr# %d\n", hidx1, hidx2, bb.hits[hidx2].fidx);
				//assert(bb.hits[hidx1].fidx == bb.hits[hidx2].fidx);
				int fr_idx = bb.hits[hidx1].fidx;
				//printf("bothing....................fr.size() = %d, fidx = %d\n", fragments.size(), fr_idx);
				fragments[fr_idx].type = 2;
				ctb += 1;

				// TODO add fr index to umiLink
				flist.push_back(fr_idx);

				
				continue;
			}

			if(bb.hits[hidx1].vlist.size() == 0 || bb.hits[hidx2].vlist.size() == 0) continue;

			fragment fr(&bb.hits[hidx1], &bb.hits[hidx2]);
			fr.type = 1;
			ctu += 1;
			fr.lpos = bb.hits[hidx1].pos;
			fr.rpos = bb.hits[hidx2].rpos;

			vector<int> v1 = decode_vlist(bb.hits[hidx1].vlist);
			vector<int> v2 = decode_vlist(bb.hits[hidx2].vlist);
			fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
			fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
			fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
			fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

			fr.b1 = true;
			if(v1.size() <= 1)
			{
				fr.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
			}

			fr.b2 = true;
                	if(v2.size() <= 1)
                	{
                        	fr.b2 = false;
                	}
                	else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
                	{
                        	if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
                	}
                	else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
                	{
                        	if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
                	}

			fragments.push_back(fr);
			bb.hits[hidx1].paired = true;
			bb.hits[hidx2].paired = true;

			int cur_fidx = fragments.size() - 1;
			flist.push_back(cur_fidx);

		}

		umiLink.push_back(flist);

	}

	//printf("total bb.hits = %lu, total fragments = %lu, paired-end fragments = %d, UMI only fragments = %d, both = %d\n", bb.hits.size(), fragments.size(), ctp, ctu, ctb);


	/*
        // print ub, hlist, umiLink
        assert(ub.size() == hlist.size());
        for(int k = 0; k < ub.size(); k++)
        {
                printf("ub: %s, hlist # %d: (  ", ub[k].c_str(), k);
                for(int kk = 0; kk < hlist[k].size(); kk++)
                {
                        printf("hit %d: ", hlist[k][kk]);
                        vector<int> v1 = decode_vlist(bb.hits[(hlist[k][kk])].vlist);
                        for(int kkk = 0; kkk < v1.size(); kkk++)
                        {
                                printf("%d ", v1[kkk]);
                        }
                        printf(";");

                }
                printf(")\n");
        }

	// print umiLink
	printf("umiLink size = %d\n", umiLink.size());
	for(int k = 0; k < umiLink.size(); k++)
	{
		printf("umiLInk %d : (", k);
		for(int kk = 0; kk < umiLink[k].size(); kk++)
		{
			printf("%d, ", umiLink[k][kk]);
		}
		printf(")\n");
	}
	*/
        

	return 0;
}

int bundle_bridge::group_fragments()
{
	if(fragments.size() == 0) return 0;

	sort(fragments.begin(), fragments.end(), compare_fragment);

	vector<fragment> ff;

	fragment fx = fragments[0];
	assert(fx.h1->vlist.size() >= 1);
	assert(fx.h2->vlist.size() >= 1);
	for(int k = 1; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		assert(fr.h1->vlist.size() >= 1);
		assert(fr.h2->vlist.size() >= 1);

		if(fx.equal(fr) == true)
		{
			fx.append(fr);
		}
		else
		{
			ff.push_back(fx);
			fx = fr;
		}
	}
	ff.push_back(fx);
	fragments = ff;

	//printf("grouped fragments = %lu\n", fragments.size());
	return 0;
}

int32_t bundle_bridge::compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v)
{
	if(v.size() == 0) return 0;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i];
		flen += regions[k].rpos - regions[k].lpos;
	}
	return flen - k1l - k2r;
}

int bundle_bridge::print(int index)
{
	printf("Bundle %d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	printf("tid = %d, #hits = %lu, #fragments = %lu, #ref-trsts = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			bb.tid, bb.hits.size(), fragments.size(), ref_trsts.size(), bb.chrm.c_str(), bb.lpos, bb.rpos, bb.strand, n0, np, nq);

	// print ref-trsts
	//for(int k = 0; k < ref_trsts.size(); k++) ref_trsts[k].write(cout);

	// print fragments 
	//for(int i = 0; i < fragments.size(); i++) fragments[i].print(i);

	/*
	// print fclusters
	for(int i = 0; i < fclusters.size(); i++) fclusters[i].print(i);
	*/

	if(verbose <= 1) return 0;

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	// print bb.hits
	for(int i = 0; i < bb.hits.size(); i++) bb.hits[i].print();

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	printf("\n");

	return 0;
}

vector<int32_t> bundle_bridge::build_accumulate_length(const vector<int> &v)
{
	int32_t x = 0;
	vector<int32_t> acc;
	acc.resize(v.size());
	for(int i = 0; i < v.size(); i++)
	{
		x += regions[v[i]].rpos - regions[v[i]].lpos;
		acc[i] = x;
	}
	return acc;
}

vector<int32_t> bundle_bridge::get_aligned_intervals(fragment &fr)
{
	vector<int32_t> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<int32_t> v = get_splices(fr);
	if(v.size() >= 1 && fr.h1->pos >= v.front()) return vv;
	if(v.size() >= 1 && fr.h2->rpos <= v.back()) return vv;

	v.insert(v.begin(), fr.h1->pos);
	v.push_back(fr.h2->rpos);
	return v;
}

vector<int32_t> bundle_bridge::get_splices(fragment &fr)
{
	vector<int32_t> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<int> v = decode_vlist(fr.paths[0].v);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() - 1; i++)
	{
		int32_t pp = regions[v[i + 0]].rpos;
		int32_t qq = regions[v[i + 1]].lpos;
		if(pp >= qq) continue;
		vv.push_back(pp);
		vv.push_back(qq);
	}
	return vv;
}
