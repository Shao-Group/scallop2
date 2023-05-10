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
	build_supplementaries();
	set_chimeric_cigar_positions(); //setting h.first_pos/second_pos etc for getting back splice positions using cigars 
	build_junctions();
	extend_junctions();

	build_regions();

	align_hits_transcripts();
	index_references();

	printf("Called from bundle_bridge\n");
	build_fragments(); //builds fragment from h1p to h2
	printf("\n");

	build_circ_fragments(); //will build fragment from h2 to h1s, added by Tasfia

	//group_fragments();

	remove_tiny_boundaries();

	bridger bdg(this);
	bdg.bridge();

	extract_circ_fragment_pairs();
	//print_circ_fragment_pairs();
	join_circ_fragment_pairs();
	print_circRNAs();

	//printf("fragments vector size after = %zu\n",fragments.size());

	/*for(int i=0;i<fragments.size();i++)
	{
		if(fragments[i].frag_type == 2)
		{
			printf("frag pi = %d\n",fragments[i].pi);
			fragments[i].print(i+1);
		}
	}*/
	return 0;
}


int bundle_bridge::build_supplementaries()
{

	int max_index = bb.hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

    vector< vector<int> > vv;
    vv.resize(max_index);

    //printf("Bundle hit size: %d\n", bb.hits.size());
    // first build index
    for(int i = 0; i < bb.hits.size(); i++)
    {
        hit &h = bb.hits[i];

        // TODO
        if((h.flag & 0x800) == 0) continue;

        //printf("%s\n",h.qname.c_str());
        //printf("%zu\n",h.qhash);

        // do not use hi; as long as qname, pos and isize are identical
        // add 0x40 and 0x80
        int k = (h.qhash % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;
        vv[k].push_back(i);
        //printf("Adding supple\n");
    }

    //printf("End of vv adding\n");

    for(int i = 0; i < bb.hits.size(); i++)
    {

        hit &h = bb.hits[i];

        //if(h.paired == true) continue; 
        //if(h.isize <= 0) continue; //commented out as this was filtering non chimeric part of an end of SRR1721290.17627808
        //if(h.vlist.size() == 0) continue;
        if((h.flag & 0x800) >= 1) continue;       // skip supplemetary

        int k = (h.qhash % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;

        for(int j = 0; j < vv[k].size(); j++)
        {
            hit &z = bb.hits[vv[k][j]];
            //if(z.hi != h.hi) continue;
            //if(z.paired == true) continue;
            //if(z.pos != h.mpos) continue;
            //if(z.isize + h.isize != 0) continue;
            //if(z.qhash != h.qhash) continue;
            if(z.qname != h.qname) continue;
            // TODO check 0x40 and 0x80 are the same
            if(((z.flag & 0x40) != (h.flag & 0x40)) || ((z.flag & 0x80) != (h.flag & 0x80))) continue;

        	h.suppl = &z;
        	break;

            //Taking the first supplementary read
        }
    }

    return 0;
}

int bundle_bridge:: set_chimeric_cigar_positions()
{
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.suppl == NULL) continue;

		/*printf("Primary hit:\n");
		h.print();
		printf("Supple hit:\n");
		h.suppl->print();


		printf("cigar prim:");
		for(int p=0;p<h.n_cigar;p++)
		{
			printf("%d%c",h.cigar_vector[p].second,h.cigar_vector[p].first);
		}
		printf("\n");

		printf("cigar supp:");
		for(int p=0;p<h.suppl->n_cigar;p++)
		{
			printf("%d%c",h.suppl->cigar_vector[p].second,h.suppl->cigar_vector[p].first);
		}
		printf("\n");*/

		
		int32_t p;
		int32_t q;

		p = h.pos;
		int match_index = 0;;

		for(int j=0;j<h.n_cigar;j++)
		{
			if(h.cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			p -= h.cigar_vector[j].second; //subtracting cigars before match index
		}
		
		q = h.suppl->pos;
		match_index = 0;;

		for(int j=0;j<h.suppl->n_cigar;j++)
		{
			if(h.suppl->cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			q -= h.suppl->cigar_vector[j].second; //subtracting cigars before match index
		}
		

		int32_t x = p;
		int32_t diff_cigar1 = 10000;
		int32_t diff_cigar2 = 10000;
		int best_pos_flag = 0;

		for(int j=0;j<h.n_cigar-1;j++)
		{
			int32_t y = q;
			pair<char, int32_t> hp_cigar1 = h.cigar_vector[j];
			pair<char, int32_t> hp_cigar2 = h.cigar_vector[j+1];

			for(int k=0;k<h.suppl->n_cigar-1;k++)
			{

				pair<char, int32_t> hs_cigar1 = h.suppl->cigar_vector[k];
				pair<char, int32_t> hs_cigar2 = h.suppl->cigar_vector[k+1];

				//printf("check %d,%c\n",hp_cigar1.second,hp_cigar1.first);

				int32_t hp_cigar1_len = hp_cigar1.second;
				int32_t hp_cigar2_len = hp_cigar2.second;
				int32_t hs_cigar1_len = hs_cigar1.second;
				int32_t hs_cigar2_len = hs_cigar2.second;

				if(((hp_cigar1.first == 'S' || hp_cigar1.first == 'H') && hp_cigar2.first == 'M') && (hs_cigar1.first == 'M' && (hs_cigar2.first == 'S' || hs_cigar2.first == 'H')))
				{
					if(hp_cigar2.first == 'M' && j+3 < h.cigar_vector.size() && h.cigar_vector[j+3].first == 'M')
					{
						//&& h.cigar_vector[j+2].first == 'N' can be I or D as well
						for(int p=j+3;p<h.cigar_vector.size();p+=2) //traverse all Ms with 1 gap in the middle ex:SMNMNM/SIMIMIM/SMDMDMD
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar2_len += h.cigar_vector[p].second;
							}
						}			
					}

					if(hs_cigar1.first == 'M' && k-2 >= 0 && h.suppl->cigar_vector[k-2].first == 'M')
					{
						//&& h.suppl->cigar_vector[k-1].first == 'N' can be I or D as well
						for(int p=k-2;p>=0;p-=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar1_len += h.suppl->cigar_vector[p].second;
							}
						}
					}
					
					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);
				
				}
				else if(((hs_cigar1.first == 'S' || hs_cigar1.first == 'H') && hs_cigar2.first == 'M') && (hp_cigar1.first == 'M' && (hp_cigar2.first == 'S' || hp_cigar2.first == 'H')))
				{
					if(hs_cigar2.first == 'M' && k+3 < h.suppl->cigar_vector.size() && h.suppl->cigar_vector[k+3].first == 'M')
					{
						//&& h.suppl->cigar_vector[k+2].first == 'N' can be I or D as well
						for(int p=k+3;p<h.suppl->cigar_vector.size();p+=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar2_len += h.suppl->cigar_vector[p].second;
							}
						}			
					}

					if(hp_cigar1.first == 'M' && j-2 >=0 && h.cigar_vector[j-2].first == 'M')
					{
						//&& h.cigar_vector[j-1].first == 'N' can be I or D as well
						for(int p=j-2;p>=0;p-=2)
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar1_len += h.cigar_vector[p].second;
							}
						}
					}

					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);					
				}

				if(diff_cigar1 < 20 && diff_cigar2 < 20) //setting diff max 20 between complementing cigars
				{
					//printf("Found best positions\n");
					best_pos_flag = 1;

					//set first,second and third positions of prim and suppl
					h.first_pos = x;
					h.second_pos = x + hp_cigar1.second;
					h.third_pos = x + hp_cigar1.second + hp_cigar2.second;

					h.suppl->first_pos = y;
					h.suppl->second_pos = y + hs_cigar1.second;
					h.suppl->third_pos = y + hs_cigar1.second + hs_cigar2.second;

					h.left_cigar = hp_cigar1.first;
					h.left_cigar_len = hp_cigar1_len;
					h.right_cigar = hp_cigar2.first;
					h.right_cigar_len = hp_cigar2_len;

					h.suppl->left_cigar = hs_cigar1.first;
					h.suppl->left_cigar_len = hs_cigar1_len;
					h.suppl->right_cigar = hs_cigar2.first;
					h.suppl->right_cigar_len = hs_cigar2_len;

					break;

				}
				y += hs_cigar1.second;

			}
			if(best_pos_flag == 1) break;

			x += hp_cigar1.second;
		}

		//printf("set_cigar p:%d-%d-%d\n",h.first_pos,h.second_pos,h.third_pos);
		//printf("set_cigar s:%d-%d-%d\n",h.suppl->first_pos,h.suppl->second_pos,h.suppl->third_pos);
	}

	return 0;
}


int bundle_bridge::build_junctions()
{
	int min_max_boundary_quality = min_mapping_quality; //building a list of all splice pos and the hit index that includes the splice pos
	map< int64_t, vector<int> > m; // map of spos against vector of bundle base indices
	for(int i = 0; i < bb.hits.size(); i++)
	{
		vector<int64_t> v = bb.hits[i].spos;
		//printf("Spos size: %d\n", v.size());

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
	//printf("spos map size = %d\n",m.size());

	junctions.clear();
	junc_map.clear();
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

		if(junc_map.find(it->first) != junc_map.end()) continue;
		else junc_map.insert(make_pair(it->first, jc.strand));
	}
	//printf("Junctions size: %d\n", junctions.size());
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

	//for(int i=0;i<junctions.size();i++)
	//{

	//	junctions[i].print();
	//	printf("\n");
	//}

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
		//remove_tiny_boundary(bb.hits[i]);
	}

	ref_phase.resize(ref_trsts.size());
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		align_transcript(m, ref_trsts[i], ref_phase[i]);
	}

	return 0;
}

int bundle_bridge::remove_tiny_boundaries()
{
	for(int i = 0; i < bb.hits.size(); i++)
	{
		remove_tiny_boundary(bb.hits[i]);
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

int bundle_bridge::remove_tiny_boundary(hit &h1)
{
	vector<int> v1 = decode_vlist(h1.vlist);
	int n1 = v1.size();
	if(n1 >= 2 && v1[n1 - 2] + 1 == v1[n1 - 1])
	{
		int k = v1[n1 - 1];
		int32_t total = regions[k].rpos - regions[k].lpos;
		int32_t flank = h1.rpos - regions[k].lpos;

		if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
		{
			vector<int> v(v1.begin(), v1.begin() + n1 - 1);
			assert(v.size() + 1 == v1.size());
			h1.vlist = encode_vlist(v);
			h1.rpos = regions[k].lpos;
		}
	}

	vector<int> v2 = decode_vlist(h1.vlist);
	int n2 = v2.size();
	if(n2 >= 2 && v2[0] + 1 == v2[1])
	{
		int k = v2[0];
		int32_t total = regions[k].rpos - regions[k].lpos;
		int32_t flank = regions[k].rpos - h1.pos;

		if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
		{
			vector<int> v(v2.begin() + 1, v2.end());
			assert(v.size() + 1 == v2.size());
			h1.vlist = encode_vlist(v);
			h1.pos = regions[k].rpos;
		}
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
	vv.resize(max_index); //max_index slots initialized to zero, here max_index is the max hash index

	// first build index
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.isize >= 0) continue;
		if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.qhash % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;

		if(strcmp(h.qname.c_str(),"simulate:311116") == 0)
		{
			printf("hash negative: %d\n",k);
		}
		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i); //vv containes hits of the same hash
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.paired == true) continue;
		if(h.isize <= 0) continue;
		if(h.vlist.size() == 0) continue;

		int k = (h.qhash % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		if(strcmp(h.qname.c_str(),"simulate:311116") == 0)
		{
			printf("hash positive: %d\n",k);
		}

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		/*if(strcmp(h.qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in hits\n");
			printf("isize: %d\n",h.isize);
			printf("vlist size: %zu\n",h.vlist.size());
		}*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];

			if(strcmp(z.qname.c_str(),"simulate:311116") == 0)
			{
				printf("simulate:311116 is in hits\n");
			}

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

		fragment fr(&bb.hits[i], &bb.hits[x]); //h2 and h1s as param or h2s and h1 as parameter

		/*if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in fragments\n");
		}*/

		fr.frag_type = 1; //this is the first set of fragment,tasfia

		//keep it
		// ===============================
		// TODO: dit for UMI
		bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index
		bb.hits[x].pi = i;
		bb.hits[i].fidx = fragments.size();//check if used somewhere
		bb.hits[x].fidx = fragments.size();//fidx is the fragment index
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
		//keep it

		//inlcude
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
	
	/*for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in build_fragments\n");
		}
	}*/

	// by shao, exit here (no UMI)
	return 0;

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
			fr.frag_type = 1;
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


int bundle_bridge::build_circ_fragments()
{
	//parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	vector<fragment> circ_fragments;
	circ_fragments.clear(); // unbridged second set of fragments, cannot use this for finding partner fragments

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in build_circ_fragments\n");
		}

		int is_compatible = 0; //1 for h1 has a suppl and compatible, 2 for h2 has a suppl and compatible

		if(fr.h1->suppl != NULL)
		{
			//printf("supple not null\n");
			//need to check compatibility from bundle.cc

			h1_supp_count++;
			hit *h1_supple = fr.h1->suppl;

			printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h1 has a supple hit.\n");
			printf("Primary: ");
			fr.h1->print();
			printf("Supple: ");
			h1_supple->print();
			printf("fr.h2: ");
			fr.h2->print();

			printf("cigar prim:");
			for(int p=0;p<fr.h1->n_cigar;p++)
			{
				printf("%d%c",fr.h1->cigar_vector[p].second,fr.h1->cigar_vector[p].first);
			}
			printf("\n");

			printf("cigar supp:");
			for(int p=0;p<fr.h1->suppl->n_cigar;p++)
			{
				printf("%d%c",fr.h1->suppl->cigar_vector[p].second,fr.h1->suppl->cigar_vector[p].first);
			}
			printf("\n");

			printf("set_cigar p:%d-%d-%d\n",fr.h1->first_pos,fr.h1->second_pos,fr.h1->third_pos);
			printf("set_cigar s:%d-%d-%d\n",fr.h1->suppl->first_pos,fr.h1->suppl->second_pos,fr.h1->suppl->third_pos);

			if(fr.h1->first_pos == 0 || fr.h1->suppl->first_pos == 0)
			{
				string combo = "special case h1s";
				printf("%s\n",combo.c_str());
				if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
				else frag2graph_freq[combo] += 1;
				continue;
			}

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(fr.h1->left_cigar == 'H' || fr.h1->left_cigar == 'S')
			{
				len_HS += fr.h1->left_cigar_len;
			}
			else if(fr.h1->right_cigar == 'H' || fr.h1->right_cigar == 'S')
			{
				len_HS += fr.h1->right_cigar_len;
			}
			if(fr.h1->suppl->left_cigar == 'H' || fr.h1->suppl->left_cigar == 'S')
			{
				len_HS += fr.h1->suppl->left_cigar_len;
			}
			else if(fr.h1->suppl->right_cigar == 'H' || fr.h1->suppl->right_cigar == 'S')
			{
				len_HS += fr.h1->suppl->right_cigar_len;
			}

			//printf("len_HS = %d\n",len_HS);

			if(abs(len_HS - 100) > 5) //here 100 is the estimated read length, replace this with any related exisiting parameter
			{
				printf("read length criteria unsatisfied h1s.\n");
				continue;
			}


			//examples give a rule - p and s has to be edges

			//this case should not occur as h1p should always be on the left of h2
			if(h1_supple->pos <= fr.h2->pos && h1_supple->pos <= fr.h1->pos && fr.h1->rpos >= h1_supple->rpos && fr.h1->rpos >= fr.h2->rpos)
			{	
				//printf("Compatible in previous definition, h1s leftmost, h1p rightmost\n");
				
				if(h1_supple->second_pos <= fr.h2->pos && h1_supple->second_pos <= fr.h1->pos && fr.h1->second_pos >= h1_supple->rpos && fr.h1->second_pos >= fr.h2->rpos)
				{
					string combo = "Compatible h1s leftmost h1p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;				
				}
				else
				{
					string combo = "Not compatible h1s leftmost h1p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;				
				}

			}

			else if(fr.h1->pos <= fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos >= fr.h2->rpos)
			{
				//printf("Compatible in previous definition, h1p leftmost, h1s rightmost\n");

				if(fr.h1->second_pos <= fr.h2->pos && fr.h1->second_pos <= h1_supple->pos && h1_supple->second_pos >= fr.h1->rpos && h1_supple->second_pos >= fr.h2->rpos)
				{
					string combo = "Compatible h1p leftmost h1s rightmost";
					printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					is_compatible = 1;
				}
				else
				{
					string combo = "Not compatible h1p leftmost h1s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;	
					continue;				
				}

			}
			else
			{
				//printf("Not compatible in previous definition\n");
				if(fr.h1->pos > fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos >= fr.h2->rpos && fr.h1->pos - fr.h2->pos <= 5)
				{
					string combo = "alignment error-h1p_pos>h2_pos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else if(fr.h1->pos <= fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos < fr.h2->rpos && fr.h2->rpos - h1_supple->rpos <= 5)
				{
					string combo = "alignment error-h2_rpos>h1s_rpos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else
				{
					string combo = "Not compatible h1s";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;
				}
			}
		}

		if(fr.h2->suppl != NULL)
		{
			//need to check compatibility from bundle.

			h2_supp_count++;
			hit *h2_supple = fr.h2->suppl;
			
			printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h2 has a supple hit.\n");
			printf("h1: ");
			fr.h1->print();
			printf("Primary: ");
			fr.h2->print();
			printf("Supple: ");
			h2_supple->print();

			printf("cigar prim:");
			for(int p=0;p<fr.h2->n_cigar;p++)
			{
				printf("%d%c",fr.h2->cigar_vector[p].second,fr.h2->cigar_vector[p].first);
			}
			printf("\n");

			printf("cigar supp:");
			for(int p=0;p<fr.h2->suppl->n_cigar;p++)
			{
				printf("%d%c",fr.h2->suppl->cigar_vector[p].second,fr.h2->suppl->cigar_vector[p].first);
			}
			printf("\n");

			printf("set_cigar p:%d-%d-%d\n",fr.h2->first_pos,fr.h2->second_pos,fr.h2->third_pos);
			printf("set_cigar s:%d-%d-%d\n",fr.h2->suppl->first_pos,fr.h2->suppl->second_pos,fr.h2->suppl->third_pos);


			if(fr.h2->first_pos == 0 || fr.h2->suppl->first_pos == 0)
			{

				string combo = "special case h2s";
				printf("%s\n",combo.c_str());
				if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
				else frag2graph_freq[combo] += 1;
				continue;
			}

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(fr.h2->left_cigar == 'H' || fr.h2->left_cigar == 'S')
			{
				len_HS += fr.h2->left_cigar_len;
			}
			else if(fr.h2->right_cigar == 'H' || fr.h2->right_cigar == 'S')
			{
				len_HS += fr.h2->right_cigar_len;
			}
			if(fr.h2->suppl->left_cigar == 'H' || fr.h2->suppl->left_cigar == 'S')
			{
				len_HS += fr.h2->suppl->left_cigar_len;
			}
			else if(fr.h2->suppl->right_cigar == 'H' || fr.h2->suppl->right_cigar == 'S')
			{
				len_HS += fr.h2->suppl->right_cigar_len;
			}

			//printf("len_HS = %d\n",len_HS);

			if(abs(len_HS - 100) > 5) //here 100 is the estimated read length, replace this with any related exisiting parameter
			{
				printf("read length criteria unsatisfied h2s.\n");
				continue;
			}

			//general rule p and s has to be edges
			if(h2_supple->pos <= fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos >= fr.h1->rpos)
			{
				//printf("Compatible in previous definition, h2s leftmost, h2p rightmost\n");

				if(h2_supple->second_pos <= fr.h1->pos && h2_supple->second_pos <= fr.h2->pos && fr.h2->second_pos >= h2_supple->rpos && fr.h2->second_pos >= fr.h1->rpos)
				{
					string combo = "Compatible h2s leftmost h2p rightmost";
					printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					is_compatible = 2;					
				}
				else
				{
					string combo = "Not compatible h2s leftmost h2p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;			
				}

			}

			//this case should not occur as h2 p should always be on the right of h1
			else if(fr.h2->pos <= fr.h1->pos && fr.h2->pos <= h2_supple->pos && h2_supple->rpos >= fr.h2->rpos && h2_supple->rpos >= fr.h1->rpos)
			{
				//printf("Compatible in previous definition, h2p leftmost, h2s rightmost\n");

				if(fr.h2->second_pos <= fr.h1->pos && fr.h2->second_pos <= h2_supple->pos && h2_supple->second_pos >= fr.h2->rpos && h2_supple->second_pos >= fr.h1->rpos)
				{
					string combo = "Compatible h2p leftmost h2s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;					
				}
				else
				{
					string combo = "Not compatible h2p leftmost h2s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;	
					continue;				
				}

			}

			else
			{
				if(h2_supple->pos <= fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos < fr.h1->rpos && fr.h1->rpos - fr.h2->rpos <= 5)
				{
					string combo = "alignment error-h1_rpos>h2p_rpos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else if(h2_supple->pos > fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos >= fr.h1->rpos && h2_supple->pos - fr.h1->pos <= 5)
				{
					string combo = "alignment error-h2s_pos>h1_pos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else
				{
					//printf("Not compatible in previous definition\n");

					string combo = "Not compatible h2s";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;
				}
			}

			//if not compatible, continued
			
			//if compatible
			//fragment fr(fr.h2->suppl, fr.h1); //h2 and h1s as param or h2s and h1 as parameter
			//fr.frag_type = 2; //this is the second set of fragment

			//do we need to check h2s-h1-h2p order if compatible checked? no as ordering already checked, see above
			//need to handle fr.h1 paired - both fr.h1 and h2 needs to paired before when we build second fragment, see below
		}

		//if compatible h1s
		if(is_compatible == 1)
		{
			//if(bb.hits[x].vlist.size() == 0) continue;
			//printf("vlist size = %zu",fr.h2->vlist.size());
			//printf("vlist size = %zu",fr.h1->suppl->vlist.size());

			//fr.h2 and fr.h1 needs to be paired by build_fragments()
			if(fr.h2->paired != true || fr.h1->paired != true) continue; //first set of fragment needs to be paired

			fr.h1->suppl->paired = true; //setting supple paired true to avoid assertion later in build hyper set
			fragment frag(fr.h2, fr.h1->suppl);
			frag.frag_type = 2;
			frag.is_compatible = 1;

			/*// ===============================
			// TODO: dit for UMI
			bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index, not used here
			bb.hits[x].pi = i;
			bb.hits[i].fidx = fragments.size();//check if used somewhere
			bb.hits[x].fidx = fragments.size();//fidx is the fragment index
			ctp += 1;
			fr.type = 0; //needed here?
			// ================================*/

			frag.type = 0;
			frag.lpos = fr.h2->pos;
			frag.rpos = fr.h1->suppl->rpos;

			vector<int> v1 = decode_vlist(fr.h2->vlist);
			vector<int> v2 = decode_vlist(fr.h1->suppl->vlist);
			frag.k1l = frag.h1->pos - regions[v1.front()].lpos;
			frag.k1r = regions[v1.back()].rpos - frag.h1->rpos; //frag.h1?
			frag.k2l = frag.h2->pos - regions[v2.front()].lpos;
			frag.k2r = regions[v2.back()].rpos - frag.h2->rpos; //frag.h2?
			//keep it

			//inlcude
			frag.b1 = true;
			if(v1.size() <= 1) 
			{
				frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + frag.h1->nm) frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + frag.h1->nm) frag.b1 = false;
			}

			frag.b2 = true;
			if(v2.size() <= 1)
			{
				frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment1 + frag.h2->nm) frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment2 + frag.h2->nm) frag.b2 = false;
			}

			frag.pi = k;
			frag.fidx = fragments.size() + circ_fragments.size();
			fr.pi = fragments.size() + circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
			fr.fidx = k;

			circ_fragments.push_back(frag);
			
			//fr.h2->paired = true;
			//fr.h1->suppl->paired = true;

		}
		else if(is_compatible == 2)
		{

			//fr.h2 and fr.h1 needs to be paired by build_fragments()
			if(fr.h2->paired != true || fr.h1->paired != true) //first set of fragment needs to be paired
			{
				printf("fr.h2 in first set not paired\n");
				fr.h2->print();
				continue;
			}							
			fr.h2->suppl->paired = true; //setting supple paired true to avoid assertion later in build hyper set
			fragment frag(fr.h2->suppl, fr.h1);
			frag.frag_type = 2;
			frag.is_compatible = 2;

			/*// ===============================
			// TODO: dit for UMI
			bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index, not used here
			bb.hits[x].pi = i;
			bb.hits[i].fidx = fragments.size();//check if used somewhere
			bb.hits[x].fidx = fragments.size();//fidx is the fragment index
			ctp += 1;
			// ================================*/

			frag.type = 0;
			frag.lpos = fr.h2->suppl->pos;
			frag.rpos = fr.h1->rpos;

			vector<int> v1 = decode_vlist(fr.h2->suppl->vlist);
			vector<int> v2 = decode_vlist(fr.h1->vlist);
			frag.k1l = frag.h1->pos - regions[v1.front()].lpos;
			frag.k1r = regions[v1.back()].rpos - frag.h1->rpos;
			frag.k2l = frag.h2->pos - regions[v2.front()].lpos;
			frag.k2r = regions[v2.back()].rpos - frag.h2->rpos;
			//keep it

			//inlcude
			frag.b1 = true;
			if(v1.size() <= 1) 
			{
				frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + frag.h1->nm) frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + frag.h1->nm) frag.b1 = false;
			}

			frag.b2 = true;
			if(v2.size() <= 1)
			{
				frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment1 + frag.h2->nm) frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment2 + frag.h2->nm) frag.b2 = false;
			}

			frag.pi = k;
			frag.fidx = fragments.size() + circ_fragments.size();
			fr.pi = fragments.size() + circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
			fr.fidx = k;

			circ_fragments.push_back(frag);	

		}
	}



	//adding the second part fragments to fragments vector
	if(circ_fragments.size() > 0)
	{
		printf("fragments vector size before = %zu\n",fragments.size());		
	}
	
	//fragments.clear();
	for(int i=0;i<circ_fragments.size();i++)
	{
		fragments.push_back(circ_fragments[i]); //pushing the circ fragments at the end of the main fragments list, 
												//does this create any problem in the main scallop2 output??
												//no as this is a separate module just for circRNA
	}

	if(circ_fragments.size() > 0)
	{
		printf("circ fragment vector size = %zu\n",circ_fragments.size());
	}

	//check fragments order before bridging, same
	/*if(circ_fragments.size()>0)
	{
		printf("fragments before bridging:\n");
		for(int i=0;i<fragments.size();i++)
		{
			printf("%d-%d,",fragments[i].h1->hid,fragments[i].h2->hid);
		}

		printf("circ fragments:\n");

		for(int i=0;i<circ_fragments.size();i++)
		{
			printf("%d-%d,",circ_fragments[i].h1->hid,circ_fragments[i].h2->hid);
		}
	}*/

	return 0;
}

int bundle_bridge::extract_circ_fragment_pairs()
{
	//check fragments order after bridging, same
	/*printf("fragments after bridging:\n");
	for(int i=0;i<fragments.size();i++)
	{
		printf("%d-%d,",fragments[i].h1->hid,fragments[i].h2->hid);
	}*/


	vector<fragment> circ_fragments;
	circ_fragments.clear();
	circ_fragment_pairs.clear();

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		if(fr.frag_type == 2) 
		{
			circ_fragments.push_back(fr); //bridging function called already
		}
	}

	if(circ_fragments.size() > 0)
	{
		printf("After bridging fragments vector size = %zu\n\n",fragments.size());		
	}


	for(int i = 0; i < fragments.size(); i++)
	{
		fragment &fr1 = fragments[i]; 

		if(fr1.frag_type == 2) continue; //filter and take original fragments of scallop2 

		for(int j=0;j<circ_fragments.size();j++)
		{
			fragment &fr2 = circ_fragments[j];

			if(strcmp(fr1.h1->qname.c_str(),fr2.h1->qname.c_str()) != 0) continue;
			if(fr1.pi == -1 || fr2.pi == -1 || fr1.fidx == -1 || fr2.fidx == -1) continue;

			if(fr1.pi == fr2.fidx && fr2.pi == fr1.fidx)
			{
				//printf("Found partner fragment\n");
				circ_fragment_pairs.push_back(pair<fragment,fragment>(fr1,fr2)); //fr1 is original fragment of scallop2 regardless of which of h1 or h2 has a supplement
			}
		}
	}

	/*if(circ_fragment_pairs.size() > 0)
	{
		printf("Printing bridged fragment pairs: size = %zu\n\n",circ_fragment_pairs.size());

		for(int i=0;i<circ_fragment_pairs.size();i++)
		{
			circ_fragment_pairs[i].first.print(i+1);
			circ_fragment_pairs[i].second.print(i+1);
			printf("\n");
		}
	}*/

	//profiling fragment hits bridged or not
	for(int i=0;i<circ_fragment_pairs.size();i++)
	{
		fragment &fr1 = circ_fragment_pairs[i].first;
		fragment &fr2 = circ_fragment_pairs[i].second;

		//printf("pairs: is_compatible = %d-%d, %c%c-%c%c\n",fr1.is_compatible, fr2.is_compatible,fr1.h1->bridged == true ? 'T' : 'F',fr1.h2->bridged == true ? 'T' : 'F',fr2.h1->bridged == true ? 'T' : 'F',fr2.h2->bridged == true ? 'T' : 'F');

		//printf("chrm = %s\n",bb.chrm.c_str());

		//bridging info for fr1 hits
		if(fr1.h1->bridged == false && fr1.h2->bridged == false)
		{
			string combo = "fr1-FF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == false && fr1.h2->bridged == true)
		{
			string combo = "fr1-FT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == true && fr1.h2->bridged == false)
		{
			string combo = "fr1-TF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == true && fr1.h2->bridged == true)
		{
			string combo = "fr1-TT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}


		/*printf("Printing first fragment:\n");
		fr1.print(i+1);
		printf("Printing second fragment:\n");
		fr2.print(i+1);*/

		//bridging info for fr2 hits
		if(fr2.h1->bridged == false && fr2.h2->bridged == false)
		{
			string combo = "fr2-FF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == false && fr2.h2->bridged == true)
		{
			string combo = "fr2-FT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == true && fr2.h2->bridged == false)
		{
			string combo = "fr2-TF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == true && fr2.h2->bridged == true)
		{
			string combo = "fr2-TT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}

		//printf("\n");

		//counting path.size = 1
		if(fr1.paths.size() == 1)
		{
			string combo = "fr1-path.size = 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		else
		{
			string combo = "fr1-path.size != 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		if(fr2.paths.size() == 1)
		{
			string combo = "fr2-path.size = 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		else
		{
			string combo = "fr2-path.size != 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;		
		}
	}

	return 0;
}

int bundle_bridge::print_circ_fragment_pairs()
{

	if(circ_fragment_pairs.size() > 0)
	{
		printf("Printing fragment pair paths:\n");
		for(int i=0;i<circ_fragment_pairs.size();i++)
		{
			fragment &fr1 = circ_fragment_pairs[i].first;
			fragment &fr2 = circ_fragment_pairs[i].second;

			if(fr1.paths.size() != 1 || fr2.paths.size() != 1) continue; //check both frags bridged
			//if(fr1.paths[0].type != 1 || fr2.paths[0].type != 1) continue; // 1: within normal range of insertsize;

			fr1.paths[0].print(i+1);
			fr2.paths[0].print(i+1);
		}
	}
	return 0;
}

int bundle_bridge::join_circ_fragment_pairs()
{
	circ_trsts.clear(); // emptying before storing joined circRNAs
	for(int i=0;i<circ_fragment_pairs.size();i++)
	{
		fragment &fr1 = circ_fragment_pairs[i].first;
		fragment &fr2 = circ_fragment_pairs[i].second;

		if(fr1.paths.size() != 1 || fr2.paths.size() != 1) continue; //not bridged
		//if(fr1.paths[0].type != 1 || fr2.paths[0].type != 1) continue; //insert size not normal

		printf("Printing separate fragments:\n");

		printf("\nchrm = %s\n",bb.chrm.c_str());

		fr1.print(i+1);
		fr2.print(i+1);

		join_circ_fragment_pair(circ_fragment_pairs[i],0,0);		
	}

	return 0;
}

int bundle_bridge::join_circ_fragment_pair(pair<fragment,fragment> &fr_pair, int ex1, int ex2)
{

	fragment &fr1 = fr_pair.first;
	fragment &fr2 = fr_pair.second;

	vector<int> v1 = decode_vlist(fr1.paths[0].v);
	vector<int> v2 = decode_vlist(fr2.paths[0].v);
	path p;

	assert(v1.size() > ex1);
	assert(v2.size() > ex2);

	if(fr2.is_compatible == 1)
	{
		vector<int>::iterator t1 = v1.end() - ex1;
		vector<int>::iterator t2 = v2.begin() + ex2;

		int x1 = v1[v1.size() - 1 - ex1];
		int x2 = v2[ex2];

		if(x1 < x2) return 0;

		/*if(strcmp(fr2.h1->qname.c_str(),"simulate:2195") == 0)
		{
			printf("simulate:2195 is in fr2.is_compatible 1\n");
		}*/

		vector<int>::iterator it = find(t2, v2.end(), x1);
		if(it == v2.end()) return 0;

		vector<int>::iterator j1, j2;
		for(j1 = t1 - 1, j2 = it; j1 >= v1.begin() && j2 >= t2; j1--, j2--)
		{
			if((*j1) != (*j2)) return 0;
		}

		p.ex1 = ex1;
		p.ex2 = ex2;
		p.v.insert(p.v.end(), v1.begin(), t1);
		p.v.insert(p.v.end(), it + 1, v2.end());
		//p.v = encode_vlist(p.v);

		printf("Printing merged path in fr2 iscomp 1:\n");
		printv(p.v);
		printf("\n\n");
	
		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr1.lpos) + "|" + tostring(fr2.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr1.lpos;
		int32_t end = fr2.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), p.v.begin(), p.v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr1.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);
		
		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
			//circ.circRNA_id = circ.circRNA_id + tostring(circ.circ_path[i]) + "|";
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}

		circ_trsts.push_back(circ);
	}
	else if(fr2.is_compatible == 2)
	{
		vector<int>::iterator t1 = v2.end() - ex1;
		vector<int>::iterator t2 = v1.begin() + ex2;

		int x1 = v2[v2.size() - 1 - ex1];
		int x2 = v1[ex2];

		if(x1 < x2) return 0;

		vector<int>::iterator it = find(t2, v1.end(), x1);
		if(it == v1.end()) return 0;

		vector<int>::iterator j1, j2;
		for(j1 = t1 - 1, j2 = it; j1 >= v2.begin() && j2 >= t2; j1--, j2--)
		{
			if((*j1) != (*j2)) return 0;
		}

		p.ex1 = ex1;
		p.ex2 = ex2;
		p.v.insert(p.v.end(), v2.begin(), t1);
		p.v.insert(p.v.end(), it + 1, v1.end());
		//p.v = encode_vlist(p.v);

		printf("Printing merged path in fr2 iscomp 2:\n");
		printv(p.v);
		printf("\n\n");
		
		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr2.lpos) + "|" + tostring(fr1.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr2.lpos;
		int32_t end = fr1.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), p.v.begin(), p.v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr1.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);

		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
			//circ.circRNA_id = circ.circRNA_id + tostring(circ.circ_path[i]) + "|";
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}		
		
		circ_trsts.push_back(circ);
	}
	else
	{
		printf("is_compatible not 1 or 2\n");
	}

	return 0;
}

char bundle_bridge::infer_circ_strand(const vector<int> &p)
{
	if(p.size() <= 1) return '.';
	int n1 = 0, n2 = 0, n3 = 0;
	for(int k = 0; k < p.size() - 1; k++)
	{
		int32_t x1 = regions[p[k + 0]].rpos;
		int32_t x2 = regions[p[k + 1]].lpos;
		if(x1 >= x2) continue;
		int64_t jc = pack(x1, x2);
		if(junc_map.find(jc) == junc_map.end()) continue;
		char c = junc_map[jc];
		if(c == '.') n1 ++;
		if(c == '+') n2 ++;
		if(c == '-') n3 ++;
	}
	if(n2 >= 1 && n1 <= 0 && n3 <= 0) return '+';
	if(n3 >= 1 && n1 <= 0 && n2 <= 0) return '-';
	printf("MIXED JUNCTION: ./+/- = %d/%d/%d\n", n1, n2, n3);
	return '.';
}

int bundle_bridge::print_circRNAs()
{
	if(circ_trsts.size() > 0)
	{
		//printf("circ_trsts size: %lu", circ_trsts.size());
		printf("Printing circRNAs from current bundle bridge:\n");

		for(int i=0;i<circ_trsts.size();i++)
		{
			circ_trsts[i].print(i+1);
		}

		printf("\n");
	}

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
