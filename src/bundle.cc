/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle.h"
#include "bundle_bridge.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"

bundle::bundle(bundle_base &b)
	: bb(b) //, br(b)
{
	// br.build();
	prepare();
}

bundle::~bundle()
{}


int bundle::prepare()
{
	compute_strand();
	build_intervals();
	// print_fmap();
	build_junctions();
	build_regions();
	build_partial_exons();

	build_partial_exon_map();
	link_partial_exons();
	return 0;
}

int bundle::build(int mode, bool revise)
{
	build_splice_graph(mode);
	if(revise == true) revise_splice_graph();
	// printf("rebuild splice graph started ....");
	rebuild_splice_graph_using_refined_hyper_set(mode);

	// refine_splice_graph();
	refine_modified_splice_graph();
	build_hyper_set();
	// printf("rebuild splice graph completed ...");
	//refine_hyper_set();
	return 0;
}

int bundle::compute_strand()
{
	if(library_type != UNSTRANDED) assert(bb.strand != '.');
	if(library_type != UNSTRANDED) return 0;

	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	if(np > nq) bb.strand = '+';
	else if(np < nq) bb.strand = '-';
	else bb.strand = '.';

	return 0;
}

int bundle::build_intervals()
{
	fmap.clear();

	// TODO: to get rid of fragments
	// simply comment out the following lines to 93
	// for(int i = 0; i < br.fragments.size(); i++)
	// {
	// 	fragment &fr = br.fragments[i];
	// 	if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
	// 	vector<int32_t> vv = br.get_aligned_intervals(fr);
	// 	if(vv.size() <= 0) continue;
	// 	assert(vv.size() % 2 == 0);

	// 	for(int k = 0; k < vv.size() / 2; k++)
	// 	{
	// 		int32_t p = vv[2 * k + 0];
	// 		int32_t q = vv[2 * k + 1];
	// 		fmap += make_pair(ROI(p, q), 1);
	// 	}
	// }

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &ht = bb.hits[i];

		// instead of using itvm, we use the splicing positions 
		// to build the interval map of this bundle

		int32_t p1 = ht.pos;
		for(int k = 0; k < ht.spos.size(); k++)
		{
			int32_t s = high32(ht.spos[k]);
			int32_t t = low32(ht.spos[k]);

			fmap += make_pair(ROI(p1, s), 1);
			// printf("qname: %s insert to fmap: %d -- %d\n",ht.qname.c_str(), p1,s);
			p1 = t;
			// print_fmap();
		}
		fmap += make_pair(ROI(p1, ht.rpos), 1);
		// printf("insert to fmap: %d -- %d\n", p1,ht.rpos);


		// TODO: comment out the following 3 lines
		// if(ht.bridged == true) continue;
		// if((ht.flag & 0x100) >= 1) continue;
		// if(br.breads.find(ht.qname) != br.breads.end()) continue;

		/*
		for(int k = 0; k < ht.itvm.size(); k++)
		{
			int32_t s = high32(ht.itvm[k]);
			int32_t t = low32(ht.itvm[k]);
			fmap += make_pair(ROI(s, t), 1);
		}
		*/
	}

	return 0;
}

int bundle::build_junctions()
{
	map<int64_t, vector<hit*>> m;		// bridged fragments
	// TODO: comment out the following paragraph
	// for(int i = 0; i < br.fragments.size(); i++)
	// {
	// 	fragment &fr = br.fragments[i];
	// 	if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
	// 	vector<int32_t> vv = br.get_splices(fr);
	// 	if(vv.size() <= 0) continue;
	// 	assert(vv.size() % 2 == 0);

	// 	for(int k = 0; k < vv.size() / 2; k++)
	// 	{
	// 		int64_t p = pack(vv[2 * k + 0], vv[2 * k + 1]);

	// 		if(m.find(p) == m.end())
	// 		{
	// 			vector<hit*> hv;
	// 			hv.push_back(fr.h1);
	// 			//hv.push_back(fr.h2);
	// 			m.insert(pair< int64_t, vector<hit*> >(p, hv));
	// 		}
	// 		else
	// 		{
	// 			m[p].push_back(fr.h1);
	// 			//m[p].push_back(fr.h2);
	// 		}
	// 	}
	// }

	for(int i = 0; i < bb.hits.size(); i++)
	{
		// TODO: comment out the follwoing 3
		// if(bb.hits[i].bridged == true) continue;
		// if((bb.hits[i].flag & 0x100) >= 1) continue;
		// if(br.breads.find(bb.hits[i].qname) != br.breads.end()) continue;

		vector<int64_t> v = bb.hits[i].spos;
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
			if(m.find(p) == m.end())
			{
				vector<hit*> hv;
				hv.push_back(&(bb.hits[i]));
				m.insert(pair< int64_t, vector<hit*> >(p, hv));
			}
			else
			{
				m[p].push_back(&(bb.hits[i]));
			}
		}
	}

	map<int64_t, vector<hit*>>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<hit*> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		int32_t p1 = high32(it->first);
		int32_t p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			nm += v[k]->nm;
			if(v[k]->xs == '.') s0++;
			if(v[k]->xs == '+') s1++;
			if(v[k]->xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", bb.chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		// TODO: be careful about the junction strand
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

int bundle::build_regions()
{
	MPI s;
	s.insert(PI(bb.lpos, START_BOUNDARY));
	s.insert(PI(bb.rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		double ave, dev, max;
		evaluate_rectangle(fmap, jc.lpos, jc.rpos, ave, dev, max);		

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	assert(pexons.size() == 0);
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		if(s.find(p.lpos) != s.end()) s.insert(PI(p.lpos, p.ltype));
		if(s.find(p.rpos) != s.end()) s.insert(PI(p.rpos, p.rtype));
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

		// printf("Build region: %d - %d\n", l, r);

		regions.push_back(region(l, r, ltype, rtype, &fmap, &(bb.imap)));
	}
	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	regional.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		for(int k = 0; k < r.pexons.size(); k++)
		{
			partial_exon &pe = r.pexons[k];
			pe.rid = i;
			pe.pid = pexons.size();
			pexons.push_back(pe);
			
			if((pe.lpos != bb.lpos || pe.rpos != bb.rpos) && pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) regional.push_back(true);
			else regional.push_back(false);
		}
	}
	return 0;
}

int bundle::build_partial_exon_map()
{
	pmap.clear();
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		pmap += make_pair(ROI(p.lpos, p.rpos), i + 1);
	}
	return 0;
}

int bundle::locate_left_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x, x + 1));
	if(it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p2 >= x);
	assert(p1 <= x);

	if(x - p1 > min_flank_length && p2 - x < min_flank_length) k++;

	if(k >= pexons.size()) return -1;
	return k;
}

int bundle::locate_right_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x - 1, x));
	if(it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p1 < x);
	assert(p2 >= x);

	if(p2 - x > min_flank_length && x - p1 <= min_flank_length) k--;
	return k;
}

vector<int> bundle::align_hit(hit &h)
{
	bool b = true;
	vector<int> sp2;
	vector<int64_t> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return sp2;

	for(int k = 0; k < v.size(); k++)
	{
		int32_t p1 = high32(v[k]);
		int32_t p2 = low32(v[k]);

		int k1 = locate_left_partial_exon(p1);
		int k2 = locate_right_partial_exon(p2);
		if(k1 < 0 || k2 < 0) b = false;
		if(b == false) break;

		for(int j = k1; j <= k2; j++) sp2.push_back(j);
	}

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}

vector<int> bundle::align_fragment(fragment &fr)
{
	bool b = true;
	vector<int> sp2;
	//if(fr.paths.size() != 1 || fr.paths[0].type != 1) return sp2;
	// vector<int32_t> v = br.get_aligned_intervals(fr);
	// assert(v.size() % 2 == 0);
	// if(v.size() == 0) return sp2;

	// for(int k = 0; k < v.size() / 2; k++)
	// {
	// 	int32_t p1 = v[2 * k + 0];
	// 	int32_t p2 = v[2 * k + 1];
	// 	int k1 = locate_left_partial_exon(p1);
	// 	int k2 = locate_right_partial_exon(p2);
	// 	if(k1 < 0 || k2 < 0) b = false;
	// 	if(b == false) break;
	// 	for(int j = k1; j <= k2; j++) sp2.push_back(j);
	// }

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}

int bundle::link_partial_exons()
{
	if(pexons.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for(int i = 0; i < pexons.size(); i++)
	{
		int32_t l = pexons[i].lpos;
		int32_t r = pexons[i].rpos;

		assert(lm.find(l) == lm.end());
		assert(rm.find(r) == rm.end());
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];
		MPI::iterator li = rm.find(b.lpos);
		MPI::iterator ri = lm.find(b.rpos);

		if(ri == lm.end()) printf("ri-lm: l=%d, r=%d\n", b.lpos, b.rpos);
		// assert(li != rm.end());
		// assert(ri != lm.end());

		if(li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else
		{
			b.lexon = b.rexon = -1;
		}
	}
	return 0;
}

int bundle::build_splice_graph(int mode)
{
	gr.clear();

	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = bb.lpos;
	vi0.rpos = bb.lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		if (length < reliability_threshold) pexons[i].rel = false;
		else pexons[i].rel = true; //TODO: add else later
		gr.add_vertex();
		if(mode == 1) gr.set_vertex_weight(i + 1, r.max < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.max);
		if(mode == 2) gr.set_vertex_weight(i + 1, r.ave < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;
		vi.regional = regional[i];
		vi.type = pexons[i].type;
		// if (length < reliability_threshold) vi.reliable = false;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = bb.rpos;
	vin.rpos = bb.rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		ei.strand = b.strand;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
	}

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if(r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].max;
			if(mode == 2 && i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].max;
			if(mode == 2 && i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}
	}

	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if(x.rpos != y.lpos) continue;

		assert(x.rpos == y.lpos);
		
		int xd = gr.out_degree(i + 1);
		int yd = gr.in_degree(i + 2);
		double wt = min_guaranteed_edge_weight;
		if(mode == 1) wt = (xd < yd) ? x.max: y.max;
		if(mode == 2) wt = (xd < yd) ? x.ave: y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		double w = (wt < min_guaranteed_edge_weight) ? min_guaranteed_edge_weight : wt;
		gr.set_edge_weight(p, w);
		edge_info ei;
		ei.weight = w;
		gr.set_edge_info(p, ei);
	}

	gr.strand = bb.strand;
	gr.chrm = bb.chrm;
	return 0;
}

int bundle::revise_splice_graph()
{
	bool b = false;
	while(true)
	{
		b = tackle_false_boundaries();
		if(b == true) continue;

		b = remove_false_boundaries();
		if(b == true) continue;

		b = remove_inner_boundaries();
		if(b == true) continue;

		b = remove_small_exons();
		if(b == true) continue;

		b = remove_intron_contamination();
		if(b == true) continue;

		b = remove_small_junctions();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		b = extend_start_boundaries();
		if(b == true) continue;

		b = extend_end_boundaries();
		if(b == true) continue;

		b = extend_boundaries();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		b = keep_surviving_edges();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		break;
	}

	//find_contamination_chain();
	refine_splice_graph();

	return 0;
}

int bundle::refine_splice_graph()
{
	while(true)
	{
		bool b = false;
		for(int i = 1; i < gr.num_vertices() - 1; i++)
		{
			if(gr.degree(i) == 0) continue;
			if(gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
			gr.clear_vertex(i);
			b = true;
		}
		if(b == false) break;
	}
	return 0;
}

int bundle::refine_modified_splice_graph()
{
	while(true)
	{
		bool b = false;
		for(int i = 1; i < new_gr.num_vertices() - 1; i++)
		{
			if(new_gr.degree(i) == 0) continue;
			if(new_gr.in_degree(i) >= 1 && new_gr.out_degree(i) >= 1) continue;
			new_gr.clear_vertex(i);
			b = true;
		}
		if(b == false) break;
	}
	return 0;
}
bool bundle::extend_start_boundaries()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		PEB p = gr.edge(0, i);
		if(p.second == true) continue;

		double wv = gr.get_vertex_weight(i);
		double we = 0;
		PEEI pei = gr.in_edges(i);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			we += gr.get_edge_weight(*it);
		}

		if(wv < we || wv < 10 * we * we + 10) continue;

		edge_descriptor ee = gr.add_edge(0, i);
		gr.set_edge_weight(ee, wv - we);
		gr.set_edge_info(ee, edge_info());

		vertex_info vi = gr.get_vertex_info(i);
		if(verbose >= 2) printf("extend start boundary: vertex = %d, wv = %.2lf, we = %.2lf, pos = %d\n", i, wv, we, vi.lpos);

		flag = true;
	}
	return flag;
}

bool bundle::extend_end_boundaries()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		PEB p = gr.edge(i, gr.num_vertices() - 1);
		if(p.second == true) continue;

		double wv = gr.get_vertex_weight(i);
		double we = 0;
		PEEI pei = gr.out_edges(i);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			we += gr.get_edge_weight(*it);
		}

		if(wv < we || wv < 10 * we * we + 10) continue;

		edge_descriptor ee = gr.add_edge(i, gr.num_vertices() - 1);
		gr.set_edge_weight(ee, wv - we);
		gr.set_edge_info(ee, edge_info());

		vertex_info vi = gr.get_vertex_info(i);
		if(verbose >= 2) printf("extend end boundary: vertex = %d, wv = %.2lf, we = %.2lf, pos = %d\n", i, wv, we, vi.rpos);

		flag = true;
	}
	return flag;
}

bool bundle::extend_boundaries()
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p = gr.get_vertex_info(t).lpos - gr.get_vertex_info(s).rpos;
		double we = gr.get_edge_weight(e);
		double ws = gr.get_vertex_weight(s);
		double wt = gr.get_vertex_weight(t);

		if(p <= 0) continue;
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		bool b = false;
		if(gr.out_degree(s) == 1 && ws >= 10.0 * we * we + 10.0) b = true;
		if(gr.in_degree(t) == 1 && wt >= 10.0 * we * we + 10.0) b = true;

		if(b == false) continue;

		if(gr.out_degree(s) == 1)
		{
			edge_descriptor ee = gr.add_edge(s, gr.num_vertices() - 1);
			gr.set_edge_weight(ee, ws);
			gr.set_edge_info(ee, edge_info());
		}
		if(gr.in_degree(t) == 1)
		{
			edge_descriptor ee = gr.add_edge(0, t);
			gr.set_edge_weight(ee, wt);
			gr.set_edge_info(ee, edge_info());
		}
		gr.remove_edge(e);

		if(verbose >= 2) printf("extend boundary: s = %d, t = %d, ws = %.1lf, we = %.1lf\n", s, t, ws, we);
		return true;
	}

	return false;
}

VE bundle::compute_maximal_edges()
{
	typedef pair<double, edge_descriptor> PDE;
	vector<PDE> ve;

	undirected_graph ug;
	edge_iterator it1, it2;
	PEEI pei;
	for(int i = 0; i < gr.num_vertices(); i++) ug.add_vertex();
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		ug.add_edge(s, t);
		ve.push_back(PDE(w, e));
	}

	vector<int> vv = ug.assign_connected_components();

	sort(ve.begin(), ve.end());

	for(int i = 1; i < ve.size(); i++) assert(ve[i - 1].first <= ve[i].first);

	VE x;
	set<int> sc;
	for(int i = ve.size() - 1; i >= 0; i--)
	{
		edge_descriptor e = ve[i].second;
		double w = gr.get_edge_weight(e);
		if(w < 1.5) break;
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices()) continue;
		int c1 = vv[s];
		int c2 = vv[t];
		assert(c1 == c2);
		if(sc.find(c1) != sc.end()) continue;
		x.push_back(e);
		sc.insert(c1);
	}
	return x;
}

bool bundle::keep_surviving_edges()
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(w < min_surviving_edge_weight) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	VE me = compute_maximal_edges();
	for(int i = 0; i < me.size(); i++)
	{
		edge_descriptor ee = me[i];
		se.insert(ee);
		sv1.insert(ee->target());
		sv2.insert(ee->source());
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		if(verbose >= 2) printf("remove edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}

bool bundle::remove_small_exons()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;

		if(p2 - p1 >= min_exon_length) continue;
		if(gr.degree(i) <= 0) continue;

		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			//if(gr.out_degree(s) <= 1) b = false;
			if(s != 0 && gr.get_vertex_info(s).rpos == p1) b = false;
			if(b == false) break;
		}
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			//if(gr.in_degree(t) <= 1) b = false;
			if(t != gr.num_vertices() - 1 && gr.get_vertex_info(t).lpos == p2) b = false;
			if(b == false) break;
		}

		if(b == false) continue;

		// only consider boundary small exons
		if(gr.edge(0, i).second == false && gr.edge(i, gr.num_vertices() - 1).second == false) continue;

		//gr.clear_vertex(i);
		if(verbose >= 2) printf("remove small exon: length = %d, pos = %d-%d\n", p2 - p1, p1, p2);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
}

bool bundle::remove_small_junctions()
{
	SE se;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;
		double wi = gr.get_vertex_weight(i);

		// compute max in-adjacent edge
		double ws = 0;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_vertex_weight(s);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos != p1) continue;
			if(w < ws) continue;
			ws = w;
		}

		// remove small in-junction
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_edge_weight(e);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == p1) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

		// compute max out-adjacent edge
		double wt = 0;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			double w = gr.get_vertex_weight(t);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos != p2) continue;
			if(w < wt) continue;
			wt = w;
		}

		// remove small in-junction
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			int t = e->target();
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos == p2) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

	}

	if(se.size() <= 0) return false;

	for(SE::iterator it = se.begin(); it != se.end(); it++)
	{
		edge_descriptor e = (*it);
		vertex_info v1 = gr.get_vertex_info(e->source());
		vertex_info v2 = gr.get_vertex_info(e->target());
		if(verbose >= 2) printf("remove small junction: length = %d, pos = %d-%d\n", v2.lpos - v1.rpos, v2.lpos, v1.rpos);
		gr.remove_edge(e);
	}

	return true;
}

bool bundle::remove_inner_boundaries()
{
	bool flag = false;
	int n = gr.num_vertices() - 1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;

		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		PEEI pei = gr.in_edges(i);
		edge_iterator it1 = pei.first, it2 = pei.second;
		edge_descriptor e1 = (*it1);

		pei = gr.out_edges(i);
		it1 = pei.first;
		it2 = pei.second;
		edge_descriptor e2 = (*it1);

		int s = e1->source();
		int t = e2->target();

		if(s != 0 && t != n) continue;
		if(s != 0 && gr.out_degree(s) == 1) continue;
		if(t != n && gr.in_degree(t) == 1) continue;

		if(vi.stddev >= 0.01) continue;

		if(verbose >= 2) printf("remove inner boundary: vertex = %d, weight = %.2lf, length = %d, pos = %d-%d\n",
				i, gr.get_vertex_weight(i), vi.length, vi.lpos, vi.rpos);

		//gr.clear_vertex(i);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);
		flag = true;
	}
	return flag;
}

bool bundle::remove_intron_contamination()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;

		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		PEEI pei = gr.in_edges(i);
		it1 = pei.first;
		edge_descriptor e1 = (*it1);
		pei = gr.out_edges(i);
		it1 = pei.first;
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		double wv = gr.get_vertex_weight(i);

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		if(gr.get_vertex_info(s).rpos != vi.lpos) continue;
		if(gr.get_vertex_info(t).lpos != vi.rpos) continue;

		PEB p = gr.edge(s, t);
		if(p.second == false) continue;

		edge_descriptor ee = p.first;
		double we = gr.get_edge_weight(ee);

		if(wv > we) continue;
		if(wv > max_intron_contamination_coverage) continue;

		if(verbose >= 2) printf("clear intron contamination %d, weight = %.2lf, length = %d, edge weight = %.2lf\n", i, wv, vi.length, we);

		//gr.clear_vertex(i);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
}

// add by Mingfu -- to use paired-end reads to remove false boundaries
bool bundle::remove_false_boundaries()
{
	// map<int, int> fb1;		// end
	// map<int, int> fb2;		// start
	// for(int i = 0; i < br.fragments.size(); i++)
	// {
	// 	fragment &fr = br.fragments[i];
	// 	if(fr.paths.size() == 1 && fr.paths[0].type == 1) continue;
	// 	//if(fr.h1->bridged == true || fr.h2->bridged == true) continue;

	// 	// only use uniquely aligned reads
	// 	//if(fr.h1->nh >= 2 || fr.h2->nh >= 2) continue;
	// 	if(br.breads.find(fr.h1->qname) != br.breads.end()) continue;

	// 	// calculate actual length
	// 	vector<int> v = align_fragment(fr);
		
	// 	if(v.size() <= 1) continue;

	// 	int32_t tlen = 0;
	// 	int32_t offset1 = (fr.lpos - pexons[v.front()].lpos);
	// 	int32_t offset2 = (pexons[v.back()].rpos - fr.rpos);
	// 	for(int i = 0; i < v.size(); i++)
	// 	{
	// 		int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
	// 		tlen += l;
	// 	}
	// 	tlen -= offset1;
	// 	tlen -= offset2;

	// 	int u1 = gr.locate_vertex(fr.h1->rpos - 1);
	// 	int u2 = gr.locate_vertex(fr.h2->pos);

	// 	if(u1 < 0 || u2 < 0) continue;
	// 	if(u1 >= u2) continue;

	// 	vertex_info v1 = gr.get_vertex_info(u1);
	// 	vertex_info v2 = gr.get_vertex_info(u2);

	// 	int types = 0;
	// 	int32_t lengths = 0;
	// 	for(int k = 0; k < fr.paths.size(); k++) types += fr.paths[k].type;
	// 	for(int k = 0; k < fr.paths.size(); k++) lengths += fr.paths[k].length;

	// 	bool use = true;
	// 	if(fr.paths.size() == 1 && types == 2 && tlen > 10000) use = false;
	// 	//if(fr.paths.size() == 1 && types == 2 && lengths <= 1.5 * insertsize_high) use = false;
	// 	//if(fr.paths.size() == 1 && types == 2 && tlen <= 1.5 * insertsize_high) use = false;
	// 	//if(fr.paths.size() == 1 && types == 2 && lengths <= 2 * tlen) use = false;

	// 	if(verbose >= 2) printf("%s: u1 = %d, %d-%d, u2 = %d, %d-%d, h1.rpos = %d, h2.lpos = %d, #bridging = %lu, types = %d, lengths = %d, tlen = %d, use = %c\n", 
	// 			fr.h1->qname.c_str(), u1, v1.lpos, v1.rpos, u2, v2.lpos, v2.rpos, fr.h1->rpos, fr.h2->pos, fr.paths.size(), types, lengths, tlen, use ? 'T' : 'F');

	// 	if(use == false) continue;

	// 	//if(gr.get_vertex_info(u1).rpos == fr.h1->rpos)
	// 	{
	// 		if(fb1.find(u1) != fb1.end()) fb1[u1]++;
	// 		else fb1.insert(make_pair(u1, 1));
	// 	}

	// 	//if(gr.get_vertex_info(u2).lpos == fr.h2->pos)
	// 	{
	// 		if(fb2.find(u2) != fb2.end()) fb2[u2]++;
	// 		else fb2.insert(make_pair(u2, 1));
	// 	}
	// }

	bool b = false;
	// for(auto &x : fb1)
	// {
	// 	PEB p = gr.edge(x.first, gr.num_vertices() - 1);
	// 	vertex_info vi = gr.get_vertex_info(x.first);
	// 	if(vi.type == EMPTY_VERTEX) continue;
	// 	if(p.second == false) continue;
	// 	double w = gr.get_vertex_weight(x.first);
	// 	double z = log(1 + w) / log(1 + x.second);
	// 	double s = log(1 + w) - log(1 + x.second);
	// 	if(s > 1.5) continue;
	// 	if(verbose >= 2) printf("detect false end boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.rpos, x.second, x.first, w, vi.type, z, s); 
	// 	//gr.remove_edge(p.first);
	// 	vi.type = EMPTY_VERTEX;
	// 	gr.set_vertex_info(x.first, vi);
	// 	b = true;
	// }

	// for(auto &x : fb2)
	// {
	// 	PEB p = gr.edge(0, x.first);
	// 	vertex_info vi = gr.get_vertex_info(x.first);
	// 	if(vi.type == EMPTY_VERTEX) continue;
	// 	if(p.second == false) continue;
	// 	double w = gr.get_vertex_weight(x.first);
	// 	double z = log(1 + w) / log(1 + x.second);
	// 	double s = log(1 + w) - log(1 + x.second);
	// 	if(s > 1.5) continue;
	// 	if(verbose >= 2) printf("detect false start boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.lpos, x.second, x.first, w, vi.type, z, s); 
	// 	//gr.remove_edge(p.first);
	// 	vi.type = EMPTY_VERTEX;
	// 	gr.set_vertex_info(x.first, vi);
	// 	b = true;
	// }
	return b;
}

bool bundle::tackle_false_boundaries()
{
	bool b = false;
	// vector<int> points(pexons.size(), 0);
	// for(int k = 0; k < br.fragments.size(); k++)
	// {
	// 	fragment &fr = br.fragments[k];

	// 	if(fr.paths.size() != 1) continue;
	// 	if(fr.paths[0].type != 2) continue;
	// 	if(br.breads.find(fr.h1->qname) != br.breads.end()) continue;

	// 	vector<int> v = align_fragment(fr);
	// 	if(v.size() <= 1) continue;

	// 	int32_t offset1 = (fr.lpos - pexons[v.front()].lpos);
	// 	int32_t offset2 = (pexons[v.back()].rpos - fr.rpos);

	// 	int32_t tlen = 0;
	// 	for(int i = 0; i < v.size(); i++)
	// 	{
	// 		int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
	// 		tlen += l;
	// 	}
	// 	tlen -= offset1;
	// 	tlen -= offset2;

	// 	// print
	// 	//fr.print(99);
	// 	if(verbose >= 2) printf("break fragment %s: total-length = %d, bridge-length = %d\n", fr.h1->qname.c_str(), tlen, fr.paths[0].length);
	// 	/*
	// 	for(int i = 0; i < v.size(); i++)
	// 	{
	// 		int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
	// 		if(i == 0) l -= offset1;
	// 		if(i == v.size() - 1) l -= offset2;
	// 		printf(" vertex %d: length = %d, region = %d-%d -> %d\n", v[i], l, pexons[v[i]].lpos, pexons[v[i]].rpos, pexons[v[i]].rpos - pexons[v[i]].lpos);
	// 	}
	// 	*/

	// 	if(tlen < insertsize_low / 2.0) continue;
	// 	if(tlen > insertsize_high * 2.0) continue;
	// 	if(tlen >= fr.paths[0].length) continue;

	// 	for(int i = 0; i < v.size() - 1; i++)
	// 	{
	// 		partial_exon &px = pexons[v[i + 0]];
	// 		partial_exon &py = pexons[v[i + 1]];
	// 		if(px.rtype == END_BOUNDARY) 
	// 		{
	// 			if(verbose >= 2) printf("break ending vertex %d, pos = %d\n", v[i], px.rpos);
	// 			points[v[i + 0]] += 1;
	// 		}
	// 		if(py.ltype == START_BOUNDARY) 
	// 		{
	// 			if(verbose >= 2) printf("break starting vertex %d, pos = %d\n", v[i + 1], py.lpos);
	// 			points[v[i + 1]] += 1;
	// 		}
	// 	}
	// }

	// for(int k = 0; k < points.size(); k++)
	// {
	// 	if(points[k] <= 0) continue;
	// 	vertex_info vi = gr.get_vertex_info(k + 1);
	// 	if(vi.type == EMPTY_VERTEX) continue;
	// 	PEB p = gr.edge(k + 1, gr.num_vertices() - 1);
	// 	if(p.second == false) continue;
	// 	double w = gr.get_vertex_weight(k + 1);
	// 	double z = log(1 + w) / log(1 + points[k]);
	// 	double s = log(1 + w) - log(1 + points[k]);
	// 	if(verbose >= 2) printf("tackle false end boundary %d with %d reads, vertex = %d, w = %.2lf, z = %.2lf, s = %.2lf\n", pexons[k].rpos, points[k], k + 1, w, z, s);
	// 	if(s > 1.5) continue;
	// 	vi.type = EMPTY_VERTEX;
	// 	gr.set_vertex_info(k + 1, vi);
	// 	b = true;
	// }

	// for(int k = 0; k < points.size(); k++)
	// {
	// 	if(points[k] <= 0) continue;
	// 	vertex_info vi = gr.get_vertex_info(k + 1);
	// 	if(vi.type == EMPTY_VERTEX) continue;
	// 	PEB p = gr.edge(0, k + 1);
	// 	if(p.second == false) continue;
	// 	double w = gr.get_vertex_weight(k + 1);
	// 	double z = log(1 + w) / log(1 + points[k]);
	// 	double s = log(1 + w) - log(1 + points[k]);
	// 	if(verbose >= 2) printf("tackle false start boundary %d with %d reads, vertex = %d, w = %.2lf, z = %.2lf, s = %.2lf\n", pexons[k].lpos, points[k], k + 1, w, z, s);
	// 	if(s > 1.5) continue;
	// 	vi.type = EMPTY_VERTEX;
	// 	gr.set_vertex_info(k + 1, vi);
	// 	b = true;
	// }

	return b;
}

int bundle::find_contamination_chain()
{
	int min_vertices = 5;
	double max_coverage = 4.0;
	int32_t max_distance = 2000;

	vector<int> chain;
	vector<string> types;
	for(int i = 1; i < pexons.size() - 1; i++)
	{
		string type = "";
		partial_exon &pe = pexons[i];
		if(pe.max > max_coverage) continue;

		if(pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) type = "island";
		if(pe.ltype == START_BOUNDARY && pe.rtype == RIGHT_SPLICE) type = "start";
		if(pe.ltype == LEFT_SPLICE && pe.rtype == RIGHT_SPLICE && gr.edge(i - 1, i + 1).second == true) type = "intron";
		if(pe.ltype == LEFT_SPLICE && pe.rtype == END_BOUNDARY) type = "end";

		if(type == "") continue;

		chain.push_back(i);
		types.push_back(type);
	}

	if(chain.size() == 0) return 0;

	int32_t pre = 0;
	for(int k = 0; k < chain.size(); k++)
	{
		partial_exon &pe = pexons[chain[k]];
		printf("chain %d, pexon = %d, type = %s, pos = %d-%d, len = %d, cov = %.2lf, dist = %d\n", k, chain[k], types[k].c_str(), pe.lpos, pe.rpos, pe.rpos - pe.lpos, pe.max, pe.lpos - pre);
		pre = pe.rpos;
	}

	int k1 = -1;
	pre = 0 - max_distance - 1;
	for(int k = 0; k < chain.size(); k++)
	{
		partial_exon &pe = pexons[chain[k]];
		int32_t dist = pe.lpos - pre;
		if(dist > max_distance)
		{
			if(k - k1 > min_vertices)
			{
				printf("delete chain: %d-%d\n", k1 + 1, k - 1);
				for(int i = k1 + 1; i < k; i++)
				{
					gr.clear_vertex(chain[k] + 1);
				}
			}
			k1 = k - 1;
		}
		pre = pe.rpos;
	}

	if((int)(chain.size()) > min_vertices + k1)
	{
		printf("delete chain: %d-%d\n", k1 + 1, (int)(chain.size()) - 1);
		for(int i = k1 + 1; i < chain.size(); i++)
		{
			gr.clear_vertex(chain[i] + 1);
		}
	}
	return 0;
}

int bundle::count_junctions() const
{
	int x = 0;
	for(int i = 0; i < junctions.size(); i++)
	{
		x += junctions[i].count;
	}
	return x;
}

int bundle::print(int index)
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

	printf("tid = %d, #bb.hits = %lu, #partial-exons = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			bb.tid, bb.hits.size(), pexons.size(), bb.chrm.c_str(), bb.lpos, bb.rpos, bb.strand, n0, np, nq);

	if(verbose <= 1) return 0;

	// print bb.hits
	//for(int i = 0; i < bb.hits.size(); i++) bb.hits[i].print();

	// print regions
	printf("Printing regions:\n");
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	printf("Printing junctions:\n");
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	// print partial exons
	printf("Printing pexons:\n");
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}

	printf("\n");

	return 0;
}

int bundle::output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const
{
	for(int i = 0; i < p.size(); i++)
	{
		string tid = gid + "." + tostring(i);
		output_transcript(fout, p[i], gid, tid);
	}
	return 0;
}

int bundle::output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const
{
	fout.precision(2);
	fout<<fixed;

	const vector<int> &v = p.v;
	double coverage = p.abd;		// number of molecular

	assert(v[0] == 0);
	assert(v[v.size() - 1] == pexons.size() + 1);
	if(v.size() < 2) return 0;

	int ss = v[1];
	int tt = v[v.size() - 2];
	int32_t ll = pexons[ss - 1].lpos;
	int32_t rr = pexons[tt - 1].rpos;

	fout<<bb.chrm.c_str()<<"\t";		// chromosome name
	fout<<algo.c_str()<<"\t";		// source
	fout<<"transcript\t";			// feature
	fout<<ll + 1<<"\t";				// left position
	fout<<rr<<"\t";					// right position
	fout<<1000<<"\t";				// score, now as abundance
	fout<<bb.strand<<"\t";				// bb.strand
	fout<<".\t";					// frame
	fout<<"gene_id \""<<gid.c_str()<<"\"; ";
	fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
	fout<<"coverage \""<<coverage<<"\";"<<endl;

	join_interval_map jmap;
	for(int k = 1; k < v.size() - 1; k++)
	{
		const partial_exon &r = pexons[v[k] - 1];
		jmap += make_pair(ROI(r.lpos, r.rpos), 1);
	}

	int cnt = 0;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		fout<<bb.chrm.c_str()<<"\t";			// chromosome name
		fout<<algo.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<lower(it->first) + 1<<"\t";	// left position
		fout<<upper(it->first)<<"\t";		// right position
		fout<<1000<<"\t";					// score
		fout<<bb.strand<<"\t";					// bb.strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gid.c_str()<<"\"; ";
		fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
		fout<<"exon_number \""<<++cnt<<"\"; ";
		fout<<"coverage \""<<coverage<<"\";"<<endl;
	}
	return 0;
}

int bundle::output_transcripts(vector<transcript> &trsts, const vector<path> &p, const string &gid) const
{
	trsts.clear();
	for(int i = 0; i < p.size(); i++)
	{
		string tid = gid + "." + tostring(i);
		transcript trst;
		output_transcript(trst, p[i], gid, tid);
		trsts.push_back(trst);
	}
	return 0;
}

int bundle::output_transcripts(gene &gn, const vector<path> &p, const string &gid) const
{
	for(int i = 0; i < p.size(); i++)
	{
		string tid = gid + "." + tostring(i);
		transcript trst;
		output_transcript(trst, p[i], gid, tid);
		gn.add_transcript(trst);
	}
	return 0;
}

int bundle::output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const
{
	trst.seqname = bb.chrm;
	trst.source = algo;
	trst.gene_id = gid;
	trst.transcript_id = tid;
	trst.coverage = p.abd;
	trst.strand = bb.strand;

	const vector<int> &v = p.v;
	join_interval_map jmap;
	for(int k = 1; k < v.size() - 1; k++)
	{
		const partial_exon &r = pexons[v[k] - 1];
		jmap += make_pair(ROI(r.lpos, r.rpos), 1);
	}
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		trst.add_exon(lower(it->first), upper(it->first));
	}
	return 0;
}

int bundle::build_hyper_set()
{
	map<vector<int>, int> m;
	map<vector<int>, map<vector<int>, int>> mm;

	// TODO: skip the followoing two paragraphs
	// for(int k = 0; k < br.fragments.size(); k++)
	// {
	// 	fragment &fr = br.fragments[k];

	// 	if(fr.type != 0) continue;	// note by Qimin, skip if not paired-end fragments

	// 	if(fr.h1->paired != true) printf("error type: %d\n", fr.type);
	// 	assert(fr.h1->paired == true);
	// 	assert(fr.h2->paired == true);

	// 	if(fr.paths.size() != 1) continue;
	// 	if(fr.paths[0].type != 1) continue;

	// 	//if(fr.h1->bridged == false) continue;
	// 	//if(fr.h2->bridged == false) continue;

	// 	vector<int> v = align_fragment(fr);
		
	// 	if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, fr.cnt));
	// 	else m[v] += fr.cnt;
	// }

	
	// note by Qimin, bridge umi-linked fragments into one single long path
	// for(int k = 0; k < br.umiLink.size(); k++)
	// {
	// 	vector<int> v;
	// 	v.clear();

	// 	int cnt = 0;

	// 	// if only one fr, no need to bridge into longer one
	// 	if(br.umiLink[k].size() == 1)
	// 	{
	// 		fragment &fr = br.fragments[(br.umiLink[k][0])];

	// 		if(fr.paths.size() != 1) continue;

	// 		// TODO: "bridged" may not be correct
	// 		// if(fr.h1->bridged == false) continue;
	// 		// if(fr.h2->bridged == false) continue;

	// 		v = align_fragment(fr);
	// 		if(fr.paths.size() != 1 || fr.paths[0].type != 1) v.clear();

	// 		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, fr.cnt));
	// 		else m[v] += fr.cnt;

	// 		continue;
	// 	}

	// 	// if multiple fr in umi-link, bridge into one single long path
	// 	for(int kk = 0; kk < br.umiLink[k].size(); kk++)
	// 	{
	// 		fragment &fr = br.fragments[(br.umiLink[k][kk])];

	// 		// if unbridge, then trucate and add to m
	// 		if(fr.paths.size() != 1 || fr.h1->bridged == false || fr.h2->bridged == false)
	// 		{
	// 			if(v.size() > 0)
	// 			{
	// 				if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
	// 				else m[v] += cnt;
	// 			}

	// 			v.clear();
	// 			cnt = 0;

	// 			continue;
	// 		}

	// 		// otherwise, add and merge cur_v to v
	// 		vector<int> cur_v = align_fragment(fr);
	// 		if(fr.paths.size() != 1 || fr.paths[0].type != 1) cur_v.clear();

	// 		if(cur_v.size()==0)
	// 		{
	// 			if(v.size() > 0)
	// 			{
	// 				if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
	// 				else m[v] += cnt;
	// 			}

	// 			v.clear();
	// 			cnt = 0;

	// 			continue;

	// 		}
	// 		cnt += fr.cnt;

	// 		v.insert(v.end(), cur_v.begin(), cur_v.end());
	// 		sort(v.begin(), v.end());
	// 		vector<int>::iterator iter = unique(v.begin(),v.end());
	// 		v.erase(iter,v.end());
	// 	}

	// 	if(v.size() > 0)
	// 	{
	// 		/*
	// 		printf("v = ");
	// 		for(int ii = 0; ii < v.size(); ii++)
	// 		{
	// 			printf("%d ", v[ii]);
	// 		}
	// 		printf("\n");
	// 		*/

	// 		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
	// 		else m[v] += cnt;
	// 	}
	// }
	
	// TODO: use starring from here
	for(int k = 0; k < bb.hits.size(); k++)
	{
		hit &h = bb.hits[k];

		// bridged used here, but maybe okay
		if(h.bridged == true) continue;

		vector<int> v1 = align_hit(h);

		// v is the list of nodes without unrealizable ones
		vector<int> v;
		for(int j=0; j < v1.size(); j++)
		{
			partial_exon &p = pexons[v1[j]];
			// p.rel = true;
			// if p is unreliable, skip it
			// otherwise add it to newv
			if(p.rel == true) v.push_back(v1[j]);
		}

		if(mm.find(v) == mm.end())
		{
			map<vector<int>, int> tm;
			tm.insert(make_pair(v1, 1));
			mm.insert(tm);
		}
		else
		{
			if(mm[v].find(v1) == mm[v].end())
			{
				mm[v].insert(make_pair(v1, 1));
			}
			else
			{
				mm[v][v1] += 1;
			}
		}
		
		//if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, 1));
		//else m[v] += 1;
		
	}

	hs.clear();
	hs2.clear();
	for(auto & it = mm.begin(); it != mm.end(); it++)
	{
		const vector<int> &v = it->first;
		int cc = 0;
		map<vector<int>, int> & tm = it->second;

		int start = hs2.nodes.size();
		for(auto &ix = tm.begin(); ix != tm.end(); ix++)
		{
			const vector<int> &v1 = ix->first;
			int c = ix->second;
			hs2.add_node_list(v1, c);

			cc += c;
		}
		int end = hs2.nodes.size();

		hs.add_node_list(v, cc);
		plink.push_back(make_pair(start, end));
	}
	printf("------------------------------------\n");

	//printf("Printing the modified hyperset for Bundle %d:\n", index);
	hs.print();
	return 0;
}

int bundle::refine_hyper_set()
{
	//MVII new_nodes;
	//for(auto &x: hs.nodes)
	//{
	//	vector<int> &v = x.first;
	//	vector<int> newv;
	//	for(int k = 0; k < v.size(); k++)
	//	{
	//		partial_exon &p = pexons[v[k] - 1];
	//		// if p is unreliable, skip it
	//		// otherwise add it to newv
	//		if(p.rel) newv.push_back(p);
	//	}
	//	new_nodes.insert(make_pair(newv, x.second));
	//}
	//hs.nodes = new_nodes;
	return 0;
}

int bundle::rebuild_splice_graph_using_refined_hyper_set(int mode)
{
	// TODO: we fill up new_gr in this function
	// splice_graph new_gr;

	new_gr.clear();

	// vertices: start, each region, end
	new_gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = bb.lpos;
	vi0.rpos = bb.lpos;
	new_gr.set_vertex_weight(0, 0);
	new_gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		new_gr.add_vertex();
		if(mode == 1) new_gr.set_vertex_weight(i + 1, r.max < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.max);
		if(mode == 2) new_gr.set_vertex_weight(i + 1, r.ave < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;
		vi.regional = regional[i];
		vi.type = pexons[i].type;
		
		new_gr.set_vertex_info(i + 1, vi);
	}

	new_gr.add_vertex();
	vertex_info vin;
	vin.lpos = bb.rpos;
	vin.rpos = bb.rpos;
	new_gr.set_vertex_weight(pexons.size() + 1, 0);
	new_gr.set_vertex_info(pexons.size() + 1, vin);

	// we use reads to add edges
	// TODO: use starring from here
	for(int k = 0; k < bb.hits.size(); k++)
	{
		hit &h = bb.hits[k];

		// skip this for now
		/*
		// bridged used here, but maybe okay
		if(h.bridged == true) continue;
		*/

		// align h to old graph gr
		vector<int> v = align_hit(h);
		
		// printf("v size = %d\n", v.size());
		// small reliable pexons are merged
		// for(int j = 1; j < (int)(v.size()-1); j++){
		// 	partial_exon &p = pexons[v[j]];
		// 	partial_exon &p_prev = pexons[v[j-1]];
		// 	partial_exon &p_next = pexons[v[j+1]];

		// 	// printf("Working on [%d, %d]\n", p.lpos, p.rpos);

		// 	// if there is an unreliable pexon in the middle of two adjacent reliable pexons
		// 	if(p.rel == false && (p_prev.rel == true && p_next.rel == true) && (p.rpos == p_next.lpos && p.lpos == p_prev.rpos)) 
		// 	{
		// 		p.rel = true;
		// 	}
		// 	// if there are multiple adjacent unreliable exons which can be meged to form a single reliable exon
		// 	else if(p.rel == false && p_next.rel == false && (p.rpos == p_next.lpos))
		// 	{
		// 		int total_len = p.rpos - p.lpos;
		// 		int nxt_idx = j+1;
		// 		partial_exon &nxt = p_next;
		// 		partial_exon &cur = p;
				
		// 		while(nxt.rel == false && (nxt.lpos == cur.rpos ))
		// 		{
		// 			total_len += (nxt.rpos - nxt.lpos);
		// 			nxt_idx++;
		// 			if(nxt_idx == v.size()) break;
		// 			cur = nxt;
		// 			nxt = pexons[v[nxt_idx]];
		// 		}

		// 		if(total_len >= reliability_threshold){
		// 			for(int jj=j; jj < nxt_idx; jj++) pexons[v[jj]].rel = true;
		// 			j = nxt_idx-1;
					
		// 		}


		// 	}

		// }



		// remove unreliable vertices
		vector<int> newv;
		for(int j = 0; j < v.size(); j++)
		{
			partial_exon &p = pexons[v[j]];
			// if p is unreliable, skip it
			// otherwise add it to newv
			if(p.rel == true) newv.push_back(v[j]);
		}

		
		// example: if newv = (0, 2, 4, 5)
		// add edges or increase weight
		//printf("newv size: %lu , old_v size = %lu\n", newv.size(), v.size() );	
		//if(newv.size() <= 0) continue;

		for(int j = 0; j < (int)(newv.size() - 1) ; j++)
		{
			// either create new edge or increase its weight
			
			//printf("j = %d,", j);
			int v1 = newv[j + 0] + 1;
			int v2 = newv[j + 1] + 1;
			PEB peb = new_gr.edge(v1, v2);
			if(peb.second == false)
			{
				edge_descriptor p = new_gr.add_edge(v1, v2);
				edge_info ei;
				ei.weight = 1;
				new_gr.set_edge_info(p, ei);
				new_gr.set_edge_weight(p, 1);
			}
			else
			{
				double w = new_gr.get_edge_weight(peb.first) + 1;
				new_gr.set_edge_weight(peb.first, w);
				edge_info ei = new_gr.get_edge_info(peb.first);
				ei.weight = w;
				new_gr.set_edge_info(peb.first, ei);
			}
		}
		//printf("\nDone making newv\n");
	}

	// do not use junctions to add edges
	/*
	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		ei.strand = b.strand;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
	}
	*/

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if(r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = new_gr.add_edge(ss, i + 1);
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].max;
			if(mode == 2 && i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;
			new_gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			new_gr.set_edge_info(p, ei);
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = new_gr.add_edge(i + 1, tt);
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].max;
			if(mode == 2 && i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;
			new_gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			new_gr.set_edge_info(p, ei);
		}
	}

	// be careful here v1 = [100, 300], v2 = [300, 400]
	// in this case, if v1,v2 is not yet connected, then add edge for them
	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if(x.rpos != y.lpos) continue;
		
		if(x.rel == false) continue;
		if(y.rel == false) continue;

		if(new_gr.edge(i + 1, i + 2).second == true) continue;
		
		// a rare case that we should add an edge

		assert(x.rpos == y.lpos);
		
		/*
		int xd = new_gr.out_degree(i + 1);
		int yd = new_gr.in_degree(i + 2);
		double wt = min_guaranteed_edge_weight;
		if(mode == 1) wt = (xd < yd) ? x.max: y.max;
		if(mode == 2) wt = (xd < yd) ? x.ave: y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;
		*/

		//double w = (wt < min_guaranteed_edge_weight) ? min_guaranteed_edge_weight : wt;
		edge_descriptor p = new_gr.add_edge(i + 1, i + 2);
		new_gr.set_edge_weight(p, 1);
		edge_info ei;
		ei.weight = 1;
		new_gr.set_edge_info(p, ei);
	}

	new_gr.strand = bb.strand;
	new_gr.chrm = bb.chrm;
	// gr = new_gr;
	return 0;
}

void bundle::print_fmap(){
	printf("fmap:\n");
	SIMI it = fmap.begin();
	while(it != fmap.end()){
		printf("%d %d: %d\n", lower(it->first), upper(it->first), it->second);
		it++;
	}

}
