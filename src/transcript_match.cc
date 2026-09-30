/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>

#include "transcript_match.h"
#include "util.h"

transcript_match::transcript_match()
{
	assigned = false;
	shared_exonic_length = 0;
	transcript_exonic_length = 0;
	shared_junctions = 0;
	transcript_junctions = 0;
}

transcript_matcher::transcript_matcher(const splice_graph &gr)
{
	build(gr);
}

transcript_matcher::transcript_matcher(const string &c, char s, const vector<PI32> &v,
		const set<int64_t> &j)
	: chrm(c), strand(s), junctions(j)
{
	build_exons(v);
}

int transcript_matcher::build(const splice_graph &gr)
{
	chrm = gr.chrm;
	strand = gr.strand;
	exons.clear();
	junctions.clear();

	int n = gr.num_vertices();
	for(int i = 1; i < n - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.lpos >= vi.rpos) continue;
		exons += make_pair(ROI(vi.lpos, vi.rpos), 1);
	}

	PEEI p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		edge_descriptor e = *it;
		int s = e->source();
		int t = e->target();
		if(s <= 0 || t >= n - 1) continue;

		vertex_info vs = gr.get_vertex_info(s);
		vertex_info vt = gr.get_vertex_info(t);
		if(vs.rpos >= vt.lpos) continue;
		junctions.insert(pack(vs.rpos, vt.lpos));
	}
	return 0;
}

int transcript_matcher::build_exons(const vector<PI32> &v)
{
	exons.clear();
	for(int i = 0; i < v.size(); i++)
	{
		if(v[i].first >= v[i].second) continue;
		exons += make_pair(ROI(v[i].first, v[i].second), 1);
	}
	return 0;
}

int transcript_matcher::compute_shared_exonic_length(const transcript &t) const
{
	int z = 0;
	for(int i = 0; i < t.exons.size(); i++)
	{
		const PI32 &p = t.exons[i];
		JIMI it = exons.find(p.first);
		if(it == exons.end()) it = exons.lower_bound(ROI(p.first, p.second));

		for(; it != exons.end(); it++)
		{
			int32_t l = lower(it->first);
			int32_t r = upper(it->first);
			if(l >= p.second) break;
			if(r <= p.first) continue;
			int32_t ll = l > p.first ? l : p.first;
			int32_t rr = r < p.second ? r : p.second;
			if(ll < rr) z += rr - ll;
		}
	}
	return z;
}

transcript_match transcript_matcher::match(const indexed_transcript &t,
		double min_exon_overlap, double min_junction_overlap) const
{
	transcript_match m;
	m.transcript_exonic_length = t.exonic_length;
	m.transcript_junctions = t.junctions.size();

	if(t.trst.seqname != chrm) return m;
	if(t.trst.strand != strand) return m;
	if(m.transcript_exonic_length <= 0) return m;

	m.shared_exonic_length = compute_shared_exonic_length(t.trst);
	for(int i = 0; i < t.junctions.size(); i++)
	{
		int64_t p = pack(t.junctions[i].first, t.junctions[i].second);
		if(junctions.find(p) != junctions.end()) m.shared_junctions++;
	}

	if(m.shared_exonic_length < min_exon_overlap * m.transcript_exonic_length) return m;
	if(m.transcript_junctions >= 1 &&
		m.shared_junctions < min_junction_overlap * m.transcript_junctions) return m;
	m.assigned = true;
	return m;
}
