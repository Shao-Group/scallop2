/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "region.h"
#include "config.h"
#include "util.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype)
	:lpos(_lpos), rpos(_rpos), ltype(_ltype), rtype(_rtype)
{
	ave = 0;
	max = 0;
	dev = 1;
}

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), mmap(_mmap), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	ave = 0;
	max = 0;
	dev = 1;
	build_join_interval_map();
	smooth_join_interval_map();
	//split_join_interval_map();
	build_partial_exons();
	build_indel_coverage();
}

int region::build_indel_coverage()
{
	if(imap == NULL) return 0;
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		PSIMI x = locate_boundary_iterators(*imap, p.lpos, p.rpos);
		SIMI first = x.first, last = x.second;
		if(first == imap->end() || last == imap->end()) continue;

		p.indel_sum_cov = compute_sum_overlap(*imap, first, last);
		double denominator = p.ave * (p.rpos - p.lpos);
		p.indel_ratio = denominator > 0 ? p.indel_sum_cov / denominator : 0;
		int32_t left = std::max(p.lpos, lower(first->first));
		int32_t right = std::min(p.rpos, upper(last->first));
		p.left_indel = left - p.lpos;
		p.right_indel = p.rpos - right;
	}
	return 0;
}

region::~region()
{}

int region::build_join_interval_map()
{
	jmap.clear();

	PSIMI pei = locate_boundary_iterators(*mmap, lpos, rpos);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == mmap->end() || rit == mmap->end()) return 0;

	SIMI it = lit;
	while(true)
	{
		//if(it->second >= 2) 
		jmap += make_pair(it->first, 1);
		if(it == rit) break;
		it++;
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}

	return 0;
}

int region::smooth_join_interval_map()
{
	int32_t gap = min_subregion_gap;
	vector<PI32> v;
	int32_t p = lpos;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 >= p);
		assert(p2 > p1);
		if(p1 - p <= gap) v.push_back(PI32(p, p1));
		p = p2;
	}

	if(p < rpos && rpos - p <= gap) v.push_back(PI32(p, rpos));

	for(int i = 0; i < v.size(); i++)
	{
		jmap += make_pair(ROI(v[i].first, v[i].second), 1);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}
	return 0;
}

bool region::empty_subregion(int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	assert(p1 >= lpos && p2 <= rpos);

	//printf(" region = [%d, %d), subregion [%d, %d), length = %d\n", lpos, rpos, p1, p2, p2 - p1);
	if(p2 - p1 < min_subregion_len) return true;

	PSIMI pei = locate_boundary_iterators(*mmap, p1, p2);
	SIMI it1 = pei.first, it2 = pei.second;
	if(it1 == mmap->end() || it2 == mmap->end()) return true;

	int32_t sum = compute_sum_overlap(*mmap, it1, it2);
	int32_t max = compute_sum_overlap(*mmap, it1, it2);
	double ratio = sum * 1.0 / (p2 - p1);
	//printf(" region = [%d, %d), subregion [%d, %d), overlap = %.2lf\n", lpos, rpos, p1, p2, ratio);
	//if(ratio < min_subregion_overlap + max_intron_contamination_coverage) return true;
	if(ratio < min_subregion_ave) return true;
	if(max < min_subregion_max) return true;

	return false;
}

int region::build_partial_exons()
{
	pexons.clear();

	if(jmap.size() == 0) return 0;

	//printf("size = %lu, size2 = %lu, [%d, %d), [%d, %d)\n", jmap.size(), distance(jmap.begin(), jmap.end()), lower(jmap.begin()->first), upper(jmap.begin()->first), lpos, rpos);

	if(lower(jmap.begin()->first) == lpos && upper(jmap.begin()->first) == rpos)
	{
		partial_exon pe(lpos, rpos, ltype, rtype);
		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev, pe.max);
		pexons.push_back(pe);
		return 0;
	}

	if(ltype == RIGHT_SPLICE && jmap.find(ROI(lpos, lpos + 1)) == jmap.end())
	{
		partial_exon pe(lpos, lpos + 1, ltype, END_BOUNDARY);
		pe.ave = min_guaranteed_edge_weight;
		pe.max = min_guaranteed_edge_weight;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 < p2);
		
		bool b = empty_subregion(p1, p2);
		//printf(" subregion [%d, %d), empty = %c\n", p1, p2, b ? 'T' : 'F');

		if(p1 == lpos && ltype == RIGHT_SPLICE) b = false;
		if(p2 == rpos && rtype == LEFT_SPLICE) b = false;

		int lt = (p1 == lpos) ? ltype : START_BOUNDARY;
		int rt = (p2 == rpos) ? rtype : END_BOUNDARY;

		partial_exon pe(p1, p2, lt, rt);
		if(b == true) pe.type = EMPTY_VERTEX;
		else pe.type = 0;

		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev, pe.max);
		pexons.push_back(pe);
	}

	if(rtype == LEFT_SPLICE && jmap.find(ROI(rpos - 1, rpos)) == jmap.end())
	{
		partial_exon pe(rpos - 1, rpos, START_BOUNDARY, rtype);
		pe.ave = min_guaranteed_edge_weight;
		pe.max = min_guaranteed_edge_weight;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	return 0;
}

int region::print(int index) const
{
	printf("region %d: length = %d, pexons = %lu, type = (%d, %d), pos = [%d, %d), cov/dev = (%.2lf, %.2lf)\n", 
			index, rpos - lpos, pexons.size(), ltype, rtype, lpos, rpos, ave, dev);

	return 0;

	int32_t lc = compute_overlap(*mmap, lpos);
	int32_t rc = compute_overlap(*mmap, rpos - 1);
	printf("region %d: partial-exons = %lu, type = (%d, %d), pos = [%d, %d), boundary coverage = (%d, %d)\n", 
			index, pexons.size(), ltype, rtype, lpos, rpos, lc, rc);

	/*
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		printf(" [%d, %d) -> %d\n", lower(it->first), upper(it->first), it->second);
	}
	*/

	return 0;
}
