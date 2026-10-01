/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "interval_map.h"

int create_split(split_interval_map &imap, int32_t p)
{
	SIMI it = imap.find(p);
	if(it == imap.end()) return 0;
	int32_t l = lower(it->first);
	int32_t r = upper(it->first);
	int32_t w = it->second;
	assert(l <= p);
	assert(r >= p);
	if(l == p || r == p) return 0;
	imap -= make_pair(ROI(l, r), w);
	imap += make_pair(ROI(l, p), w);
	imap += make_pair(ROI(p, r), w);
	return 0;
}

int compute_overlap(const split_interval_map &imap, int32_t p)
{
	SIMI it = imap.find(p);
	if(it == imap.end()) return 0;
	return it->second;
}

SIMI locate_right_iterator(const split_interval_map &imap, int32_t x)
{
	return imap.upper_bound(ROI(x - 1, x));
}

SIMI locate_left_iterator(const split_interval_map &imap, int32_t x)
{
	SIMI it = imap.lower_bound(ROI(x - 1, x));
	if(it == imap.end() && it == imap.begin()) return it;
	if(it == imap.end()) it--;

	while(upper(it->first) > x)
	{
		if(it == imap.begin()) return imap.end();
		it--;
	}
	return it;
}

PSIMI locate_boundary_iterators(const split_interval_map &imap, int32_t x, int32_t y)
{
	SIMI lit, rit;
	lit = locate_right_iterator(imap, x);
	if(lit == imap.end() || upper(lit->first) > y) lit = imap.end();

	rit = locate_left_iterator(imap, y);
	if(rit == imap.end() || lower(rit->first) < x) rit = imap.end();

	if(lit == imap.end()) assert(rit == imap.end());
	if(rit == imap.end() && lit != imap.end()) 
	{
		//printf("x = %d, y = %d, lit = [%d, %d)\n", x, y, lower(lit->first), upper(lit->first));
		assert(lit == imap.end());
	}

	return PSIMI(lit, rit); 
}

int32_t compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int32_t x = it->second;
		if(x > s) s = x;
	}

	if(q != imap.end())
	{
		int32_t x = q->second;
		if(x > s) s = x;
	}

	return s;
}

int32_t compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int l = lower(it->first);
		int u = upper(it->first);
		assert(u > l);

		//printf(" AA add [%d, %d) : %d\n", lower(it->first), upper(it->first), it->second);

		s += (u - l) * it->second;
	}
	if(q != imap.end())
	{
		//printf(" BB add [%d, %d) : %d\n", lower(q->first), upper(q->first), q->second);
		s += (upper(q->first) - lower(q->first)) * q->second;
	}
	return s;
}

int evaluate_rectangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev, double &max)
{
	ave = 0;
	dev = 1.0;

	PSIMI pei = locate_boundary_iterators(imap, ll, rr);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == imap.end()) return 0;
	if(rit == imap.end()) return 0;

	ave = 1.0 * compute_sum_overlap(imap, lit, rit) / (rr - ll);
	//printf("compute average %d-%d = %.2lf\n", ll, rr, ave);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave) * (it->second - ave) * (upper(it->first) - lower(it->first));

		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	//if(dev < 1.0) dev = 1.0;

	max = 1.0 * compute_max_overlap(imap, lit, rit);
	return 0;
}
