/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __INTERVAL_MAP_H__
#define __INTERVAL_MAP_H__

// boost::interval map
#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"

#include <vector>

using namespace boost;
using namespace std;

typedef icl::right_open_interval<int32_t> ROI;
typedef std::set<int> SI;
typedef icl::interval<int32_t>::type interval32;

// join interval map
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> join_interval_map;
typedef join_interval_map::const_iterator JIMI;
typedef pair<JIMI, JIMI> PJIMI;

// split interval map
typedef icl::split_interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map;
typedef split_interval_map::const_iterator SIMI;
typedef pair<SIMI, SIMI> PSIMI;

// join interval map with associated sets
typedef icl::interval_map<int32_t, SI> interval_set_map;
typedef interval_set_map::const_iterator ISMI;
typedef pair<ISMI, ISMI> PISMI;

// if p is inside an interval, split this interval into 2
int create_split(split_interval_map &imap, int32_t p);

// return the overlap at position p
int compute_overlap(const split_interval_map &imap, int32_t p);

// find the leftmost iterator whose upper posistion <= x
SIMI locate_right_iterator(const split_interval_map &imap, int32_t x);
ISMI locate_right_iterator(const interval_set_map &ism, int32_t x);

// find the rightmost interval whose lower position >= x
SIMI locate_left_iterator(const split_interval_map &imap, int32_t x);
ISMI locate_left_iterator(const interval_set_map &ism, int32_t x);

// locate boundary iterators
PSIMI locate_boundary_iterators(const split_interval_map &imap, int32_t x, int32_t y);
PISMI locate_boundary_iterators(const interval_set_map &ism, int32_t x, int32_t y);

// return the sum of the lengths of intervals from p to q (include q)
int compute_coverage(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the maximum overlap of the intervals from p to q (include q)
int compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the sum of the overlap of the intervals from p to q (include q)
int compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

set<int> get_overlapped_set(const interval_set_map &ism, int32_t x, int32_t y);

// evaluate a region
int evaluate_rectangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev, double &max);
int evaluate_rectangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev);
int evaluate_triangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev);

// print
int print_interval_set_map(const interval_set_map &ism);
int print_split_interval_map(const split_interval_map &ism);

// testing
int test_split_interval_map();
int test_interval_set_map();

#endif
