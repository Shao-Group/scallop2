/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __TRANSCRIPT_MATCH_H__
#define __TRANSCRIPT_MATCH_H__

#include <set>
#include <vector>

#include "splice_graph.h"
#include "transcript_index.h"

using namespace std;

class transcript_match
{
public:
	transcript_match();

public:
	bool assigned;
	int shared_exonic_length;
	int transcript_exonic_length;
	int shared_junctions;
	int transcript_junctions;
};

class transcript_matcher
{
public:
	transcript_matcher(const splice_graph &gr);
	transcript_matcher(const string &chrm, char strand, const vector<PI32> &exons,
			const set<int64_t> &junctions);

private:
	string chrm;
	char strand;
	join_interval_map exons;
	set<int64_t> junctions;

public:
	transcript_match match(const indexed_transcript &t, double min_exon_overlap,
			double min_junction_overlap) const;

private:
	int build(const splice_graph &gr);
	int build_exons(const vector<PI32> &v);
	int compute_shared_exonic_length(const transcript &t) const;
};

#endif
