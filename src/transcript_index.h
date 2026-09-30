/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __TRANSCRIPT_INDEX_H__
#define __TRANSCRIPT_INDEX_H__

#include <map>
#include <set>
#include <string>
#include <vector>

#include "genome.h"
#include "interval_map.h"

using namespace std;

typedef icl::split_interval_map<int32_t, set<int>, icl::partial_absorber,
		less, icl::inplace_plus, icl::inter_section, ROI> transcript_interval_map;
typedef transcript_interval_map::const_iterator TIMI;

class indexed_transcript
{
public:
	indexed_transcript(const transcript &t);

public:
	transcript trst;
	PI32 bounds;
	vector<PI32> junctions;
	int exonic_length;
};

class transcript_index
{
public:
	transcript_index();
	transcript_index(const genome &gm);

public:
	vector<indexed_transcript> transcripts;
	map<pair<string, char>, transcript_interval_map> index;

public:
	int build(const genome &gm);
	vector<int> query(const string &chrm, char strand, int32_t lpos, int32_t rpos) const;
	const indexed_transcript& get(int k) const;
	int size() const;
};

#endif
