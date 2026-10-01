/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __TRANSCRIPT_MATCH_H__
#define __TRANSCRIPT_MATCH_H__

#include <set>
#include <vector>

#include "bundle_base.h"
#include "transcript_index.h"

using namespace std;

enum transcript_match_failure
{
	TRANSCRIPT_MATCHED,
	TRANSCRIPT_CHROMOSOME_MISMATCH,
	TRANSCRIPT_STRAND_MISMATCH,
	TRANSCRIPT_INVALID_EXONIC_LENGTH,
	TRANSCRIPT_INSUFFICIENT_EXON_OVERLAP,
	TRANSCRIPT_INSUFFICIENT_SPLICING_POSITION_OVERLAP
};

class transcript_match
{
public:
	transcript_match();
	const char* failure_reason() const;

public:
	bool assigned;
	transcript_match_failure failure;
	int shared_exonic_length;
	int transcript_exonic_length;
	int shared_splicing_positions;
	int transcript_splicing_positions;
};

class transcript_matcher
{
public:
	transcript_matcher(const bundle_base &bb);
	transcript_matcher(const string &chrm, char strand, const vector<PI32> &exons,
			const set<int32_t> &splicing_positions);

private:
	string chrm;
	char strand;
	join_interval_map exons;
	set<int32_t> splicing_positions;

public:
	transcript_match match(const indexed_transcript &t, double min_exon_overlap,
			double min_splicing_position_overlap) const;
	const string& chromosome() const;
	char get_strand() const;

private:
	int build(const bundle_base &bb);
	int build_exons(const split_interval_map &mmap);
	int build_exons(const vector<PI32> &v);
	int compute_shared_exonic_length(const transcript &t) const;
};

#endif
