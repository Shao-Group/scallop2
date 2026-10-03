/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PARTIAL_EXON_H__
#define __PARTIAL_EXON_H__

#include <stdint.h>
#include <vector>
#include <string>

using namespace std;

class partial_exon
{
public:
	partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype);

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary

	int rid;						// parental region id
	int pid;						// index in the parental pexons
	int type;						// label
	double ave;						// average abundance
	double max;						// maximum abundance
	double dev;						// standard-deviation of abundance
	int newly_added_length;			// bases supplied only by assigned transcripts
	double indel_sum_cov;			// summed insertion/deletion coverage
	double indel_ratio;				// indel coverage divided by exon coverage
	int left_indel;					// distance to the nearest indel from the left
	int right_indel;					// distance to the nearest indel from the right

public:
};

#endif
