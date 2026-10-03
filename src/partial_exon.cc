/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "partial_exon.h"
#include "util.h"
#include <cstdio>

partial_exon::partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype)
	: lpos(_lpos), rpos(_rpos), ltype(_ltype), rtype(_rtype)
{
	type = 0;
	rid = -1;
	pid = -1;
	newly_added_length = 0;
	ave = 0;
	max = 0;
	dev = 1;
	indel_sum_cov = 0;
	indel_ratio = 0;
	left_indel = -1;
	right_indel = -1;
}
