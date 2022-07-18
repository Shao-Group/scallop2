/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include "junction.h"
#include "config.h"

junction::junction()
{}

junction::junction(int64_t _p)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	junc_type = 1;
	strand = '.';
	lexon = -1;
	rexon = -1;
	lregion = -1;
	rregion = -1;
	nm = 0;
}

junction::junction(int64_t _p, int _c)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	junc_type = 1;
	strand = '.';
	lexon = -1;
	rexon = -1;
	lregion = -1;
	rregion = -1;
	nm = 0;
}

junction::junction(const junction &sp)
{
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	junc_type = sp.junc_type;
	lexon = sp.lexon;
	rexon = sp.rexon;
	strand = sp.strand;
	lregion = sp.lregion;
	rregion = sp.rregion;
	nm = sp.nm;
}

bool junction::operator<(const junction &x) const
{
	if(lpos <= x.lpos) return true;
	else return false;
}

int junction::print(const string &chrm, int index) const
{
	printf("junction %d: type = %d, region = %s:%d-%d, region = %d -> %d, pexon = %d -> %d, length = %d, count = %d, strand = %c, nm = %d\n", 
			index, junc_type, chrm.c_str(), lpos, rpos, lregion, rregion, lexon, rexon, rpos - lpos, count, strand, nm);

	return 0;
}

bool junction_cmp_length(const junction &x, const junction &y)
{
	int32_t p1 = x.rpos - x.lpos;
	int32_t p2 = y.rpos - y.lpos;
	if(p1 < p2) return true;
	else return false;
}
