/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BRIDGE_H__
#define __BRIDGE_H__

#include <stdint.h>
#include <string>

using namespace std;

class junction
{
public:
	junction();
	junction(int64_t _p);
	junction(int64_t _p, int32_t _c);
	junction(const junction &p);

	bool operator<(const junction &x) const;

public:
	// add a flag to mark if this junction is a normal one or BSJ; 
	// int flag;
	int32_t lpos;		// left position [left, right)
	int32_t rpos;		// right position
	int count;			// number of hits having this splice junction
	int junc_type;		//junc_type = 1 if normal junction, 2 if BSJ.
	char strand;		// strandness of this junction
	int nm;				// total mismatch
	char boundary_match; //L if junc lpos matches hits with H/S rpos, R if junc rpos matches hits with H/s lpos

	int lexon;			// pexon index corresponds to lpos
	int rexon;			// pexon index corresponds to rpos
	int lregion;		// region index corresponds to lpos
	int rregion;		// region index corresponds to rpos

public:
	int print(const string &chrm, int index) const;
};

bool junction_cmp_length(const junction &x, const junction &y);

#endif
