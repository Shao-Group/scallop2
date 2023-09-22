/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FRAGMENT_H__
#define __FRAGMENT_H__

#include <vector>
#include <stdint.h>

#include "hit.h"
#include "path.h"

using namespace std;

class fragment
{
public:
	fragment(hit *x1, hit *x2);

public:
	hit* h1;			// list of first mate
	hit* h2;			// list of second mate

	int frag_type;		//1 if normal fragment from h1p to h2 or h1 to h2p, 2 if from h2 to h1s or h2s to h1
	int is_compatible; //1: comaptible 3 segments case h1 has a supple, 2: compatible 3 segments case h2 has a supple
	int pi; 			//partner fragment index for circ RNA building
	int fidx;			//own fragments index
	int fake_hit_index; //used to keep track of fragment for which fake hit is created, stores index of fake_hit
	bool HS_frag;        //used to keep track whether this is a non chimeric frag with HS on both sides

	int cnt;			// count of the equal hits
	int32_t lpos;		// equals to hits[k1].pos
	int32_t rpos;		// equals to hits[k2].rpos
	int32_t k1l;		// k1-left-outside length
	int32_t k1r;		// k1-right-outside length
	int32_t k2l;		// k2-left-outside length
	int32_t k2r;		// k2-right-outside length
	bool b1;			// whether left mate can be shorten
	bool b2;			// whether right mate can be shorten
	vector<path> paths;	// possible connecting paths

	int type;		// 0: paird-end, 1: UMI-paird 2: both
	int next;		// index for next fragments in UMI-linked

public:
	bool equal(const fragment &f) const;
	int append(const fragment &f);
	int print(int index);
	int set_paired(bool b);
	int set_bridged(bool b);
	int clear();
};

bool compare_fragment(const fragment &f1, const fragment &f2);

#endif
