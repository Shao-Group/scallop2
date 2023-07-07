/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __PATH_H__
#define __PATH_H__

#include <vector>
#include <stdint.h>
#include "region.h"

using namespace std;

class path
{
public:
	path();
	~path();

public:
	bool operator< (const path& p) const;

public:
	int type;			// 1: within normal range of insertsize; 2: outside normal range of insertsize
	int ex1;
	int ex2;
	vector<int> v;
	vector<int32_t> acc;
	int32_t length;
	int fcindex;
	double abd;
	double prlen;
	double score;
	double reads;

	int nf;				// 1: empty = true, non full length path; 0: empty = false, full length path

	vector<region> path_regions;
	vector<region> merged_regions;
	vector<pair<int32_t,int32_t>> junc_regions;

public:
	int clear();
	int print(int index) const;
	int print_bridge(int index) const;
	vector<int> index(int n) const;
};

bool compare_path_abundance(const path &p1, const path &p2);
bool compare_path_vertices(const path &p1, const path &p2);
bool compare_path_score(const path &p1, const path &p2);

#endif
