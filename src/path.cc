/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "path.h"
#include "util.h"

#include <cassert>
#include <cstdio>

path::path()
{
	v.clear();
	type = 0;
	abd = 0;
	reads = 0;
	score = 0;
	length = 0;
	fcindex = -1;
	ex1 = 0;
	ex2 = 0;
}

path::~path()
{}

bool path::operator<(const path& p) const
{
	if(p.length < length) return true;
	else return false;
}

int path::clear()
{
	type = 0;
	ex1 = 0;
	ex2 = 0;
	v.clear();
	abd = 0;
	reads = 0;
	score = 0;
	length = 0;
	fcindex = -1;
	return 0;
}

int path::print(int index) const
{
	if(v.size() == 0) return 0;
	printf("path %d: type = %d, abundance = %.2lf, length = %d, reads = %.2lf, vertices = ", index, type, abd, length, reads);
	for(int i = 0; i < v.size() - 1; i++)
	{
		printf("%d, ", v[i]);
	}
	printf("%d\n", v[v.size() - 1]);
	return 0;
}

int path::print_bridge(int index) const
{
	printf("path %d: type = %d, length = %d, score = %.2lf, vertices = ( ", index, type, length, score);
	printv(v);
	printf(")\n");
	return 0;
}

vector<int> path::index(int n) const
{
	vector<int> vv;
	vv.resize(n, -1);
	for(int i = 1; i < v.size(); i++)
	{
		int s = v[i - 1];
		int t = v[i];
		assert(s >= 0 && s < n);
		assert(t >= 0 && t < n);
		vv[s] = t;
	}
	return vv;
}

bool compare_path_abundance(const path &p1, const path &p2)
{
	if(p1.abd > p2.abd) return true;
	else return false;
}

bool compare_path_vertices(const path &p1, const path &p2)
{
	for(int k = 0; k < p1.v.size() && k < p2.v.size(); k++)
	{
		if(p1.v[k] < p2.v[k]) return true;
		if(p1.v[k] > p2.v[k]) return false;
	}
	if(p1.v.size() < p2.v.size()) return true;
	if(p1.v.size() > p2.v.size()) return false;
	return false;
}

bool compare_path_score(const path &p1, const path &p2)
{
	if(p1.score > p2.score) return true;
	else return false;
}
