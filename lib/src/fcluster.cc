/*
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "fcluster.h"
#include "util.h"
#include <cstdio>

int fcluster::clear()
{
	fset.clear();
	v0.clear();
	v1.clear();
	v2.clear();
	phase.clear();
	count.clear();
	return 0;
}

int fcluster::add_phase(const vector<int> &v)
{
	for(int k = 0; k < phase.size(); k++)
	{
		if(phase[k] == v)
		{
			count[k]++;
			return 0;
		}
	}
	phase.push_back(v);
	count.push_back(1);
	return 0;
}

const vector<int> & fcluster::get_vlist() const
{
	if(type == 0) return v0;
	if(type == 1) return v1;
	if(type == 2) return v2;
	assert(false);
}


int fcluster::print(int index) const
{
	printf("fcluster %d: type = %d, #fragments = %lu, #phase = %lu, ", index, type, fset.size(), phase.size());

	printf("  v1 = ( ");
	printv(v1);
	printf("), v2 = ( ");
	printv(v2);
	printf(")\n");

	for(int k = 0; k < phase.size(); k++)
	{
		printf("  count = %d, phase %d = (", count[k], k);
		printv(phase[k]);
		printf(")\n");
	}

	return 0;
}

bool compare_fcluster(const fcluster &fx, const fcluster &fy)
{
	const vector<int> &vx = fx.get_vlist();
	const vector<int> &vy = fy.get_vlist();

	for(int k = 0; k < vx.size() && k < vy.size(); k++)
	{
		if(vx[k] < vy[k]) return true;
		if(vx[k] > vy[k]) return false;
	}

	if(vx.size() < vy.size()) return true;
	if(vx.size() > vy.size()) return false;

	if(fx.type < fy.type) return true;
	if(fx.type > fy.type) return false;

	return false;
}

bool compare_fcluster_v1_v2(const fcluster &fx, const fcluster &fy)
{
	for(int k = 0; k < fx.v1.size() && k < fy.v1.size(); k++)
	{
		if(fx.v1[k] < fy.v1[k]) return true;
		if(fx.v1[k] > fy.v1[k]) return false;
	}

	if(fx.v1.size() < fy.v1.size()) return true;
	if(fx.v1.size() > fy.v1.size()) return false;

	// NOTE: reverse sorting for v2 (not for v1)
	for(int k = 0; k < fx.v2.size() && k < fy.v2.size(); k++)
	{
		if(fx.v2[k] > fy.v2[k]) return true;
		if(fx.v2[k] < fy.v2[k]) return false;
	}

	if(fx.v2.size() > fy.v2.size()) return true;
	if(fx.v2.size() < fy.v2.size()) return false;

	return false;
}
