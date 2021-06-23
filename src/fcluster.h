/*
Part of Coral
(c) 2019 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FRAGMENT_CLUSTER_H__
#define __FRAGMENT_CLUSTER_H__

#include <vector>
#include "path.h"
#include "fragment.h"

using namespace std;

class fcluster
{
public:
	vector<fragment*> fset;			// set of fragments in this cluster
	int type;						// for multiple uses
	vector<int> v0;					// vlist for closed fragments
	vector<int> v1;					// vlist for mate1
	vector<int> v2;					// vlist for mate2
	vector< vector<int> > phase;	// possible phasing w.r.t. reference
	vector<int> count;				// phase count

public:
	int clear();
	int print(int k) const;
	int add_phase(const vector<int> &v);
	const vector<int> & get_vlist() const;
};

bool compare_fcluster(const fcluster &fx, const fcluster &fy);
bool compare_fcluster_v1_v2(const fcluster &fx, const fcluster &fy);

#endif
