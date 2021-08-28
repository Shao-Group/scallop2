/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE_BRIDGE_H__
#define __BUNDLE_BRIDGE_H__

#include "bundle_base.h"
#include "junction.h"
#include "region.h"
#include "fragment.h"
#include "transcript.h"

using namespace std;

class bundle_bridge
{
public:
	bundle_bridge(bundle_base &bb);
	virtual ~bundle_bridge();

public:
	bundle_base &bb;					// input bundle base
	set<string> breads;					// bridged reads
	vector<fragment> fragments;			// to-be-filled fragments
	vector<junction> junctions;			// splice junctions
	vector<region> regions;				// pexons
	vector<transcript> ref_trsts;		// overlaped genes in reference
	vector< vector<int> > ref_phase;	// phasing paths for ref transcripts
	vector< vector<PI> > ref_index;		// the set of trsts that contain each region

	vector< vector<int> > umiLink;			// umi linked list: fragments index

public:
	int build();
	int print(int index);
	int32_t compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v);
	vector<int32_t> build_accumulate_length(const vector<int> &v);
	vector<int32_t> get_aligned_intervals(fragment &fr);
	vector<int32_t> get_splices(fragment &fr);

public:
	int build_junctions();
	int extend_junctions();
	int build_regions();

	int remove_tiny_boundary();
	int build_fragments();
	int group_fragments();

	int align_hits_transcripts();
	int align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &v);
	int align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &v);
	int index_references();
	int locate_region(int32_t x);
};

#endif
