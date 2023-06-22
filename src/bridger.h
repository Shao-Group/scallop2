/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BRIDGER_H__
#define __BRIDGER_H__

#include "bundle_bridge.h"
#include "fcluster.h"

using namespace std;

class entry
{
public:
	vector<int> stack;
	int32_t length;
	int trace1;
	int trace2;

public:
	int print();
};

bool entry_compare(const entry &x, const entry &y);

class bridger
{
public:
	bridger(bundle_bridge *b);

public:
	bundle_bridge *bd;				// parent bundle
	vector<path> pnodes;			// path nodes (not used)
	vector< map<int, int> > jsetx;	// junction graph (out) 
	vector< map<int, int> > jsety;	// junction graph (in)
	vector< map<int, int> > psetx;	// path graph (out) (not used)
	vector< map<int, int> > psety;	// path graph (in) (not used)
	int max_pnode_length;			// kmer size
	int32_t length_median;
	int32_t length_low;	//DISTRBN OF FRGAMNET length 0
	int32_t length_high; //DISTRBN OF FRGAMNET length 10000


public:
	int bridge();
	int bridge_clip(int32_t p1, int32_t p2, circular_transcript &circ);
	int pick_bridge_path();
	int print();

public:
	int bridge_overlapped_fragments();
	int bridge_overlapped_fragment(fragment &fr, int ex1, int ex2);

	int bridge_phased_fragments(vector<fcluster> &fclusters);
	int phase_cluster(fcluster &fc);
	int bridge_phased_cluster(fcluster &fc);

	int remove_tiny_boundary();

	int build_junction_graph();
	int bridge_hard_fragments(vector<fcluster> &open);
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int evaluate_bridging_path(const vector<int> &pb);
	int determine_overlap(const vector<int> &vx, const vector<int> &vy, PI &p);
	int determine_overlap1(const vector<int> &vx, const vector<int> &vy, PI &p);
	bool determine_identical(const vector<int> &vx, const vector<int> &vy, int x1, int x2, int y1, int y2);

	int build_overlap_index();
	int bridge_tough_fragments();
	int dynamic_programming(int k1, int k2, vector<int> &trace, vector< vector<int> > &table_cov, vector<int32_t> &table_len);
	int compare_stack(const vector<int> &x, const vector<int> &y);
	vector<int> update_stack(const vector<int> &v, int s);

	vector<int> trace_back(int k1, int k2, const vector<int> &trace);
	vector<int> get_bridge(const vector<int> &vv, const vector<int> &v1, const vector<int> &v2);
	int32_t get_extended_length1(int k2, int p1, int p2);
	int32_t get_extended_length2(int k1, int p1, int p2);
	vector<int> get_suffix(const vector<int> &v);
	vector<int> get_prefix(const vector<int> &v);

	int cluster_open_fragments(vector<fcluster> &fclusters);
	int build_path_nodes();
	int build_path_nodes(int max_len);
	int build_path_nodes(int low, int high);
	int build_path_nodes(map<vector<int>, int> &m, const vector<int> &v, int cnt);
	int add_consecutive_path_nodes();
	int adjust_path_score(path &p);

	int set_thresholds();
	int update_length();
	int filter_paths();
	int get_paired_fragments();
	vector<int> get_bridged_fragments_type();
};

bool compare_fragment_v1(fragment *f1, fragment *f2);
bool compare_fragment_v2(fragment *f1, fragment *f2);
bool compare_fragment_v3(fragment *f1, fragment *f2);
bool compare_fragment_v3_flank(fragment *f1, fragment *f2);
bool compare_fragment_path(fragment *f1, fragment *f2);
bool check_suffix(const vector<int> &vx, const vector<int> &vy);

#endif
