/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "interval_map.h"
#include "bundle_base.h"
#include "bundle_bridge.h"
#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "path.h"
#include "gene.h"
#include "transcript.h"
#include "berth.h"

using namespace std;

class tss_tes
{
public:
	tss_tes(int type = 0);
	tss_tes(int type, int32_t pos, int weight_berth, int weight_sg);
	void calculate_junction_cnt(vector<junction> &sorted_junctions_start, vector<junction> &sorted_junctions_end);
	void calculate_clip_length(vector<hit> &hits);
	void calculate_anchor_features(vector<hit> &hits);
public:
	int type; // 0: TSS, 1: TES
	int32_t pos; // position
	int weight_berth; // weight from berth
	int weight_sg; // weight from splice graph
	int read_density; // count of read starting/ending in the neighborhood
	float leading_clip_length; // average leading soft clip length of reads starting/ending in the neighborhood
	float trailing_clip_length; // average trailing soft clip length of reads starting/ending in the neighborhood
	int junction_start_cnt; // count of junction starting in the neighborhood
	int junction_end_cnt; // count of junction ending in the neighborhood
	int junction_cross_cnt; // count of junction crossing the neighborhood
	int left_anchor_cnt; // count of anchors in leftclip of reads in the neighborhood
	int right_anchor_cnt; // count of anchors in rightclip of reads in the neighborhood
	int left_anchor_padding_mean; // average anchor padding length of leftclip of reads starting/ending in the neighborhood
	int right_anchor_padding_mean; // average anchor padding length of rightclip reads starting/ending in the neighborhood
	// int anchor_padding_stddev; // stddev of anchor padding length of reads starting/ending in the neighborhood
};

class bundle
{
public:
	bundle(bundle_base &bb);
	bundle();
	virtual ~bundle();

public:
	bundle_base &bb;				// input bundle base	
	berth       bth;					

	// TODO: we aim for comment out the br
	// bundle_bridge br;				// contains fragments
	split_interval_map fmap;		// matched interval map
	vector<junction> junctions;		// splice junctions
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	vector<bool> regional;			// if a pe is regional
	split_interval_map pmap;		// partial exon map
	splice_graph gr;				// splice graph
	splice_graph new_gr;				// splice graph
	hyper_set hs;					// hyper set w/o unreliable vertices
	// hyper_set hs2;					// hyper set w/ unreliaable vertices
	// hyper_set hs_majority;          // hyper set with majority voting	
	vector<pair<int32_t,int>> tss_list_sg;		// TSS as (pos, cnt)
	vector<pair<int32_t,int>> tes_list_sg;		// TES as (pos, cnt)
	vector<tss_tes> tss_merged;		// merged TSS
	vector<tss_tes> tes_merged;		// merged TES
	
	vector<int32_t> left_anchorpos_list;
	vector<int32_t> right_anchorpos_list;

public:
	virtual int build(int mode, bool revise);
	int output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const;	
	int output_transcripts(gene &gn, const vector<path> &p, const string &gid) const;	
	int output_transcripts(vector<transcript> &trsts, const vector<path> &p, const string &gid) const;	
	int output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const;	
	int output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const;	
	int count_junctions() const;
	int print(int index);

public:
	int prepare();

	// check and init
	int check_left_ascending();
	int check_right_ascending();
	int compute_strand();

	// splice graph
	int build_intervals();
	int build_junctions();
	int correct_junctions();
	int build_regions();
	int build_partial_exons();
	int link_partial_exons();
	int build_splice_graph(int mode);
	int rebuild_splice_graph_using_refined_hyper_set(int mode);
	int build_partial_exon_map();
	int locate_left_partial_exon(int32_t x);
	int locate_right_partial_exon(int32_t x);
	vector<int> align_hit(hit &h);
	vector<int> align_fragment(fragment &f);
	void print_fmap(); // Debug

	// TSS-TES
	int build_tss_tes();
	void write_tss_tes();
	void write_tss_tes_features();
	int merge_tss_tes();
	int build_anchors();

	// revise splice graph
	VE compute_maximal_edges();
	int revise_splice_graph();
	int refine_splice_graph();
	int refine_modified_splice_graph();
	bool keep_surviving_edges();
	bool extend_boundaries();
	bool extend_start_boundaries();
	bool extend_end_boundaries();
	bool remove_small_junctions();
	bool remove_small_exons();
	bool remove_inner_boundaries();
	bool remove_intron_contamination();
	bool remove_false_boundaries();
	bool tackle_false_boundaries();
	int find_contamination_chain();

	// hyper set
	int build_hyper_set();
	int refine_hyper_set();
	int build_majority_hyper_set();
};



#endif
