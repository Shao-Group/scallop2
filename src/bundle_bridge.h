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
#include "circular_transcript.h"
#include "htslib/faidx.h"

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
	vector<fragment> circ_fragments;	// to-be-filled fragments
	faidx_t *fai;						//pointer to fetch fasta seq from region

	vector<pair<fragment,fragment>> circ_fragment_pairs;	//bridged fragment pairs for circular RNA
	vector<circular_transcript> circ_trsts; //a vector of circular transcripts class objs, with duplicates
	vector<circular_transcript> circ_trsts_HS; ////a vector of circular transcripts from all possible H/S reads, with duplicates

	vector<string> HS_both_side_reads; //for statistics of RO reads from CIRI-full
	vector<string> chimeric_reads; //for statistics of RO reads from CIRI-full
	int RO_count;

	vector<junction> junctions;			// splice junctions
	map<int64_t, char> junc_map;		// map junction to strandness
	vector<region> regions;				// pexons
	vector<partial_exon> pexons;		// partial exons
	vector<transcript> ref_trsts;		// overlaped genes in reference
	vector< vector<int> > ref_phase;	// phasing paths for ref transcripts
	vector< vector<PI> > ref_index;		// the set of trsts that contain each region

	vector< vector<int> > umiLink;		// umi linked list: fragments index

public:
	int build(map <string, int> RO_reads_map, faidx_t *_fai);
	int print(int index);
	int32_t compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v);
	vector<int32_t> build_accumulate_length(const vector<int> &v);
	vector<int32_t> get_aligned_intervals(fragment &fr);
	vector<int32_t> get_splices(fragment &fr);

public:
	int set_hits_RO_parameter(map <string, int> RO_reads_map);
	int set_chimeric_cigar_positions();
	int build_supplementaries();
	int build_junctions();
	int extend_junctions();
	int build_regions();
	int build_partial_exons();
	int build_fragments();
	int get_frags_with_HS_on_both_sides();
	int get_RO_frags_with_HS();
	string get_fasta_seq(int32_t pos1, int32_t pos2);
	int min(int x, int y, int z);
	int get_edit_distance(string s, string t);
	int get_more_chimeric();
	int fix_alignment_boundaries();

	int build_circ_fragments();
	int extract_all_non_supple_HS_hits();
	int extract_nonsupple_HS_hits();
	int extract_RO_circRNA();
	int extract_circ_fragment_pairs();
	int print_circ_fragment_pairs();
	int join_circ_fragment_pairs();
	int join_circ_fragment_pair(pair<fragment,fragment> &fr_pair, int ex1, int ex2);
	int print_circRNAs();
	char infer_circ_strand(const vector<int> &p);
	
	int group_fragments();

	int align_hits_transcripts();
	int align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &v);
	int align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &v);
	int remove_tiny_boundaries();
	int remove_tiny_boundary(hit &h);
	int set_fragment_lengths();
	int set_fragment_length(fragment &fr);
	int index_references();
	int locate_region(int32_t x);
};

#endif
