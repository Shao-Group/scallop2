/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define EMPTY -1
#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

#define EMPTY_VERTEX -9

#define TRANSCRIPT_COUNT_ADD_COVERAGE_ADD 1
#define TRANSCRIPT_COUNT_ADD_COVERAGE_NUL 2
#define TRANSCRIPT_COUNT_ADD_COVERAGE_MAX 3
#define TRANSCRIPT_COUNT_ADD_COVERAGE_MIN 4

//// parameters

//for circRNA
extern int read_length;
extern int min_pathscore;
extern int max_softclip_to_junction_gap;
extern int min_junction_count;
extern double min_junction_count_ratio;
extern int max_circ_vsize;
extern int max_exon_length;
extern int same_chain_circ_end_diff;
extern int alignment_boundary_error;
extern double max_fset_score;
extern int min_soft_clip_len; 

// for bam file and reads
extern int min_flank_length;
extern int max_num_cigar;
extern int max_edit_distance;
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern int min_num_splices_in_bundle;
extern uint32_t min_mapping_quality;
extern int32_t min_splice_boundary_hits;
extern bool uniquely_mapped_only;
extern bool use_second_alignment;
extern int library_type;

// for preview
extern bool preview_only;
extern int max_preview_reads;
extern int max_preview_spliced_reads;
extern int min_preview_spliced_reads;
extern double preview_infer_ratio;
extern double insertsize_ave;
extern double insertsize_std;
extern int insertsize_median;
extern int insertsize_low;
extern int insertsize_high;
extern double insertsize_low_percentile;
extern double insertsize_high_percentile;
extern vector<double> insertsize_profile;

// for bridging
extern double min_bridging_score;
extern int max_num_path_nodes;
extern int dp_solution_size;
extern int dp_stack_size;
extern bool use_overlap_scoring;
extern int32_t max_clustering_flank; 
extern int32_t flank_tiny_length;
extern double flank_tiny_ratio;

// for identifying subgraphs
extern int32_t min_subregion_gap;
extern int32_t min_subregion_len;
extern int32_t min_subregion_max;
extern double min_subregion_ave;
extern double min_guaranteed_edge_weight;

// for filtering transcripts
extern double min_transcript_numreads;
extern double min_transcript_coverage;
extern double min_single_exon_coverage;
extern double min_transcript_coverage_ratio; 
extern int min_transcript_length_base;
extern int min_transcript_length_increase;
extern int min_exon_length;
extern int max_num_exons;

// input and output
extern string algo;
extern string input_file;
extern string fasta_file;
extern string ref_file;
extern string ref_file1;
extern string ref_file2;
extern string output_file;
extern string output_file1;
extern string output_circ_file;
extern string cirifull_file;

// for controling
extern int batch_bundle_size;
extern int verbose;
extern string version;

//for extracting BSJ
extern int BSJ_threshold;
extern map<string, int> frag2graph_freq; //adding data structure to keep frequency profile of fragment to splice graph align cases
extern int h1_supp_count; //keeps count of the number of h1 hits in ffragment with supple
extern int h2_supp_count; //keeps count of the number of h2 hits in ffragment with supple
extern map<string, int> circ_frag_bridged_freq; //keeps count of fragment pair bridged cases, TT,TF,FT,FF

// parse arguments
int print_command_line(int argc, const char ** argv);
int parse_arguments(int argc, const char ** argv);
int print_parameters();
int print_copyright();
int print_help();

#endif
