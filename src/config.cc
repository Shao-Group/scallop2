/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//// parameters


//for circRNA
int read_length = 100;
int min_pathscore = 1; //remove paths <= this of higher exists in bridger
int max_softclip_to_junction_gap = 100000; //discard junction for checking seq match with softclip if > this 
int min_junction_count = 10; // filter junction
double min_junction_count_ratio = 0.05; //filter junction
int max_circ_vsize = 15; //discard circRNA if vsize > this
int max_single_exon_length = 2000; //discard single exon circRNA if indivudal exon length > this
int max_multi_exon_length = 1000; //discard multi exon circRNA if indivudal exon length > this
int same_chain_circ_end_diff = 50; //allow this much diff for merging circRNA with same intron chain
int alignment_boundary_error = 5; //allow this much error when fixing alignment boundary
double max_fset_score = 1.5; //bottleneck should be close to number of reads, define close as this
int min_soft_clip_len = 10; // for new chimeric reads soft_len needs to be greater than this to be considered
double min_jaccard = 0.5; //new chimeric similarity threshold

// for bam file and reads
int min_flank_length = 3;
int max_num_cigar = 1000;
int max_edit_distance = 10;
int32_t min_bundle_gap = 100;		
int min_num_hits_in_bundle = 5;	
int min_num_splices_in_bundle = 15;	// not used; accept bundle if #hits with splices is at least this number
uint32_t min_mapping_quality = 1;
int32_t min_splice_boundary_hits = 1;
bool use_second_alignment = false; //change if needed
bool uniquely_mapped_only = false;
int library_type = EMPTY;

// for preview
int max_preview_reads = 2000000;
int max_preview_spliced_reads = 50000;
int min_preview_spliced_reads = 10000;
double preview_infer_ratio = 0.85;
bool preview_only = false;
double insertsize_ave = 300;
double insertsize_std = 50;
int insertsize_median = -1;
int insertsize_low = -1;
int insertsize_high = -1;
double insertsize_low_percentile = 0.005;
double insertsize_high_percentile = 0.998;
vector<double> insertsize_profile;

// for bridging
double min_bridging_score = 0.5;
int max_num_path_nodes = 10000;
int dp_solution_size = 10;
int dp_stack_size = 5;
bool use_overlap_scoring = false;
int32_t max_clustering_flank = 30;
int32_t flank_tiny_length = 15;
double flank_tiny_ratio = 0.4;

// for identifying subgraphs
int32_t min_subregion_gap = 10; //previously = 3
int32_t min_subregion_len = 15;
int32_t min_subregion_max = 3;
double min_subregion_ave = 1.5;
double min_guaranteed_edge_weight = 0.01;

// for selecting paths
double min_transcript_coverage = 1.5;
double min_single_exon_coverage = 20;
double min_transcript_numreads = 10;
int min_transcript_length_base = 150;
int min_transcript_length_increase = 50;
int min_exon_length = 20;
int max_num_exons = 1000;

// input and output
//sinitializing as empty string for verifying absence
string algo = "scallop2";
string input_file = "";
string fasta_file = "";
string ref_file = "";
string ref_file1 = "";
string ref_file2 = "";
string output_file = "";
string output_file1 = "";
string output_circ_file = "";
string cirifull_file = "";

// for controling
int batch_bundle_size = 100;
int verbose = 1;
string version = "v1.1.2";


//for extracting BSJ
int BSJ_threshold = 5;
map<string, int> frag2graph_freq; //adding data structure to keep frequency profile of fragment to splice graph align cases
map<string, int> circ_frag_bridged_freq; //keeps count of fragment pair bridged cases, TT,TF,FT,FF
int h1_supp_count = 0; //keeps count of the number of h1 hits in ffragment with supple
int h2_supp_count = 0; //keeps count of the number of h2 hits in ffragment with supple

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		// necessary ones
		if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-f")
		{
			output_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-c")
		{
			output_circ_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-ro")
		{
			cirifull_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-fa")
		{
			fasta_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--transcript_fragments")
		{
			output_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r")
		{
			ref_file = string(argv[i + 1]);
			i++;
		}

		// internal use
		else if(string(argv[i]) == "-a")
		{
			algo = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r1")
		{
			ref_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r2")
		{
			ref_file2 = string(argv[i + 1]);
			i++;
		}

		// user specified
		else if(string(argv[i]) == "--version")
		{
			printf("%s\n", version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_copyright();
			print_help();
			printf("\n");
			exit(0);
		}
		else if(string(argv[i]) == "--read_length")
		{
			read_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_jaccard")
		{
			min_jaccard = stod(argv[i + 1]);
			printf("min_jaccard=%lf\n",min_jaccard);
			i++;
		}
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_cigar")
		{
			max_num_cigar = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_edit_distance")
		{
			max_edit_distance = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bundle_gap")
		{
			min_bundle_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_hits_in_bundle")
		{
			min_num_hits_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_splices_in_bundle")
		{
			min_num_splices_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_mapping_quality")
		{
			min_mapping_quality = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_boundary_hits")
		{
			min_splice_boundary_hits = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_preview_spliced_reads")
		{
			max_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_preview_spliced_reads")
		{
			min_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview")
		{
			preview_only = true;
		}
		else if(string(argv[i]) == "--max_preview_reads")
		{
			max_preview_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview_infer_ratio")
		{
			preview_infer_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_gap")
		{
			min_subregion_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_len")
		{
			min_subregion_len = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_ave")
		{
			min_subregion_ave = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_max")
		{
			min_subregion_max = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
			if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
		}
		else if(string(argv[i]) == "--min_single_exon_coverage")
		{
			min_single_exon_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_numreads")
		{
			min_transcript_numreads = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_base")
		{
			min_transcript_length_base = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_increase")
		{
			min_transcript_length_increase = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_exon_length")
		{
			min_exon_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_exons")
		{
			max_num_exons = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--library_type")
		{
			string s(argv[i + 1]);
			if(s == "empty") library_type = EMPTY;
			if(s == "unstranded") library_type = UNSTRANDED;
			if(s == "first") library_type = FR_FIRST;
			if(s == "second") library_type = FR_SECOND;
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			if(s == "true") use_second_alignment = true;
			else use_second_alignment = false;
			i++;
		}
		else if(string(argv[i]) == "--uniquely_mapped_only")
		{
			string s(argv[i + 1]);
			if(s == "true") uniquely_mapped_only = true;
			else uniquely_mapped_only = false;
			i++;
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--batch_bundle_size")
		{
			batch_bundle_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bridging_score")
		{
			min_bridging_score = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--dp_solution_size")
		{
			dp_solution_size = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--dp_stack_size")
		{
			dp_stack_size = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_clustering_flank")
		{
			max_clustering_flank = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--flank_tiny_length")
		{
			flank_tiny_length = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--flank_tiny_ratio")
		{
			flank_tiny_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_median")
		{
			insertsize_median = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_low")
		{
			insertsize_low = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--insertsize_high")
		{
			insertsize_high = atof(argv[i + 1]);
			i++;
		}

	}

	/*
	if(min_surviving_edge_weight < 0.1 + min_transcript_coverage) 
	{
		min_surviving_edge_weight = 0.1 + min_transcript_coverage;
		if(min_surviving_edge_weight > 10.0) min_surviving_edge_weight = 10.0;
	}
	*/

	// verify arguments
	if(input_file == "")
	{
		printf("error: input-file is missing.\n");
		exit(0);
	}

	if(output_file == "" && preview_only == false)
	{
		printf("error: output-file is missing.\n");
		exit(0);
	}

	/*if(output_circ_file == "" && preview_only == false)
	{
		printf("error: output-file for circRNA is missing.\n");
		exit(0);
	}*/

	return 0;
}

int print_parameters()
{
	printf("parameters:\n");

	// for bam file and reads
	printf("min_flank_length = %d\n", min_flank_length);
	printf("max_num_cigar = %d\n", max_num_cigar);
	printf("max_edit_distance = %d\n", max_edit_distance);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);

	// for preview
	printf("preview_only = %c\n", preview_only ? 'T' : 'F');
	printf("max_preview_reads = %d\n", max_preview_reads);
	printf("max_preview_spliced_reads = %d\n", max_preview_spliced_reads);
	printf("min_preview_spliced_reads = %d\n", min_preview_spliced_reads);
	printf("preview_infer_ratio = %.3lf\n", preview_infer_ratio);

	// for identifying subgraphs
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_subregion_len = %d\n", min_subregion_len);
	printf("min_subregion_max = %d\n", min_subregion_max);
	printf("min_subregion_ave = %.2lf\n", min_subregion_ave);

	// for input and output
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("ref_file = %s\n", ref_file.c_str());
	printf("ref_file1 = %s\n", ref_file1.c_str());
	printf("ref_file2 = %s\n", ref_file2.c_str());
	printf("output_file = %s\n", output_file.c_str());
	printf("output_file1 = %s\n", output_file1.c_str());
	printf("output_circ_file = %s\n", output_circ_file.c_str());

	// for controling
	printf("library_type = %d\n", library_type);
	printf("use_second_alignment = %c\n", use_second_alignment ? 'T' : 'F');
	printf("uniquely_mapped_only = %c\n", uniquely_mapped_only ? 'T' : 'F');
	printf("verbose = %d\n", verbose);
	printf("batch_bundle_size = %d\n", batch_bundle_size);

	printf("\n");

	return 0;
}

int print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}

int print_help()
{
	printf("\n");
	printf("Usage: scallop2 -i <bam-file> -o <gtf-file> [options]\n");
	printf("\n");
	printf("Options:\n");
	printf(" %-42s  %s\n", "--help",  "print usage of Scallop and exit");
	printf(" %-42s  %s\n", "--version",  "print current version of Scallop and exit");
	printf(" %-42s  %s\n", "--preview",  "determine fragment-length-range and library-type and exit");
	printf(" %-42s  %s\n", "--verbose <0, 1, 2>",  "0: quiet; 1: one line for each graph; 2: with details, default: 1");
	printf(" %-42s  %s\n", "-f/--transcript_fragments <filename>",  "file to which the assembled non-full-length transcripts will be written to");
	printf(" %-42s  %s\n", "--library_type <first, second, unstranded>",  "library type of the sample, default: unstranded");
	printf(" %-42s  %s\n", "--min_transcript_coverage <float>",  "minimum coverage required for a multi-exon transcript, default: 1.5");
	printf(" %-42s  %s\n", "--min_single_exon_coverage <float>",  "minimum coverage required for a single-exon transcript, default: 20");
	printf(" %-42s  %s\n", "--min_transcript_length_increase <integer>",  "default: 50");
	printf(" %-42s  %s\n", "--min_transcript_length_base <integer>",  "default: 150, minimum length of a transcript would be");
	printf(" %-42s  %s\n", "",  "--min_transcript_length_base + --min_transcript_length_increase * num-of-exons");
	printf(" %-42s  %s\n", "--min_mapping_quality <integer>",  "ignore reads with mapping quality less than this value, default: 1");
	printf(" %-42s  %s\n", "--max_num_cigar <integer>",  "ignore reads with CIGAR size larger than this value, default: 1000");
	printf(" %-42s  %s\n", "--min_bundle_gap <integer>",  "minimum distances required to start a new bundle, default: 100");
	printf(" %-42s  %s\n", "--min_num_hits_in_bundle <integer>",  "minimum number of reads required in a gene locus, default: 5");
	printf(" %-42s  %s\n", "--min_flank_length <integer>",  "minimum match length in each side for a spliced read, default: 3");

	return 0;
}

int print_copyright()
{
	printf("Scallop2 %s (c) 2021 Qimin Zhang, Qian Shi, and Mingfu Shao, The Pennsylvania State University\n", version.c_str());
	return 0;
}

