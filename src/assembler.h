/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>
#include "bundle_base.h"
#include "transcript.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "transcript_set.h"
#include "reference.h"

#include "circular_transcript.h"
#include "RO_read.h"
#include "htslib/faidx.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

private:
	samFile *sfn;
	faidx_t *fai;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	reference ref;
	bundle_base bb1;		// +
	bundle_base bb2;		// -
	vector<bundle_base> pool;

	int hid;
	int index;
	bool terminate;
	int qcnt;
	double qlen;
	vector<transcript> trsts;
	vector<transcript> non_full_trsts;

	vector<circular_transcript> circular_trsts; //a vector of circular transcripts class objs from all bundles
	map <string, pair<circular_transcript, int>> circ_trst_map; // a map of distinct circ trsts with circRNA_id as key and the corresponding circRNA object
	map <string, pair<circular_transcript, int>> circ_trst_merged_map; // map with circRNAs having few bp diff ends but same intron chains merged

	vector<circular_transcript> circular_trsts_HS;///a vector of circular transcripts class objs from all HS reads from all bundles
	vector<string> HS_both_side_reads; //for statistics of RO reads from CIRI-full
	vector<string> chimeric_reads; //for statistics of RO reads from CIRI-full
	vector<RO_read> RO_reads; //list of RO reads from CIRI-full simu_ro2_info.list
	map <string, int> RO_reads_map; //map of RO read name concatenated with chrm
	int RO_count;

public:
	int assemble();

private:
	int process(int n);
	int remove_duplicate_circ_trsts();
	vector<string> split_str(string str, string delimiter);
	int print_circular_trsts();
	int assemble(const splice_graph &gr, const hyper_set &hs, transcript_set &ts1, transcript_set &ts2);
	int assign_RPKM();
	int write_RO_info();
	int write();
	int write_circular_boundaries();
	int write_circular();
	int read_cirifull_file();
	int split(const std::string &s, char delim, std::vector<std::string> &elems);
	int compare(splice_graph &gr, const string &ref, const string &tex = "");
	bool determine_regional_graph(splice_graph &gr);
};

#endif
