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

#include "circular_transcript.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
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

public:
	int assemble();

private:
	int process(int n);
	int print_circular_trsts();
	int assemble(const splice_graph &gr, const hyper_set &hs, transcript_set &ts1, transcript_set &ts2);
	int assign_RPKM();
	int write();
	int write_circular();
	int compare(splice_graph &gr, const string &ref, const string &tex = "");
	bool determine_regional_graph(splice_graph &gr);
};

#endif
