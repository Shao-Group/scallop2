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
#include "genome.h"
#include "transcript_index.h"
#include "feature_writer.h"

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

	genome gm;			// for gtf genome (-b)
	transcript_index tridx;		// annotation transcript index
	feature_writer features;		// streaming single-sample feature output
	vector<int64_t> transcript_candidate_bundle_counts;
	vector<int64_t> transcript_assigned_bundle_counts;

	int hid;
	int index;

public:
	int assemble();

private:
	int process_gnn(int n);
	int write_transcript_bundle_counts() const;
};

#endif
