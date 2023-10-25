/*
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __GTF_CIRCULAR_TRANSCRIPT_H__
#define __GTF_CIRCULAR_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "../../src/region.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

class circular_transcript
{
public:
    circular_transcript();
    circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path);
    circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path, int32_t junc_reads, int32_t non_junc_reads);
	int write(ostream &fout, double cov2 = -1, int count = -1) const;
    int print(int id);
    ~circular_transcript();
public:

	string circRNA_id;
    string seqname; //chromosome id
	string source;
	string feature;
	string gene_id; //is it needed for circRNA? this infor comes from splice graphs of each gene, difficult to extract for circRNA
	string transcript_id;
	string gene_type;
	string transcript_type;
	int32_t start;
	int32_t end;
	double score;
	char strand; //is it needed for circRNA? this infor comes from splice graphs of each gene, difficult to extract for circRNA
	int frame;
	int coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

	//storing chimeric and bridging features
	bool fake_supple;
	int32_t supple_len; //len of supple read
	double path_score; //bridging paths core
	int path_type; //1/2 for ref path, 3/4 for read path 

    vector<int> circ_path;
	vector<region> circ_path_regions;
	vector<region> merged_regions;
    
    int32_t junc_reads;
    int32_t non_junc_reads;
	size_t bundle_size;
};

#endif