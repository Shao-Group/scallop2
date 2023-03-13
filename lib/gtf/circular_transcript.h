
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

	//string chrm_id;
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
	double coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

    vector<int> circ_path;
	vector<region> circ_path_regions;

    //string circRNA_ID;
	//int32_t start;
	//int32_t end;
    
    int32_t junc_reads;
    int32_t non_junc_reads;
};

#endif