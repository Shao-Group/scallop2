
#ifndef __GTF_CIRCULAR_TRANSCRIPT_H__
#define __GTF_CIRCULAR_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace std;

class circular_transcript
{
public:
    circular_transcript();
    circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path);
    circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path, int32_t junc_reads, int32_t non_junc_reads);
    int print(int id);
    ~circular_transcript();
public:
    string circRNA_ID;
    string chrm_id;
	int32_t start;
	int32_t end;
    vector<int> circ_path;

    int32_t junc_reads;
    int32_t non_junc_reads;
};

#endif