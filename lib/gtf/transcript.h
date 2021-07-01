/*
Part of aletsch
(c) 2020 by  Mingfu Shao, The Pennsylvania State University.
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GTF_TRANSCRIPT_H__
#define __GTF_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "item.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

class transcript
{
public:
	transcript(const item &ie);
	transcript();
	~transcript();

public:
	bool operator< (const transcript &t) const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;
	string transcript_id;
	string gene_type;
	string transcript_type;
	int32_t start;
	int32_t end;
	double score;
	char strand;
	int frame;
	double coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

	vector<PI32> exons;

public:
	int add_exon(int s, int t);
	int add_exon(const item &e);
	int assign_RPKM(double factor);
	int sort();
	int clear();
	int shrink();
	int assign(const item &e);
	int length() const;
	PI32 get_bounds() const;
	PI32 get_first_intron() const;
	vector<PI32> get_intron_chain() const;
	size_t get_intron_chain_hashing() const;
	bool intron_chain_match(const transcript &t) const;
	int intron_chain_compare(const transcript &t) const;
	bool equal1(const transcript &t, double single_exon_overlap) const;
	int compare1(const transcript &t, double single_exon_overlap) const;
	int extend_bounds(const transcript &t);
	string label() const;
	int write(ostream &fout, double cov2 = -1, int count = -1) const;
};

#endif
