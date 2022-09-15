/*
(c) 2022 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __TRANSFORMER_H__
#define __TRANSFORMER_H__

#include <fstream>
#include <string>
#include "bundle_base.h"
#include "transcript.h"
#include "genome.h"

using namespace std;

class transformer
{
public:
	transformer();
	~transformer();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	bundle_base bb1;		// +
	bundle_base bb2;		// -
	genome gm;				// genome
	vector<bundle_base> pool;
	vector<transcript> trsts;
	map<string, int> tmap;

	int hid;
	int index;
	int qcnt;
	double qlen;

public:
	int transform();

private:
	int process(int n);
	int process(bundle_base &bb);
	int build_reference();
	int write();
};

#endif
