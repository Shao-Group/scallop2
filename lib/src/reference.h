/*
Part of Coral
(c) 2019 by  Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __REFERENCE_H__
#define __REFERENCE_H__

#include <map>

#include "genome.h"
#include "interval_map.h"

using namespace std;

class reference : public genome
{
public:
	reference(const string &file);

public:
	map<string, interval_set_map> isms0;		// indexed by chrm name
	map<string, interval_set_map> isms1;		// indexed by chrm name
	map<string, interval_set_map> isms2;		// indexed by chrm name

public:
	int build_interval_set_map();
	int add_interval_set(map<string, interval_set_map> &isms, int k);
	int print();
	vector<transcript> get_overlapped_transcripts(string chrm, char c, int32_t x, int32_t y) const;
	vector<transcript> get_overlapped_transcripts(const map<string, interval_set_map> &isms, string chrm, int32_t x, int32_t y) const;
};

#endif
