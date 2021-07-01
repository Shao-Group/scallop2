/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#ifndef __TRANSCRIPT_SET_H__
#define __TRANSCRIPT_SET_H__

#include <map>
#include <vector>
#include <string>
#include "transcript.h"

using namespace std;

class trans_item
{
public:
	trans_item();
	trans_item(const transcript &t, int count, int sid);

public:
	transcript trst;
	int count;
	map<int, double> samples;

public:
	int merge(const trans_item &ti, int mode);
};

int merge_sorted_trans_items(vector<trans_item> &vx, const vector<trans_item> &vy, int mode1, int mode2, double single_exon_overlap);

class transcript_set
{
public:
	transcript_set(const string &c, double single_exon_overlap);
	transcript_set(const transcript &t, int count, int sid, double single_exon_overlap);

public:
	string chrm;
	map<size_t, vector<trans_item>> mt;
	double single_exon_overlap;

public:
	int add(const transcript &t, int count, int sid, int mode1, int mode2);
	int add(const transcript_set &ts, int mode1, int mode2);
	int increase_count(int count);
	int filter(int min_count);
	int print() const;
	pair<bool, trans_item> query(const transcript &t) const;
	vector<transcript> get_transcripts(int min_count) const;
	vector<transcript> get_transcripts(int min1, int min2) const;
};

#endif
