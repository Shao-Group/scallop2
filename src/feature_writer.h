/* Single-sample, Aletsch-compatible splice-graph feature serialization. */
#ifndef __FEATURE_WRITER_H__
#define __FEATURE_WRITER_H__

#include <fstream>
#include <string>
#include "bundle.h"

using namespace std;

class feature_writer
{
public:
	feature_writer(const string &prefix, const string &sample);
	bool good() const;
	int write_bundle(bundle &bd, int bundle_index);

private:
	string sample;
	ofstream nodes;
	ofstream edges;
	ofstream phasing;
	ofstream paths;

	static string csv(const string &s);
	static string quoted(const string &s);
	static string sequence(const vector<int> &v, int offset);
};

#endif
