/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __SCALLOP3B_H__
#define __SCALLOP3B_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop3
{
public:
	scallop3();
	scallop3(const splice_graph &gr, const hyper_set &hs);
	virtual ~scallop3();

public:
	int assemble();

public:
	splice_graph gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	hyper_set hs;						// hyper edges

    //DP entries
    vector< vector< vector<int> > > top_paths; //top_paths[v]: top-3 paths until v
    vector< vector<int> > top_btn; //top_btn[v]: bottleneck of top-3 paths until v
    int topnum;

    //output
	vector<path> paths;					// predicted paths
	vector<transcript> trsts;			// predicted transcripts

private:
    int calculate_critical_edges(set<int> &critical_edge);
    int calculate_max_bottleneck_path();
    int print_phasing_path(const vector<int> &phasing_edge);
    int edge_within_start_exon(const vector<int>& p, set<int>& start_edges);
    int edge_within_end_exon(const vector<int>& p, set<int>& end_edges);


};

#endif
