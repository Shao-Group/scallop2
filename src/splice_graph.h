/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

#include "directed_graph.h"
#include "vertex_info.h"
#include "edge_info.h"

#include <map>
#include <cassert>

#define SMIN 0.00001

using namespace std;

typedef map<edge_descriptor, edge_info> MEIF;
typedef pair<edge_descriptor, edge_info> PEIF;

class splice_graph : public directed_graph
{
public:
	splice_graph();
	splice_graph(const splice_graph &gr);
	virtual ~splice_graph();

public:
	string chrm;
	string gid;
	char strand;

	vector<double> vwrt;
	vector<vertex_info> vinf;
	MED ewrt;
	MEIF einf;

public:
	// get and set properties
	double get_vertex_weight(int v) const;
	double get_edge_weight(edge_base *e) const;
	vertex_info get_vertex_info(int v) const;
	edge_info get_edge_info(edge_base *e) const;

	int set_vertex_weight(int v, double w);
	int set_vertex_info(int v, const vertex_info &vi);
	int set_edge_weight(edge_base *e, double w);
	int set_edge_info(edge_base *e, const edge_info &ei);

	edge_descriptor max_out_edge(int v);
	edge_descriptor max_in_edge(int v);

	// modify the splice_graph
	int clear();
	int copy(const splice_graph &gr, MEE &x2y, MEE &y2x);

	// read and write splice graph
	int build(const string &file);
	int write(const string &file) const;
	int locate_vertex(int32_t p);
	int locate_vertex(int32_t p, int a, int b);

	// draw and print
	int draw(const string &file);
	int draw(const string &file, const MIS &mis, const MES &mes, double len);
	int draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp);
	int print_nontrivial_vertices();
	int print_weights();
	int print();
};

#endif
