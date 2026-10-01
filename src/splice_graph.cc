/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "splice_graph.h"
#include "util.h"
#include "config.h"
#include "interval_map.h"
#include <sstream>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <cstdlib>

using namespace std;

splice_graph::splice_graph()
{}

splice_graph::splice_graph(const splice_graph &gr)
{
	chrm = gr.chrm;
	gid = gr.gid;
	strand = gr.strand;

	MEE x2y;
	MEE y2x;
	copy(gr, x2y, y2x);
}

int splice_graph::copy(const splice_graph &gr, MEE &x2y, MEE &y2x)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
		set_vertex_weight(i, gr.get_vertex_weight(i));
		set_vertex_info(i, gr.get_vertex_info(i));
	}

	PEEI p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		edge_descriptor e = add_edge((*it)->source(), (*it)->target());
		set_edge_weight(e, gr.get_edge_weight(*it));
		set_edge_info(e, gr.get_edge_info(*it));

		assert(e != NULL);
		assert(ewrt.find(e) != ewrt.end());
		assert(einf.find(e) != einf.end());
		assert(x2y.find(*it) == x2y.end());
		assert(y2x.find(e) == y2x.end());

		x2y.insert(PEE(*it, e));
		y2x.insert(PEE(e, *it));
	}

	return 0;
}

int splice_graph::clear()
{
	directed_graph::clear();
	vwrt.clear();
	vinf.clear();
	ewrt.clear();
	einf.clear();
	return 0;
}

splice_graph::~splice_graph()
{}

double splice_graph::get_vertex_weight(int v) const
{
	assert(v >= 0 && v < vwrt.size());
	return vwrt[v];
}

vertex_info splice_graph::get_vertex_info(int v) const
{
	assert(v >= 0 && v < vinf.size());
	return vinf[v];
}

double splice_graph::get_edge_weight(edge_base *e) const
{
	MED::const_iterator it = ewrt.find(e);
	assert(it != ewrt.end());
	return it->second;
}

edge_info splice_graph::get_edge_info(edge_base *e) const
{
	MEIF::const_iterator it = einf.find(e);
	assert(it != einf.end());
	return it->second;
}

int splice_graph::set_vertex_weight(int v, double w) 
{
	assert(v >= 0 && v < vv.size());
	if(vwrt.size() != vv.size()) vwrt.resize(vv.size());
	vwrt[v] = w;
	return 0;
}

int splice_graph::set_vertex_info(int v, const vertex_info &vi) 
{
	assert(v >= 0 && v < vv.size());
	if(vinf.size() != vv.size()) vinf.resize(vv.size());
	vinf[v] = vi;
	return 0;
}

int splice_graph::set_edge_weight(edge_base* e, double w) 
{
	if(ewrt.find(e) != ewrt.end()) ewrt[e] = w;
	else ewrt.insert(PED(e, w));
	return 0;
}

int splice_graph::set_edge_info(edge_base* e, const edge_info &ei) 
{
	if(einf.find(e) != einf.end()) einf[e] = ei;
	else einf.insert(PEIF(e, ei));
	return 0;
}

edge_descriptor splice_graph::max_out_edge(int v)
{
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	PEEI pei;
	double ww = 0;
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

edge_descriptor splice_graph::max_in_edge(int v)
{
	PEEI pei;
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	double ww = 0;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

int splice_graph::build(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}

	char line[10240];
	// get the number of vertices
	fin.getline(line, 10240, '\n');	
	int n = atoi(line);

	for(int i = 0; i < n; i++)
	{
		char name[10240];
		double weight;
		vertex_info vi;
		fin.getline(line, 10240, '\n');	
		stringstream sstr(line);
		sstr>>name>>weight>>vi.length;

		add_vertex();
		set_vertex_weight(i, weight);
		set_vertex_info(i, vi);
	}

	while(fin.getline(line, 10240, '\n'))
	{
		int x, y;
		double weight;
		edge_info ei;
		stringstream sstr(line);
		sstr>>x>>y>>weight>>ei.length;

		assert(x != y);
		assert(x >= 0 && x < num_vertices());
		assert(y >= 0 && y < num_vertices());

		edge_descriptor p = add_edge(x, y);
		set_edge_weight(p, weight);
		set_edge_info(p, ei);
	}

	fin.close();
	return 0;
}

int splice_graph::write(const string &file) const
{
	ofstream fin(file.c_str());
	if(fin.fail()) 
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}
	
	fin<<fixed;
	fin.precision(2);
	int n = num_vertices();
	
	fin<<n<<endl;
	for(int i = 0; i < n; i++)
	{
		string name = "scallop2";
		double weight = get_vertex_weight(i);
		vertex_info vi = get_vertex_info(i);
		fin<<name.c_str()<<" "<<weight<<" "<<vi.length<<endl;
	}

	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source(); 
		int t = (*it1)->target();
		double weight = get_edge_weight(*it1);
		edge_info ei = get_edge_info(*it1);
		fin<<s<<" "<<t<<" "<<weight<<" "<<ei.length<<endl;
	}
	fin.close();
	return 0;
}

int splice_graph::locate_vertex(int32_t p)
{
	return locate_vertex(p, 0, num_vertices());
}

int splice_graph::locate_vertex(int32_t p, int a, int b)
{
	if(a >= b) return -1;
	int m = (a + b) / 2;
	assert(m >= 0 && m < num_vertices());
	const vertex_info &v = get_vertex_info(m);
	if(p >= v.lpos && p < v.rpos) return m;
	if(p < v.lpos) return locate_vertex(p, a, m);
	return locate_vertex(p, m + 1, b);
}

int splice_graph::draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp)
{
	return directed_graph::draw(file, mis, mes, len, tp);
}

int splice_graph::draw(const string &file, const MIS &mis, const MES &mes, double len)
{
	return directed_graph::draw(file, mis, mes, len);
}

int splice_graph::draw(const string &file)
{
	MIS mis;
	char buf[10240];

	for(int i = 0; i < num_vertices(); i++)
	{
		double w = get_vertex_weight(i);
		vertex_info vi = get_vertex_info(i);
		int ll = vi.lpos % 100000;
		int rr = vi.rpos % 100000;
		sprintf(buf, "%.1lf:%d-%d", w, ll, rr);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		sprintf(buf, "%.1lf", w);
		mes.insert(PES(*it1, buf));
	}
	draw(file, mis, mes, 4.5);
	return 0;
}

int splice_graph::print_nontrivial_vertices()
{
	int k = 0;
	for(int i = 1; i < num_vertices() - 1; i++)
	{
		if(in_degree(i) <= 1) continue;
		if(out_degree(i) <= 1) continue;
		printf("nontrivial vertex %d length = %d\n", k++, get_vertex_info(i).length);
	}
	return 0;
}

int splice_graph::print()
{
	for(int i = 0; i < num_vertices(); i++)
	{
		//if(degree(i) <= 1) continue;
		vertex_info vi = get_vertex_info(i);
		edge_iterator it1, it2;
		PEEI pei;
		printf("vertex %d, range = [%d, %d), length = %d\n", i, vi.lpos, vi.rpos, vi.rpos - vi.lpos);
		printf(" in-vertices ="); 
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			printf(" %d, ", (*it1)->source());
		}
		printf("\n out-vertices = ");
		for(pei = out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			printf(" %d, ", (*it1)->target());
		}
		printf("\n");
	}
	return 0;
}

int splice_graph::print_weights()
{
	for(int i = 0; i < num_vertices(); i++)
	{
		//if(degree(i) <= 1) continue;
		vertex_info vi = get_vertex_info(i);
		edge_iterator it1, it2;
		PEEI pei;
		printf("vertex %d, range = [%d, %d), length = %d\n", i, vi.lpos, vi.rpos, vi.rpos - vi.lpos);
	}

	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p1 = get_vertex_info(s).rpos;
		int32_t p2 = get_vertex_info(t).lpos;
		double w1 = get_edge_weight(e);
		double w2 = get_edge_info(e).weight;
		printf("edge (%d, %d) pos = %d-%d length = %d weight = (%.2lf, %.2lf)\n", s, t, p1, p2, p2 - p1 + 1, w1, w2);
	}
	return 0;
}
