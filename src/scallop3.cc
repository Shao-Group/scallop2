/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "scallop3.h"
#include "config.h"

#include <cstdio>
#include <iostream>
#include <climits>
#include <cfloat>
#include <algorithm>

scallop3::scallop3()
{}

scallop3::scallop3(const splice_graph &g, const hyper_set &h)
	: gr(g), hs(h)
{
	gr.get_edge_indices(i2e, e2i);
	hs.build(gr, e2i);

    topnum = 3;
    top_paths.resize(gr.num_vertices(), vector< vector<int> >(topnum, vector<int>()));
    top_btn.resize(gr.num_vertices(), vector<int>(topnum, -1));
}

scallop3::~scallop3()
{
}

int scallop3::assemble()
{
	
    /*gr.print_weights();
	hs.print();

    cout << "\n#Phasing path: " << hs.edges.size() << endl;
    for(int i = 0; i < hs.edges.size(); i++)
    {
        printf("%d: weight = %d, list = ",i, hs.ecnts[i]);
        print_phasing_path(hs.edges[i]);
		printf("\n");

    }
    printf("\n");*/

	printf("process splice graph %s, vertices = %lu, edges = %lu, phasing paths = %lu\n", gr.gid.c_str(), gr.num_vertices(), gr.num_edges(), hs.edges.size());

    hs.add_edge_not_phased(gr.num_edges());

    //DP with maximum bottleneck
    int t = gr.num_vertices()-1;
    for(int v = 1; v <= t; v++)
    {
        set< vector<int> > new_paths;
        for(auto p = hs.edges.begin(); p != hs.edges.end(); p++)
        {
            int backv = i2e[p->back()]->target();
            if(backv != v) continue;
            //printf("-----Current phasing: ");
            //print_phasing_path(*p);
            for(int i = 0; i < topnum; i++)
            {
                int frontv = i2e[p->front()]->source();
                vector<int> new_path = top_paths[frontv][i];
                if(frontv != 0 && new_path.size() <= 0) continue;
                if(frontv == 0 && i>0) continue;

                new_path.insert(new_path.end(), p->begin(), p->end());
                new_paths.insert(new_path);
                //printf("New path: ");
                //print_phasing_path(new_path);

            }
        }

        for(auto np = new_paths.begin(); np != new_paths.end(); np++)
        {
            int btn = hs.get_compatible_bottleneck(*np);
            for(int i = 0; i < topnum; i++)
            {
                int b = top_btn[v][i];
                if(btn > b)
			    {
				    for(int j = topnum-1; j>i; j--)
                    {
                        top_btn[v][j] = top_btn[v][j-1];
                        top_paths[v][j] = top_paths[v][j-1];
                    }
                    top_btn[v][i] =  btn;
				    top_paths[v][i] = *np;
                    break;
			    }
            }
        }

    }

    paths.clear();
    for(int i = 0; i < 1; i++)
    {
        if(top_btn[t][i] <= 0) continue;
        path p;
        vector<int> &H = top_paths[t][i];
        for(auto it = H.begin(); it != H.end(); it++)
        {
            int s = i2e[*it]->source();
            int t = i2e[*it]->target();
            if(it == H.begin())
                p.v.push_back(s);
            p.v.push_back(t);
        }
        p.abd = top_btn[t][i];
        paths.push_back(p);
    }

    trsts.clear();
	gr.output_transcripts(trsts, paths);

	if(verbose >= 0) 
	{
		for(int i = 0; i < paths.size(); i++) paths[i].print(i);
		printf("finish assemble bundle %s\n\n", gr.gid.c_str());
	}

	return 0;
}

int scallop3::print_phasing_path(const vector<int> &phasing_edge)
{
    set<int> vpath;
    vector<int>::const_iterator it;
    for(it = phasing_edge.begin(); it != phasing_edge.end(); it++)
    {
        assert(*it < i2e.size());
        edge_descriptor e = i2e[*it];
        int s = e->source();
        int t = e->target();
        printf("%d(%d, %d), ", *it, s, t);
        vpath.insert(s);
        vpath.insert(t);
    }
    printf("\n");
    set<int>::const_iterator vit;

	//assert(phasing_edge.size() + 1 == vpath.size());
    if(phasing_edge.size()+1 != vpath.size()) printf("print phasing path wrong.\n");
    printf("p = ");
    for(vit=vpath.begin(); vit != vpath.end(); vit++)
    {
        cout << *vit << ' ';
    }
    printf("\n");
    return 0;
}
