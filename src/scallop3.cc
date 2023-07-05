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
    //top_paths.resize(gr.num_vertices(), vector< vector<int> >(topnum, vector<int>()));
    //top_btn.resize(gr.num_vertices(), vector<int>(topnum, -1));
}

scallop3::~scallop3()
{
}

int scallop3::assemble()
{
	
    gr.print_weights();
	//hs.print();

    printf("process splice graph %s, vertices = %lu, edges = %lu, phasing paths = %lu\n", gr.gid.c_str(), gr.num_vertices(), gr.num_edges(), hs.edges.size());

    set<int> critical_edge;
    //calculate_critical_edges(critical_edge);
    hs.add_edge_not_phased(gr.num_edges(), critical_edge);

    cout << "\n#Phasing path: " << hs.edges.size() << endl;
    for(int i = 0; i < hs.edges.size(); i++)
    {
        printf("%d: weight = %d, list = ",i, hs.ecnts[i]);
        print_phasing_path(hs.edges[i]);
		printf("\n");

    }
    printf("\n");

    int round = 1, r = 0;
    int t = gr.num_vertices()-1;
    map<vector<int>, int> top_rounds;
    while(r < round)
    {
        round--;

        top_paths.clear();
        top_btn.clear();
        top_paths.resize(gr.num_vertices(), vector< vector<int> >(topnum, vector<int>()));
        top_btn.resize(gr.num_vertices(), vector<int>(topnum, -1));

        calculate_max_bottleneck_path();
        /*for(int i = 0; i < 1; i++)
        {
            if(top_rounds.find(top_paths[t][i]) == top_rounds.end() && top_btn[t][i]>0)
            {
                top_rounds[top_paths[t][i]] = top_btn[t][i];
                hs.update_edge_count(top_paths[t][i], top_btn[t][i]);

            }
        }

        printf("\nAfter round %d----Phasing path:\n", r);
        for(int i = 0; i < hs.edges.size(); i++)
        {
            if(hs.edges[i].size() > 1) continue;
            printf("%d: weight = %d, list = ",i, hs.ecnts[i]);
            print_phasing_path(hs.edges[i]);
            printf("\n");
        } 
        printf("\n");
        if(top_btn[t][0] <= 0) break;
        if(topnum>1 && top_btn[t][1]<0) break;*/

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

    /*for(auto it1 = top_rounds.begin(); it1 != top_rounds.end(); it1++)
    {
        path p;
        const vector<int> &H = it1->first;
        for(auto it2 = H.begin(); it2 != H.end(); it2++)
        {
            int s = i2e[*it2]->source();
            int t = i2e[*it2]->target();
            if(it2 == H.begin())
                p.v.push_back(s);
            p.v.push_back(t);
        }
        p.abd = it1->second;
        paths.push_back(p);
    }*/

    trsts.clear();
	gr.output_transcripts(trsts, paths);

	if(verbose >= 0) 
	{
		for(int i = 0; i < paths.size(); i++) paths[i].print(i);
		printf("finish assemble bundle %s\n\n", gr.gid.c_str());
	}

    return 0;
}


int scallop3::calculate_max_bottleneck_path()
{
    //DP with maximum bottleneck
    int t = gr.num_vertices()-1;
    for(int v = 1; v <= t; v++)
    {
        set< vector<int> > new_paths;

        /*PEEI pei = gr.in_edges(v);
		for(auto it = pei.first; it != pei.second; it++)
		{
			edge_descriptor e = (*it);
            printf("-----Current edge: (%d, %d)\n", e->source(), e->target());

            for(int i = 0; i < topnum; i++)
            {
                int frontv = e->source();
                vector<int> new_path = top_paths[frontv][i];
                if(frontv != 0 && new_path.size() <= 0) continue;
                if(frontv == 0 && i>0) continue;

                new_path.push_back(e2i[e]);
                new_paths.insert(new_path);
                //printf("New path: ");
                //print_phasing_path(new_path);
            }
        }*/
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
            printf("New path: ");
            print_phasing_path(*np);

            set<int> edge_btn_ignored;
            edge_within_start_exon(*np, edge_btn_ignored);
            if(v == t) edge_within_end_exon(*np, edge_btn_ignored);
           
            PI btnP = hs.get_compatible_bottleneck(*np, edge_btn_ignored);
            //PI btnP = hs.get_compatible_bottleneck(*np);
            int btn = btnP.second;
            int btn_edge = btnP.first;
            if(btn_edge == -1)
                printf("Currently no btn.\n");
            else
                printf("Edge(%d, %d) is btn = %d\n", i2e[btn_edge]->source(), i2e[btn_edge]->target(), btn);
            for(int i = 0; i < topnum; i++)
            {
                int b = top_btn[v][i];
                if(btn > b)
			    {
				    for(int j = topnum-1; j>i; j--)
                    {
                        if(top_paths[v][j-1].size() == 0)continue;
                        top_btn[v][j] = top_btn[v][j-1];
                        top_paths[v][j] = top_paths[v][j-1];
                    }
                    top_btn[v][i] =  btn;
				    top_paths[v][i] = *np;
                    break;
			    }
            }
        }
        printf("Top %d for vertex %d:\n", topnum, v);
        for(int i = 0; i < topnum; i++)
        {
            print_phasing_path(top_paths[v][i]);
            printf("Btn = %d\n", top_btn[v][i]);
        }
        printf("\n");

    }

	return 0;
}

int scallop3::calculate_critical_edges(set<int> &critical_edge)
{
    for(int v = 0; v < gr.num_vertices(); v++)
    {
        int critical = true;
        if(gr.out_degree(v) == 1)
        {
            edge_descriptor e = *(gr.out_edges(v).first);
            int u = e->target();
            if(gr.in_degree(u) == 1)
            {
                PEEI pei = gr.edges();
                for(auto it = pei.first; it != pei.second; it++)
                {
                    edge_descriptor e2 = *it;
                    int s = e2->source(), t = e2->target();
                    if(s<=v && t>=u && e!=e2)
                    {
                        critical = false;
                        break;
                    }
                }
                if(critical)
                {
                    critical_edge.insert(e2i[e]);
                    printf("Critical: (%d, %d)\n", v, u);
                }
            }
        }
    }
    if(critical_edge.size() == gr.num_edges()) critical_edge.clear();
    return 0;
}

int scallop3::print_phasing_path(const vector<int> &phasing_edge)
{
    if(phasing_edge.size() <= 0) return 0;
    set<int> vpath;
    vector<int>::const_iterator it;
    for(it = phasing_edge.begin(); it != phasing_edge.end(); it++)
    {
        assert(*it < i2e.size());
        edge_descriptor e = i2e[*it];
        int s = e->source();
        int t = e->target();
        //printf("%d(%d, %d), ", *it, s, t);
        vpath.insert(s);
        vpath.insert(t);
    }
    //printf("\n");
    set<int>::const_iterator vit;

	assert(phasing_edge.size() + 1 == vpath.size());
    //if(phasing_edge.size()+1 != vpath.size()) printf("print phasing path wrong.\n");
    //printf("p = ");
    for(vit=vpath.begin(); vit != vpath.end(); vit++)
    {
        cout << *vit << ' ';
    }
    printf("\n");
    return 0;
}

//compute continuous edges within the start exon of the path
int scallop3::edge_within_start_exon(const vector<int>& p, set<int>& output)
{
    //assert(i2e[p[0]]->source() == 0);
    set<int> start_edges;
    bool longer = false;
    printf("Start exon: ");
    for(int i = 0; i < p.size(); i++)
    {
        edge_descriptor e = i2e[p[i]];
        int s = e->source();
        int t = e->target();
        int32_t p1 = gr.get_vertex_info(s).rpos;
        int32_t p2 = gr.get_vertex_info(t).lpos;
        int32_t eLen = p2-p1+1;
        if(gr.in_degree(s)>1) longer=true;
        if(s == 0 || eLen <= 1 || t == gr.num_vertices()-1)
        {
            start_edges.insert(p[i]);
            printf("(%d -> %d) ", s, t);
        }
        else
        {
            if(!longer)
                output.insert(start_edges.begin(), start_edges.end());
            break;
        }
    }
    printf("\n");
    return 0;
}

//compute continuous edges within the ending exon of the path
int scallop3::edge_within_end_exon(const vector<int>& p, set<int>& output)
{
    set<int> end_edges;
    bool longer = false;
    printf("Ending exon: ");
    for(int i = p.size()-1; i >= 0; i--)
    {
        edge_descriptor e = i2e[p[i]];
        int s = e->source();
        int t = e->target();
        int32_t p1 = gr.get_vertex_info(s).rpos;
        int32_t p2 = gr.get_vertex_info(t).lpos;
        int32_t eLen = p2-p1+1;
        if(gr.out_degree(t)>1) longer = true;
        if(t == gr.num_vertices()-1 || eLen <= 1 || s == 0)
        {
            end_edges.insert(p[i]);
            printf("(%d <- %d) ", t, s);
        }
        else
        {
            if(!longer)
                output.insert(end_edges.begin(), end_edges.end());
            break;
        }
    }
    printf("\n");
    return 0;
}
