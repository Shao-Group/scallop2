/*
Part of Coral
(c) 2019 by  Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "reference.h"

reference::reference(const string &file)
	: genome(file)
{
	build_interval_set_map();
}

int reference::build_interval_set_map()
{
	for(int k = 0; k < genes.size(); k++)
	{
		char c = genes[k].get_strand();
		if(c == '.') add_interval_set(isms0, k);
		if(c == '+') add_interval_set(isms1, k);
		if(c == '-') add_interval_set(isms2, k);
	}
	return 0;
}

int reference::add_interval_set(map<string, interval_set_map> &isms, int k)
{
	if(genes[k].transcripts.size() <= 0) return 0;

	PI32 p = genes[k].get_bounds();

	if(p.first >= p.second) return 0;
	interval32 x(p.first, p.second);

	set<int> s;
	s.insert(k);

	string chrm = genes[k].get_seqname();

	printf("shao: add gene %d-%d to chrm %s, gene-id = %s\n", p.first, p.second, chrm.c_str(), genes[k].get_gene_id().c_str());

	if(isms.find(chrm) == isms.end())
	{
		interval_set_map ism;
		ism += make_pair(x, s);
		isms.insert(pair<string, interval_set_map>(chrm, ism));
	}
	else
	{
		isms[chrm] += make_pair(x, s);
	}
	return 0;
}

vector<transcript> reference::get_overlapped_transcripts(string chrm, char c, int32_t x, int32_t y) const
{
	vector<transcript> v;
	if(c == '.') return get_overlapped_transcripts(isms0, chrm, x, y);
	if(c == '+') return get_overlapped_transcripts(isms1, chrm, x, y);
	if(c == '-') return get_overlapped_transcripts(isms2, chrm, x, y);
	return v;
}

vector<transcript> reference::get_overlapped_transcripts(const map<string, interval_set_map> &isms, string chrm, int32_t x, int32_t y) const
{
	vector<transcript> v;
	map<string, interval_set_map>::const_iterator it = isms.find(chrm);
	if(it == isms.end()) return v;
	set<int> s = get_overlapped_set(it->second, x, y);

	printf("shao: query %d-%d of chrm %s and found %lu genes\n", x, y, chrm.c_str(), s.size());

	for(set<int>::iterator x = s.begin(); x != s.end(); x++)
	{
		v.insert(v.end(), genes[*x].transcripts.begin(), genes[*x].transcripts.end());
	}
	return v;
}

int reference::print()
{
	printf("refernece annotation contains %lu genes\n", genes.size());
	for(int i = 0; i < genes.size(); i++)
	{
		printf("gene %d with %lu transcripts\n", i, genes[i].transcripts.size());
	}
	for(map<string, interval_set_map>::iterator it = isms0.begin(); it != isms0.end(); it++)
	{
		printf("chromosomes %s with strand = .\n", it->first.c_str());
		print_interval_set_map(it->second);
	}
	for(map<string, interval_set_map>::iterator it = isms1.begin(); it != isms1.end(); it++)
	{
		printf("chromosomes %s with strand = +\n", it->first.c_str());
		print_interval_set_map(it->second);
	}
	for(map<string, interval_set_map>::iterator it = isms2.begin(); it != isms2.end(); it++)
	{
		printf("chromosomes %s with strand = -\n", it->first.c_str());
		print_interval_set_map(it->second);
	}
	return 0;
}
