#include "feature_writer.h"

#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

feature_writer::feature_writer(const string &prefix, const string &s)
	: sample(s), nodes((prefix + ".node.csv").c_str()),
	  edges((prefix + ".edge.csv").c_str()),
	  phasing((prefix + ".phasing.csv").c_str()),
	  paths((prefix + ".path.label.csv").c_str())
{
	nodes << "chr,graph_id,node_id,start_pos,end_pos,extPathSupport,weight,length,maxcov,stddev,indel_sum_cov,indel_ratio,left_indel,right_indel,sample\n";
	edges << "chr,graph_id,source,target,start_pos,end_pos,extPathSupport,weight,length,sample\n";
	phasing << "chr,graph_id,path_id,node_sequence,count,sample\n";
	paths << "chr,graph_id,path_id,node_sequence,splice_source,splice_target,abundance,label,sample\n";
}

bool feature_writer::good() const
{
	return nodes.good() && edges.good() && phasing.good() && paths.good();
}

string feature_writer::csv(const string &s)
{
	if(s.find_first_of(",\"\r\n") == string::npos) return s;
	return quoted(s);
}

string feature_writer::quoted(const string &s)
{
	string z = "\"";
	for(size_t i = 0; i < s.size(); i++)
	{
		if(s[i] == '"') z += "\"\"";
		else z += s[i];
	}
	return z + "\"";
}

string feature_writer::sequence(const vector<int> &v, int offset)
{
	ostringstream ss;
	for(size_t i = 0; i < v.size(); i++)
	{
		if(i > 0) ss << ",";
		ss << v[i] + offset;
	}
	return ss.str();
}

int feature_writer::write_bundle(bundle &bd, int bundle_index)
{
	const splice_graph &gr = bd.gr;
	ostringstream id;
	if(gr.chrm.compare(0, 3, "chr") != 0) id << "chr";
	id << gr.chrm << ".instance." << bundle_index << "." << sample;
	string graph_id = id.str();

	set<int> supported_nodes;
	set<pair<int, int> > supported_edges;
	vector<vector<int> > transcript_paths(bd.assigned_transcripts.size());
	for(size_t i = 0; i < bd.assigned_transcripts.size(); i++)
	{
		if(!bd.get_transcript_path(bd.assigned_transcripts[i], transcript_paths[i])) return 1;
		const vector<int> &v = transcript_paths[i];
		for(size_t j = 0; j < v.size(); j++) supported_nodes.insert(v[j]);
		for(size_t j = 1; j < v.size(); j++) supported_edges.insert(make_pair(v[j - 1], v[j]));
	}

	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		const partial_exon &pe = bd.pexons[i - 1];
		nodes << csv(gr.chrm) << "," << csv(graph_id) << "," << i - 1 << ","
			<< vi.lpos + 1 << "," << vi.rpos << "," << fixed << setprecision(1)
			<< (supported_nodes.count(i) ? 1.0 : 0.0) << "," << setprecision(4)
			<< gr.get_vertex_weight(i) << "," << vi.rpos - vi.lpos << ","
			<< pe.max << "," << pe.dev << "," << pe.indel_sum_cov << "," << pe.indel_ratio << ","
			<< pe.left_indel << "," << pe.right_indel << "," << csv(sample) << "\n";
	}

	PEEI all_edges = gr.edges();
	for(edge_iterator it = all_edges.first; it != all_edges.second; ++it)
	{
		edge_descriptor e = *it;
		int s = e->source(), t = e->target();
		if(s == 0 || t == gr.num_vertices() - 1) continue;
		vertex_info sv = gr.get_vertex_info(s), tv = gr.get_vertex_info(t);
		edges << csv(gr.chrm) << "," << csv(graph_id) << "," << s - 1 << "," << t - 1 << ","
			<< sv.rpos << "," << tv.lpos + 1 << ","
			<< fixed << setprecision(1) << (supported_edges.count(make_pair(s, t)) ? 1.0 : 0.0) << ","
			<< setprecision(4) << gr.get_edge_weight(e) << ","
			<< tv.lpos - sv.rpos + 1 << "," << csv(sample) << "\n";
	}

	MVII selected_phases = bd.hs.nodes;
	if(gr.num_vertices() > 100 && selected_phases.size() > static_cast<size_t>(gr.num_vertices()))
	{
		size_t max_length = 0;
		int max_count = 0;
		for(MVII::const_iterator it = selected_phases.begin(); it != selected_phases.end(); ++it)
		{
			max_length = std::max(max_length, it->first.size());
			max_count = std::max(max_count, it->second);
		}
		typedef pair<double, PVII> scored_phase;
		vector<scored_phase> scored;
		for(MVII::const_iterator it = selected_phases.begin(); it != selected_phases.end(); ++it)
		{
			if(it->first.size() <= 2) continue;
			double score = 1.0 * it->first.size() / max_length + 1.0 * it->second / max_count;
			scored.push_back(make_pair(score, *it));
		}
		sort(scored.begin(), scored.end(), [](const scored_phase &a, const scored_phase &b)
		{
			if(a.first != b.first) return a.first > b.first;
			return a.second.first < b.second.first;
		});
		selected_phases.clear();
		for(size_t i = 0; i < scored.size() && i < static_cast<size_t>(gr.num_vertices()); i++)
			selected_phases.insert(scored[i].second);
	}

	int phase_id = 0;
	for(MVII::const_iterator it = selected_phases.begin(); it != selected_phases.end(); ++it)
	{
		if(it->first.size() <= 2) continue;
		phasing << csv(gr.chrm) << "," << csv(graph_id) << ","
			<< csv(graph_id + ".phasing." + to_string(phase_id++)) << ","
			<< quoted(sequence(it->first, -1)) << "," << it->second << "," << csv(sample) << "\n";
	}

	for(size_t i = 0; i < bd.assigned_transcripts.size(); i++)
	{
		const transcript &tr = bd.assigned_transcripts[i];
		const vector<int> &v = transcript_paths[i];
		vector<int> splice_sources, splice_targets;
		for(size_t j = 1; j < v.size(); j++)
		{
			if(gr.get_vertex_info(v[j - 1]).rpos != gr.get_vertex_info(v[j]).lpos)
			{
				splice_sources.push_back(v[j - 1]);
				splice_targets.push_back(v[j]);
			}
		}
		paths << csv(gr.chrm) << "," << csv(graph_id) << ","
			<< csv(graph_id + "." + tr.transcript_id) << "," << quoted(sequence(v, -1)) << ","
			<< quoted(sequence(splice_sources, -1)) << "," << quoted(sequence(splice_targets, -1)) << ","
			// No truth annotation is available from -i/-b alone; -1 means unlabeled.
			<< tr.coverage << ",-1," << csv(sample) << "\n";
	}

	return good() ? 0 : 1;
}
