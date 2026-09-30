/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>

#include "transcript_index.h"

indexed_transcript::indexed_transcript(const transcript &t)
	: trst(t)
{
	bounds = trst.get_bounds();
	junctions = trst.get_intron_chain();
	exonic_length = trst.length();
}

transcript_index::transcript_index()
{}

transcript_index::transcript_index(const genome &gm)
{
	build(gm);
}

int transcript_index::build(const genome &gm)
{
	transcripts.clear();
	index.clear();

	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			if(t.exons.size() == 0) continue;

			indexed_transcript z(t);
			if(z.bounds.first < 0 || z.bounds.first >= z.bounds.second) continue;
			if(z.exonic_length <= 0) continue;

			int id = transcripts.size();
			transcripts.push_back(z);

			set<int> s;
			s.insert(id);
			pair<string, char> key(t.seqname, t.strand);
			index[key] += make_pair(ROI(z.bounds.first, z.bounds.second), s);
		}
	}
	return 0;
}

vector<int> transcript_index::query(const string &chrm, char strand, int32_t lpos, int32_t rpos) const
{
	vector<int> v;
	if(lpos >= rpos) return v;

	pair<string, char> key(chrm, strand);
	map<pair<string, char>, transcript_interval_map>::const_iterator mi = index.find(key);
	if(mi == index.end()) return v;

	const transcript_interval_map &imap = mi->second;
	TIMI it = imap.find(lpos);
	if(it == imap.end()) it = imap.lower_bound(ROI(lpos, rpos));

	set<int> s;
	for(; it != imap.end(); it++)
	{
		if(lower(it->first) >= rpos) break;
		if(upper(it->first) <= lpos) continue;
		s.insert(it->second.begin(), it->second.end());
	}

	v.insert(v.end(), s.begin(), s.end());
	return v;
}

const indexed_transcript& transcript_index::get(int k) const
{
	assert(k >= 0 && k < transcripts.size());
	return transcripts[k];
}

int transcript_index::size() const
{
	return transcripts.size();
}
