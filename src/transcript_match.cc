/*
Part of Scallop2
(c) 2026 by  Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>

#include "transcript_match.h"
#include "util.h"

transcript_match::transcript_match()
{
	assigned = false;
	failure = TRANSCRIPT_MATCHED;
	shared_exonic_length = 0;
	transcript_exonic_length = 0;
	shared_splicing_positions = 0;
	transcript_splicing_positions = 0;
}

const char* transcript_match::failure_reason() const
{
	switch(failure)
	{
	case TRANSCRIPT_MATCHED: return "matched";
	case TRANSCRIPT_CHROMOSOME_MISMATCH: return "chromosome-mismatch";
	case TRANSCRIPT_STRAND_MISMATCH: return "strand-mismatch";
	case TRANSCRIPT_INVALID_EXONIC_LENGTH: return "invalid-exonic-length";
	case TRANSCRIPT_INSUFFICIENT_EXON_OVERLAP: return "insufficient-exon-overlap";
	case TRANSCRIPT_INSUFFICIENT_SPLICING_POSITION_OVERLAP: return "insufficient-splicing-position-overlap";
	}
	return "unknown";
}

transcript_matcher::transcript_matcher(const bundle_base &bb)
{
	build(bb);
}

transcript_matcher::transcript_matcher(const string &c, char s, const vector<PI32> &v,
		const set<int32_t> &p)
	: chrm(c), strand(s), splicing_positions(p)
{
	build_exons(v);
}

int transcript_matcher::build(const bundle_base &bb)
{
	chrm = bb.chrm;
	strand = bb.strand;
	if(strand == '.')
	{
		int np = 0, nq = 0;
		for(int i = 0; i < bb.hits.size(); i++)
		{
			if(bb.hits[i].xs == '+') np++;
			if(bb.hits[i].xs == '-') nq++;
		}
		if(np > nq) strand = '+';
		else if(np < nq) strand = '-';
	}

	build_exons(bb.mmap);
	splicing_positions.clear();
	for(int i = 0; i < bb.hits.size(); i++)
	{
		const vector<int64_t> &v = bb.hits[i].spos;
		for(int k = 0; k < v.size(); k++)
		{
			splicing_positions.insert(high32(v[k]));
			splicing_positions.insert(low32(v[k]));
		}
	}
	return 0;
}

int transcript_matcher::build_exons(const split_interval_map &mmap)
{
	exons.clear();
	for(SIMI it = mmap.begin(); it != mmap.end(); it++)
	{
		exons += make_pair(it->first, 1);
	}
	return 0;
}

int transcript_matcher::build_exons(const vector<PI32> &v)
{
	exons.clear();
	for(int i = 0; i < v.size(); i++)
	{
		if(v[i].first >= v[i].second) continue;
		exons += make_pair(ROI(v[i].first, v[i].second), 1);
	}
	return 0;
}

int transcript_matcher::compute_shared_exonic_length(const transcript &t) const
{
	int z = 0;
	for(int i = 0; i < t.exons.size(); i++)
	{
		const PI32 &p = t.exons[i];
		JIMI it = exons.find(p.first);
		if(it == exons.end()) it = exons.lower_bound(ROI(p.first, p.second));

		for(; it != exons.end(); it++)
		{
			int32_t l = lower(it->first);
			int32_t r = upper(it->first);
			if(l >= p.second) break;
			if(r <= p.first) continue;
			int32_t ll = l > p.first ? l : p.first;
			int32_t rr = r < p.second ? r : p.second;
			if(ll < rr) z += rr - ll;
		}
	}
	return z;
}

const string& transcript_matcher::chromosome() const
{
	return chrm;
}

char transcript_matcher::get_strand() const
{
	return strand;
}

transcript_match transcript_matcher::match(const indexed_transcript &t,
		double min_exon_overlap, double min_splicing_position_overlap) const
{
	transcript_match m;
	m.transcript_exonic_length = t.exonic_length;
	set<int32_t> transcript_positions;
	for(int i = 0; i < t.junctions.size(); i++)
	{
		transcript_positions.insert(t.junctions[i].first);
		transcript_positions.insert(t.junctions[i].second);
	}
	m.transcript_splicing_positions = transcript_positions.size();

	if(t.trst.seqname != chrm)
	{
		m.failure = TRANSCRIPT_CHROMOSOME_MISMATCH;
		return m;
	}
	if(t.trst.strand != strand)
	{
		m.failure = TRANSCRIPT_STRAND_MISMATCH;
		return m;
	}
	if(m.transcript_exonic_length <= 0)
	{
		m.failure = TRANSCRIPT_INVALID_EXONIC_LENGTH;
		return m;
	}

	m.shared_exonic_length = compute_shared_exonic_length(t.trst);
	for(set<int32_t>::const_iterator it = transcript_positions.begin();
			it != transcript_positions.end(); it++)
	{
		if(splicing_positions.find(*it) != splicing_positions.end())
			m.shared_splicing_positions++;
	}

	if(m.shared_exonic_length < min_exon_overlap * m.transcript_exonic_length)
	{
		m.failure = TRANSCRIPT_INSUFFICIENT_EXON_OVERLAP;
		return m;
	}
	if(m.transcript_splicing_positions >= 1 &&
		m.shared_splicing_positions <
		min_splicing_position_overlap * m.transcript_splicing_positions)
	{
		m.failure = TRANSCRIPT_INSUFFICIENT_SPLICING_POSITION_OVERLAP;
		return m;
	}
	m.assigned = true;
	return m;
}
