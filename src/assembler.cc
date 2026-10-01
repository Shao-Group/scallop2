/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "genome.h"
#include "assembler.h"
#include "bundle.h"
#include "transcript_match.h"

assembler::assembler()
	: gm(gtf_file), tridx(gm)
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;

	if(verbose >= 1) printf("loaded %d annotation transcripts\n", tridx.size());
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int assembler::assemble()
{
	while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if(p.tid < 0) continue;
		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;							// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, hid++);
		ht.set_tags(b1t);
		ht.set_strand();
		//ht.print();

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;

		//if(p.tid > 1) break;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
		}

		// process
		if(process_gnn(batch_bundle_size) != 0) return 1;

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.' && ht.spos.size() <= 0) bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	if(process_gnn(0) != 0) return 1;
	
	return 0;
}

int assembler::process_gnn(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];

		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);
		bb.chrm = string(buf);

		transcript_matcher matcher(bb);
		vector<int> candidates = tridx.query(matcher.chromosome(), matcher.get_strand(),
				bb.lpos, bb.rpos);
		vector<int> assigned_transcripts;
		vector<transcript_match> matches;
		vector<int> unassigned_transcripts;
		vector<transcript_match> unassigned_matches;

		for(int k = 0; k < candidates.size(); k++)
		{
			int id = candidates[k];
			transcript_match m = matcher.match(tridx.get(id), min_bundle_transcript_exon_overlap,
					min_bundle_transcript_splicing_position_overlap);
			if(m.assigned == false)
			{
				if(verbose >= 2)
				{
					unassigned_transcripts.push_back(id);
					unassigned_matches.push_back(m);
				}
				continue;
			}
			assigned_transcripts.push_back(id);
			matches.push_back(m);
		}

		vector<transcript> assigned_models;
		for(int k = 0; k < assigned_transcripts.size(); k++)
		{
			assigned_models.push_back(tridx.get(assigned_transcripts[k]).trst);
		}
		bundle bd(bb, assigned_models);
		if(bd.build(1, true) != 0)
		{
			printf("error: failed to preserve assigned transcripts in bundle splice graph.\n");
			return 1;
		}

		int bundle_index = index++;
		if(verbose >= 1)
		{
			printf("bundle %d: candidate-transcripts = %lu, assigned-transcripts = %lu\n",
					bundle_index, candidates.size(), assigned_transcripts.size());
		}
		if(verbose >= 2)
		{
			for(int k = 0; k < assigned_transcripts.size(); k++)
			{
				const indexed_transcript &t = tridx.get(assigned_transcripts[k]);
				const transcript_match &m = matches[k];
				printf("bundle %d: transcript = %s, gene = %s, exonic-overlap = %d/%d, splicing-position-overlap = %d/%d\n",
						bundle_index, t.trst.transcript_id.c_str(), t.trst.gene_id.c_str(),
						m.shared_exonic_length, m.transcript_exonic_length,
						m.shared_splicing_positions, m.transcript_splicing_positions);
			}
			for(int k = 0; k < unassigned_transcripts.size(); k++)
			{
				int id = unassigned_transcripts[k];
				const indexed_transcript &t = tridx.get(id);
				const transcript_match &m = unassigned_matches[k];
				double exon_fraction = m.transcript_exonic_length > 0 ?
					1.0 * m.shared_exonic_length / m.transcript_exonic_length : 0;
				double splicing_position_fraction = m.transcript_splicing_positions > 0 ?
					1.0 * m.shared_splicing_positions / m.transcript_splicing_positions : 1;

				printf("bundle %d: unassigned-transcript = %s, transcript-index = %d, gene = %s, region = %s:%d-%d, strand = %c, reason = %s, exonic-overlap = %d/%d (%.6f; required %.6f), splicing-position-overlap = %d/%d (%.6f; required %.6f), exons = ",
						bundle_index, t.trst.transcript_id.c_str(), id, t.trst.gene_id.c_str(),
						t.trst.seqname.c_str(), t.bounds.first, t.bounds.second, t.trst.strand,
						m.failure_reason(), m.shared_exonic_length, m.transcript_exonic_length,
						exon_fraction, min_bundle_transcript_exon_overlap,
						m.shared_splicing_positions, m.transcript_splicing_positions,
						splicing_position_fraction, min_bundle_transcript_splicing_position_overlap);
				for(int j = 0; j < t.trst.exons.size(); j++)
				{
					if(j >= 1) printf(",");
					printf("[%d,%d)", t.trst.exons[j].first, t.trst.exons[j].second);
				}
				printf("\n");
			}
		}
	}
	pool.clear();
	return 0;
}
