/*
(c) 2022 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "gtf.h"
#include "genome.h"
#include "transformer.h"

transformer::transformer()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	qlen = 0;
	qcnt = 0;

	if(ref_file != "")
	{
		gm.read(ref_file);
		for(auto x: gm.g2i)
		{
			printf("gene %s has %lu transcripts\n", x.first.c_str(), gm.genes[x.second].transcripts.size());
		}
	}
}

transformer::~transformer()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int transformer::transform()
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

		qlen += ht.qlen;
		qcnt += 1;

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
		process(batch_bundle_size);

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	write();
	return 0;
}

int transformer::process(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];

		/*
		// calculate the number of hits with splices
		int splices = 0;
		for(int k = 0; k < bb.hits.size(); k++)
		{
			if(bb.hits[k].spos.size() >= 1) splices++;
		}
		if(bb.hits.size() < min_num_hits_in_bundle && splices < min_num_splices_in_bundle) continue;
		//printf("bundle %d has %lu reads, %d reads have splices\n", i, bb.hits.size(), splices);
		*/

		int cnt1 = 0;
		int cnt2 = 0;
		for(int k = 0; k < bb.hits.size(); k++)
		{
			//counts += (1 + bb.hits[k].spos.size());
			if(bb.hits[k].spos.size() >= 1) cnt1 ++;
			else cnt2++;
		}

		if(cnt1 + cnt2 < min_num_hits_in_bundle) continue;
		//if(cnt1 < 5 && cnt1 * 2 + cnt2 < min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);
		bb.chrm = string(buf);

		process(bb);
	}
	return 0;
}

int transformer::process(bundle_base &bb)
{
	return 0;
}

int transformer::write()
{
	ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;
	fout.close();
	return 0;
}
