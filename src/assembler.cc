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
#include "gtf.h"
#include "genome.h"
#include "assembler.h"
#include "bundle.h"
#include "scallop.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
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
		if(terminate == true) return 0;

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

		if(p.tid > 1) break;

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

	assign_RPKM();

	filter ft(trsts);
	ft.merge_single_exon_transcripts();
	trsts = ft.trs;

	filter ft1(non_full_trsts);
	ft1.merge_single_exon_transcripts();
	non_full_trsts = ft1.trs;

	write();
	
	return 0;
}

int assembler::process(int n)
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

		if(bb.hits.size() < min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);
		bb.chrm = string(buf);

		transcript_set ts1(bb.chrm, 0.9);		// full-length set
		transcript_set ts2(bb.chrm, 0.9);		// non-full-length set

		bundle bd(bb);

		bd.build(1, true);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);

		bd.build(2, true);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);

		/*
		bd.build(1, false);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);

		bd.build(2, false);
		bd.print(index++);
		assemble(bd.gr, bd.hs, ts1, ts2);
		*/

		int sdup = assemble_duplicates / 1 + 1;
		int mdup = assemble_duplicates / 2 + 0;

		vector<transcript> gv1 = ts1.get_transcripts(sdup, mdup);
		vector<transcript> gv2 = ts2.get_transcripts(sdup, mdup);

		for(int k = 0; k < gv1.size(); k++)
		{
			if(gv1[k].exons.size() >= 2) gv1[k].coverage /= (1.0 * assemble_duplicates);
		}
		for(int k = 0; k < gv2.size(); k++) 
		{
			if(gv2[k].exons.size() >= 2) gv2[k].coverage /= (1.0 * assemble_duplicates);
		}

		filter ft1(gv1);
		ft1.filter_length_coverage();
		ft1.remove_nested_transcripts();
		if(ft1.trs.size() >= 1) trsts.insert(trsts.end(), ft1.trs.begin(), ft1.trs.end());

		filter ft2(gv2);
		ft2.filter_length_coverage();
		ft2.remove_nested_transcripts();
		if(ft2.trs.size() >= 1) non_full_trsts.insert(non_full_trsts.end(), ft2.trs.begin(), ft2.trs.end());
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, transcript_set &ts1, transcript_set &ts2)
{
	super_graph sg(gr0, hs0);
	sg.build();

	/*
	vector<transcript> gv;
	vector<transcript> gv1;
	*/

	for(int k = 0; k < sg.subs.size(); k++)
	{
		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		if(determine_regional_graph(gr) == true) continue;
		if(gr.num_edges() <= 0) continue;

		for(int r = 0; r < assemble_duplicates; r++)
		{
			string gid = "gene." + tostring(index) + "." + tostring(k) + "." + tostring(r);
			gr.gid = gid;
			scallop sc(gr, hs, r == 0 ? false : true);
			sc.assemble();

			if(verbose >= 2)
			{
				printf("assembly with r = %d, total %lu transcripts:\n", r, sc.trsts.size());
				for(int i = 0; i < sc.trsts.size(); i++) sc.trsts[i].write(cout);
			}

			for(int i = 0; i < sc.trsts.size(); i++)
			{
				ts1.add(sc.trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(int i = 0; i < sc.non_full_trsts.size(); i++)
			{
				ts2.add(sc.non_full_trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}

			/*
			filter ft(sc.trsts);
			//ft.join_single_exon_transcripts();
			ft.filter_length_coverage();
			if(ft.trs.size() >= 1) gv.insert(gv.end(), ft.trs.begin(), ft.trs.end());

			if(verbose >= 2)
			{
			printf("transcripts after filtering:\n");
			for(int i = 0; i < ft.trs.size(); i++) ft.trs[i].write(cout);

			printf("non full length transcripts:\n");
			for(int i = 0; i < sc.non_full_trsts.size(); i++) sc.non_full_trsts[i].write(cout);
			}

			filter ft1(sc.non_full_trsts);
			//ft.join_single_exon_transcripts();
			ft1.filter_length_coverage();
			if(ft1.trs.size() >= 1) gv1.insert(gv1.end(), ft1.trs.begin(), ft1.trs.end());

			if(verbose >= 2)
			{
			printf("non full transcripts after filtering:\n");
			for(int i = 0; i < ft1.trs.size(); i++) ft1.trs[i].write(cout);
			}
			*/
		}
	}

	return 0;
}

bool assembler::determine_regional_graph(splice_graph &gr)
{
	bool all_regional = true;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.get_vertex_info(i).regional == false) all_regional = false;
		if(all_regional == false) break;
	}
	return all_regional;
}

int assembler::assign_RPKM()
{
	double factor = 1e9 / qlen;
	for(int i = 0; i < trsts.size(); i++)
	{
		trsts[i].assign_RPKM(factor);
	}
	return 0;
}

int assembler::write()
{
	ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;
	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
	}
	fout.close();

        ofstream fout1(output_file1.c_str());
        if(fout1.fail()) return 0;
        for(int i = 0; i < non_full_trsts.size(); i++)
        {
                transcript &t = non_full_trsts[i];
                t.write(fout1);
        }
        fout1.close();

	return 0;
}

int assembler::compare(splice_graph &gr, const string &file, const string &texfile)
{
	if(file == "") return 0;

	genome g(file);
	if(g.genes.size() <= 0) return 0;

	gtf gg(g.genes[0]);

	splice_graph gt;
	gg.build_splice_graph(gt);

	sgraph_compare sgc(gt, gr);
	sgc.compare(texfile);

	return 0;
}
