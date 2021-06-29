/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "previewer.h"
#include "config.h"
#include "bundle_bridge.h"
#include "bridger.h"

int previewer::open_file()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	return 0;
}

int previewer::close_file()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	return 0;
}

int previewer::preview()
{
	if(library_type == EMPTY)
	{
		open_file();
		solve_strandness();
		close_file();
	}

	if(insertsize_median < 0)
	{
		open_file();
		solve_insertsize();
		close_file();
	}

	return 0;
}

int previewer::solve_strandness()
{
	int total = 0;
	int single = 0;
	int paired = 0;

	int first = 0;
	int second = 0;
	vector<int> sp1;
	vector<int> sp2;

	int hid = 0;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(total >= max_preview_reads) break;
		if(sp1.size() >= max_preview_spliced_reads && sp2.size() >= max_preview_spliced_reads) break;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// qstrandary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		total++;

		hit ht(b1t, hid++);
		ht.set_tags(b1t);

		if((ht.flag & 0x1) >= 1) paired ++;
		if((ht.flag & 0x1) <= 0) single ++;

		if(ht.xs == '.') continue;
		if(ht.xs == '+' && sp1.size() >= max_preview_spliced_reads) continue;
		if(ht.xs == '-' && sp2.size() >= max_preview_spliced_reads) continue;

		// predicted strand
		char xs = '.';

		// for paired read
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '+';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '+';

		// for single read
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) <= 0) xs = '-';
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) >= 1) xs = '+';

		if(xs == '+' && xs == ht.xs) sp1.push_back(1);
		if(xs == '-' && xs == ht.xs) sp2.push_back(1);
		if(xs == '+' && xs != ht.xs) sp1.push_back(2);
		if(xs == '-' && xs != ht.xs) sp2.push_back(2);
	}

	int sp = sp1.size() < sp2.size() ? sp1.size() : sp2.size();

	for(int k = 0; k < sp; k++)
	{
		if(sp1[k] == 1) first++;
		if(sp2[k] == 1) first++;
		if(sp1[k] == 2) second++;
		if(sp2[k] == 2) second++;
	}

	vector<string> vv;
	vv.push_back("empty");
	vv.push_back("unstranded");
	vv.push_back("first");
	vv.push_back("second");

	int s1 = UNSTRANDED;
	if(sp >= min_preview_spliced_reads && first > preview_infer_ratio * 2.0 * sp) s1 = FR_FIRST;
	if(sp >= min_preview_spliced_reads && second > preview_infer_ratio * 2.0 * sp) s1 = FR_SECOND;

	if(verbose >= 1)
	{
		printf("preview strandness: sampled reads = %d, single = %d, paired = %d, first = %d, second = %d, inferred = %s, given = %s\n",
			total, single, paired, first, second, vv[s1 + 1].c_str(), vv[library_type + 1].c_str());
	}

	if(library_type == EMPTY) library_type = s1;

	return 0;
}

int previewer::solve_insertsize()
{
	map<int32_t, int> m;
	bundle_base bb1;
	bundle_base bb2;
	bb1.strand = '+';
	bb2.strand = '-';
	int cnt = 0;
	int hid = 0;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1) continue;										// secondary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, hid++);
		ht.set_tags(b1t);
		ht.set_strand();

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			cnt += process_bundle(bb1, m);
			bb1.clear();
			bb1.strand = '+';
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			cnt += process_bundle(bb2, m);
			bb2.clear();
			bb2.strand = '-';
		}

		//if(cnt >= 1000000) break;
		if(cnt >= 500000) break;

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

	int total = 0;
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		total += it->second;
	}

	if(total < 100) // single-cell data may have low number of paired-end reads
	//if(total < 10000)
	{
		printf("not enough paired-end reads to create the profile (%d collected)\n", total);
		exit(0);
	}

	int n = 0;
	insertsize_ave = 0;
	double sx2 = 0;
	insertsize_low = -1;
	insertsize_high = -1;
	insertsize_median = -1;
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		n += it->second;
		if(n >= 0.5 * total && insertsize_median < 0) insertsize_median = it->first;
		insertsize_ave += it->second * it->first;
		sx2 += it->second * it->first * it->first;
		if(insertsize_low == -1 && n >= 0.005 * total) insertsize_low = it->first;
		if(insertsize_high == -1 && n >= 0.98 * total) insertsize_high = it->first;
		if(n >= 0.998 * total) break;
	}
	
	insertsize_ave = insertsize_ave * 1.0 / n;
	insertsize_std = sqrt((sx2 - n * insertsize_ave * insertsize_ave) * 1.0 / n);

	insertsize_profile.assign(insertsize_high, 1);
	n = insertsize_high;
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		if(it->first >= insertsize_high) continue;
		insertsize_profile[it->first] += it->second;
		n += it->second;
	}

	if(verbose >= 1)
	{
		printf("preview insertsize: sampled reads = %d, isize = %.2lf +/- %.2lf, median = %d, low = %d, high = %d\n", 
				total, insertsize_ave, insertsize_std, insertsize_median, insertsize_low, insertsize_high);
	}

	for(int i = 0; i < insertsize_profile.size(); i++)
	{
		insertsize_profile[i] = insertsize_profile[i] * 1.0 / n;
		//printf("insertsize_profile %d %.8lf\n", i, insertsize_profile[i]);
	}

	// further relax bounds of insertsize
	//insertsize_low = insertsize_low / 1.15;
	//insertsize_high = insertsize_high * 1.15;

	return 0;
}

int previewer::process_bundle(bundle_base &bb, map<int32_t, int>& m)
{
	if(bb.hits.size() < 50) return 0;
	if(bb.tid < 0) return 0;

	int cnt = 0;

	bundle_bridge br(bb);

	br.build_junctions();
	br.extend_junctions();
	br.build_regions();

	br.align_hits_transcripts();
	br.index_references();

	br.build_fragments();
	//br.group_fragments();

	bridger bdg(&br);
	bdg.bridge_overlapped_fragments();

	for(int k = 0; k < br.fragments.size(); k++)
	{
		fragment& fr = br.fragments[k];
		if(fr.paths.size() != 1) continue;

		int32_t len = fr.paths[0].length;

		/*
		fr.print(k);
		printf("fragment: len = %d, v = (", len);
		printv(decode_vlist(fr.paths[0].v));
		printf(")\n");
		*/

		cnt++;

		if(m.find(len) != m.end()) m[len]++;
		else m.insert(pair<int, int>(len, 1));

		if(cnt >= 1000) return cnt;
	}
	return cnt;
}
