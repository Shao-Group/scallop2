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
	build_reference();
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	qlen = 0;
	qcnt = 0;
}


int transformer::build_reference()
{
	if(ref_file == "") return 0;

	gm.read(ref_file);
	//for(auto x: gm.g2i) printf("gene %s has %lu transcripts\n", x.first.c_str(), gm.genes[x.second].transcripts.size());

	trsts = gm.collect_transcripts();
	for(int k = 0; k < trsts.size(); k++)
	{
		transcript &t = trsts[k];
		assert(tmap.find(t.transcript_id) == tmap.end());
		tmap.insert(make_pair(t.transcript_id, k));
	}
	return 0;
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

	pool.clear();
	return 0;
}

string transformer::get_transcript_id(string s)
{
	string delim = ":";

	auto start = 0U;
	auto end = s.find(delim);
	for(int k = 0; k < 2; k++)
	{
		start = end + delim.length();
		end = s.find(delim, start);
	}
	//printf("name = %s\n", s.substr(start, end - start).c_str());
	return s.substr(start, end - start);
}

int transformer::process(bundle_base &bb)
{
	map<string, int> ttm;
	for(int k = 0; k < bb.hits.size(); k++)
	{
		hit &h = bb.hits[k];
		string tid = get_transcript_id(h.qname);

		if(ttm.find(tid) == ttm.end()) ttm.insert(make_pair(tid, 1));
		else ttm[tid]++;
	}

	map<string, int> ggm;
	for(auto x: ttm)
	{
		if(tmap.find(x.first) == tmap.end()) printf("transcript name = %s count = %d\n", x.first.c_str(), x.second);
		assert(tmap.find(x.first) != tmap.end());
		string gid = trsts[tmap[x.first]].gene_id;
		//assert(gm.g2i.find(gid) != gm.g2i.end());
		//int gi = gm.g2i[gid];
		if(ggm.find(gid) == ggm.end()) ggm.insert(make_pair(gid, 1));
		else ggm[gid]++;
	}

	int num = 0;
	for(auto &x: ggm)
	{
		assert(gm.g2i.find(x.first) != gm.g2i.end());
		num += gm.genes[gm.g2i[x.first]].transcripts.size();
	}

	// filtering
	if(ggm.size() != 1 || ttm.size() != num) return 0;

	printf("## instance bundle with %lu reads in %lu transcripts and in %lu genes (total %d transcripts)\n", 
			bb.hits.size(), ttm.size(), ggm.size(), num);

	map<vector<int64_t>, int> rmap;
	for(int k = 0; k < bb.hits.size(); k++)
	{
		hit &h = bb.hits[k];
		vector<int64_t> v;
		h.get_aligned_intervals(v);
		if(v.size() <= 0) continue;
		if(rmap.find(v) == rmap.end()) rmap.insert(make_pair(v, 1));
		else rmap[v]++;
	}

	for(auto &x: rmap)
	{
		const vector<int64_t> &v = x.first;
		printf("%d, ", x.second);
		for(int i = 0; i < v.size(); i++)
		{
			printf("%d-%d, ", low32(v[i]), high32(v[i]));
		}
		printf("\n");
	}

	printf("# ground-truth transcripts\n"); 

	for(auto x: ttm)
	{
		assert(tmap.find(x.first) != tmap.end());
		transcript &t = trsts[tmap[x.first]];
		for(int k = 0; k < t.exons.size(); k++)
		{
			printf("%d-%d, ", t.exons[k].first, t.exons[k].second);
		}
		printf("\n");
	}

	printf("\n");

	return 0;
}

int transformer::write()
{
	ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;
	fout.close();
	return 0;
}
