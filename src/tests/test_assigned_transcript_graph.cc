#include <cassert>
#include <cmath>
#include <cstdio>

#include "bundle.h"
#include "config.h"
#include "util.h"

transcript make_assigned_transcript()
{
	transcript t;
	t.seqname = "chr1";
	t.strand = '+';
	t.gene_id = "g1";
	t.transcript_id = "assigned";
	t.exons.push_back(PI32(80, 150));
	t.exons.push_back(PI32(200, 250));
	t.exons.push_back(PI32(300, 340));
	return t;
}

transcript make_nested_transcript()
{
	transcript t;
	t.seqname = "chr1";
	t.strand = '+';
	t.gene_id = "g1";
	t.transcript_id = "nested";
	t.exons.push_back(PI32(100, 150));
	t.exons.push_back(PI32(200, 225));
	return t;
}

void test_add_only_missing_intervals()
{
	bundle_base bb;
	bb.chrm = "chr1";
	bb.strand = '+';
	bb.lpos = 90;
	bb.rpos = 160;
	bundle bd(bb);
	bd.fmap += make_pair(ROI(100, 120), 5);
	bd.fmap += make_pair(ROI(130, 140), 3);

	transcript t1;
	t1.exons.push_back(PI32(90, 135));
	transcript t2;
	t2.exons.push_back(PI32(110, 150));
	bd.assigned_transcripts.push_back(t1);
	bd.assigned_transcripts.push_back(t2);
	bd.add_assigned_transcript_intervals();

	assert(compute_overlap(bd.fmap, 105) == 5);
	assert(compute_overlap(bd.fmap, 135) == 3);
	assert(compute_overlap(bd.fmap, 95) == 1);
	assert(compute_overlap(bd.fmap, 125) == 1);
	assert(compute_overlap(bd.fmap, 145) == 1);
	assert(compute_overlap(bd.newly_added_intervals, 105) == 0);
	assert(compute_overlap(bd.newly_added_intervals, 95) == 1);
	assert(compute_overlap(bd.newly_added_intervals, 125) == 1);
	assert(compute_overlap(bd.newly_added_intervals, 145) == 1);

	int missing_length = 0;
	for(SIMI it = bd.newly_added_intervals.begin(); it != bd.newly_added_intervals.end(); it++)
		missing_length += upper(it->first) - lower(it->first);
	assert(missing_length == 30);

	bd.regions.clear();
	bd.build_regions();
	bd.build_partial_exons();
	int partial_exon_missing_length = 0;
	bool found_partially_added_exon = false;
	for(int i = 0; i < bd.pexons.size(); i++)
	{
		int length = bd.pexons[i].rpos - bd.pexons[i].lpos;
		partial_exon_missing_length += bd.pexons[i].newly_added_length;
		if(bd.pexons[i].newly_added_length > 0 && bd.pexons[i].newly_added_length < length)
			found_partially_added_exon = true;
	}
	assert(partial_exon_missing_length == 30);
	assert(found_partially_added_exon);
	bd.build_splice_graph(1);
	for(int i = 0; i < bd.pexons.size(); i++)
		assert(bd.gr.get_vertex_info(i + 1).newly_added_length == bd.pexons[i].newly_added_length);
}

int main()
{
	library_type = FR_SECOND;
	verbose = 0;
	test_add_only_missing_intervals();

	bundle_base bb;
	bb.chrm = "chr1";
	bb.strand = '+';
	bb.lpos = 100;
	bb.rpos = 320;
	bb.mmap += make_pair(ROI(100, 150), 5);
	bb.mmap += make_pair(ROI(200, 225), 5);
	bb.imap += make_pair(ROI(110, 115), 2);

	transcript t = make_assigned_transcript();
	transcript nested = make_nested_transcript();
	vector<transcript> assigned;
	assigned.push_back(t);
	assigned.push_back(nested);
	bundle bd(bb, assigned);
	assert(bd.assigned_transcripts.size() == 2);
	assert(bb.lpos == 80 && bb.rpos == 340);

	set<int64_t> observed;
	for(int i = 0; i < bd.junctions.size(); i++)
		observed.insert(pack(bd.junctions[i].lpos, bd.junctions[i].rpos));
	assert(observed.find(pack(150, 200)) != observed.end());
	assert(observed.find(pack(250, 300)) != observed.end());
	int shared_junction_count = 0;
	for(int i = 0; i < bd.junctions.size(); i++)
		if(bd.junctions[i].lpos == 150 && bd.junctions[i].rpos == 200)
			shared_junction_count++;
	assert(shared_junction_count == 1);

	assert(bd.regions.size() >= 1);
	assert(bd.regions.front().lpos == 80);
	assert(bd.regions.back().rpos == 340);
	bool found_nested_start = false;
	bool found_nested_end = false;
	for(int i = 0; i < bd.regions.size(); i++)
	{
		if(bd.regions[i].lpos == 100 || bd.regions[i].rpos == 100)
			found_nested_start = true;
		if(bd.regions[i].lpos == 225 || bd.regions[i].rpos == 225)
			found_nested_end = true;
	}
	assert(found_nested_start && found_nested_end);
	bool found_indel_features = false;
	for(int i = 0; i < bd.pexons.size(); i++)
	{
		const partial_exon &pe = bd.pexons[i];
		if(pe.lpos != 100 || pe.rpos != 150) continue;
		assert(pe.indel_sum_cov == 10);
		assert(fabs(pe.indel_ratio - 0.20) < 1e-9);
		assert(pe.left_indel == 10);
		assert(pe.right_indel == 35);
		found_indel_features = true;
	}
	assert(found_indel_features);
	assert(bd.build(1, true) == 0);
	assert(bd.has_transcript_path(t));
	assert(bd.has_transcript_path(nested));
	int newly_added_length = 0;
	for(int i = 0; i < bd.pexons.size(); i++)
	{
		newly_added_length += bd.pexons[i].newly_added_length;
		assert(bd.gr.get_vertex_info(i + 1).newly_added_length == bd.pexons[i].newly_added_length);
	}
	assert(newly_added_length == 160);

	printf("assigned transcript graph tests passed\n");
	return 0;
}
