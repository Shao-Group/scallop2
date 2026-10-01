#include <cassert>
#include <cstdio>

#include "transcript_index.h"
#include "transcript_match.h"
#include "util.h"

transcript build_transcript(const string &chrm, char strand, const string &gid,
		const string &tid, const vector<PI32> &exons)
{
	transcript t;
	t.seqname = chrm;
	t.strand = strand;
	t.gene_id = gid;
	t.transcript_id = tid;
	t.exons = exons;
	return t;
}

int test_index()
{
	genome gm;
	gene g;

	vector<PI32> e0;
	e0.push_back(PI32(10, 20));
	e0.push_back(PI32(30, 40));
	g.add_transcript(build_transcript("chr1", '+', "g1", "t0", e0));

	vector<PI32> e1;
	e1.push_back(PI32(35, 50));
	g.add_transcript(build_transcript("chr1", '+', "g1", "t1", e1));

	vector<PI32> e2;
	e2.push_back(PI32(10, 50));
	g.add_transcript(build_transcript("chr1", '-', "g1", "t2", e2));
	gm.add_gene(g);

	gene h;
	vector<PI32> e3;
	e3.push_back(PI32(10, 50));
	h.add_transcript(build_transcript("chr2", '+', "g2", "t3", e3));
	gm.add_gene(h);

	transcript_index ti(gm);
	assert(ti.size() == 4);
	assert(ti.query("chr1", '+', 0, 10).size() == 0);
	assert(ti.query("chr1", '+', 0, 11) == vector<int>(1, 0));

	vector<int> v = ti.query("chr1", '+', 36, 37);
	assert(v.size() == 2 && v[0] == 0 && v[1] == 1);
	v = ti.query("chr1", '+', 40, 41);
	assert(v.size() == 1 && v[0] == 1);
	v = ti.query("chr1", '-', 36, 37);
	assert(v.size() == 1 && v[0] == 2);
	v = ti.query("chr2", '+', 36, 37);
	assert(v.size() == 1 && v[0] == 3);
	return 0;
}

int test_matcher()
{
	vector<PI32> te;
	te.push_back(PI32(10, 20));
	te.push_back(PI32(30, 40));
	te.push_back(PI32(50, 60));
	te.push_back(PI32(70, 80));
	indexed_transcript t(build_transcript("chr1", '+', "g1", "t0", te));

	vector<PI32> ge;
	ge.push_back(PI32(10, 20));
	ge.push_back(PI32(30, 40));
	set<int32_t> gp;
	gp.insert(20);
	gp.insert(30);
	gp.insert(40);
	gp.insert(50);
	transcript_matcher m1("chr1", '+', ge, gp);
	transcript_match x = m1.match(t, 0.5, 0.5);
	assert(x.assigned == true);
	assert(x.shared_exonic_length == 20 && x.transcript_exonic_length == 40);
	assert(x.shared_splicing_positions == 4 && x.transcript_splicing_positions == 6);

	gp.erase(40);
	gp.erase(50);
	transcript_matcher m2("chr1", '+', ge, gp);
	transcript_match x2 = m2.match(t, 0.5, 0.5);
	assert(x2.assigned == false);
	assert(x2.failure == TRANSCRIPT_INSUFFICIENT_SPLICING_POSITION_OVERLAP);

	ge[1] = PI32(30, 39);
	gp.insert(40);
	gp.insert(50);
	transcript_matcher m3("chr1", '+', ge, gp);
	transcript_match x3 = m3.match(t, 0.5, 0.5);
	assert(x3.assigned == false);
	assert(x3.failure == TRANSCRIPT_INSUFFICIENT_EXON_OVERLAP);

	transcript_matcher m4("chr1", '-', ge, gp);
	transcript_match x4 = m4.match(t, 0.5, 0.5);
	assert(x4.assigned == false);
	assert(x4.failure == TRANSCRIPT_STRAND_MISMATCH);

	vector<PI32> se;
	se.push_back(PI32(100, 110));
	indexed_transcript single(build_transcript("chr1", '+', "g1", "single", se));
	vector<PI32> half;
	half.push_back(PI32(100, 105));
	set<int32_t> empty;
	transcript_matcher m5("chr1", '+', half, empty);
	assert(m5.match(single, 0.5, 0.5).assigned == true);

	half[0] = PI32(110, 120);
	transcript_matcher m6("chr1", '+', half, empty);
	transcript_match x6 = m6.match(single, 0.5, 0.5);
	assert(x6.assigned == false);
	assert(x6.failure == TRANSCRIPT_INSUFFICIENT_EXON_OVERLAP);
	return 0;
}

int test_bundle_matcher()
{
	bundle_base bb;
	bb.chrm = "chr1";
	bb.strand = '+';
	bb.lpos = 10;
	bb.rpos = 40;
	bb.mmap += make_pair(ROI(10, 15), 1);
	bb.mmap += make_pair(ROI(15, 20), 2);
	bb.mmap += make_pair(ROI(30, 40), 1);

	vector<PI32> te;
	te.push_back(PI32(10, 20));
	te.push_back(PI32(30, 40));
	indexed_transcript t(build_transcript("chr1", '+', "g1", "t0", te));
	transcript_matcher matcher(bb);
	transcript_match m = matcher.match(t, 0.5, 0.0);
	assert(matcher.chromosome() == "chr1" && matcher.get_strand() == '+');
	assert(m.assigned == true);
	assert(m.shared_exonic_length == 20);
	assert(m.shared_splicing_positions == 0 && m.transcript_splicing_positions == 2);
	return 0;
}

int main()
{
	test_index();
	test_matcher();
	test_bundle_matcher();
	printf("transcript assignment tests passed\n");
	return 0;
}
