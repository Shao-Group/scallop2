/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <fstream>

#include "hit.h"
#include "config.h"
#include "util.h"
#include "aligner.h"

/*
hit::hit(int32_t p)
{
	bam1_core_t::pos = p;
	strand = '.';
	xs = '.';
	ts = '.';
	hi = -1;
	nh = -1;
	nm = 0;
	qlen = 0;
	cigar = NULL;
}
*/

hit& hit::operator=(const hit &h)
{
	bam1_core_t::operator=(h);
	hid = h.hid;
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nh = h.nh;
	nm = h.nm;

	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
	itvc1 = h.itvc1;
	itvc2 = h.itvc2;
	left_anchor_padding = h.left_anchor_padding;
	right_anchor_padding = h.right_anchor_padding;		

	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;

	return *this;
}

hit::hit(const hit &h) 
	:bam1_core_t(h)
{
	hid = h.hid;
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nh = h.nh;
	nm = h.nm;

	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
	itvc1 = h.itvc1;
	itvc2 = h.itvc2;
	left_anchor_padding = h.left_anchor_padding;
	right_anchor_padding = h.right_anchor_padding;		

	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;
}

hit::hit(bam1_t *b, int id) 
	:bam1_core_t(b->core), hid(id)
{
	// fetch query name
	qname = get_qname(b);
	qhash = string_hash(qname);
	paired = false;
	bridged = false;
	next = NULL;

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));
	qlen = (int32_t)bam_cigar2qlen(n_cigar, bam_get_cigar(b));

	// get cigar
	assert(n_cigar <= max_num_cigar);
	assert(n_cigar >= 1);
	uint32_t * cigar = bam_get_cigar(b);

	// build splice positions
	spos.clear();
	int32_t p = pos;
	int32_t q = 0;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			q += bam_cigar_oplen(cigar[k]);

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;
		// if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		// if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;

		// consider ALL splice positions
		//if(bam_cigar_oplen(cigar[k-1]) < min_flank_length) continue;
		//if(bam_cigar_oplen(cigar[k+1]) < min_flank_length) continue;

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		spos.push_back(pack(s, p));
	}

	// open for scallop+coral
	itvm.clear();
	itvi.clear();
	itvd.clear();
	p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
		{
			p += bam_cigar_oplen(cigar[k]);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvm.push_back(pack(s, p));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			itvi.push_back(pack(p - 1, p + 1));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvd.push_back(pack(s, p));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
		{
			assert(bam_cigar_op(cigar[k]) != BAM_CHARD_CLIP);
			assert (k == 0 || k == n_cigar - 1);
			if (k == 0)	
			{
				itvc1 = {p - bam_cigar_oplen(cigar[k]), p};
			}
			else 
			{
				itvc2 = {p, p + bam_cigar_oplen(cigar[k])};
			}
		}
	}
	
	set_anchors(b);
	//printf("call regular constructor\n");
}

int hit::set_anchors(bam1_t *b)
{
	cout << "hitid : " << qname << endl ;
	vector<string> sc_info;
	sc_info.push_back(qname);
 	left_anchor_padding  = -1;
	right_anchor_padding = -1;

	int anchor_start_nm = anchor_nm_threshold >= 0? anchor_nm_threshold : 10;
	int anchor_end_nm   = anchor_nm_threshold >= 0? anchor_nm_threshold : 10;

	if (berth_mode == 0) return 0;

	// whether SEQ in bam is reverse complemente
	bool seqrev;
	if((flag & 0x10) >= 1) seqrev = true;
	if((flag & 0x10) <= 0) seqrev = false;
	
	// whether second in pair	// TODO: what if seq attachment to 1st strand but sequence 2nd strand
	bool second_in_pair;
	if((flag & 0x1) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) second_in_pair = true;

	// left clipped sequence
	if ((anchor_start != "" && !seqrev) || (anchor_end != "" && seqrev))
	{
		int ql1 = itvc1.first;
		int ql2 = itvc1.second;
		int seqlen = ql2 - ql1;
		cout << ql1 << ", " << ql2 << ", " << qlen <<  "," << pos << endl;
		assert (seqlen >= 0);	

		// get left clip seq
		string leftclipseq(seqlen, 'N');	
		uint8_t *seq_ptr = bam_get_seq (b);
		for (int i = 0; i < seqlen; i++)
		{
			leftclipseq[i] = seq_nt16_str[bam_seqi(seq_ptr, i)];
		}
		cout << "leftclipseq: " << leftclipseq << endl; //CLEAN:
		sc_info.push_back(leftclipseq);

		// get left anchor position
		int anchorpos = -1;
		if (seqrev)	
		{
			pair<int, int> pos_pair = subseq_pos(anchor_end, revcomp(leftclipseq), anchor_end_nm);
			left_anchor_padding = pos_pair.first >= 0 ? pos_pair.first : -1;
		}
		else
		{
			pair<int, int> pos_pair = subseq_pos(anchor_start, leftclipseq, anchor_start_nm);
			left_anchor_padding = pos_pair.second >= 0 ? seqlen - pos_pair.second : -1;
		}

		// left_anchor_padding = anchorpos >= 0? seqlen - anchorpos: -1;
		cout << "left anchor padding: " << left_anchor_padding << endl;
		sc_info.push_back(to_string(left_anchor_padding));
	}

	// right clipped sequence
	if ((anchor_end != "" && !seqrev) || (anchor_start != "" && seqrev)) 
	{
		int ql1 = itvc2.first;
		int ql2 = itvc2.second;
		int seqlen = ql2 - ql1;
		cout << ql1 << ", " << ql2 << ", " << qlen <<  "," << pos << endl;
		assert (seqlen >= 0);	

		// get right clip seq
		string rightclipseq(seqlen, 'N');	
		uint8_t *seq_ptr = bam_get_seq (b);
		for (int i = 0; i < seqlen; i++)
		{	
			// assert( i + ql1 - pos < qlen);
			rightclipseq[i] = seq_nt16_str[bam_seqi(seq_ptr, qlen - seqlen + i)];
			// rightclipseq[i] = seq_nt16_str[bam_seqi(seq_ptr, i + ql1)
// -			];
		}
		cout << "rightclipseq: " << rightclipseq << endl; //CLEAN:
		sc_info.push_back(rightclipseq);
		// get right anchor position
		int anchorpos = -1;
		if(seqrev) 
		{
			pair<int, int> pos_pair = subseq_pos(anchor_start, revcomp(rightclipseq), anchor_start_nm);
			right_anchor_padding = pos_pair.second >= 0 ? pos_pair.second : -1;
		}
		else
		{
			pair<int, int> pos_pair = subseq_pos(anchor_end, rightclipseq, anchor_end_nm);
			right_anchor_padding = pos_pair.first >= 0 ? seqlen - pos_pair.first : -1;
		}
	
		// right_anchor_padding = anchorpos >= 0? anchorpos: -1;
		cout << "right anchor padding: " << right_anchor_padding << endl;
		sc_info.push_back(to_string(right_anchor_padding));
	}

	std::ofstream anchor_file(anchor_file_name, std::ios::app);
	for(int i=0; i<sc_info.size(); i++) anchor_file << sc_info[i] << "\t";
	anchor_file << endl;
	anchor_file.close();

	return 0;
}

bool hit::is_anchor_satisfactory(int side, int padding_max = 20)
{
	bool seqrev;
	if((flag & 0x10) >= 1) seqrev = true;
	if((flag & 0x10) <= 0) seqrev = false;

	assert (side == 0 || side == 1);

	if (side == 0) // left
	{
		if ((anchor_start != "" && !seqrev) || (anchor_end != "" && seqrev))
		{
			if (left_anchor_padding < 0) return false; 
			if (left_anchor_padding > padding_max) return false;
			return true;
		}
		
	}

	if (side == 1) // right
	{
		if ((anchor_end != "" && !seqrev) || (anchor_start != "" && seqrev)) 
		{
			if (right_anchor_padding < 0) return false; 
			if (right_anchor_padding > padding_max) return false;
			return true;
		}
	}
	
	// anchor seq not provided
	return true;
}

int hit::get_aligned_intervals(vector<int64_t> &v) const
{
	v.clear();
	int32_t p1 = pos;
	for(int k = 0; k < spos.size(); k++)
	{
		int32_t p2 = high32(spos[k]);
		v.push_back(pack(p1, p2));
		p1 = low32(spos[k]);
	}
	v.push_back(pack(p1, rpos));
	return 0;
}

string hit::get_qname(bam1_t *b)
{
	char buf[1024];
	char *q = bam_get_qname(b);
	int l = strlen(q);
	memcpy(buf, q, l);
	buf[l] = '\0';
	return string(buf);
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts");
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);
	if(p0 && (*p0) == 'a') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS");
	if(p1 && (*p1) == 'A') xs = bam_aux2A(p1);
	if(p1 && (*p1) == 'a') xs = bam_aux2A(p1);

	if(xs == '.' && ts != '.')
	{
		// convert ts to xs
		if((flag & 0x10) >= 1 && ts == '+') xs = '-';
		if((flag & 0x10) >= 1 && ts == '-') xs = '+';
		if((flag & 0x10) <= 0 && ts == '+') xs = '+';
		if((flag & 0x10) <= 0 && ts == '-') xs = '-';
	}

	hi = -1;
	uint8_t *p2 = bam_aux_get(b, "HI");
	if(p2 && (*p2) == 'C') hi = bam_aux2i(p2);
	if(p2 && (*p2) == 'c') hi = bam_aux2i(p2);

	nh = -1;
	uint8_t *p3 = bam_aux_get(b, "NH");
	if(p3 && (*p3) == 'C') nh = bam_aux2i(p3);
	if(p3 && (*p3) == 'c') nh = bam_aux2i(p3);

	nm = 0;
	uint8_t *p4 = bam_aux_get(b, "nM");
	if(p4 && (*p4) == 'C') nm = bam_aux2i(p4);
	if(p4 && (*p4) == 'c') nm = bam_aux2i(p4);

	uint8_t *p5 = bam_aux_get(b, "NM");
	if(p5 && (*p5) == 'C') nm = bam_aux2i(p5);
	if(p5 && (*p5) == 'c') nm = bam_aux2i(p5);

	// set umi
        umi = "";
        uint8_t *p6 = bam_aux_get(b, "UB");
        if(p6 && (*p6) == 'H') umi = bam_aux2Z(p6);
	if(p6 && (*p6) == 'Z') umi = bam_aux2Z(p6);

	/*	
	// TODO: check if UB = UX
	uint8_t *p7 = bam_aux_get(b, "UX");
	string raw_umi;
	if(p7&& (*p7) == 'Z') raw_umi = bam_aux2Z(p7);
	if(raw_umi != "")
	{
		if(umi != raw_umi) 
		{
			printf("raw_umi: %s, umi = %s\n", raw_umi.c_str(), umi.c_str());
			umi = "";
		}
	}
	*/
	

        //printf("qname: %s, umi = %s\n", qname.c_str(), umi.c_str());

	return 0;
}

int hit::set_concordance()
{
	bool concordant = false;
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// F1R2
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// R1F2
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// F2R1
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// R2F1
	return 0;
}

int hit::set_strand()
{
	strand = '.';
	
	if(library_type == FR_FIRST && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
	}

	if(library_type == FR_SECOND && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
	}

	if(library_type == FR_FIRST && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '-';
		if((flag & 0x10) >= 1) strand = '+';
	}

	if(library_type == FR_SECOND && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '+';
		if((flag & 0x10) >= 1) strand = '-';
	}

	return 0;
}

bool hit::operator<(const hit &h) const
{
	if(qname < h.qname) return true;
	if(qname > h.qname) return false;
	if(hi != -1 && h.hi != -1 && hi < h.hi) return true;
	if(hi != -1 && h.hi != -1 && hi > h.hi) return false;
	return (pos < h.pos);
}

int hit::print() const
{
	// print basic information
	printf("Hit %s: hid = %d, [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %lu, qlen = %d, hi = %d, nh = %d, umi = %s, bridged = %c\n", 
			qname.c_str(), hid, pos, rpos, mpos, flag, qual, strand, xs, ts, isize, qlen, hi, nh, umi.c_str(), bridged ? 'T' : 'F');

	/*
	printf(" start position (%d - )\n", pos);
	for(int i = 0; i < spos.size(); i++)
	{
		int64_t p = spos[i];
		int32_t p1 = high32(p);
		int32_t p2 = low32(p);
		printf(" splice position (%d - %d)\n", p1, p2);
	}
	printf(" end position (%d - )\n", rpos);
	*/

	return 0;
}

vector<int> encode_vlist(const vector<int> &v)
{
	vector<int> vv;
	if(v.size() <= 0) return vv;

	int p = v[0];
	int k = 1;
	for(int i = 1; i < v.size(); i++)
	{
		if(v[i] == v[i - 1] + 1)
		{
			k++;
		}
		else
		{
			assert(k >= 1);
			vv.push_back(p);
			vv.push_back(k);
			p = v[i];
			k = 1;
		}
	}
	vv.push_back(p);
	vv.push_back(k);

	/*
	printf("encode: (");
	printv(v);
	printf(") -> (");
	printv(vv);
	printf(")\n");
	*/
	return vv;
}

vector<int> decode_vlist(const vector<int> &v)
{
	vector<int> vv;
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() / 2; i++)
	{
		int p = v[i * 2 + 0];
		int k = v[i * 2 + 1];
		for(int j = p; j < p + k; j++)
		{
			vv.push_back(j);
		}
	}
	return vv;
}
