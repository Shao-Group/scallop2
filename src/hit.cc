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

#include "hit.h"
#include "config.h"
#include "util.h"

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

hit::hit()
{
	flag = 0;
	bin = 0;
	qual = 0;
	l_extranul = 0;

	hid = 0;
	qname = "";
	qhash = 0;
	paired = false;
	bridged = false;
	next = NULL;
	suppl = NULL;

	supple_pos = 0;
	is_reverse_overlap = false;
	is_fake = false;
	fake_hit_index = -1;
	soft_clip_side = 0;
	
	left_cigar = '.';					// S=soft clip, H=hard clip, M=match, .=default
	right_cigar = '.';					// S=soft clip, H=hard clip, M=match, .=default
	left_cigar_len = 0;
	right_cigar_len = 0;
	first_pos = 0;						//.H.M. the three dots are the 1st, 2nd, and 3rd pos respectively
	second_pos = 0;
	third_pos = 0;

	l_qseq = 0;
	seq = "";
	soft_left_clip_seqs.clear();
	soft_right_clip_seqs.clear();

	rpos = 0;
	qlen = 0;
	nh = 0;
	hi = 0;
	nm = 0;
	sa = "";

	cigar_vector.clear();
	spos.clear();
	vlist.clear();
	itvm.clear();

	itvi.clear();
	itvd.clear();
	
	concordant = true;
	paired = false;
	bridged = false;
	strand = '.';
	xs = '.';
	ts = '.';

	umi = "";
	pi = 0;
	fidx = 0;
}

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
	sa = h.sa;
	supple_pos = h.supple_pos;
	suppl = h.suppl;
	is_reverse_overlap = h.is_reverse_overlap;
	is_fake = h.is_fake;
	fake_hit_index = h.fake_hit_index;
	soft_clip_side = h.soft_clip_side;

	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;

	vlist = h.vlist;
	concordant = h.concordant;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;
	pi = h.pi;
	fidx = h.fidx;

	cigar_vector = h.cigar_vector;
	left_cigar = h.left_cigar;					// S=soft clip, H=hard clip, M=match, .=default
	right_cigar = h.right_cigar;					// S=soft clip, H=hard clip, M=match, .=default
	left_cigar_len = h.left_cigar_len;
	right_cigar_len = h.right_cigar_len;
	first_pos = h.first_pos;						//.H.M. the three dots are the 1st, 2nd, and 3rd pos respectively
	second_pos = h.second_pos;
	third_pos = h.third_pos;

	l_qseq = h.l_qseq;
	seq = h.seq;
	soft_left_clip_seqs = h.soft_left_clip_seqs;
	soft_right_clip_seqs = h.soft_right_clip_seqs;

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
	sa = h.sa;
	supple_pos = h.supple_pos;
	is_reverse_overlap = h.is_reverse_overlap;
	is_fake = h.is_fake;
	fake_hit_index = h.fake_hit_index;
	soft_clip_side = h.soft_clip_side;
	suppl = h.suppl;

	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;

	vlist = h.vlist;
	concordant = h.concordant;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;
	pi = h.pi;
	fidx = h.fidx;

	cigar_vector = h.cigar_vector;
	left_cigar = h.left_cigar;					// S=soft clip, H=hard clip, M=match, .=default
	right_cigar = h.right_cigar;					// S=soft clip, H=hard clip, M=match, .=default
	left_cigar_len = h.left_cigar_len;
	right_cigar_len = h.right_cigar_len;
	first_pos = h.first_pos;						//.H.M. the three dots are the 1st, 2nd, and 3rd pos respectively
	second_pos = h.second_pos;
	third_pos = h.third_pos;

	l_qseq = h.l_qseq;
	seq = h.seq;
	soft_left_clip_seqs = h.soft_left_clip_seqs;
	soft_right_clip_seqs = h.soft_right_clip_seqs;
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
	suppl = NULL;

	supple_pos = 0;
	is_reverse_overlap = false;
	is_fake = false;
	fake_hit_index = -1;
	soft_clip_side = 0;
	
	left_cigar = '.';					// S=soft clip, H=hard clip, M=match, .=default
	right_cigar = '.';					// S=soft clip, H=hard clip, M=match, .=default
	left_cigar_len = 0;
	right_cigar_len = 0;
	first_pos = 0;						//.H.M. the three dots are the 1st, 2nd, and 3rd pos respectively
	second_pos = 0;
	third_pos = 0;

	cigar_vector.clear();
	spos.clear();
	vlist.clear();
	itvm.clear();
	itvi.clear();
	itvd.clear();

	l_qseq = 0;
	seq = "";
	soft_left_clip_seqs.clear();
	soft_right_clip_seqs.clear();

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));
	qlen = (int32_t)bam_cigar2qlen(n_cigar, bam_get_cigar(b));
	//printf("rpos = %d, qlen = %d\n", rpos, qlen);

	/*if(strcmp(qname.c_str(),"SRR1721290.17627808") == 0)
	{
		print();
	}*/

	//print();

	// get cigar
	assert(n_cigar <= max_num_cigar);
	assert(n_cigar >= 1);
	uint32_t * cigar = bam_get_cigar(b); //commented by Tasfia

	// build splice positions
	spos.clear();
	int32_t p = pos;
	int32_t q = 0;

	//printf("qname: %s, n_cigar size: %d\n",qname.c_str(),n_cigar);

    for(int k = 0; k < n_cigar; k++)
	{
		/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference/**/

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);
		//printf("p=%d,",p);	

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			q += bam_cigar_oplen(cigar[k]);
		//printf("q=%d,",q);

		//printf("cigar size = %d\n",bam_cigar_op(cigar[0]));

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue; //BAM_CREF_SKIP junction
		if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue; //match
		if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue; //match

		// consider ALL splice positions
		//if(bam_cigar_oplen(cigar[k-1]) < min_flank_length) continue;
		//if(bam_cigar_oplen(cigar[k+1]) < min_flank_length) continue;
		//printf("\nPassed conditions of spos\n");

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		//printf("s=%d\n",s);
		spos.push_back(pack(s, p));

	}

	cigar_vector.clear();

	for(int k = 0; k < n_cigar; k++)
	{
		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			cigar_vector.push_back(pair<char, int32_t>('M',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('S',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('H',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CREF_SKIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('N',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			cigar_vector.push_back(pair<char, int32_t>('I',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			cigar_vector.push_back(pair<char, int32_t>('D',bam_cigar_oplen(cigar[k])));
		}
		else
		{
			cigar_vector.push_back(pair<char, int32_t>('.',0));
		}
	}
	
	//putting a placeholder to avoid core dumped when cigar_vector[0] is accessed
	if(cigar_vector.size() == 0)
	{
		cigar_vector.push_back(pair<char, int32_t>('.',0));
	}

	//assign booleans to see if left splice position H/S and right M or vie versa and stor their lengths
	/*if(n_cigar == 2)
	{
		if(bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP && bam_cigar_op(cigar[1]) == BAM_CMATCH)
		{
			//printf("First case\n");
			left_cigar = 'S';
			right_cigar = 'M';
			first_pos = pos - bam_cigar_oplen(cigar[0]);
			second_pos = pos;
			third_pos = rpos;

			//if(hid == 11789) printf("SM positions: %d-%d-%d\n", first_pos,second_pos,third_pos);
		}
		else if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP && bam_cigar_op(cigar[1]) == BAM_CMATCH)
		{
			//printf("Second case\n");
			left_cigar = 'H';
			right_cigar = 'M';
			first_pos = pos - bam_cigar_oplen(cigar[0]);
			second_pos = pos;
			third_pos = rpos;

			//if(hid == 11789) printf("HM positions: %d-%d-%d\n", first_pos,second_pos,third_pos);
		}
		else if(bam_cigar_op(cigar[0]) == BAM_CMATCH && bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP)
		{
			//printf("Third case\n");
			left_cigar = 'M';
			right_cigar = 'S';
			first_pos = pos;
			second_pos = rpos;
			third_pos = rpos + bam_cigar_oplen(cigar[1]);;

			//if(hid == 20562) printf("MS positions: %d-%d-%d\n", first_pos,second_pos,third_pos);
		}
		else if(bam_cigar_op(cigar[0]) == BAM_CMATCH && bam_cigar_op(cigar[1]) == BAM_CHARD_CLIP)
		{
			//printf("Fourth case\n");
			left_cigar = 'M';
			right_cigar = 'H';
			first_pos = pos;
			second_pos = rpos;
			third_pos = rpos + bam_cigar_oplen(cigar[1]);;

			//if(hid == 11789) printf("MH positions: %d-%d-%d\n", first_pos,second_pos,third_pos);
		}
	}*/


	//printf("spos size: %d\n", spos.size());

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
	}

	//printf("call regular constructor\n");
	//printf("end ...................................\n");
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

string hit::convert_to_IUPAC(vector<int> code)
{
	string seq = "";
	for(int i=0;i<code.size();i++)
	{
		if(code[i] == 1)
		{
			seq = seq + 'A';
		}
		else if(code[i] == 2)
		{
			seq = seq + 'C';
		}
		else if(code[i] == 4)
		{
			seq = seq + 'G';
		}
		else if(code[i] == 8)
		{
			seq = seq + 'T';
		}
		else if(code[i] == 15)
		{
			seq = seq + 'N';
		}
	}

	return seq;
}

int hit::set_seq(bam1_t *b)
{
	uint32_t seq_len = b->core.l_qseq;
	l_qseq = seq_len;
	uint8_t *q = bam_get_seq(b); 
	vector<int> code;

	for(int i=0;i<seq_len;i++)
	{
		code.push_back(bam_seqi(q,i)); //gets nucleotide id
	}

	string hit_seq = convert_to_IUPAC(code);
	seq = hit_seq;

	/*printf("%s test print seq:\n",qname.c_str());
	print_cigar();

	if((flag & 0x10) >= 1)
	{
		printf("rev comp\n");
	}
	printf("seq: %s\n",seq.c_str());*/

	set_soft_clip_seq_combo();

	return 0;
}

int hit::set_soft_clip_seq_combo()
{
	if(cigar_vector[0].first == 'S')
	{
		int32_t len = cigar_vector[0].second;
		//index 0, extract start len bp
		string str0 = "";
		for(int i=0;i<len;i++)
		{
			str0 = str0 + seq[i];
		}
		soft_left_clip_seqs.push_back(str0);

		//index 1, extract start len bp rev comp
		soft_left_clip_seqs.push_back(get_reverse_complement(soft_left_clip_seqs[0]));
	}
	if(cigar_vector[cigar_vector.size()-1].first == 'S')
	{
		int32_t len = cigar_vector[cigar_vector.size()-1].second;

		//index 0, extract end len bp
		string str2 = "";
		for(int i=seq.size()-len;i<seq.size();i++)
		{
			str2 = str2 + seq[i];
		}
		soft_right_clip_seqs.push_back(str2);

		//index 1,extract end len bp rev comp
		soft_right_clip_seqs.push_back(get_reverse_complement(soft_right_clip_seqs[0]));
	}

	/*if(strcmp(qname.c_str(),"simulate:122")== 0)
	{
		printf("cigar of simulate:122:\n");
		for(int i=0;i<cigar_vector.size();i++)
		{
			printf("%d%c",cigar_vector[i].second,cigar_vector[i].first);
		}
	}*/

	/*printf("Printing four combos:\n");
	printf("start: %s\n",soft_clip_seqs[0].c_str());
	printf("start RC: %s\n",soft_clip_seqs[1].c_str());
	printf("end: %s\n",soft_clip_seqs[2].c_str());
	printf("end RC: %s\n",soft_clip_seqs[3].c_str());*/

	return 0;
}

string hit::get_reverse_complement(string str)
{
	string out = "";
	for(int i=str.size()-1;i>=0;i--)
	{
		if(str[i] == 'A')
		{
			out = out + 'T';
		}
		else if(str[i] == 'T')
		{
			out = out + 'A';
		}
		else if(str[i] == 'C')
		{
			out = out + 'G';
		}
		else if(str[i] == 'G')
		{
			out = out + 'C';
		}
		else if(str[i] == 'N')
		{
			out = out + str[i];
		}
	}
	return out;
}


string hit::get_complement(string str)
{
	string out = "";
	for(int i=0;i<str.size();i++)
	{
		if(str[i] == 'A')
		{
			out = out + 'T';
		}
		else if(str[i] == 'T')
		{
			out = out + 'A';
		}
		else if(str[i] == 'C')
		{
			out = out + 'G';
		}
		else if(str[i] == 'G')
		{
			out = out + 'C';
		}
		else if(str[i] == 'N')
		{
			out = out + str[i];
		}
	}
	return out;
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts"); //ts used by minimap2
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);
	if(p0 && (*p0) == 'a') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS"); //used by star and hisat
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

	sa = "";
	uint8_t *p7 = bam_aux_get(b, "SA"); //sa tag has the supple pos and cigar of curr hit, ex: SA:Z:chr1,14068602,+,51M99H,255,0;
	if(p7 && (*p7) == 'Z') sa = bam_aux2Z(p7);
	if(p7 && (*p7) == 'z') sa = bam_aux2Z(p7);

	//printf("umi tag:%s\n",umi.c_str());
	//printf("sa tag:%s\n",sa.c_str());

	string sa_tag = sa.c_str();

	vector<string> results;
	stringstream  ss(sa_tag);
	string str;
	while (getline(ss, str, ',')) {
		results.push_back(str);
	}

	if(results.size() > 0)
	{
		/*for(int i=0;i<results.size();i++)
		{
			printf("%s ",results[i].c_str());
		}
		printf("\n");*/

		supple_pos = atoi(results[1].c_str());
		//printf("supple_pos = %d\n",supple_pos);
	}


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
	
	if(library_type == FR_FIRST && ((flag & 0x1) >= 1)) //data strand specific if FR_FIRST/FR_SECOND
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
	printf("Hit %s: hid = %d, [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, qlen = %d, hi = %d, nh = %d, umi = %s, sa = %s, supple_pos = %d, bridged = %c, paired = %c, is_fake = %c, fake_hit_index = %d, spos_size = %zu\n", qname.c_str(), hid, pos, rpos, mpos, flag, qual, strand, xs, ts, isize, qlen, hi, nh, umi.c_str(), sa.c_str(), supple_pos, bridged ? 'T' : 'F',paired ? 'T' : 'F', is_fake ? 'T' : 'F', fake_hit_index, spos.size());
	//printf("vlist size %lu\n",vlist.size());
	for(int i=0;i<vlist.size();i++)
	{
		printf("%d,",vlist[i]);
	}
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

int hit::print_cigar()
{
	printf("%s cigar:",qname.c_str());
	for(int j=0;j<cigar_vector.size();j++)
	{
		printf("%d%c",cigar_vector[j].second,cigar_vector[j].first);
	}
	printf("\n");
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
	vv.clear();
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
