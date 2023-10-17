/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral package
(c) 2019 by Mingfu Shao, The Pennsylvania State University
Part of Scallop2                                               
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/
                                   
#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "htslib/sam.h"
#include "config.h"

using namespace std;


/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;
*/

struct hit : public bam1_core_t
{
public:
	//hit(int32_t p);
	hit();
	hit(bam1_t *b, int id);
	hit(const hit &h);
	bool operator<(const hit &h) const;
	hit& operator=(const hit &h);

public:
	int hid;								// hit-id

	char left_cigar;						// S=soft clip, H=hard clip, M=match, I-insertion ... etc
	char right_cigar;						// S=soft clip, H=hard clip, M=match, I-insertion ... etc

	int32_t left_cigar_len;					// store len of left cigar
	int32_t right_cigar_len;				// store len of right cigar

	int32_t first_pos;						//.H.M. the three dots are the 1st, 2nd, and 3rd pos respectively
	int32_t second_pos;
	int32_t third_pos;
	vector<pair<char, int32_t>> cigar_vector; 	//stores all cigars of a hit with length


	//vector<uint32_t> cigar_positions;		// stores putative back splice positions
	vector<int64_t> spos;					// splice positions
	vector<int> vlist;						// list of spanned vertices in the junction graph
	int32_t rpos;							// right position mapped to reference [pos, rpos)
	int32_t qlen;							// read length
	int32_t nh;								// NH aux in sam
	int32_t hi;								// HI aux in sam
	int32_t nm;								// NM aux in sam
	string sa;								// SA aux in sam
	int32_t supple_pos;						// stores position of supple from SA tag
	bool is_reverse_overlap;				// whether this is a RO read
	size_t qhash;							// hash code for qname
	hit *suppl;								// supplementary hit
	bool is_fake;							// whether this is a fake hit
	int fake_hit_index; 					//used to keep track of fragment for which fake hit is created, stores index of partner fragment
	int soft_clip_side;						//used to keep track of whether the fake hit comes from a soft left clip (1) or soft right clip (2)

	// scallop+coral
	vector<int64_t> itvm;					// matched interval
	vector<int64_t> itvi;					// insert interval
	vector<int64_t> itvd;					// delete interval

	bool concordant;						// whether it is concordant
	bool paired;							// whether this hit has been paired
	bool bridged;							// whether this hit has been bridged 
	char strand;							// strandness
	char xs;								// XS aux in sam
	char ts;								// ts tag used in minimap2
	string qname;							// query name
	uint32_t l_qseq;						// length of sequence information of the read except H clip
	string seq;								// sequence information of the read except H clip
	vector<string> soft_left_clip_seqs;		// index 0:start seq,1:start seq rev comp
	vector<string> soft_right_clip_seqs;	// index 0:ending seq,1:ending seq rev comp
	hit *next;								// next hit that is equivalent with current one

	// UMI
	string umi;
	int pi;							// paired hits index
	int fidx;						// its fragments index

public:
 	static string get_qname(bam1_t *b);
	string get_reverse_complement(string str);
	string get_complement(string str);
	string convert_to_IUPAC(vector<int> code);
	int set_soft_clip_seq_combo();
	int set_seq(bam1_t *b);
	int set_tags(bam1_t *b);
	int set_strand();
	int set_concordance();
	int get_aligned_intervals(vector<int64_t> &v) const;
	int print() const;
	int print_cigar();
};

vector<int> encode_vlist(const vector<int> &v);
vector<int> decode_vlist(const vector<int> &v);

//inline bool hit_compare_by_name(const hit &x, const hit &y);

#endif
