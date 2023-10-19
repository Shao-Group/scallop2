/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle_bridge.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "bridger.h"

bundle_bridge::bundle_bridge(bundle_base &b, reference &r)
	: bb(b), ref(r)
{
	circ_trsts.clear(); // emptying before storing circRNAs
	circ_trsts_HS.clear();
	circ_fragments.clear();
	RO_count = 0;

	compute_strand();
	ref_trsts = ref.get_overlapped_transcripts(bb.chrm, bb.strand, bb.lpos, bb.rpos);
	//build(RO_reads_map, fai);
}

bundle_bridge::bundle_bridge(bundle_base &b, reference &r, map <string, int> RO_reads_map, faidx_t *fai)
	: bb(b), ref(r)
{
	circ_trsts.clear(); // emptying before storing circRNAs
	circ_trsts_HS.clear();
	circ_fragments.clear();
	RO_count = 0;

	compute_strand();
	ref_trsts = ref.get_overlapped_transcripts(bb.chrm, bb.strand, bb.lpos, bb.rpos);
	build(RO_reads_map, fai);
}

bundle_bridge::~bundle_bridge()
{}

int bundle_bridge::compute_strand()
{
	if(library_type != UNSTRANDED) assert(bb.strand != '.');
	if(library_type != UNSTRANDED) return 0;

	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	// a new way of defining strandness
	// resulting to more "." strands
	if(np > 0 && nq == 0) bb.strand = '+';
	else if(nq > 0 && np == 0) bb.strand = '-';
	else bb.strand = '.';

	/*
	if(np > nq) bb.strand = '+';
	else if(np < nq) bb.strand = '-';
	else bb.strand = '.';
	*/
	
	return 0;
}

int bundle_bridge::build(map <string, int> RO_reads_map, faidx_t *_fai)
{
	fai = _fai;

	/*if(fai != NULL)
	{
		printf("extracting fasta seq from region:\n");
		int32_t seqlen;
		char* seq = faidx_fetch_seq(fai, "1", 39511735, 39511804, &seqlen);
		if(seq != NULL && seqlen > 0)
		{
			printf("seqlen = %d, seq = %s\n",seqlen,seq);
		}
	}*/

	set_hits_RO_parameter(RO_reads_map);

	build_supplementaries();
	set_chimeric_cigar_positions(); //setting h.first_pos/second_pos etc for getting back splice positions using cigars 
	build_junctions();
	extend_junctions();

	build_regions();
	build_partial_exons();

	align_hits_transcripts();
	index_references();

	//printf("Called from bundle_bridge\n");
	build_fragments(); //builds fragment from h1p to h2
	//printf("\n");

	fix_alignment_boundaries();
	build_circ_fragments(); //will build fragment from h2 to h1s, added by Tasfia

	//group_fragments();

	//extract_all_non_supple_HS_hits();
	//extract_nonsupple_HS_hits();

	remove_tiny_boundaries();
	set_fragment_lengths();

	//RO reads statistics
	//get_frags_with_HS_on_both_sides();

	//create circ fragments from RO reads with H/S on both sides uisng ciri-full
	//get_RO_frags_with_HS();

	//create circ fragments from frags with H/S on both sides using our data instead of ciri-full
	//get_frags_with_HS_from_data();

	//find more chimeric reads from soft clip reads
	get_more_chimeric();

	//create vlist of fake hits
	align_fake_hits();

	//create fragments from fake hits
	create_fake_fragments();

	// call remove_tiny_boundary for the newly 
	// created (fake) hits here
	remove_tiny_boundaries_fake_hits();

	set_circ_fragment_lengths();

	bridger bdg(this);
	bdg.bridge_normal_fragments();
	bdg.bridge_circ_fragments();

	//extract_RO_circRNA();
	//extract_HS_frags_circRNA();

	extract_circ_fragment_pairs();
	//print_circ_fragment_pairs();
	join_circ_fragment_pairs(bdg.length_high);
	print_circRNAs();

	//printf("fragments vector size after = %zu\n",fragments.size());

	/*for(int i=0;i<fragments.size();i++)
	{
		if(fragments[i].frag_type == 2)
		{
			printf("frag pi = %d\n",fragments[i].pi);
			fragments[i].print(i+1);
		}
	}*/
	return 0;
}

int bundle_bridge::set_hits_RO_parameter(map <string, int> RO_reads_map)
{
	for(int i=0;i<bb.hits.size();i++)
	{
		hit *curr_hit = &bb.hits[i];
		string hash = "";
		hash = hash + bb.chrm + ":" + curr_hit->qname;

		if(RO_reads_map.find(hash) != RO_reads_map.end())
		{
			//hit is a RO read
			curr_hit->is_reverse_overlap = true;
		}
	}

	int cnt = 0;
	for(int i=0;i<bb.hits.size();i++)
	{
		if(bb.hits[i].is_reverse_overlap == true)
		{
			cnt++;
		}
	}
	//printf("# RO reads in bundle %d\n",cnt);
	RO_count = cnt;

	return 0;
}

int bundle_bridge::get_frags_with_HS_on_both_sides()
{
	HS_both_side_reads.clear();
	chimeric_reads.clear();

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		if(fr.h1->pos <= fr.h2->pos && (fr.h1->cigar_vector[0].first == 'S' || fr.h1->cigar_vector[0].first == 'H') && (fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'S' || fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'H'))
		{
			string chrm_id = bb.chrm.c_str();
			string name = bb.chrm+":"+fr.h1->qname;
			HS_both_side_reads.push_back(name);
		}
		else if(fr.h1->pos > fr.h2->pos && (fr.h2->cigar_vector[0].first == 'S' || fr.h2->cigar_vector[0].first == 'H') && (fr.h1->cigar_vector[fr.h1->cigar_vector.size()-1].first == 'S' || fr.h1->cigar_vector[fr.h1->cigar_vector.size()-1].first == 'H'))
		{
			string chrm_id = bb.chrm.c_str();
			string name = bb.chrm+":"+fr.h1->qname;
			HS_both_side_reads.push_back(name);
		}

		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			string chrm_id = bb.chrm.c_str();
			string name = bb.chrm+":"+fr.h1->qname;
			chimeric_reads.push_back(name);
		}
	}

	return 0;
}

int bundle_bridge::get_frags_with_HS_from_data()
{
	int junc_range = 10;
	for(int k = 0; k < fragments.size(); k++)
	{
		fragment fr = fragments[k];

		//has a supple
		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			continue;
		}

		//is a supple
		if((fr.h1->flag & 0x800) >= 1 || (fr.h2->flag & 0x800) >= 1) continue;

		int left_boundary_flag = 0;
		int right_boundary_flag = 0;

		if(fr.h1->pos > fr.h2->pos) continue;
		if(fr.h1->cigar_vector[0].first != 'S' || fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first != 'S') continue;
		if(fr.h1->cigar_vector[0].second < 5 && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second < 5) continue;

		printf("HS paired hit case 1: pos %d, rpos %d\n",fr.lpos,fr.rpos);
		printf("chrm %s\n",bb.chrm.c_str());
		printf("Hit 1: ");
		fr.h1->print();

		//checking if reads junction matches left boundary
		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];

			//if(jc.rpos <= fr.h1->pos+junc_range && jc.rpos >= fr.h1->pos-junc_range)
			if(jc.rpos == fr.h1->pos)
			{
				printf("reads junction present left %d\n",jc.rpos);
				left_boundary_flag = 1;
				break;
			}
		}

		//checking if ref junction matches left boundary
		int temp_flag = 0;
		for(int t=0;t<ref_trsts.size();t++)
		{
			transcript trst = ref_trsts[t];
			vector<PI32> chain = trst.get_intron_chain();

			for(int k=0;k<chain.size();k++)
			{
				//if(chain[k].second <= fr.h1->pos+junc_range && chain[k].second >= fr.h1->pos-junc_range)
				if(chain[k].second == fr.h1->pos)
				{
					printf("ref junction present left %d\n",chain[k].second);
					temp_flag = 1;
					left_boundary_flag = 1;
					break;
				}
			}

			if(temp_flag == 1)
			{
				break;
			}
		}

		//checking if pexon boundary matches left boundary if ref not given
		if(ref_file == "")
		{
			for(int p=0;p<pexons.size();p++)
			{
				if(pexons[p].lpos == fr.h1->pos && pexons[p].ltype == START_BOUNDARY)
				{
					left_boundary_flag = 1;
					break;
				}
			}
		}

		printf("Hit 2: ");
		fr.h2->print();

		//checking if reads junction matches right boundary
		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];

			//if(jc.lpos <= fr.h2->rpos+junc_range && jc.lpos >= fr.h2->rpos-junc_range)
			if(jc.lpos == fr.h2->rpos)
			{
				printf("reads junction present right %d\n",jc.lpos);
				right_boundary_flag = 1;
				break;
			}
		}

		//checking if ref junction matches right boundary
		temp_flag = 0;
		for(int t=0;t<ref_trsts.size();t++)
		{
			transcript trst = ref_trsts[t];
			vector<PI32> chain = trst.get_intron_chain();

			for(int k=0;k<chain.size();k++)
			{
				//if(chain[k].first <= fr.h2->rpos+junc_range && chain[k].first >= fr.h2->rpos-junc_range)
				if(chain[k].first == fr.h2->rpos)
				{
					printf("ref junction present right %d\n",chain[k].first);
					right_boundary_flag = 1;
					temp_flag = 1;
					break;
				}
			}

			if(temp_flag == 1)
			{
				break;
			}
		}

		//checking if pexon boundary matches left boundary if ref not given
		if(ref_file == "")
		{
			for(int p=0;p<pexons.size();p++)
			{
				if(pexons[p].rpos == fr.h2->rpos && pexons[p].rtype == END_BOUNDARY)
				{
					right_boundary_flag = 1;
					break;
				}
			}
		}

		printf("HS frags: left_boundary_flag = %d, right_boundary_flag = %d\n\n",left_boundary_flag,right_boundary_flag);

		/*if(left_boundary_flag == 1 && right_boundary_flag == 1)
		{
			fr.HS_frag = true;
			circ_fragments.push_back(fr);
		}*/
		
		if(left_boundary_flag != 1 || right_boundary_flag != 1) continue;

		int seq_match_left_hit = 0;
		int seq_match_right_hit = 0;
		int max_read_to_junction_gap = 100000;

		int32_t soft_len = fr.h1->cigar_vector[0].second;
		if(soft_len < 10) continue;

		//create a hash map for the kmers in soft clip region and pass that for similarity testing
		int kmer_length = 10;
		map<string,int> kmer_map;
		kmer_map.clear();

		assert(fr.h1->soft_left_clip_seqs.size() > 0);
		string s = fr.h1->soft_left_clip_seqs[0];
		//printf("s size %lu\n",s.size());
		
		for(int i=0;i<=s.size()-kmer_length;i++)
		{
			string kmer = s.substr(i,kmer_length);

			if(kmer_map.find(kmer) == kmer_map.end())
			{
				kmer_map.insert(pair<string,int>(kmer,1));
			}
		}

		int32_t prev_pos2 = 0;
		int jc_multiple = 0;

		//match left hit soft seq with junctions seq
		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];
			if(jc.lpos <= fr.h2->rpos || jc.lpos <= fr.h1->rpos) continue;
			if(abs(jc.lpos-fr.h2->rpos) > max_read_to_junction_gap) continue;

			int32_t pos1 = jc.lpos-soft_len+1;
			int32_t pos2 = jc.lpos;

			string junc_seq = get_fasta_seq(pos1,pos2);

			int edit_match = 0;
			assert(junc_seq.size() == fr.h1->soft_left_clip_seqs[0].size());

			bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
			if(is_similar == false) continue;

			int edit = get_edit_distance(junc_seq,fr.h1->soft_left_clip_seqs[0]);

			if(edit <= floor(soft_len/10))
			{
				edit_match = 1;
				printf("RO recover editmatch case 1: junc seq pos1=%d, pos2=%d, junc_seqlen = %lu\n",pos1,pos2,junc_seq.size());
			}

			if(edit_match == 1)
			{
				if(prev_pos2 != 0 && pos2 != prev_pos2)
				{
					printf("pos2 = %d, prev_pos2 = %d\n",pos2,prev_pos2);
					jc_multiple = 1;
					break;
				}

				prev_pos2 = pos2;
			}
		}

		if(jc_multiple == 1) 
		{
			printf("Ro recover jc_multiple = %d\n",jc_multiple);
			continue;
		}

		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];

			if(jc.lpos <= fr.h2->rpos || jc.lpos <= fr.h1->rpos) continue;
			if(abs(jc.lpos-fr.h2->rpos) > max_read_to_junction_gap) continue;

			int32_t pos1 = jc.lpos-soft_len+1;
			int32_t pos2 = jc.lpos;

			string junc_seq = get_fasta_seq(pos1,pos2);

			assert(junc_seq.size() == fr.h1->soft_left_clip_seqs[0].size());

			bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
			if(is_similar == false) continue;

			int edit = get_edit_distance(junc_seq,fr.h1->soft_left_clip_seqs[0]);

			//if(edit == 0 || edit == 1)
			if(edit <= floor(soft_len/10))
			{
				printf("both sides S clip left: chrm=%s, read=%s, read_pos=%d\n",bb.chrm.c_str(),fr.h1->qname.c_str(),fr.h1->pos);
				if((fr.h1->flag & 0x10) >= 1)
				{
					//printf("rev comp 0x10 = on\n");
				}
				else
				{
					//printf("rev comp 0x10 = off\n");
				}
				if((fr.h1->flag & 0x4) >= 1)
				{
					//printf("seg unmapped 0x4 = on\n");
				}
				else
				{
					//printf("seg unmapped 0x4 = off\n");
				}
				printf("read seq combo index=0, combo_seq=%s, edit=%d\n",fr.h1->soft_left_clip_seqs[0].c_str(),edit);
				printf("junction lpos = %d, rpos = %d\n",jc.lpos,jc.rpos);
				printf("junc seq pos1=%d, pos2=%d, junc_seqlen = %lu, junc_seq=%s\n",pos1,pos2,junc_seq.size(),junc_seq.c_str());

				seq_match_left_hit = 1;
				break;
			}
		}

		soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;
		if(soft_len < 10) continue;

		//create a hash map for the kmers in soft clip region and pass that for similarity testing
		kmer_map.clear();

		assert(fr.h2->soft_right_clip_seqs.size() > 0);
		s = fr.h2->soft_right_clip_seqs[0];
		printf("s size %lu\n",s.size());

		for(int i=0;i<=s.size()-kmer_length;i++)
		{
			string kmer = s.substr(i,kmer_length);

			if(kmer_map.find(kmer) == kmer_map.end())
			{
				kmer_map.insert(pair<string,int>(kmer,1));
			}
		}

		//discard if seq match with multiple junction
		int32_t prev_pos1 = 0;
		jc_multiple = 0;

		//match right hit soft seq with junctions seq
		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];

			if(jc.rpos >= fr.h2->pos || jc.rpos >= fr.h1->pos) continue;
			if(abs(fr.h1->pos-jc.rpos) > max_read_to_junction_gap) continue;

			int32_t pos1 = jc.rpos;
			int32_t pos2 = jc.rpos+soft_len-1;

			string junc_seq = get_fasta_seq(pos1,pos2);

			int edit_match = 0;
			assert(junc_seq.size() == fr.h2->soft_right_clip_seqs[0].size());

			bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
			if(is_similar == false) continue;

			int edit = get_edit_distance(junc_seq,fr.h2->soft_right_clip_seqs[0]);
			
			if(edit <= floor(soft_len/10))
			{
				edit_match = 1;
				printf("RO recover editmatch case 2: junc seq pos1=%d, pos2=%d, junc_seqlen = %lu\n",pos1,pos2,junc_seq.size());
			}

			if(edit_match == 1)
			{
				if(prev_pos1 != 0 && pos1 != prev_pos1)
				{
					printf("pos1 = %d, prev_pos1 = %d\n",pos1,prev_pos1);
					jc_multiple = 1;
					break;
				}

				prev_pos1 = pos1;
			}
		}

		if(jc_multiple == 1) 
		{
			printf("RO recover jc_multiple = %d\n",jc_multiple);
			continue;
		}

		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];

			if(jc.rpos >= fr.h2->pos || jc.rpos >= fr.h1->pos) continue;
			if(abs(fr.h1->pos-jc.rpos) > max_read_to_junction_gap) continue;

			int32_t pos1 = jc.rpos;
			int32_t pos2 = jc.rpos+soft_len-1;

			string junc_seq = get_fasta_seq(pos1,pos2);

			assert(junc_seq.size() == fr.h2->soft_right_clip_seqs[0].size());

			bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
			if(is_similar == false) continue;

			int edit = get_edit_distance(junc_seq,fr.h2->soft_right_clip_seqs[0]);
			
			//if(edit == 0 || edit == 1)
			if(edit <= floor(soft_len/10))
			{
				printf("both sides S clip right: chrm=%s, read=%s, read_pos=%d\n",bb.chrm.c_str(),fr.h2->qname.c_str(),fr.h2->pos);
				if((fr.h2->flag & 0x10) >= 1)
				{
					//printf("rev comp 0x10 = on\n");
				}
				else
				{
					//printf("rev comp 0x10 = off\n");
				}
				if((fr.h2->flag & 0x4) >= 1)
				{
					//printf("seg unmapped 0x4 = on\n");
				}
				else
				{
					//printf("seg unmapped 0x4 = off\n");
				}
				printf("read seq combo index=0, combo_seq=%s, edit=%d\n",fr.h2->soft_right_clip_seqs[0].c_str(),edit);
				printf("junction lpos = %d, rpos = %d\n",jc.lpos,jc.rpos);
				printf("junc seq pos1=%d, pos2=%d, junc_seqlen = %lu, junc_seq=%s\n",pos1,pos2,junc_seq.size(),junc_seq.c_str());

				seq_match_right_hit = 1;
				break;
			}
		}

		if(seq_match_left_hit == 1 && seq_match_right_hit == 1)
		{
			fr.HS_frag = true;
			circ_fragments.push_back(fr);
		}
	}

	return 0;
}

int bundle_bridge::get_RO_frags_with_HS()
{
	int junc_range = 10;

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment fr = fragments[k];

		//has a supple
		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			continue;
		}

		//is a supple
		if((fr.h1->flag & 0x800) >= 1 || (fr.h2->flag & 0x800) >= 1) continue;

		//not a RO frag
		if(fr.h1->is_reverse_overlap == false && fr.h2->is_reverse_overlap == false)
		{
			continue;
		}

		int left_boundary_flag = 0;
		int right_boundary_flag = 0;

		if(fr.h1->pos <= fr.h2->pos && (fr.h1->cigar_vector[0].first == 'S' || fr.h1->cigar_vector[0].first == 'H') && (fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'S' || fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'H'))
		{
			printf("RO paired hit case 1: pos %d, rpos %d\n",fr.lpos,fr.rpos);
			printf("chrm %s\n",bb.chrm.c_str());
			printf("Hit 1: ");
			fr.h1->print();

			//checking if reads junction matches left boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				//if(jc.rpos <= fr.h1->pos+junc_range && jc.rpos >= fr.h1->pos-junc_range)
				if(jc.rpos == fr.h1->pos)
				{
					printf("reads junction present left %d\n",jc.rpos);
					left_boundary_flag = 1;
					break;
				}
			}

			//checking if ref junction matches left boundary
			int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					//if(chain[k].second <= fr.h1->pos+junc_range && chain[k].second >= fr.h1->pos-junc_range)
					if(chain[k].second == fr.h1->pos)
					{
						printf("ref junction present left %d\n",chain[k].second);
						temp_flag = 1;
						left_boundary_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}

			printf("Hit 2: ");
			fr.h2->print();

			//checking if reads junction matches right boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				//if(jc.lpos <= fr.h2->rpos+junc_range && jc.lpos >= fr.h2->rpos-junc_range)
				if(jc.lpos == fr.h2->rpos)
				{
					printf("reads junction present right %d\n",jc.lpos);
					right_boundary_flag = 1;
					break;
				}
			}

			//checking if ref junction matches right boundary
			temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					//if(chain[k].first <= fr.h2->rpos+junc_range && chain[k].first >= fr.h2->rpos-junc_range)
					if(chain[k].first == fr.h2->rpos)
					{
						printf("ref junction present right %d\n",chain[k].first);
						right_boundary_flag = 1;
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}

			printf("RO frags: left_boundary_flag = %d, right_boundary_flag = %d\n\n",left_boundary_flag,right_boundary_flag);
		
		}

		/*else if(fr.h1->pos > fr.h2->pos && (fr.h2->cigar_vector[0].first == 'S' || fr.h2->cigar_vector[0].first == 'H') && (fr.h1->cigar_vector[fr.h1->cigar_vector.size()-1].first == 'S' || fr.h1->cigar_vector[fr.h1->cigar_vector.size()-1].first == 'H'))
		{
			printf("RO paired hit case 2: pos %d, rpos %d\n",fr.lpos,fr.rpos);
			printf("chrm %s\n",bb.chrm.c_str());
			printf("Hit 1: ");
			fr.h2->print();

			//checking if reads junction matches left boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				//if(jc.rpos <= fr.h2->pos+junc_range && jc.rpos >= fr.h2->pos-junc_range)
				if(jc.rpos == fr.h2->pos)
				{
					printf("reads junction present left %d\n",jc.rpos);
					left_boundary_flag = 1;
					break;
				}
			}

			//checking if ref junction matches left boundary
			int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					//if(chain[k].second <= fr.h2->pos+junc_range && chain[k].second >= fr.h2->pos-junc_range)
					if(chain[k].second == fr.h2->pos)
					{
						printf("ref junction present left %d\n",chain[k].second);
						left_boundary_flag = 1;
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}

			printf("Hit 2: ");
			fr.h1->print();

			//checking if reads junction matches right boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				//if(jc.lpos <= fr.h1->rpos+junc_range && jc.lpos >= fr.h1->rpos-junc_range)
				if(jc.lpos == fr.h1->rpos)
				{
					printf("reads junction present right %d\n",jc.lpos);
					right_boundary_flag = 1;
					break;
				}
			}

			//checking if ref junction matches right boundary
			temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					//if(chain[k].first <= fr.h1->rpos+junc_range && chain[k].first >= fr.h1->rpos-junc_range)
					if(chain[k].first == fr.h1->rpos)
					{
						printf("ref junction present right %d\n",chain[k].first);
						right_boundary_flag = 1;
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}

			printf("RO frags: left_boundary_flag = %d, right_boundary_flag = %d\n\n",left_boundary_flag,right_boundary_flag);
			
		}*/

		if(left_boundary_flag == 1 && right_boundary_flag == 1)
		{
			circ_fragments.push_back(fr);
		}
	}

	return 0;
}

string bundle_bridge::get_fasta_seq(int32_t pos1, int32_t pos2)
{
	string out = "";
	if(fai != NULL)
	{
		//printf("extracting fasta seq from region:\n");
		int32_t seqlen;
		char* seq = faidx_fetch_seq(fai, bb.chrm.c_str(), pos1, pos2, &seqlen);
		if(seq != NULL && seqlen > 0)
		{
			//printf("seqlen = %d, seq = %s\n",seqlen,seq);
			for(int i=0;i<seqlen;i++)
			{
				out = out + seq[i];
			}
		}
	}
	return out;
}


int bundle_bridge::min(int x, int y, int z) { return std::min(std::min(x, y), z); }

int bundle_bridge::get_edit_distance(string s, string t)
{
	int n = s.size();
	int m =t.size();
	int** dp = new int*[n+1];
	for(int i=0;i<=n;i++)
	{
		dp[i] = new int[m+1];
		for(int j=0;j<=m;j++)
		{
			dp[i][j]=0;
			if(i==0)dp[i][j]=j;
			else if(j==0)dp[i][j] = i;
		}
	}
	s = " " + s;
	t = " " + t;
	for(int i =1;i<=n;i++)
	{
		for(int j=1;j<=m;j++)
		{
			if(s[i] !=t[j])
			{
				dp[i][j] = 1+min(dp[i-1][j],dp[i][j-1],dp[i-1][j-1]);
			}
			else
			{
				dp[i][j] = dp[i-1][j-1];
			}
		}
	}

	int res = dp[n][m];
	for(int i=0;i<=n;i++)
	{
		delete [] dp[i];
	}
	delete [] dp;

	return res;
}

double bundle_bridge::get_Jaccard(int kmer_length, map<string,int> kmer_map, string t)
{
	int match_count = 0;

	for(int i=0;i<=t.size()-kmer_length;i++)
	{
		string kmer = t.substr(i,kmer_length);

		if(kmer_map.find(kmer) != kmer_map.end())
		{
			match_count++;
		}
	}

	int kmer_number = t.size()-kmer_length+1;
	double jaccard = 0;
	jaccard = (double) match_count/ (double) (kmer_number + kmer_number - match_count);
	//printf("match count = %d, kmer_number = %d, jaccard %lf\n",match_count,kmer_number,jaccard);
	return jaccard;
	
}

bool bundle_bridge::are_strings_similar(int kmer_length, map<string,int> kmer_map, string t)
{
	int match_count = 0;
	//printf("t size = %lu\n",t.size());

	for(int i=0;i<=t.size()-kmer_length;i++)
	{
		string kmer = t.substr(i,kmer_length);

		if(kmer_map.find(kmer) != kmer_map.end())
		{
			match_count++;
		}
	}
	//printf("end of similar\n");

	if(match_count >= floor(t.size()/15)) return true;

	return false;
}

int bundle_bridge::get_more_chimeric()
{
	int max_read_to_junction_gap = 100000;
	int min_soft_clip_len = 15;

	map<string, pair<int32_t, int32_t>> left_soft; //key:pos and seq, val junc pos pair
	map<string, pair<int32_t, int32_t>> right_soft;

	left_soft.clear();
	right_soft.clear();
	
	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		//has a supple
		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			continue;
		}

		//is a supple
		if((fr.h1->flag & 0x800) >= 1 || (fr.h2->flag & 0x800) >= 1) continue;

		//is a RO frag
		/*if(fr.h1->is_reverse_overlap == true || fr.h2->is_reverse_overlap == true)
		{
			continue;
		}*/

		//if h2 pos is to the left of h1 pos
		if(fr.h1->pos > fr.h2->pos) continue;

		//if none of the sides have soft clip, continue
		if(fr.h1->cigar_vector[0].first != 'S' && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first != 'S') continue;

		int soft_clip_side = 0; //if left soft kept, this will be 1, if right soft kept, this will be 2

		//check which side has soft clip and set flag accordingly, if both sides have soft clip, keep larger
		if(fr.h1->cigar_vector[0].first == 'S' && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first != 'S')
		{
			soft_clip_side = 1;
		}
		else if(fr.h1->cigar_vector[0].first != 'S' && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'S')
		{
			soft_clip_side = 2;
		}
		else if(fr.h1->cigar_vector[0].first == 'S' && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].first == 'S')
		{
			int32_t left_len = fr.h1->cigar_vector[0].second;
			int32_t right_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;

			if(left_len > right_len) soft_clip_side = 1;
			else if(left_len < right_len) soft_clip_side = 2;
			else //if both side len same, skip for now
			{
				if(left_len >= min_soft_clip_len)
				{
					printf("both side soft clips same length in new chimeric, chrm=%s, fr.qname=%s, fr.h1.pos=%d\n",bb.chrm.c_str(),fr.h1->qname.c_str(),fr.h1->pos);
				}
				continue;
			}
		}

		//discard if soft len < min_soft_clip_len
		if(soft_clip_side == 1 && fr.h1->cigar_vector[0].second < min_soft_clip_len) continue;
		if(soft_clip_side == 2 && fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second < min_soft_clip_len) continue;

		bool exists = false;

		//check if already exists in left_soft
		if(soft_clip_side == 1)
		{
			int32_t soft_len = fr.h1->cigar_vector[0].second + fr.h1->tiny_boundary;
			assert(fr.h1->soft_left_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h1->cigar_vector[0].second; i<soft_len;i++)
			{
				if(i<0 || i>fr.h1->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h1->cigar_vector[0].second;
					break;
				}
				tiny = tiny + fr.h1->seq[i];
			}

			string s = fr.h1->soft_left_clip_seqs[0] + tiny;
			printf("soft_len:%d\n",soft_len);
			printf("tiny:%s\n clip:%s, s=%s \n seq=%s\n",tiny.c_str(),fr.h1->soft_left_clip_seqs[0].c_str(),s.c_str(),fr.h1->seq.c_str());

			string hash = to_string(fr.h1->pos) + "|" + s;
			if(left_soft.find(hash) != left_soft.end())
			{
				if(left_soft[hash].first != -1 && left_soft[hash].second != -1)
				{
					//present in map and junc exitsts, create supple
					printf("soft left clip:exists, fake created\n");
					create_fake_supple(k,fr,soft_len,left_soft[hash].first,left_soft[hash].second,soft_clip_side);
				}
				exists = true;
			}
		}
		//check if already exists in right_soft
		else if(soft_clip_side == 2)
		{
			int32_t soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second + fr.h2->tiny_boundary;
			assert(fr.h2->soft_right_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h2->seq.size()-soft_len; i<fr.h2->seq.size()-fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;i++)
			{
				if(i<0 || i>fr.h2->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;
					break;
				}
				tiny = tiny + fr.h2->seq[i];
			}

			string s = tiny + fr.h2->soft_right_clip_seqs[0];
			printf("soft_len:%d\n",soft_len);
			printf("tiny:%s\n clip:%s, s=%s \n seq=%s\n",tiny.c_str(),fr.h2->soft_right_clip_seqs[0].c_str(),s.c_str(),fr.h2->seq.c_str());

			string hash = to_string(fr.h2->rpos) + "|" + s;
			if(right_soft.find(hash) != right_soft.end())
			{
				if(right_soft[hash].first != -1 && right_soft[hash].second != -1)
				{	
					//present in map and junc exists, create supple
					printf("soft right clip:exists, fake created\n");
					create_fake_supple(k,fr,soft_len,right_soft[hash].first,right_soft[hash].second,soft_clip_side);
				}
				exists = true;
			}
		}

		if(exists == true) continue;

		int left_boundary_match = 0;
		int right_boundary_match = 0;

		int pexon_range = 5;
		
		if(soft_clip_side == 1)
		{
			//check if soft clip end matches a read junction
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.rpos == fr.h1->pos)
				{
					//printf("jc.rpos = %d\n",jc.rpos);
					left_boundary_match = 1;
					break;
				}
			}

			//check if soft clip end matches a ref junction
			/*int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int p=0;p<chain.size();p++)
				{
					// assert(chain[p].first < chain[p].second);
					// if(chain[p].first <= bb.lpos) continue;
					// if(chain[p].second >= bb.rpos) continue;

					if(chain[p].second == fr.h1->pos)
					{
						left_boundary_match = 1;
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/

			//check if soft clip end matches pexon boundary of ref not given
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					//if(pexons[p].lpos <= fr.h1->pos+pexon_range && pexons[p].lpos >= fr.h1->pos-pexon_range && pexons[p].ltype == START_BOUNDARY)
					if(pexons[p].lpos == fr.h1->pos && pexons[p].ltype == START_BOUNDARY)
					{
						left_boundary_match = 1;
						break;
					}
				}
			}
		}
		else if(soft_clip_side == 2)
		{
			//check if soft clip end matches a read junction
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.lpos == fr.h2->rpos)
				{
					//printf("jc.rpos = %d\n",jc.rpos);
					right_boundary_match = 1;
					break;
				}
			}

			//check if soft clip end matches a ref junction
			/*int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int p=0;p<chain.size();p++)
				{
					// assert(chain[p].first < chain[p].second);
					// if(chain[p].first <= bb.lpos) continue;
					// if(chain[p].second >= bb.rpos) continue;

					if(chain[p].first == fr.h2->rpos)
					{
						right_boundary_match = 1;
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/

			//check if soft clip end matches pexon boundary of ref not given
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					//if(pexons[p].rpos <= fr.h2->rpos+pexon_range && pexons[p].rpos >= fr.h2->rpos-pexon_range && pexons[p].rtype == END_BOUNDARY)
					if(pexons[p].rpos == fr.h2->rpos && pexons[p].rtype == END_BOUNDARY)
					{
						right_boundary_match = 1;
						break;
					}
				}
			}
		}

		if(left_boundary_match == 0 && right_boundary_match == 0) continue;

		if(soft_clip_side == 1 && left_boundary_match == 0) //add to map that this frag left soft clip is invalid
		{
			int32_t soft_len = fr.h1->cigar_vector[0].second + fr.h1->tiny_boundary;
			assert(fr.h1->soft_left_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h1->cigar_vector[0].second; i<soft_len;i++)
			{
				if(i<0 || i>fr.h1->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h1->cigar_vector[0].second;
					break;
				}
				tiny = tiny + fr.h1->seq[i];
			}

			string s = fr.h1->soft_left_clip_seqs[0] + tiny;
			string hash = to_string(fr.h1->pos) + "|" + s;
			left_soft[hash] = pair<int32_t,int32_t> (-1,-1);
		}
		else if(soft_clip_side == 2 && right_boundary_match == 0) //add to map that this frag right soft clip is invalid
		{
			int32_t soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second + fr.h2->tiny_boundary;
			assert(fr.h2->soft_right_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h2->seq.size()-soft_len; i<fr.h2->seq.size()-fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;i++)
			{
				if(i<0 || i>fr.h2->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;
					break;
				}
				tiny = tiny + fr.h2->seq[i];
			}

			string s = tiny + fr.h2->soft_right_clip_seqs[0];
			string hash = to_string(fr.h2->rpos) + "|" + s;
			right_soft[hash] = pair<int32_t,int32_t> (-1,-1);
		}

		//printf("more chimeric instances:\n");
		if(soft_clip_side == 1  && left_boundary_match == 1)
		{
			int32_t soft_len = fr.h1->cigar_vector[0].second + fr.h1->tiny_boundary;

			//create a hash map for the kmers in soft clip region and pass that for similarity testing
			int kmer_length = 10;
			map<string,int> kmer_map;
			kmer_map.clear();

			assert(fr.h1->soft_left_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h1->cigar_vector[0].second; i<soft_len;i++)
			{
				if(i<0 || i>fr.h1->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h1->cigar_vector[0].second;
					break;
				}
				tiny = tiny + fr.h1->seq[i];
			}

			string s = fr.h1->soft_left_clip_seqs[0] + tiny;
			printf("soft_len:%d\n",soft_len);
			printf("tiny:%s\n clip:%s, s=%s \n seq=%s\n",tiny.c_str(),fr.h1->soft_left_clip_seqs[0].c_str(),s.c_str(),fr.h1->seq.c_str());

			string hash = to_string(fr.h1->pos) + "|" + s;

			for(int i=0;i<=s.size()-kmer_length;i++)
			{
				string kmer = s.substr(i,kmer_length);

				if(kmer_map.find(kmer) == kmer_map.end())
				{
					kmer_map.insert(pair<string,int>(kmer,1));
				}
			}

			//discard if seq match with multiple junction
			int32_t prev_pos2 = 0;
			int jc_multiple = 0;

			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.lpos <= fr.h2->rpos || jc.lpos <= fr.h1->rpos) continue;
				if(abs(jc.lpos-fr.h2->rpos) > max_read_to_junction_gap) continue;

				int32_t pos1 = jc.lpos-soft_len+1;
				int32_t pos2 = jc.lpos;

				string junc_seq = get_fasta_seq(pos1,pos2);

				int edit_match = 0;

				assert(junc_seq.size() == s.size());

				//bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
				//if(is_similar == false) continue;

				//int edit = get_edit_distance(junc_seq,fr.h1->soft_left_clip_seqs[0]);

				double similarity = get_Jaccard(kmer_length,kmer_map,junc_seq);
				/*if(strcmp(fr.h1->qname.c_str(),"ST-E00299:245:HKTJJALXX:6:2107:30563:57952")==0)
				{
					printf("ST-E00299:245:HKTJJALXX:6:2107:30563:57952 soft len=%d, edit=%d, sim==%lf\n",soft_len,edit,similarity);
				}*/
				
				if(similarity > 0.6)
				//if(edit <= floor(soft_len/10))
				{
					edit_match = 1;
					printf("editmatch case 1: junc seq pos1=%d, pos2=%d, junc_seqlen = %lu\n",pos1,pos2,junc_seq.size());
				}
				
				if(edit_match == 1)
				{
					if(prev_pos2 != 0 && pos2 != prev_pos2)
					{
						printf("pos2 = %d, prev_pos2 = %d\n",pos2,prev_pos2);
						jc_multiple = 1;
						break;
					}

					prev_pos2 = pos2;
				}
			}

			if(jc_multiple == 1) 
			{
				printf("jc_multiple = %d, frag_name=%s\n",jc_multiple,fr.h1->qname.c_str());
				left_soft[hash] = pair<int32_t,int32_t> (-1,-1);
				continue;
			}

			int jc_flag = 0;
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.lpos <= fr.h2->rpos || jc.lpos <= fr.h1->rpos) continue;
				if(abs(jc.lpos-fr.h2->rpos) > max_read_to_junction_gap) continue;

				int32_t pos1 = jc.lpos-soft_len+1;
				int32_t pos2 = jc.lpos;

				string junc_seq = get_fasta_seq(pos1,pos2);

				assert(junc_seq.size() == s.size());

				//bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
				//if(is_similar == false) continue;

				//int edit = get_edit_distance(junc_seq,fr.h1->soft_left_clip_seqs[0]);

				double similarity = get_Jaccard(kmer_length,kmer_map,junc_seq);
				
				if(similarity > 0.6)
				//if(edit <= floor(soft_len/10))
				{
					printf("soft left clip: combo index=0, chrm=%s, read=%s, read_pos=%d, combo_seq=%s, similarity=%lf\n",bb.chrm.c_str(),fr.h1->qname.c_str(),fr.h1->pos,s.c_str(),similarity);
					printf("junction lpos = %d, rpos = %d\n",jc.lpos,jc.rpos);
					printf("junc seq pos1=%d, pos2=%d, junc_seqlen = %lu, junc_seq=%s\n",pos1,pos2,junc_seq.size(),junc_seq.c_str());
					create_fake_supple(k,fr,soft_len,pos1,pos2,soft_clip_side);
					jc_flag = 1;
					left_soft[hash] = pair<int32_t,int32_t> (pos1,pos2);
					break;
				}
			}
		}
		else if(soft_clip_side == 2 && right_boundary_match == 1)
		{
			int32_t soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second + fr.h2->tiny_boundary;

			//create a hash map for the kmers in soft clip region and pass that for similarity testing
			int kmer_length = 10;
			map<string,int> kmer_map;
			kmer_map.clear();

			assert(fr.h2->soft_right_clip_seqs.size() > 0);

			//extract tiny seq
			string tiny = "";
			for(int i=fr.h2->seq.size()-soft_len; i<fr.h2->seq.size()-fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;i++)
			{
				if(i<0 || i>fr.h2->seq.size()-1)
				{
					tiny = "";
					soft_len = fr.h2->cigar_vector[fr.h2->cigar_vector.size()-1].second;
					break;
				}
				tiny = tiny + fr.h2->seq[i];
			}

			string s = tiny + fr.h2->soft_right_clip_seqs[0];
			printf("soft_len:%d\n",soft_len);
			printf("tiny:%s\n clip:%s, s=%s \n seq=%s\n",tiny.c_str(),fr.h2->soft_right_clip_seqs[0].c_str(),s.c_str(),fr.h2->seq.c_str());

			string hash = to_string(fr.h2->rpos) + "|" + s;

			for(int i=0;i<=s.size()-kmer_length;i++)
			{
				string kmer = s.substr(i,kmer_length);

				if(kmer_map.find(kmer) == kmer_map.end())
				{
					kmer_map.insert(pair<string,int>(kmer,1));
				}
			}

			//discard if seq match with multiple junction
			int32_t prev_pos1 = 0;
			int jc_multiple = 0;

			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.rpos >= fr.h2->pos || jc.rpos >= fr.h1->pos) continue;
				if(abs(fr.h1->pos-jc.rpos) > max_read_to_junction_gap) continue;

				int32_t pos1 = jc.rpos;
				int32_t pos2 = jc.rpos+soft_len-1;

				string junc_seq = get_fasta_seq(pos1,pos2);

				int edit_match = 0;

				assert(junc_seq.size() == s.size());

				//bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
				//if(is_similar == false) continue;
				
				//int edit = get_edit_distance(junc_seq,fr.h2->soft_right_clip_seqs[0]);

				double similarity = get_Jaccard(kmer_length,kmer_map,junc_seq);

				if(similarity > 0.6)
				//if(edit <= floor(soft_len/10))
				{
					edit_match = 1;
					printf("editmatch case 2: junc seq pos1=%d, pos2=%d, junc_seqlen = %lu\n",pos1,pos2,junc_seq.size());
				}
				
				if(edit_match == 1)
				{
					if(prev_pos1 != 0 && pos1 != prev_pos1)
					{
						printf("pos1 = %d, prev_pos1 = %d\n",pos1,prev_pos1);
						jc_multiple = 1;
						break;
					}

					prev_pos1 = pos1;
				}
			}

			if(jc_multiple == 1) 
			{
				printf("jc_multiple = %d, frag_name=%s\n",jc_multiple,fr.h1->qname.c_str());
				right_soft[hash] = pair<int32_t,int32_t> (-1,-1);
				continue;
			}

			int jc_flag = 0;
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.rpos >= fr.h2->pos || jc.rpos >= fr.h1->pos) continue;
				if(abs(fr.h1->pos-jc.rpos) > max_read_to_junction_gap) continue;

				int32_t pos1 = jc.rpos;
				int32_t pos2 = jc.rpos+soft_len-1;

				string junc_seq = get_fasta_seq(pos1,pos2);

				assert(junc_seq.size() == s.size());

				//bool is_similar = are_strings_similar(kmer_length,kmer_map,junc_seq);
				//if(is_similar == false) continue;
				
				//int edit = get_edit_distance(junc_seq,fr.h2->soft_right_clip_seqs[0]);

				double similarity = get_Jaccard(kmer_length,kmer_map,junc_seq);

				if(similarity > 0.6)
				//if(edit <= floor(soft_len/10))
				{
					//printf("read seq combo index=%d, combo_seq=%s, edit=%d\n",i,fr.h2->soft_right_clip_seqs[i].c_str(),edit);
					printf("soft right clip: combo index=0, chrm=%s, read=%s, read_pos=%d, combo_seq=%s, similarity=%lf\n",bb.chrm.c_str(),fr.h2->qname.c_str(),fr.h2->pos,s.c_str(),similarity);
					printf("junction lpos = %d, rpos = %d\n",jc.lpos,jc.rpos);
					printf("junc seq pos1=%d, pos2=%d, junc_seqlen = %lu, junc_seq=%s\n",pos1,pos2,junc_seq.size(),junc_seq.c_str());
					create_fake_supple(k,fr,soft_len,pos1,pos2,soft_clip_side);
					jc_flag = 1;
					right_soft[hash] = pair<int32_t,int32_t> (pos1,pos2);
					break;
				}
			}
		}
	}
	return 0;
}

int bundle_bridge::create_fake_supple(int fr_index, fragment &fr, int32_t soft_len, int32_t pos1, int32_t pos2, int soft_clip_side)
{
	hit new_hit;

	if(bb.fake_hits.size() == 0 && bb.hits.size() > 0)
	{
		new_hit.hid = bb.hits[bb.hits.size()-1].hid + 1;
	}
	else if(bb.fake_hits.size() > 0)
	{
		new_hit.hid = bb.fake_hits[bb.fake_hits.size()-1].hid + 1;
	}

	//hit id

	new_hit.is_fake = true;
	new_hit.soft_clip_side = soft_clip_side;
	new_hit.tid = fr.h1->tid;
	new_hit.qname = fr.h1->qname;
	new_hit.qhash = fr.h1->qhash;

	new_hit.pos = pos1;
	new_hit.rpos = pos2;
	new_hit.n_cigar = 1;
	new_hit.qlen = soft_len;

	new_hit.cigar_vector.push_back(pair<char, int32_t>('M',soft_len)); //100M
	new_hit.itvm.push_back(pack(pos1, pos2)); //check

	//check vlist 

	//printf("Printing fake hit:\n");
	//new_hit.print();

	if(soft_clip_side == 1)
	{
		//printf("cond1\n");
		new_hit.mpos = fr.h1->mpos;
		new_hit.mtid = fr.h1->mtid;
		new_hit.isize = -(pos2-fr.h1->pos+1); //isize includes soft clip?
		new_hit.strand = fr.h1->strand;
	}
	else if(soft_clip_side == 2)
	{
		//printf("cond2\n");
		new_hit.mpos = fr.h2->mpos;
		new_hit.mtid = fr.h2->mtid;
		new_hit.isize = (fr.h2->rpos-pos1+1); //isize includes soft clip?
		new_hit.strand = fr.h2->strand;
	}

	//set fake_hit_index for fr and new_hit
	fr.fake_hit_index = bb.fake_hits.size();
	new_hit.fake_hit_index = fr_index;

	//push to bb.fake_hits
	bb.fake_hits.push_back(new_hit);

	return 0;
}


int bundle_bridge::build_supplementaries()
{

	int max_index = bb.hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

    vector< vector<int> > vv;
    vv.resize(max_index);

    //printf("Bundle hit size: %d\n", bb.hits.size());
    // first build index
    for(int i = 0; i < bb.hits.size(); i++)
    {
        hit &h = bb.hits[i];

		if(strcmp(h.qname.c_str(),"simulate:432512") == 0)
		{
			printf("Found read simulate:432512\n");
		}

        // TODO
        if((h.flag & 0x800) == 0) continue;

        //printf("%s\n",h.qname.c_str());
        //printf("%zu\n",h.qhash);

        // do not use hi; as long as qname, pos and isize are identical
        // add 0x40 and 0x80
        int k = (h.qhash % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;

		if(strcmp(h.qname.c_str(),"simulate:432512") == 0)
		{
			printf("Found read supple simulate:432512, k=%d\n",k);
		}
        vv[k].push_back(i);
        //printf("Adding supple\n");
    }

    //printf("End of vv adding\n");

    for(int i = 0; i < bb.hits.size(); i++)
    {
        hit &h = bb.hits[i];

        //if(h.paired == true) continue; 
        //if(h.isize <= 0) continue; //commented out as this was filtering non chimeric part of an end of SRR1721290.17627808
        //if(h.vlist.size() == 0) continue;
        if((h.flag & 0x800) >= 1) continue;       // skip supplemetary

        int k = (h.qhash % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;

        for(int j = 0; j < vv[k].size(); j++)
        {
            hit &z = bb.hits[vv[k][j]];

			if(strcmp(h.qname.c_str(),"simulate:432512") == 0)
			{
				printf("Found read z simulate:432512, k=%d\n",k);
				z.print();
				printf("Read h corr to read z:\n");
				h.print();
			}
			
            //if(z.hi != h.hi) continue;
            //if(z.paired == true) continue;
            //if(z.pos != h.mpos) continue;
            //if(z.isize + h.isize != 0) continue;
            //if(z.qhash != h.qhash) continue;
            if(z.qname != h.qname) continue;

            // TODO check 0x40 and 0x80 are the same
            if(((z.flag & 0x40) != (h.flag & 0x40)) || ((z.flag & 0x80) != (h.flag & 0x80))) continue;

        	h.suppl = &z;
        	break;

            //Taking the first supplementary read
        }
    }

	//printf("End of bundle bridge extract supple.\n");
    return 0;
}

int bundle_bridge:: set_chimeric_cigar_positions()
{
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.suppl == NULL) continue;

		/*printf("Primary hit:\n");
		h.print();
		printf("Supple hit:\n");
		h.suppl->print();


		printf("cigar prim:");
		for(int p=0;p<h.n_cigar;p++)
		{
			printf("%d%c",h.cigar_vector[p].second,h.cigar_vector[p].first);
		}
		printf("\n");

		printf("cigar supp:");
		for(int p=0;p<h.suppl->n_cigar;p++)
		{
			printf("%d%c",h.suppl->cigar_vector[p].second,h.suppl->cigar_vector[p].first);
		}
		printf("\n");*/

		
		int32_t p;
		int32_t q;

		p = h.pos;
		int match_index = 0;;

		for(int j=0;j<h.n_cigar;j++)
		{
			if(h.cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			p -= h.cigar_vector[j].second; //subtracting cigars before match index
		}
		
		q = h.suppl->pos;
		match_index = 0;;

		for(int j=0;j<h.suppl->n_cigar;j++)
		{
			if(h.suppl->cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			q -= h.suppl->cigar_vector[j].second; //subtracting cigars before match index
		}
		

		int32_t x = p;
		int32_t diff_cigar1 = 10000;
		int32_t diff_cigar2 = 10000;
		int best_pos_flag = 0;

		for(int j=0;j<h.n_cigar-1;j++)
		{
			int32_t y = q;
			pair<char, int32_t> hp_cigar1 = h.cigar_vector[j];
			pair<char, int32_t> hp_cigar2 = h.cigar_vector[j+1];

			for(int k=0;k<h.suppl->n_cigar-1;k++)
			{

				pair<char, int32_t> hs_cigar1 = h.suppl->cigar_vector[k];
				pair<char, int32_t> hs_cigar2 = h.suppl->cigar_vector[k+1];

				//printf("check %d,%c\n",hp_cigar1.second,hp_cigar1.first);

				int32_t hp_cigar1_len = hp_cigar1.second;
				int32_t hp_cigar2_len = hp_cigar2.second;
				int32_t hs_cigar1_len = hs_cigar1.second;
				int32_t hs_cigar2_len = hs_cigar2.second;

				if(((hp_cigar1.first == 'S' || hp_cigar1.first == 'H') && hp_cigar2.first == 'M') && (hs_cigar1.first == 'M' && (hs_cigar2.first == 'S' || hs_cigar2.first == 'H')))
				{
					if(hp_cigar2.first == 'M' && j+3 < h.cigar_vector.size() && h.cigar_vector[j+3].first == 'M')
					{
						//&& h.cigar_vector[j+2].first == 'N' can be I or D as well
						for(int p=j+3;p<h.cigar_vector.size();p+=2) //traverse all Ms with 1 gap in the middle ex:SMNMNM/SIMIMIM/SMDMDMD
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar2_len += h.cigar_vector[p].second;
							}
						}			
					}

					if(hs_cigar1.first == 'M' && k-2 >= 0 && h.suppl->cigar_vector[k-2].first == 'M')
					{
						//&& h.suppl->cigar_vector[k-1].first == 'N' can be I or D as well
						for(int p=k-2;p>=0;p-=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar1_len += h.suppl->cigar_vector[p].second;
							}
						}
					}
					
					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);
				
				}
				else if(((hs_cigar1.first == 'S' || hs_cigar1.first == 'H') && hs_cigar2.first == 'M') && (hp_cigar1.first == 'M' && (hp_cigar2.first == 'S' || hp_cigar2.first == 'H')))
				{
					if(hs_cigar2.first == 'M' && k+3 < h.suppl->cigar_vector.size() && h.suppl->cigar_vector[k+3].first == 'M')
					{
						//&& h.suppl->cigar_vector[k+2].first == 'N' can be I or D as well
						for(int p=k+3;p<h.suppl->cigar_vector.size();p+=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar2_len += h.suppl->cigar_vector[p].second;
							}
						}			
					}

					if(hp_cigar1.first == 'M' && j-2 >=0 && h.cigar_vector[j-2].first == 'M')
					{
						//&& h.cigar_vector[j-1].first == 'N' can be I or D as well
						for(int p=j-2;p>=0;p-=2)
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar1_len += h.cigar_vector[p].second;
							}
						}
					}

					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);					
				}

				if(diff_cigar1 < 20 && diff_cigar2 < 20) //setting diff max 20 between complementing cigars
				{
					//printf("Found best positions\n");
					best_pos_flag = 1;

					//set first,second and third positions of prim and suppl
					h.first_pos = x;
					h.second_pos = x + hp_cigar1.second;
					h.third_pos = x + hp_cigar1.second + hp_cigar2.second;

					h.suppl->first_pos = y;
					h.suppl->second_pos = y + hs_cigar1.second;
					h.suppl->third_pos = y + hs_cigar1.second + hs_cigar2.second;

					h.left_cigar = hp_cigar1.first;
					h.left_cigar_len = hp_cigar1_len;
					h.right_cigar = hp_cigar2.first;
					h.right_cigar_len = hp_cigar2_len;

					h.suppl->left_cigar = hs_cigar1.first;
					h.suppl->left_cigar_len = hs_cigar1_len;
					h.suppl->right_cigar = hs_cigar2.first;
					h.suppl->right_cigar_len = hs_cigar2_len;

					break;

				}
				y += hs_cigar1.second;

			}
			if(best_pos_flag == 1) break;

			x += hp_cigar1.second;
		}

		//printf("set_cigar p:%d-%d-%d\n",h.first_pos,h.second_pos,h.third_pos);
		//printf("set_cigar s:%d-%d-%d\n",h.suppl->first_pos,h.suppl->second_pos,h.suppl->third_pos);
	}

	return 0;
}


int bundle_bridge::build_junctions()
{
	int min_max_boundary_quality = min_mapping_quality; //building a list of all splice pos and the hit index that includes the splice pos
	map< int64_t, vector<int> > m; // map of spos against vector of hits indices
	for(int i = 0; i < bb.hits.size(); i++)
	{
		vector<int64_t> v = bb.hits[i].spos;
		//printf("Spos size: %d\n", v.size());

		if(v.size() == 0) continue;

		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i); 
			}
		}
	}
	//printf("spos map size = %d\n",m.size());

	junctions.clear();
	junc_map.clear();
	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		int32_t p1 = high32(it->first);
		int32_t p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bb.hits[v[k]];
			nm += h.nm;
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", bb.chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(it->first, v.size());
		jc.nm = nm;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);

		if(junc_map.find(it->first) != junc_map.end()) continue;
		else junc_map.insert(make_pair(it->first, jc.strand));
	}

	//printf("old Junctions size: %lu\n", junctions.size());

	// do some filtering here?
	// let M be the maximum count among all junctions

	int max_count = 0;
	for(int j=0;j<junctions.size();j++)
	{
		junction jc = junctions[j];
		if(jc.count > max_count)
		{
			max_count = junctions[j].count;
		}
	}

	//printf("junction max_count = %d\n",max_count);

	// keep a junction J if: 
	// either J.count >= ratio * M, say ratio = 0.01, 
	// or J.count >= a fixed threshold, say 10 

	double ratio = 0.05;
	filtered_junctions.clear();

	for(int j=0;j<junctions.size();j++)
	{
		junction jc = junctions[j];
		
		if(jc.count >= ratio*max_count || jc.count >= 10) //discard if both < 10 and < 0.01*max_count
		{
			filtered_junctions.push_back(jc);
		}
	}

	junctions.clear();
	for(int j=0;j<filtered_junctions.size();j++)
	{
		junction jc = filtered_junctions[j];
		junctions.push_back(jc);
	}

	/*if(junctions.size() != filtered_junctions.size())
	{
		printf("junctions size: %lu\n",junctions.size());
		printf("filtered junctions size: %lu\n",filtered_junctions.size());
	}*/

	//printf("new Junctions size: %lu\n", junctions.size());

	return 0;
}

int bundle_bridge::extend_junctions()
{
	map< int64_t, vector<int> > m;
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		vector<PI32> v = ref_trsts[i].get_intron_chain();
		for(int k = 0; k < v.size(); k++)
		{
			assert(v[k].first < v[k].second);
			// TODO, TODO
			//if(v[k].first <= bb.lpos) continue;
			//if(v[k].second >= bb.rpos) continue;
			int64_t p = pack(v[k].first, v[k].second);

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		for(int k = 0; k < v.size(); k++)
		{
			char c = ref_trsts[v[k]].strand;
			if(c == '.') s0++;
			if(c == '+') s1++;
			if(c == '-') s2++;
		}

		junction jc(it->first, 0 - v.size());
		jc.nm = 0;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);
		filtered_junctions.push_back(jc);

		// triggle to split intervals
		//printf("EXTEND junction: %d-%d\n", low32(it->first), high32(it->first));
		bb.mmap += make_pair(ROI(high32(it->first), high32(it->first) + 1), 1);
		bb.mmap += make_pair(ROI(high32(it->first), high32(it->first) + 1), -1);
		bb.mmap += make_pair(ROI(low32(it->first), low32(it->first) + 1), 1);
		bb.mmap += make_pair(ROI(low32(it->first), low32(it->first) + 1), -1);
	}

	//for(int i=0;i<junctions.size();i++)
	//{

	//	junctions[i].print();
	//	printf("\n");
	//}

	return 0;
}


int bundle_bridge::build_regions()
{
	MPI s;
	s.insert(PI(bb.lpos, START_BOUNDARY));
	s.insert(PI(bb.rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	for(int k = 0; k < v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second; 
		int rtype = v[k + 1].second; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		//region rr(l, r, ltype, rtype);
		region rr(l, r, ltype, rtype, &(bb.mmap), &(bb.imap));
		evaluate_rectangle(bb.mmap, l, r, rr.ave, rr.dev, rr.max);
		regions.push_back(rr);
	}

	return 0;
}

int bundle_bridge::build_partial_exons()
{
	pexons.clear();
	//regional.clear();
	//printf("regions size = %lu\n",regions.size());
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		//printf("r.pexons size = %lu\n",r.pexons.size());
		for(int k = 0; k < r.pexons.size(); k++)
		{
			partial_exon &pe = r.pexons[k];
			pe.rid = i;
			pe.pid = pexons.size();
			//if(pe.rpos == 12248702) printf("rpos is 12248702\n");
			pexons.push_back(pe);
			//if((pe.lpos != bb.lpos || pe.rpos != bb.rpos) && pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) regional.push_back(true);
			//else regional.push_back(false);
		}
	}

	for(int i=0;i<pexons.size();i++)
	{
		//pexons[i].print(i+1);
	}
	return 0;
}

int bundle_bridge::align_hits_transcripts()
{
	map<int32_t, int> m;
	for(int k = 0; k < regions.size(); k++)
	{
		if(k >= 1) assert(regions[k - 1].rpos == regions[k].lpos);
		m.insert(pair<int32_t, int>(regions[k].lpos, k));
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		align_hit(m, bb.hits[i], bb.hits[i].vlist);
		bb.hits[i].vlist = encode_vlist(bb.hits[i].vlist);
		//remove_tiny_boundary(bb.hits[i]);
	}

	ref_phase.resize(ref_trsts.size());
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		align_transcript(m, ref_trsts[i], ref_phase[i]);
	}

	return 0;
}

int bundle_bridge::align_fake_hits()
{
	map<int32_t, int> m;
	for(int k = 0; k < regions.size(); k++)
	{
		if(k >= 1) assert(regions[k - 1].rpos == regions[k].lpos);
		m.insert(pair<int32_t, int>(regions[k].lpos, k));
	}

	for(int i = 0; i < bb.fake_hits.size(); i++)
	{
		align_hit(m, bb.fake_hits[i], bb.fake_hits[i].vlist);
		bb.fake_hits[i].vlist = encode_vlist(bb.fake_hits[i].vlist);
	}

	/*if(bb.fake_hits.size() > 0)
	{
		printf("Fake hits print\n");
		for(int i = 0; i < bb.fake_hits.size(); i++)
		{
			bb.fake_hits[i].print();
		}
	}*/

	return 0;
}

int bundle_bridge::create_fake_fragments()
{
	/*for(int k=0; k<fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		if(fr.fake_hit_index == -1) continue;

		for(int i=0;i<bb.fake_hits.size();i++)
		{
			hit &z = bb.fake_hits[i];

			if(z.vlist.size() == 0) continue;
			if(z.fake_hit_index == -1) continue;
			if(fr.h1->qname != z.qname) continue;

			if(fr.fake_hit_index == i && z.fake_hit_index == k)
			{
				if(z.soft_clip_side == 1)
				{
					printf("create_fake_frag: soft_clip_side 1\n");
					fr.h1->suppl = &z;
					fr.h1->supple_pos = z.pos;
					fr.h1->suppl->paired = true;
					fragment frag(fr.h2, fr.h1->suppl);
					frag.frag_type = 2;
					frag.is_compatible = 1;
					frag.type = 0;
					frag.pi = k;
					frag.fidx = circ_fragments.size();
					fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
					fr.fidx = k;
					circ_fragments.push_back(frag);
				}
				else if(z.soft_clip_side == 2)
				{
					printf("create_fake_frag: soft_clip_side 2\n");
					fr.h2->suppl = &z;
					fr.h2->supple_pos = z.pos;
					fr.h2->suppl->paired = true;
					fragment frag(fr.h2->suppl, fr.h1);
					frag.frag_type = 2;
					frag.is_compatible = 2;
					frag.type = 0;
					frag.pi = k;
					frag.fidx = circ_fragments.size();
					fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
					fr.fidx = k;
					circ_fragments.push_back(frag);	
				}
			}
		}
	}*/

	//making above extraction more efficient

	map<string, int> circ_map;
	for(int j=0;j<bb.fake_hits.size();j++)
	{
		hit &z = bb.fake_hits[j];
		if(circ_map.find(z.qname) == circ_map.end())
		{
			circ_map.insert(make_pair(z.qname, j));
		}
	}

	for(int k=0; k<fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		if(fr.fake_hit_index == -1) continue;

		if(circ_map.find(fr.h1->qname) == circ_map.end()) continue;

		int j = circ_map[fr.h1->qname];
		hit &z = bb.fake_hits[j];

		if(z.vlist.size() == 0) continue;
		if(z.fake_hit_index == -1) continue;

		assert(fr.h1->qname == z.qname);

		if(fr.fake_hit_index == j && z.fake_hit_index == k)
		{
			if(z.soft_clip_side == 1)
			{
				fr.h1->suppl = &z;
				fr.h1->supple_pos = z.pos;
				fr.h1->suppl->paired = true;
				fragment frag(fr.h2, fr.h1->suppl);
				frag.frag_type = 2;
				frag.is_compatible = 1;
				frag.type = 0;
				frag.pi = k;
				frag.fidx = circ_fragments.size();
				fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
				fr.fidx = k;
				circ_fragments.push_back(frag);
			}
			else if(z.soft_clip_side == 2)
			{
				fr.h2->suppl = &z;
				fr.h2->supple_pos = z.pos;
				fr.h2->suppl->paired = true;
				fragment frag(fr.h2->suppl, fr.h1);
				frag.frag_type = 2;
				frag.is_compatible = 2;
				frag.type = 0;
				frag.pi = k;
				frag.fidx = circ_fragments.size();
				fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
				fr.fidx = k;
				circ_fragments.push_back(frag);	
			}
		}
	}

	return 0;
}

int bundle_bridge::remove_tiny_boundaries_fake_hits()
{
	for(int i = 0; i < bb.fake_hits.size(); i++)
	{
		remove_tiny_boundary(bb.fake_hits[i]);
	}
	return 0;
}

int bundle_bridge::remove_tiny_boundaries()
{
	for(int i = 0; i < bb.hits.size(); i++)
	{
		remove_tiny_boundary(bb.hits[i]);
	}
	return 0;
}

int bundle_bridge::set_fragment_lengths()
{
	for(int k = 0; k < fragments.size(); k++)
	{
		set_fragment_length(fragments[k]);
	}
	return 0;
}

int bundle_bridge::set_circ_fragment_lengths()
{
	for(int k = 0; k < circ_fragments.size(); k++)
	{
		set_fragment_length(circ_fragments[k]);
	}
	return 0;
}

int bundle_bridge::set_fragment_length(fragment &fr)
{
	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	fr.lpos = fr.h1->pos;
	fr.rpos = fr.h2->rpos;
	
	if(fr.h1->vlist.size() % 2 != 0)
	{
		printf("debug vlist\n");
		printf("size %lu\n",fr.h1->vlist.size());	
		return 0;
	}
	if(fr.h2->vlist.size() % 2 != 0)
	{
		printf("debug vlist\n");
		printf("size %lu\n",fr.h2->vlist.size());
		return 0;
	}

	vector<int> v1 = decode_vlist(fr.h1->vlist);
	vector<int> v2 = decode_vlist(fr.h2->vlist);
	fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
	fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
	fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
	fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

	fr.b1 = true;
	if(v1.size() <= 1) 
	{
		fr.b1 = false;
	}
	else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
	{
		if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
	}
	else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
	{
		if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
	}

	fr.b2 = true;
	if(v2.size() <= 1)
	{
		fr.b2 = false;
	}
	else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
	{
		if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
	}
	else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
	{
		if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
	}

	return 0;
}

int bundle_bridge::align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &vv)
{
	vv.clear();
	vector<int64_t> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return 0;

	vector<PI> sp;
	sp.resize(v.size());

	int32_t p1 = high32(v.front());
	int32_t p2 = low32(v.back());

	sp[0].first = locate_region(p1);
	for(int k = 1; k < v.size(); k++)
	{
		p1 = high32(v[k]);

		map<int32_t, int>::const_iterator it = m.find(p1);

		//assert(it != m.end());
		if(it == m.end()) return 0;
		sp[k].first = it->second;
	}

	sp[sp.size() - 1].second = locate_region(p2 - 1);
	for(int k = 0; k < v.size() - 1; k++)
	{
		p2 = low32(v[k]);
		map<int32_t, int>::const_iterator it = m.find(p2);
		
		//assert(it != m.end());
		if(it == m.end()) return 0;
		sp[k].second = it->second - 1; 
	}

	for(int k = 0; k < sp.size(); k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > 0) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle_bridge::align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &vv)
{
	vv.clear();
	int k1 = -1;
	int k2 = -1;
	for(int k = 0; k < t.exons.size(); k++)
	{
		if(t.exons[k].second > bb.lpos)
		{
			k1 = k;
			break;
		}
	}
	for(int k = t.exons.size() - 1; k >= 0; k--)
	{
		if(t.exons[k].first < bb.rpos)
		{
			k2 = k;
			break;
		}
	}

	if(k1 > k2) return 0;
	if(k1 == -1 || k2 == -1) return 0;

	vector<PI> sp;
	sp.resize(k2 + 1);

	int32_t p1 = t.exons[k1].first > bb.lpos ? t.exons[k1].first : bb.lpos;
	int32_t p2 = t.exons[k2].second < bb.rpos ? t.exons[k2].second : bb.rpos;

	sp[k1].first = locate_region(p1);
	for(int k = k1 + 1; k <= k2; k++)
	{
		p1 = t.exons[k].first;
		map<int32_t, int>::const_iterator it = m.find(p1);
		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[k2].second = locate_region(p2 - 1);
	for(int k = k1; k < k2; k++)
	{
		p2 = t.exons[k].second;
		map<int32_t, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = k1; k <= k2; k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > k1) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle_bridge::remove_tiny_boundary(hit &h1)
{
	vector<int> v1 = decode_vlist(h1.vlist);
	int n1 = v1.size();
	if(n1 >= 2 && v1[n1 - 2] + 1 == v1[n1 - 1])
	{
		int k = v1[n1 - 1];
		int32_t total = regions[k].rpos - regions[k].lpos;
		int32_t flank = h1.rpos - regions[k].lpos;

		if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
		{
			vector<int> v(v1.begin(), v1.begin() + n1 - 1);
			assert(v.size() + 1 == v1.size());
			h1.vlist = encode_vlist(v);
			h1.tiny_boundary = abs(h1.rpos-regions[k].lpos);
			h1.rpos = regions[k].lpos;
		}
	}

	vector<int> v2 = decode_vlist(h1.vlist);
	int n2 = v2.size();
	if(n2 >= 2 && v2[0] + 1 == v2[1])
	{
		int k = v2[0];
		int32_t total = regions[k].rpos - regions[k].lpos;
		int32_t flank = regions[k].rpos - h1.pos;

		if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
		{
			vector<int> v(v2.begin() + 1, v2.end());
			assert(v.size() + 1 == v2.size());
			h1.vlist = encode_vlist(v);
			h1.tiny_boundary = abs(regions[k].rpos-h1.pos);
			h1.pos = regions[k].rpos;
		}
	}
	return 0;
}

int bundle_bridge::index_references()
{
	ref_index.clear();
	ref_index.resize(regions.size());
	for(int k = 0; k < ref_phase.size(); k++)
	{
		vector<int> &v = ref_phase[k];
		for(int j = 0; j < v.size(); j++)
		{
			int x = v[j];
			ref_index[x].push_back(PI(k, j));
		}
	}
	return 0;
}

int bundle_bridge::locate_region(int32_t x)
{
	if(regions.size() == 0) return -1;

	int k1 = 0;
	int k2 = regions.size();
	while(k1 < k2)
	{
		int m = (k1 + k2) / 2;
		region &r = regions[m];
		if(x >= r.lpos && x < r.rpos) return m;
		else if(x < r.lpos) k2 = m;
		else k1 = m;
	}
	return -1;
}


int bundle_bridge::build_fragments()
{

	int ctp = 0;// count fragments number from paired-end reads
	int ctu = 0;// count fragments number from UMI
	int ctb = 0;// count fragments number for both

	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	fragments.clear();
	if(bb.hits.size() == 0) return 0;

	int max_index = bb.hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index); //max_index slots initialized to zero, here max_index is the max hash index

	// first build index
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.isize >= 0) continue;
		if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.qhash % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;

		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i); //vv containes hits of the same hash
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];

		if(h.paired == true) continue;
		if(h.isize <= 0) continue;
		if(h.vlist.size() == 0) continue;

		int k = (h.qhash % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		/*if(strcmp(h.qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in hits\n");
			printf("isize: %d\n",h.isize);
			printf("vlist size: %zu\n",h.vlist.size());
		}*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];

			/*if(strcmp(z.qname.c_str(),"simulate:311116") == 0)
			{
				printf("simulate:311116 is in hits\n");
			}*/

			//if(z.hi != h.hi) continue;
			if(z.paired == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;

			x = vv[k][j];
			break;
		}

		/*
		SI si(h.qname, h.hi);
		MSI::iterator it = vv[k].find(si);
		if(it == vv[k].end()) continue;
		int x = it->second;
		*/

		//printf("HIT: i = %d, x = %d, bb.hits[i].vlist = %lu | ", i, x, bb.hits[i].vlist.size(), bb.hits[i].qname.c_str()); bb.hits[i].print();

		if(x == -1) continue;
		if(bb.hits[x].vlist.size() == 0) continue;

		fragment fr(&bb.hits[i], &bb.hits[x]); //h2 and h1s as param or h2s and h1 as parameter

		/*if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in fragments\n");
		}*/

		fr.frag_type = 1; //this is the first set of fragment,tasfia

		//keep it
		// ===============================
		// TODO: dit for UMI
		bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index
		bb.hits[x].pi = i;
		bb.hits[i].fidx = fragments.size();//check if used somewhere
		bb.hits[x].fidx = fragments.size();//fidx is the fragment index
		ctp += 1;
		fr.type = 0; 

		/*
		// shao: moved to a separate function
		// and do it later (after remove-tiny-boundary)

		fr.lpos = h.pos;
		fr.rpos = bb.hits[x].rpos;

		vector<int> v1 = decode_vlist(bb.hits[i].vlist);
		vector<int> v2 = decode_vlist(bb.hits[x].vlist);
		fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
		fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
		fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
		fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;
		//keep it

		//inlcude
		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
		}
		*/

		fragments.push_back(fr);

		bb.hits[i].paired = true;
		bb.hits[x].paired = true;

	}

	//printf("total bb.hits = %lu, total fragments = %lu\n", bb.hits.size(), fragments.size());
	
	/*for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in build_fragments\n");
		}
	}*/

	// by shao, exit here (no UMI)
	return 0;
}

int bundle_bridge::fix_alignment_boundaries()
{
	int end_error = 5;
	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		if(fr.h1->suppl != NULL)
		{
			hit *h1_supple = fr.h1->suppl;

			/*printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h1 has a supple hit.\n");
			printf("Primary: ");
			fr.h1->print();
			printf("Supple: ");
			h1_supple->print();
			printf("fr.h2: ");
			fr.h2->print();*/

			if(fr.h1->pos > fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos >= fr.h2->rpos && fr.h1->pos - fr.h2->pos <= end_error)
			{
				string combo = "alignment error-h1p_pos>h2_pos";
				//printf("%s\n",combo.c_str());

				int32_t diff = fr.h1->pos - fr.h2->pos;

				int first_M_len = 0;
				for(int i=0;i<fr.h2->cigar_vector.size();i++)
				{
					if(fr.h2->cigar_vector[i].first == 'M')
					{
						first_M_len = fr.h2->cigar_vector[i].second;
						break;
					}
				}

				if(first_M_len > end_error)
				{
					fr.h2->pos = fr.h2->pos + diff;
				}
			}
			else if(fr.h1->pos <= fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos < fr.h2->rpos && fr.h2->rpos - h1_supple->rpos <= end_error)
			{
				string combo = "alignment error-h2_rpos>h1s_rpos";
				//printf("%s\n",combo.c_str());

				int32_t diff = fr.h2->rpos - h1_supple->rpos;

				int last_M_len = 0;
				for(int i=0;i<fr.h2->cigar_vector.size();i++)
				{
					if(fr.h2->cigar_vector[i].first == 'M')
					{
						last_M_len = fr.h2->cigar_vector[i].second;
					}
				}

				if(last_M_len > end_error)
				{
					fr.h2->rpos = fr.h2->rpos - diff;
				}
			}
		}

		if(fr.h2->suppl != NULL)
		{
			hit *h2_supple = fr.h2->suppl;
			
			/*printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h2 has a supple hit.\n");
			printf("h1: ");
			fr.h1->print();
			printf("Primary: ");
			fr.h2->print();
			printf("Supple: ");
			h2_supple->print();*/

			if(h2_supple->pos <= fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos < fr.h1->rpos && fr.h1->rpos - fr.h2->rpos <= end_error)
			{
				string combo = "alignment error-h1_rpos>h2p_rpos";
				//printf("%s\n",combo.c_str());

				int32_t diff = fr.h1->rpos - fr.h2->rpos;

				int last_M_len = 0;
				for(int i=0;i<fr.h1->cigar_vector.size();i++)
				{
					if(fr.h1->cigar_vector[i].first == 'M')
					{
						last_M_len = fr.h1->cigar_vector[i].second;
					}
				}

				if(last_M_len > end_error)
				{
					fr.h1->rpos = fr.h1->rpos - diff;
				}
			}
			else if(h2_supple->pos > fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos >= fr.h1->rpos && h2_supple->pos - fr.h1->pos <= end_error)
			{
				string combo = "alignment error-h2s_pos>h1_pos";
				//printf("%s\n",combo.c_str());

				int32_t diff = h2_supple->pos - fr.h1->pos;

				int first_M_len = 0;
				for(int i=0;i<fr.h1->cigar_vector.size();i++)
				{
					if(fr.h1->cigar_vector[i].first == 'M')
					{
						first_M_len = fr.h1->cigar_vector[i].second;
						break;
					}
				}

				if(first_M_len > end_error)
				{
					fr.h1->pos = fr.h1->pos + diff;
				}
			}
		}
	}

	return 0;
}

int bundle_bridge::build_circ_fragments()
{
	//parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	//vector<fragment> circ_fragments;
	//circ_fragments.clear(); // unbridged second set of fragments, cannot use this for finding partner fragments

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		/*if(strcmp(fr.h1->qname.c_str(),"simulate:311116") == 0)
		{
			printf("simulate:311116 is in build_circ_fragments\n");
		}*/

		int is_compatible = 0; //1 for h1 has a suppl and compatible, 2 for h2 has a suppl and compatible

		if(fr.h1->suppl != NULL)
		{
			//printf("supple not null\n");
			//need to check compatibility from bundle.cc

			h1_supp_count++;
			hit *h1_supple = fr.h1->suppl;

			printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h1 has a supple hit.\n");
			printf("Primary: ");
			fr.h1->print();
			printf("Supple: ");
			h1_supple->print();
			printf("fr.h2: ");
			fr.h2->print();

			printf("cigar prim:");
			for(int p=0;p<fr.h1->n_cigar;p++)
			{
				printf("%d%c",fr.h1->cigar_vector[p].second,fr.h1->cigar_vector[p].first);
			}
			printf("\n");

			printf("cigar supp:");
			for(int p=0;p<fr.h1->suppl->n_cigar;p++)
			{
				printf("%d%c",fr.h1->suppl->cigar_vector[p].second,fr.h1->suppl->cigar_vector[p].first);
			}
			printf("\n");

			printf("set_cigar p:%d-%d-%d\n",fr.h1->first_pos,fr.h1->second_pos,fr.h1->third_pos);
			printf("set_cigar s:%d-%d-%d\n",fr.h1->suppl->first_pos,fr.h1->suppl->second_pos,fr.h1->suppl->third_pos);

			if(fr.h1->first_pos == 0 || fr.h1->suppl->first_pos == 0)
			{
				string combo = "special case h1s";
				printf("%s\n",combo.c_str());
				if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
				else frag2graph_freq[combo] += 1;
				continue;
			}

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(fr.h1->left_cigar == 'H' || fr.h1->left_cigar == 'S')
			{
				len_HS += fr.h1->left_cigar_len;
			}
			else if(fr.h1->right_cigar == 'H' || fr.h1->right_cigar == 'S')
			{
				len_HS += fr.h1->right_cigar_len;
			}
			if(fr.h1->suppl->left_cigar == 'H' || fr.h1->suppl->left_cigar == 'S')
			{
				len_HS += fr.h1->suppl->left_cigar_len;
			}
			else if(fr.h1->suppl->right_cigar == 'H' || fr.h1->suppl->right_cigar == 'S')
			{
				len_HS += fr.h1->suppl->right_cigar_len;
			}

			//printf("len_HS = %d\n",len_HS);

			if(abs(len_HS - read_length) > 5) //here 100 is the estimated read length, replace this with any related exisiting parameter
			{
				printf("read length criteria unsatisfied h1s.\n");
				continue;
			}


			//examples give a rule - p and s has to be edges

			//this case should not occur as h1p should always be on the left of h2
			if(h1_supple->pos <= fr.h2->pos && h1_supple->pos <= fr.h1->pos && fr.h1->rpos >= h1_supple->rpos && fr.h1->rpos >= fr.h2->rpos)
			{	
				//printf("Compatible in previous definition, h1s leftmost, h1p rightmost\n");
				
				if(h1_supple->second_pos <= fr.h2->pos && h1_supple->second_pos <= fr.h1->pos && fr.h1->second_pos >= h1_supple->rpos && fr.h1->second_pos >= fr.h2->rpos)
				{
					string combo = "Compatible h1s leftmost h1p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;				
				}
				else
				{
					string combo = "Not compatible h1s leftmost h1p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;				
				}

			}

			else if(fr.h1->pos <= fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos >= fr.h2->rpos)
			{
				//printf("Compatible in previous definition, h1p leftmost, h1s rightmost\n");

				if(fr.h1->second_pos <= fr.h2->pos && fr.h1->second_pos <= h1_supple->pos && h1_supple->second_pos >= fr.h1->rpos && h1_supple->second_pos >= fr.h2->rpos)
				{
					string combo = "Compatible h1p leftmost h1s rightmost";
					printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					is_compatible = 1;
				}
				else
				{
					string combo = "Not compatible h1p leftmost h1s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;	
					continue;				
				}

			}
			else
			{
				//printf("Not compatible in previous definition\n");
				if(fr.h1->pos > fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos >= fr.h2->rpos && fr.h1->pos - fr.h2->pos <= 5)
				{
					string combo = "alignment error-h1p_pos>h2_pos";
					//printf("%s\n",combo.c_str());
					//fr.print(0);

					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else if(fr.h1->pos <= fr.h2->pos && fr.h1->pos <= h1_supple->pos && h1_supple->rpos >= fr.h1->rpos && h1_supple->rpos < fr.h2->rpos && fr.h2->rpos - h1_supple->rpos <= 5)
				{
					string combo = "alignment error-h2_rpos>h1s_rpos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else
				{
					string combo = "Not compatible h1s";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;
				}
			}
		}

		if(fr.h2->suppl != NULL)
		{
			//need to check compatibility from bundle.

			h2_supp_count++;
			hit *h2_supple = fr.h2->suppl;
			
			printf("\nchrm = %s\n",bb.chrm.c_str());
			printf("fr.h2 has a supple hit.\n");
			printf("h1: ");
			fr.h1->print();
			printf("Primary: ");
			fr.h2->print();
			printf("Supple: ");
			h2_supple->print();

			printf("cigar prim:");
			for(int p=0;p<fr.h2->n_cigar;p++)
			{
				printf("%d%c",fr.h2->cigar_vector[p].second,fr.h2->cigar_vector[p].first);
			}
			printf("\n");

			printf("cigar supp:");
			for(int p=0;p<fr.h2->suppl->n_cigar;p++)
			{
				printf("%d%c",fr.h2->suppl->cigar_vector[p].second,fr.h2->suppl->cigar_vector[p].first);
			}
			printf("\n");

			printf("set_cigar p:%d-%d-%d\n",fr.h2->first_pos,fr.h2->second_pos,fr.h2->third_pos);
			printf("set_cigar s:%d-%d-%d\n",fr.h2->suppl->first_pos,fr.h2->suppl->second_pos,fr.h2->suppl->third_pos);


			if(fr.h2->first_pos == 0 || fr.h2->suppl->first_pos == 0)
			{

				string combo = "special case h2s";
				printf("%s\n",combo.c_str());
				if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
				else frag2graph_freq[combo] += 1;
				continue;
			}

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(fr.h2->left_cigar == 'H' || fr.h2->left_cigar == 'S')
			{
				len_HS += fr.h2->left_cigar_len;
			}
			else if(fr.h2->right_cigar == 'H' || fr.h2->right_cigar == 'S')
			{
				len_HS += fr.h2->right_cigar_len;
			}
			if(fr.h2->suppl->left_cigar == 'H' || fr.h2->suppl->left_cigar == 'S')
			{
				len_HS += fr.h2->suppl->left_cigar_len;
			}
			else if(fr.h2->suppl->right_cigar == 'H' || fr.h2->suppl->right_cigar == 'S')
			{
				len_HS += fr.h2->suppl->right_cigar_len;
			}

			//printf("len_HS = %d\n",len_HS);

			if(abs(len_HS - read_length) > 5) //here 100 is the estimated read length, replace this with any related exisiting parameter
			{
				printf("read length criteria unsatisfied h2s.\n");
				continue;
			}

			//general rule p and s has to be edges
			if(h2_supple->pos <= fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos >= fr.h1->rpos)
			{
				//printf("Compatible in previous definition, h2s leftmost, h2p rightmost\n");

				if(h2_supple->second_pos <= fr.h1->pos && h2_supple->second_pos <= fr.h2->pos && fr.h2->second_pos >= h2_supple->rpos && fr.h2->second_pos >= fr.h1->rpos)
				{
					string combo = "Compatible h2s leftmost h2p rightmost";
					printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					is_compatible = 2;					
				}
				else
				{
					string combo = "Not compatible h2s leftmost h2p rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;			
				}

			}

			//this case should not occur as h2 p should always be on the right of h1
			else if(fr.h2->pos <= fr.h1->pos && fr.h2->pos <= h2_supple->pos && h2_supple->rpos >= fr.h2->rpos && h2_supple->rpos >= fr.h1->rpos)
			{
				//printf("Compatible in previous definition, h2p leftmost, h2s rightmost\n");

				if(fr.h2->second_pos <= fr.h1->pos && fr.h2->second_pos <= h2_supple->pos && h2_supple->second_pos >= fr.h2->rpos && h2_supple->second_pos >= fr.h1->rpos)
				{
					string combo = "Compatible h2p leftmost h2s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;					
				}
				else
				{
					string combo = "Not compatible h2p leftmost h2s rightmost";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;	
					continue;				
				}

			}

			else
			{
				if(h2_supple->pos <= fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos < fr.h1->rpos && fr.h1->rpos - fr.h2->rpos <= 5)
				{
					string combo = "alignment error-h1_rpos>h2p_rpos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else if(h2_supple->pos > fr.h1->pos && h2_supple->pos <= fr.h2->pos && fr.h2->rpos >= h2_supple->rpos && fr.h2->rpos >= fr.h1->rpos && h2_supple->pos - fr.h1->pos <= 5)
				{
					string combo = "alignment error-h2s_pos>h1_pos";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
				}
				else
				{
					//printf("Not compatible in previous definition\n");

					string combo = "Not compatible h2s";
					//printf("%s\n",combo.c_str());
					if(frag2graph_freq.find(combo) == frag2graph_freq.end()) frag2graph_freq.insert(pair<string, int>(combo, 1));
					else frag2graph_freq[combo] += 1;
					continue;
				}
			}

			//if not compatible, continued
			
			//if compatible
			//fragment fr(fr.h2->suppl, fr.h1); //h2 and h1s as param or h2s and h1 as parameter
			//fr.frag_type = 2; //this is the second set of fragment

			//do we need to check h2s-h1-h2p order if compatible checked? no as ordering already checked, see above
			//need to handle fr.h1 paired - both fr.h1 and h2 needs to paired before when we build second fragment, see below
		}

		//if compatible h1s
		if(is_compatible == 1)
		{
			//if(bb.hits[x].vlist.size() == 0) continue;
			//printf("vlist size = %zu",fr.h2->vlist.size());
			//printf("vlist size = %zu",fr.h1->suppl->vlist.size());

			//fr.h2 and fr.h1 needs to be paired by build_fragments()
			if(fr.h2->paired != true || fr.h1->paired != true) continue; //first set of fragment needs to be paired
			if(fr.h1->suppl->vlist.size() == 0) continue;

			fr.h1->suppl->paired = true; //setting supple paired true to avoid assertion later in build hyper set
			fragment frag(fr.h2, fr.h1->suppl);
			frag.frag_type = 2;
			frag.is_compatible = 1;

			/*// ===============================
			// TODO: dit for UMI
			bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index, not used here
			bb.hits[x].pi = i;
			bb.hits[i].fidx = fragments.size();//check if used somewhere
			bb.hits[x].fidx = fragments.size();//fidx is the fragment index
			ctp += 1;
			fr.type = 0; //needed here?
			// ================================*/

			frag.type = 0;
			
			/*
			frag.lpos = fr.h2->pos;
			frag.rpos = fr.h1->suppl->rpos;

			vector<int> v1 = decode_vlist(fr.h2->vlist);
			vector<int> v2 = decode_vlist(fr.h1->suppl->vlist);
			frag.k1l = frag.h1->pos - regions[v1.front()].lpos;
			frag.k1r = regions[v1.back()].rpos - frag.h1->rpos; //frag.h1?
			frag.k2l = frag.h2->pos - regions[v2.front()].lpos;
			frag.k2r = regions[v2.back()].rpos - frag.h2->rpos; //frag.h2?
			//keep it

			//inlcude
			frag.b1 = true;
			if(v1.size() <= 1) 
			{
				frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + frag.h1->nm) frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + frag.h1->nm) frag.b1 = false;
			}

			frag.b2 = true;
			if(v2.size() <= 1)
			{
				frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment1 + frag.h2->nm) frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment2 + frag.h2->nm) frag.b2 = false;
			}
			*/

			frag.pi = k;
			frag.fidx = circ_fragments.size();
			fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
			fr.fidx = k;

			circ_fragments.push_back(frag);
			
			//fr.h2->paired = true;
			//fr.h1->suppl->paired = true;

		}
		else if(is_compatible == 2)
		{

			//fr.h2 and fr.h1 needs to be paired by build_fragments()
			if(fr.h2->paired != true || fr.h1->paired != true) //first set of fragment needs to be paired
			{
				printf("fr.h2 in first set not paired\n");
				fr.h2->print();
				continue;
			}
			if(fr.h2->suppl->vlist.size() == 0) continue;

			fr.h2->suppl->paired = true; //setting supple paired true to avoid assertion later in build hyper set
			fragment frag(fr.h2->suppl, fr.h1);
			frag.frag_type = 2;
			frag.is_compatible = 2;

			/*// ===============================
			// TODO: dit for UMI
			bb.hits[i].pi = x; //index of other hit of a fragment stored, i-x are fragment pair indices, partner index, not used here
			bb.hits[x].pi = i;
			bb.hits[i].fidx = fragments.size();//check if used somewhere
			bb.hits[x].fidx = fragments.size();//fidx is the fragment index
			ctp += 1;
			// ================================*/

			frag.type = 0;

			/*
			frag.lpos = fr.h2->suppl->pos;
			frag.rpos = fr.h1->rpos;

			vector<int> v1 = decode_vlist(fr.h2->suppl->vlist);
			vector<int> v2 = decode_vlist(fr.h1->vlist);
			frag.k1l = frag.h1->pos - regions[v1.front()].lpos;
			frag.k1r = regions[v1.back()].rpos - frag.h1->rpos;
			frag.k2l = frag.h2->pos - regions[v2.front()].lpos;
			frag.k2r = regions[v2.back()].rpos - frag.h2->rpos;
			//keep it

			//inlcude
			frag.b1 = true;
			if(v1.size() <= 1) 
			{
				frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + frag.h1->nm) frag.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(frag.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + frag.h1->nm) frag.b1 = false;
			}

			frag.b2 = true;
			if(v2.size() <= 1)
			{
				frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment1 + frag.h2->nm) frag.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
			{
				if(regions[v2.front()].rpos - frag.h2->pos > max_misalignment2 + frag.h2->nm) frag.b2 = false;
			}
			*/

			frag.pi = k;
			frag.fidx = circ_fragments.size();
			fr.pi = circ_fragments.size(); //pi not set for all fragments, only those that have second frags, check pi before using if it is -1
			fr.fidx = k;

			circ_fragments.push_back(frag);	

		}
	}



	//adding the second part fragments to fragments vector
	if(circ_fragments.size() > 0)
	{
		printf("fragments vector size = %zu\n",fragments.size());		
	}
	
	//fragments.clear();
	/*for(int i=0;i<circ_fragments.size();i++)
	{
		fragments.push_back(circ_fragments[i]); //pushing the circ fragments at the end of the main fragments list, 
												//does this create any problem in the main scallop2 output??
												//no as this is a separate module just for circRNA
	}*/

	if(circ_fragments.size() > 0)
	{
		printf("circ fragment vector size = %zu\n",circ_fragments.size());
	}

	//check fragments order before bridging, same
	/*if(circ_fragments.size()>0)
	{
		printf("fragments before bridging:\n");
		for(int i=0;i<fragments.size();i++)
		{
			printf("%d-%d,",fragments[i].h1->hid,fragments[i].h2->hid);
		}

		printf("circ fragments:\n");

		for(int i=0;i<circ_fragments.size();i++)
		{
			printf("%d-%d,",circ_fragments[i].h1->hid,circ_fragments[i].h2->hid);
		}
	}*/

	return 0;
}

int bundle_bridge::extract_all_non_supple_HS_hits()
{
	vector<hit> left_boundary_hits;
	vector<hit> right_boundary_hits;

	left_boundary_hits.clear();
	right_boundary_hits.clear();

	for(int i = 0; i < bb.hits.size(); i++)
    {
		hit h = bb.hits[i];
        
        if((h.flag & 0x800) >= 1) continue; //is a supple hit
		if(h.suppl != NULL) continue; //has a supple hit

		if(h.cigar_vector[0].first == 'S' || h.cigar_vector[0].first == 'H')
		{
			left_boundary_hits.push_back(h);
		}

		if(h.cigar_vector[h.cigar_vector.size()-1].first == 'S' || h.cigar_vector[h.cigar_vector.size()-1].first == 'H')
		{
			right_boundary_hits.push_back(h);
		}
	}

	for(int i=0;i<left_boundary_hits.size();i++)
	{
		hit h1 = left_boundary_hits[i];
		for(int j=0;j<right_boundary_hits.size();j++)
		{
			hit h2 = right_boundary_hits[j];

			circular_transcript circ;
			circ.seqname = bb.chrm.c_str();
			circ.start = h1.pos;
			circ.end = h2.rpos;

			circ_trsts_HS.push_back(circ);
		}
	}

	//printf("size of circ_trsts_HS:%lu\n",circ_trsts_HS.size());
	return 0;
}

int bundle_bridge::extract_nonsupple_HS_hits()
{
	map <int32_t,pair<junction, vector<hit>>> junc_HS_map; //key junc lpos/rpos, stores junction and set of hits supporting the junction
	junc_HS_map.clear();

	for(int i = 0; i < bb.hits.size(); i++)
    {
		hit h = bb.hits[i];
        
        if((h.flag & 0x800) >= 1) continue; //is a supple hit
		if(h.suppl != NULL) continue; //has a supple hit

		for(int j=0;j<junctions.size();j++)
		{
			junction jc = junctions[j];
		
			if(h.cigar_vector[0].first == 'S' || h.cigar_vector[0].first == 'H')
			{
				if(h.pos >= jc.rpos - 4 && h.pos <= jc.rpos + 4)
				{
					if(junc_HS_map.find(jc.rpos) != junc_HS_map.end())
					{
						junc_HS_map[jc.rpos].second.push_back(h);
					}
					else
					{
						vector<hit> temp;
						temp.push_back(h);
						jc.boundary_match = 'R';
						junc_HS_map.insert(pair<int32_t,pair<junction,vector<hit>>>(jc.rpos,pair<junction,vector<hit>>(jc,temp)));
					}
				}
			}
			if(h.cigar_vector[h.cigar_vector.size()-1].first == 'S' || h.cigar_vector[h.cigar_vector.size()-1].first == 'H')
			{
				if(h.rpos >= jc.lpos - 4 && h.rpos <= jc.lpos + 4)
				{
					if(junc_HS_map.find(jc.lpos) != junc_HS_map.end())
					{
						junc_HS_map[jc.lpos].second.push_back(h);
					}
					else
					{
						vector<hit> temp;
						temp.push_back(h);
						jc.boundary_match = 'L';
						junc_HS_map.insert(pair<int32_t,pair<junction,vector<hit>>>(jc.lpos,pair<junction,vector<hit>>(jc,temp)));
					}
				}				
			}
		}
	}

	if(junc_HS_map.size() > 0)
	{
		printf("chrm=%s\n",bb.chrm.c_str());
		printf("junc_HS_map size = %lu\n",junc_HS_map.size());

		map<int32_t, pair<junction, vector<hit>>>::iterator itn;
		for(itn = junc_HS_map.begin(); itn != junc_HS_map.end(); itn++)
		{ 
			printf("key = %d, boundary = %c, count = %lu\n",itn->first,itn->second.first.boundary_match,itn->second.second.size());
			for(int i=0;i<itn->second.second.size();i++)
			{
				hit h = itn->second.second[i];
				printf("%s cigar:",h.qname.c_str());
				for(int j=0;j<h.cigar_vector.size();j++)
				{
					printf("%d%c",h.cigar_vector[j].second,h.cigar_vector[j].first);
				}
				printf("\n");
			}
		}
	}

	junction max_jc_lpos, max_jc_rpos;
	vector<hit> support_jc_lpos, support_jc_rpos;
	support_jc_lpos.clear();
	support_jc_rpos.clear();
	int cnt_lpos = 0, cnt_rpos = 0;

	//extract jc boundary lpos and rpos with support >= 5
	map<int32_t, pair<junction, vector<hit>>>::iterator itn;
	for(itn = junc_HS_map.begin(); itn != junc_HS_map.end(); itn++)
	{
		junction jc = itn->second.first;
		int support = itn->second.second.size();
		vector<hit> support_vec = itn->second.second;

		//printf("inside map check, boundary_match = %c\n",jc.boundary_match);

		if(jc.boundary_match == 'L')
		{
			if(support > cnt_lpos)
			{
				//printf("inside find max lpos\n");
				cnt_lpos = support;
				max_jc_lpos = jc;
				support_jc_lpos.clear();
				support_jc_lpos.insert(support_jc_lpos.begin(),support_vec.begin(),support_vec.end());
			}
		}

		if(jc.boundary_match == 'R')
		{
			if(support > cnt_rpos)
			{
				//printf("inside find max rpos\n");
				cnt_rpos = support;
				max_jc_rpos = jc;
				support_jc_rpos.clear();
				support_jc_rpos.insert(support_jc_rpos.begin(),support_vec.begin(),support_vec.end());
			}
		}
	}

	if(cnt_lpos >= 4 && cnt_rpos >= 4)
	{
		printf("chrm=%s\n",bb.chrm.c_str());
		printf("Left junction: %d, support: %d\n", max_jc_lpos.lpos, cnt_lpos); //left junc matching with read rpos
		printf("Right junction: %d, support: %d\n", max_jc_rpos.rpos, cnt_rpos); //right junc matching with read lpos

		bridger brdg(this);
		circular_transcript circ;

		brdg.bridge_clip(max_jc_rpos.rpos, max_jc_lpos.lpos, circ);

		if(circ.start != 0 && circ.end != 0)
		{
			string chrm_id = bb.chrm.c_str();
			string circRNA_id = "chrm" + chrm_id + ":" + tostring(circ.start) + "|" + tostring(circ.end) + "|";
			char strand = bb.strand;

			circ.circRNA_id = circRNA_id;
			circ.seqname = chrm_id;
			circ.source = "scallop2_HS";
			circ.feature = "circRNA";
			circ.gene_id = "gene";
			circ.strand = infer_circ_strand(circ.circ_path);

			for(int k=0;k<support_jc_rpos.size();k++)
			{
				circ.transcript_id = circ.transcript_id + support_jc_rpos[k].qname.c_str() + "|";
			}
			for(int k=0;k<support_jc_lpos.size();k++)
			{
				circ.transcript_id = circ.transcript_id + support_jc_lpos[k].qname.c_str() + "|";
			}

			for(int i=0;i<circ.merged_regions.size();i++)
			{
				region r = circ.merged_regions[i];
				circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
			}

			//circ.print(0);
			circ_trsts.push_back(circ);
		}
	}

	return 0;
}

int bundle_bridge::extract_HS_frags_circRNA()
{
	for(int j=0;j<circ_fragments.size();j++)
	{
		fragment &fr = circ_fragments[j];

		//has a supple
		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			continue;
		}

		//is a supple
		if((fr.h1->flag & 0x800) >= 1 || (fr.h2->flag & 0x800) >= 1) continue;

		if(fr.HS_frag == false) continue;

		if(fr.paths.size() != 1) continue;

		//printf("Some HS circRNA\n");

		vector<int> v = decode_vlist(fr.paths[0].v);

		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr.lpos) + "|" + tostring(fr.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr.lpos;
		int32_t end = fr.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), v.begin(), v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2_HS";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);
		
		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}

		circ_trsts.push_back(circ); 
	}

	return 0;
}

int bundle_bridge::extract_RO_circRNA()
{
	for(int j=0;j<circ_fragments.size();j++)
	{
		fragment &fr = circ_fragments[j];

		//has a supple
		if(fr.h1->suppl != NULL || fr.h2->suppl != NULL)
		{
			continue;
		}

		//is a supple
		if((fr.h1->flag & 0x800) >= 1 || (fr.h2->flag & 0x800) >= 1) continue;

		//not RO frag
		if(fr.h1->is_reverse_overlap == false && fr.h2->is_reverse_overlap == false)
		{
			continue;
		}

		if(fr.paths.size() != 1) continue;

		vector<int> v = decode_vlist(fr.paths[0].v);

		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr.lpos) + "|" + tostring(fr.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr.lpos;
		int32_t end = fr.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), v.begin(), v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2_RO";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);
		
		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}

		circ_trsts.push_back(circ); 
	}
	return 0;
}


int bundle_bridge::extract_circ_fragment_pairs()
{
	circ_fragment_pairs.clear();

	if(circ_fragments.size() > 0)
	{
		printf("After bridging fragments vector size = %zu\n\n",fragments.size());		
	}

	/*for(int i = 0; i < fragments.size(); i++)
	{
		fragment &fr1 = fragments[i]; 

		if(fr1.frag_type == 2) continue; //filter and take original fragments of scallop2 

		for(int j=0;j<circ_fragments.size();j++)
		{
			fragment &fr2 = circ_fragments[j];

			if(strcmp(fr1.h1->qname.c_str(),fr2.h1->qname.c_str()) != 0) continue;
			if(fr1.pi == -1 || fr2.pi == -1 || fr1.fidx == -1 || fr2.fidx == -1) continue;

			if(fr1.pi == fr2.fidx && fr2.pi == fr1.fidx)
			{
				//printf("Found partner fragment\n");
				circ_fragment_pairs.push_back(pair<fragment,fragment>(fr1,fr2)); //fr1 is original fragment of scallop2 regardless of which of h1 or h2 has a supplement
			}
		}
	}*/

	// above quadratic implementation seems slow down the program
	// an index might improve (see below)
	map<string, int> circ_map;
	for(int j=0;j<circ_fragments.size();j++)
	{
		fragment &fr2 = circ_fragments[j];
		if(circ_map.find(fr2.h1->qname) == circ_map.end())
		{
			circ_map.insert(make_pair(fr2.h1->qname, j));
		}
	}
	for(int i = 0; i < fragments.size(); i++)
	{
		fragment &fr1 = fragments[i]; 
		if(circ_map.find(fr1.h1->qname) == circ_map.end()) continue;

		int j = circ_map[fr1.h1->qname];
		// TODO: check frag, pi, fidx, etc
		// TODO: add pair to circ_fragment_pairs

		fragment &fr2 = circ_fragments[j]; 
		if(fr1.pi == -1 || fr2.pi == -1 || fr1.fidx == -1 || fr2.fidx == -1) continue;

		if(fr1.pi == fr2.fidx && fr2.pi == fr1.fidx)
		{
			//printf("Found partner fragment\n");
			circ_fragment_pairs.push_back(pair<fragment,fragment>(fr1,fr2)); //fr1 is original fragment of scallop2 regardless of which of h1 or h2 has a supplement
		}

	}

	//printf("Printing bridged fragment pairs: size = %zu\n\n",circ_fragment_pairs.size());

	/*if(circ_fragment_pairs.size() > 0)
	{
		printf("Printing bridged fragment pairs: size = %zu\n\n",circ_fragment_pairs.size());

		for(int i=0;i<circ_fragment_pairs.size();i++)
		{
			circ_fragment_pairs[i].first.print(i+1);
			circ_fragment_pairs[i].second.print(i+1);
			printf("\n");
		}
	}*/

	//profiling fragment hits bridged or not
	for(int i=0;i<circ_fragment_pairs.size();i++)
	{
		fragment &fr1 = circ_fragment_pairs[i].first;
		fragment &fr2 = circ_fragment_pairs[i].second;

		//printf("pairs: is_compatible = %d-%d, %c%c-%c%c\n",fr1.is_compatible, fr2.is_compatible,fr1.h1->bridged == true ? 'T' : 'F',fr1.h2->bridged == true ? 'T' : 'F',fr2.h1->bridged == true ? 'T' : 'F',fr2.h2->bridged == true ? 'T' : 'F');

		//printf("chrm = %s\n",bb.chrm.c_str());

		//bridging info for fr1 hits
		if(fr1.h1->bridged == false && fr1.h2->bridged == false)
		{
			string combo = "fr1-FF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == false && fr1.h2->bridged == true)
		{
			string combo = "fr1-FT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == true && fr1.h2->bridged == false)
		{
			string combo = "fr1-TF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr1.h1->bridged == true && fr1.h2->bridged == true)
		{
			string combo = "fr1-TT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}


		/*printf("Printing first fragment:\n");
		fr1.print(i+1);
		printf("Printing second fragment:\n");
		fr2.print(i+1);*/

		//bridging info for fr2 hits
		if(fr2.h1->bridged == false && fr2.h2->bridged == false)
		{
			string combo = "fr2-FF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == false && fr2.h2->bridged == true)
		{
			string combo = "fr2-FT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == true && fr2.h2->bridged == false)
		{
			string combo = "fr2-TF";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}
		else if(fr2.h1->bridged == true && fr2.h2->bridged == true)
		{
			string combo = "fr2-TT";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;
		}

		//printf("\n");

		//counting path.size = 1
		if(fr1.paths.size() == 1)
		{
			string combo = "fr1-path.size = 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		else
		{
			string combo = "fr1-path.size != 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		if(fr2.paths.size() == 1)
		{
			string combo = "fr2-path.size = 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;			
		}
		else
		{
			string combo = "fr2-path.size != 1";
			//printf("%s\n",combo.c_str());
			if(circ_frag_bridged_freq.find(combo) == circ_frag_bridged_freq.end()) circ_frag_bridged_freq.insert(pair<string, int>(combo, 1));
			else circ_frag_bridged_freq[combo] += 1;		
		}
	}

	return 0;
}

int bundle_bridge::print_circ_fragment_pairs()
{

	if(circ_fragment_pairs.size() > 0)
	{
		printf("Printing fragment pair paths:\n");
		for(int i=0;i<circ_fragment_pairs.size();i++)
		{
			fragment &fr1 = circ_fragment_pairs[i].first;
			fragment &fr2 = circ_fragment_pairs[i].second;

			if(fr1.paths.size() != 1 || fr2.paths.size() != 1) continue; //check both frags bridged
			//if(fr1.paths[0].type != 1 || fr2.paths[0].type != 1) continue; // 1: within normal range of insertsize;

			fr1.print(i+1);
			fr2.print(i+1);
			fr1.paths[0].print(i+1);
			fr2.paths[0].print(i+1);
		}
	}
	return 0;
}

int bundle_bridge::join_circ_fragment_pairs(int32_t length_high)
{
	for(int i=0;i<circ_fragment_pairs.size();i++)
	{
		fragment &fr1 = circ_fragment_pairs[i].first;
		fragment &fr2 = circ_fragment_pairs[i].second;

		printf("\nPrinting separate fragments: chrm = %s",bb.chrm.c_str());

		fr1.print(i+1);
		fr2.print(i+1);

		if(fr1.paths.size() != 1)
		{
			printf("Not valid: fr1 paths size = %lu, read_name=%s, hit1_pos=%d, hit1_rpos=%d, hit2_pos=%d, hit2_rpos=%d\n",fr1.paths.size(),fr1.h1->qname.c_str(),fr1.h1->pos,fr1.h1->rpos,fr1.h2->pos,fr1.h2->rpos); 
		}

		if(fr2.paths.size() != 1)
		{
			if((fr2.h1->flag & 0x800) >= 1 || (fr2.h2->flag & 0x800) >= 1)// print only not valid chimeric frags
			{
				printf("Not valid chimeric: fr2 paths size = %lu, read_name=%s, hit1_pos=%d, hit1_rpos=%d, hit2_pos=%d, hit2_rpos=%d\n",fr2.paths.size(),fr2.h1->qname.c_str(),fr2.h1->pos,fr2.h1->rpos,fr2.h2->pos,fr2.h2->rpos);
			}
			if(fr2.h1->is_fake == true || fr2.h2->is_fake == true)// print only not valid fake frags
			{
				printf("Not valid fake: fr2 paths size = %lu, read_name=%s, hit1_pos=%d, hit1_rpos=%d, hit2_pos=%d, hit2_rpos=%d\n",fr2.paths.size(),fr2.h1->qname.c_str(),fr2.h1->pos,fr2.h1->rpos,fr2.h2->pos,fr2.h2->rpos);
			}
		}

		if(fr1.paths.size() != 1 || fr2.paths.size() != 1)//not bridged
		{
			continue; 
		}
		
		if(fr1.paths[0].length > length_high && fr2.paths[0].length > length_high)
		{
			printf("Not valid: both fragment length > length high\n");
			printf("length_high = %d, fr1.path[0].length = %d, fr2.path[0].length = %d\n",length_high,fr1.paths[0].length,fr2.paths[0].length);
			continue;
		}

		int left_boundary_flag = 0;
		int right_boundary_flag = 0;
		int pexon_left_flag = 0;
		int pexon_right_flag = 0;
		int pexon_range = 5;
		int bundle_range = 5;
		int junc_range = 5;

		char junc_match = '.';
		char ref_match = '.';

		//printf("pexons size = %lu\n",pexons.size());

		if(fr2.is_compatible == 1)
		{

			//checking if reads junction matches left boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.rpos == fr1.lpos)
				//if(jc.rpos <= fr1.lpos+junc_range && jc.rpos >= fr1.lpos-junc_range)
				{
					//printf("jc.rpos = %d\n",jc.rpos);
					left_boundary_flag = 1;
					junc_match = 'L';
					break;
				}
			}

			//checking if ref junction matches left boundary
			/*int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					// assert(chain[k].first < chain[k].second);
					// if(chain[k].first <= bb.lpos) continue;
					// if(chain[k].second >= bb.rpos) continue;

					if(chain[k].second == fr1.lpos)
					//if(chain[k].second <= fr1.lpos+junc_range && chain[k].second >= fr1.lpos-junc_range)
					{
						left_boundary_flag = 1;
						ref_match = 'L';
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/

			//checking if pexon matches left boundary when ref file is not given
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					if(pexons[p].lpos <= fr1.lpos+pexon_range && pexons[p].lpos >= fr1.lpos-pexon_range && pexons[p].ltype == START_BOUNDARY)
					//if(pexons[p].lpos == fr1.lpos && pexons[p].ltype == START_BOUNDARY)
					{
						//left_boundary_flag = 1;
						pexon_left_flag = 1;
						break;
					}
				}
			}

			//checking if reads junction matches right boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];
				
				if(jc.lpos == fr2.rpos)
				//if(jc.lpos <= fr2.rpos+junc_range && jc.lpos >= fr2.rpos-junc_range)
				{
					//printf("jc.lpos = %d\n",jc.lpos);
					right_boundary_flag = 1;
					junc_match = 'R';
					break;
				}
			}
			
			//checking if ref junction matches right boundary
			/*temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					// assert(chain[k].first < chain[k].second);
					// if(chain[k].first <= bb.lpos) continue;
					// if(chain[k].second >= bb.rpos) continue;

					if(chain[k].first == fr2.rpos)
					//if(chain[k].first <= fr2.rpos+junc_range && chain[k].first >= fr2.rpos-junc_range)
					{
						right_boundary_flag = 1;
						ref_match = 'R';
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/
		
			//checking if pexon matches right boundary when ref file is not given
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					if(pexons[p].rpos <= fr2.rpos+pexon_range && pexons[p].rpos >= fr2.rpos-pexon_range && pexons[p].rtype == END_BOUNDARY)
					//if(pexons[p].rpos == fr2.rpos && pexons[p].rtype == END_BOUNDARY)
					{
						//right_boundary_flag = 1;
						pexon_right_flag = 1;
						break;
					}
				}
			}

			if((left_boundary_flag == 1 || pexon_left_flag == 1) && (right_boundary_flag == 1 || pexon_right_flag == 1))
			{
				printf("Found a case with junc comp 1\n");
				printf("valid: left_boundary_flag = %d, right_boundary_flag = %d, circ left = %d, circ right = %d, bundle left = %d, bundle right = %d\n",left_boundary_flag, right_boundary_flag, fr1.lpos, fr2.rpos, bb.lpos, bb.rpos);
				join_circ_fragment_pair(circ_fragment_pairs[i],0,0,left_boundary_flag,right_boundary_flag);
			}
			else
			{
				printf("Not valid: left_boundary_flag = %d, right_boundary_flag = %d, circ left = %d, circ right = %d, bundle left = %d, bundle right = %d\n",left_boundary_flag, right_boundary_flag, fr1.lpos, fr2.rpos, bb.lpos, bb.rpos);
			}
		}
		else if(fr2.is_compatible == 2)
		{
			//checking if reads junction matches left boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.rpos == fr2.lpos)
				//if(jc.rpos <= fr2.lpos+junc_range && jc.rpos >= fr2.lpos-junc_range)
				{
					//printf("jc.rpos = %d\n",jc.rpos);
					left_boundary_flag = 1;
					junc_match = 'L';
					break;
				}
			}

			//checking if ref junction matches left boundary
			/*int temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					// assert(chain[k].first < chain[k].second);
					// if(chain[k].first <= bb.lpos) continue;
					// if(chain[k].second >= bb.rpos) continue;

					if(chain[k].second == fr2.lpos)
					//if(chain[k].second <= fr2.lpos+junc_range && chain[k].second >= fr2.lpos-junc_range)
					{
						left_boundary_flag = 1;
						ref_match = 'L';
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/

			//checking if pexon matches left boundary
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					if(pexons[p].lpos <= fr2.lpos+pexon_range && pexons[p].lpos >= fr2.lpos-pexon_range && pexons[p].ltype == START_BOUNDARY)
					//if(pexons[p].lpos == fr2.lpos && pexons[p].ltype == START_BOUNDARY)
					{
						//left_boundary_flag = 1;
						pexon_left_flag = 1;
						break;
					}
				}
			}

			//checking if reads junction matches right boundary
			for(int j=0;j<junctions.size();j++)
			{
				junction jc = junctions[j];

				if(jc.lpos == fr1.rpos)
				//if(jc.lpos <= fr1.rpos+junc_range && jc.lpos >= fr1.rpos-junc_range)
				{
					//printf("jc.lpos = %d\n",jc.lpos);
					right_boundary_flag = 1;
					junc_match = 'R';
					break;
				}
			}

			//checking if ref junction matches right boundary
			/*temp_flag = 0;
			for(int t=0;t<ref_trsts.size();t++)
			{
				transcript trst = ref_trsts[t];
				vector<PI32> chain = trst.get_intron_chain();

				for(int k=0;k<chain.size();k++)
				{
					// assert(chain[k].first < chain[k].second);
					// if(chain[k].first <= bb.lpos) continue;
					// if(chain[k].second >= bb.rpos) continue;

					if(chain[k].first == fr1.rpos)
					//if(chain[k].first <= fr1.rpos+junc_range && chain[k].first >= fr1.rpos-junc_range)
					{
						right_boundary_flag = 1;
						ref_match = 'R';
						temp_flag = 1;
						break;
					}
				}

				if(temp_flag == 1)
				{
					break;
				}
			}*/

			//checking if pexon matches right boundary
			if(ref_file == "")
			{
				for(int p=0;p<pexons.size();p++)
				{
					if(pexons[p].rpos <= fr1.rpos+pexon_range && pexons[p].rpos >= fr1.rpos-pexon_range && pexons[p].rtype == END_BOUNDARY)
					//if(pexons[p].rpos == fr1.rpos && pexons[p].rtype == END_BOUNDARY)
					{
						//right_boundary_flag = 1;
						pexon_right_flag = 1;
						break;
					}
				}
			}

			if((left_boundary_flag == 1 || pexon_left_flag == 1) && (right_boundary_flag == 1 || pexon_right_flag == 1))
			{
				printf("Found a case with junc comp 2\n");
				printf("valid: left_boundary_flag = %d, right_boundary_flag = %d, circ left = %d, circ right = %d, bundle left = %d, bundle right = %d\n",left_boundary_flag, right_boundary_flag, fr2.lpos, fr1.rpos, bb.lpos, bb.rpos);
				join_circ_fragment_pair(circ_fragment_pairs[i],0,0,left_boundary_flag,right_boundary_flag);
			}
			else
			{
				printf("Not valid: left_boundary_flag = %d, right_boundary_flag = %d, circ left = %d, circ right = %d, bundle left = %d, bundle right = %d\n",left_boundary_flag, right_boundary_flag, fr2.lpos, fr1.rpos, bb.lpos, bb.rpos);
			}
		}
	}

	return 0;
}

int bundle_bridge::join_circ_fragment_pair(pair<fragment,fragment> &fr_pair, int ex1, int ex2, int left_boundary_flag, int right_boundary_flag)
{

	fragment &fr1 = fr_pair.first;
	fragment &fr2 = fr_pair.second;

	vector<int> v1 = decode_vlist(fr1.paths[0].v);
	vector<int> v2 = decode_vlist(fr2.paths[0].v);
	path p;

	assert(v1.size() > ex1);
	assert(v2.size() > ex2);

	if(fr2.is_compatible == 1)
	{
		vector<int>::iterator t1 = v1.end() - ex1;
		vector<int>::iterator t2 = v2.begin() + ex2;

		int x1 = v1[v1.size() - 1 - ex1];
		int x2 = v2[ex2];

		if(x1 < x2) return 0;

		/*if(strcmp(fr2.h1->qname.c_str(),"simulate:2195") == 0)
		{
			printf("simulate:2195 is in fr2.is_compatible 1\n");
		}*/

		vector<int>::iterator it = find(t2, v2.end(), x1);
		if(it == v2.end()) return 0;

		vector<int>::iterator j1, j2;
		for(j1 = t1 - 1, j2 = it; j1 >= v1.begin() && j2 >= t2; j1--, j2--)
		{
			if((*j1) != (*j2)) return 0;
		}

		p.ex1 = ex1;
		p.ex2 = ex2;
		p.v.insert(p.v.end(), v1.begin(), t1);
		p.v.insert(p.v.end(), it + 1, v2.end());
		//p.v = encode_vlist(p.v);

		printf("Printing merged path in fr2 iscomp 1:\n");
		printv(p.v);
		printf("\n\n");
	
		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr1.lpos) + "|" + tostring(fr2.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr1.lpos;
		int32_t end = fr2.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), p.v.begin(), p.v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2";
		if(fr2.h1->is_fake == true || fr2.h2->is_fake == true) circ.source = "scallop2_MC";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr1.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);
		
		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
			//circ.circRNA_id = circ.circRNA_id + tostring(circ.circ_path[i]) + "|";
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}

		//return if single exon circRNA and any one side junction flag not 1
		if(circ.merged_regions.size() == 1 && (left_boundary_flag != 1 || right_boundary_flag != 1))
		{
			return 0;
		}
		circ_trsts.push_back(circ);
	}
	else if(fr2.is_compatible == 2)
	{
		vector<int>::iterator t1 = v2.end() - ex1;
		vector<int>::iterator t2 = v1.begin() + ex2;

		int x1 = v2[v2.size() - 1 - ex1];
		int x2 = v1[ex2];

		if(x1 < x2) return 0;

		vector<int>::iterator it = find(t2, v1.end(), x1);
		if(it == v1.end()) return 0;

		vector<int>::iterator j1, j2;
		for(j1 = t1 - 1, j2 = it; j1 >= v2.begin() && j2 >= t2; j1--, j2--)
		{
			if((*j1) != (*j2)) return 0;
		}

		p.ex1 = ex1;
		p.ex2 = ex2;
		p.v.insert(p.v.end(), v2.begin(), t1);
		p.v.insert(p.v.end(), it + 1, v1.end());
		//p.v = encode_vlist(p.v);

		printf("Printing merged path in fr2 iscomp 2:\n");
		printv(p.v);
		printf("\n\n");
		
		string chrm_id = bb.chrm.c_str();
		string circRNA_id = "chrm" + chrm_id + ":" + tostring(fr2.lpos) + "|" + tostring(fr1.rpos) + "|";
		//printf("circularRNA = %s\n",circRNA_id.c_str());
		
		char strand = bb.strand;
		int32_t start = fr2.lpos;
		int32_t end = fr1.rpos;
    	vector<int> circ_path;
		circ_path.insert(circ_path.begin(), p.v.begin(), p.v.end());

		circular_transcript circ;
		circ.circRNA_id = circRNA_id;
		circ.seqname = chrm_id;
		circ.source = "scallop2";
		if(fr2.h1->is_fake == true || fr2.h2->is_fake == true) circ.source = "scallop2_MC";
		circ.feature = "circRNA";
		circ.gene_id = "gene"; //later change this to bundle id
		circ.transcript_id = fr1.h1->qname; //use hit qname, same for all hits in fragments
		circ.start = start;
		circ.end = end;
		circ.circ_path.insert(circ.circ_path.begin(),circ_path.begin(),circ_path.end());
		circ.strand = infer_circ_strand(circ.circ_path);

		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(regions[circ.circ_path[i]]);
			//circ.circRNA_id = circ.circRNA_id + tostring(circ.circ_path[i]) + "|";
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}

		if(circ.merged_regions.size() == 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[0].rpos = circ.end; 
		}
		else if(circ.merged_regions.size() > 1)
		{
			circ.merged_regions[0].lpos = circ.start;
			circ.merged_regions[circ.merged_regions.size()-1].rpos = circ.end; 			
		}

		for(int i=0;i<circ.merged_regions.size();i++)
    	{
			region r = circ.merged_regions[i];
			circ.circRNA_id = circ.circRNA_id + tostring(r.lpos) + "|" + tostring(r.rpos) + "|";
		}		
		
		//return if single exon circRNA and any one side junction flag not 1
		if(circ.merged_regions.size() == 1 && (left_boundary_flag != 1 || right_boundary_flag != 1))
		{
			return 0;
		}
		circ_trsts.push_back(circ);
	}
	else
	{
		printf("is_compatible not 1 or 2\n");
	}

	return 0;
}

char bundle_bridge::infer_circ_strand(const vector<int> &p)
{
	if(p.size() <= 1) return '.';
	int n1 = 0, n2 = 0, n3 = 0;
	for(int k = 0; k < p.size() - 1; k++)
	{
		int32_t x1 = regions[p[k + 0]].rpos;
		int32_t x2 = regions[p[k + 1]].lpos;
		if(x1 >= x2) continue;
		int64_t jc = pack(x1, x2);
		if(junc_map.find(jc) == junc_map.end()) continue;
		char c = junc_map[jc];
		if(c == '.') n1 ++;
		if(c == '+') n2 ++;
		if(c == '-') n3 ++;
	}
	if(n2 >= 1 && n1 <= 0 && n3 <= 0) return '+';
	if(n3 >= 1 && n1 <= 0 && n2 <= 0) return '-';
	printf("MIXED JUNCTION: ./+/- = %d/%d/%d\n", n1, n2, n3);
	return '.';
}

int bundle_bridge::print_circRNAs()
{
	if(circ_trsts.size() > 0)
	{
		//printf("circ_trsts size: %lu", circ_trsts.size());
		printf("Printing circRNAs from current bundle bridge:\n");

		for(int i=0;i<circ_trsts.size();i++)
		{
			circ_trsts[i].print(i+1);
		}

		printf("\n");
	}

	return 0;
}

int bundle_bridge::group_fragments()
{
	if(fragments.size() == 0) return 0;

	sort(fragments.begin(), fragments.end(), compare_fragment);

	vector<fragment> ff;

	fragment fx = fragments[0];
	assert(fx.h1->vlist.size() >= 1);
	assert(fx.h2->vlist.size() >= 1);
	for(int k = 1; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		assert(fr.h1->vlist.size() >= 1);
		assert(fr.h2->vlist.size() >= 1);

		if(fx.equal(fr) == true)
		{
			fx.append(fr);
		}
		else
		{
			ff.push_back(fx);
			fx = fr;
		}
	}
	ff.push_back(fx);
	fragments = ff;

	//printf("grouped fragments = %lu\n", fragments.size());
	return 0;
}

int32_t bundle_bridge::compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v)
{
	if(v.size() == 0) return 0;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i];
		flen += regions[k].rpos - regions[k].lpos;
	}
	return flen - k1l - k2r;
}

int bundle_bridge::print(int index)
{
	printf("Bundle %d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	printf("tid = %d, #hits = %lu, #fragments = %lu, #ref-trsts = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			bb.tid, bb.hits.size(), fragments.size(), ref_trsts.size(), bb.chrm.c_str(), bb.lpos, bb.rpos, bb.strand, n0, np, nq);

	// print ref-trsts
	//for(int k = 0; k < ref_trsts.size(); k++) ref_trsts[k].write(cout);

	// print fragments 
	//for(int i = 0; i < fragments.size(); i++) fragments[i].print(i);

	/*
	// print fclusters
	for(int i = 0; i < fclusters.size(); i++) fclusters[i].print(i);
	*/

	if(verbose <= 1) return 0;

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	// print bb.hits
	for(int i = 0; i < bb.hits.size(); i++) bb.hits[i].print();

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	printf("\n");

	return 0;
}

vector<int32_t> bundle_bridge::build_accumulate_length(const vector<int> &v)
{
	int32_t x = 0;
	vector<int32_t> acc;
	acc.resize(v.size());
	for(int i = 0; i < v.size(); i++)
	{
		x += regions[v[i]].rpos - regions[v[i]].lpos;
		acc[i] = x;
	}
	return acc;
}

vector<int32_t> bundle_bridge::get_aligned_intervals(fragment &fr)
{
	vector<int32_t> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<int32_t> v = get_splices(fr);
	if(v.size() >= 1 && fr.h1->pos >= v.front()) return vv;
	if(v.size() >= 1 && fr.h2->rpos <= v.back()) return vv;

	v.insert(v.begin(), fr.h1->pos);
	v.push_back(fr.h2->rpos);
	return v;
}

vector<int32_t> bundle_bridge::get_splices(fragment &fr)
{
	vector<int32_t> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<int> v = decode_vlist(fr.paths[0].v);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() - 1; i++)
	{
		int32_t pp = regions[v[i + 0]].rpos;
		int32_t qq = regions[v[i + 1]].lpos;
		if(pp >= qq) continue;
		vv.push_back(pp);
		vv.push_back(qq);
	}
	return vv;
}
