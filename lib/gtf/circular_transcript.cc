/*
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "circular_transcript.h"
#include "util.h"

circular_transcript::circular_transcript()
{
    circRNA_id = "";
    seqname = "";
    source = "";
    feature = "";
    gene_id = "";
    transcript_id = "";
    gene_type = "";
    transcript_type = "";
	start = 0;
	end = 0;
    score = 0;
    strand = '.';
    frame = -1;
	coverage = 0;
    covratio = 0;
	RPKM = 0;
    FPKM = 0;
	TPM = 0;

    fake_supple = false;
    supple_len = 0;
    path_score = 0;
    path_type = 0;
    bundle_size = 0;
    ref_trsts_size = 0;
    exon_count = 0;
	total_exon_length = 0;
	min_exon_length = 0;
	max_exon_length = 0;
	avg_exon_length = 0;
	path_count_1 = 0;
	path_count_2 = 0;
	path_count_3 = 0;
	path_count_4 = 0;
    candidate_path_count = 0;
    fake_count = 0;

    junc_reads = 0;
    non_junc_reads = 0;
    circ_path.clear();
    circ_path_regions.clear();
    merged_regions.clear();
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path)
{
    circRNA_id = circRNA_ID;
    seqname = chrm_id;
    start = start;
    end = end;
    circ_path.clear();
    this->circ_path.insert(this->circ_path.begin(), circ_path.begin(),circ_path.end());
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path, int32_t junc_reads, int32_t non_junc_reads)
{
    circRNA_id = circRNA_ID;
    seqname = chrm_id;
    start = start;
    end = end;
    circ_path.clear();
    this->circ_path.insert(this->circ_path.begin(), circ_path.begin(),circ_path.end());

    junc_reads = junc_reads;
    non_junc_reads = non_junc_reads;
}

int circular_transcript::write(ostream &fout, double cov2, int count) const
{
    fout.precision(4);
	fout<<fixed;
    
    //fout<<"chrm_id="<<chrm_id.c_str()<<"\t";
    fout<<seqname.c_str()<<"\t";
    fout<<source.c_str()<<"\t";
    fout<<feature.c_str()<<"\t";
    fout<<start + 1<<"\t";
    fout<<end<<"\t";
	fout<<score<<"\t";							// score, for now as zero
	fout<<strand<<"\t";							// strand
	fout<<".\t";								            // frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
    fout<<"cov \""<<coverage<<"\";"<<endl;
    
    /*fout<<"path vertices= ( ";
    for(int i=0;i<circ_path.size();i++)
    {
        fout<<circ_path[i]<<" "; 
    }
    fout<<")";
    fout<<endl;*/

    /*join_interval_map jmap;
	for(int k = 0; k < circ_path_regions.size() ; k++)
	{
		int32_t p1 = circ_path_regions[k].lpos;
		int32_t p2 = circ_path_regions[k].rpos;
		jmap += make_pair(ROI(p1, p2), 1);
	}

	vector<region> merged_regions;
    merged_regions.clear();
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
        region r(lower(it->first), upper(it->first), '.', '.');
        merged_regions.push_back(r);
	}*/

    int cnt = 0;
    for(int i=0;i<merged_regions.size();i++)
    {
        fout<<seqname.c_str()<<"\t";
        fout<<source.c_str()<<"\t";
        fout<<"exon\t";

        /*if(merged_regions.size() == 1)
        {
            fout<<start + 1<<"\t";          
            fout<<end<<"\t";
        }
        else
        {
            if(i == 0)
            {
                fout<<start + 1<<"\t";           
                fout<<merged_regions[i].rpos<<"\t";             
            }
            if(i == circ_path_regions.size()-1)
            {
                fout<<merged_regions[i].lpos + 1<<"\t";           //is +1 needed here? yes
                fout<<end<<"\t";
            }
            if(i > 0 && i < circ_path_regions.size()-1)
            {
                fout<<merged_regions[i].lpos + 1<<"\t";           //is +1 needed here? yes
                fout<<merged_regions[i].rpos<<"\t";
            }
        }*/

        fout<<merged_regions[i].lpos + 1<<"\t";           
        fout<<merged_regions[i].rpos<<"\t";

        fout<<score<<"\t";							// score, for now as zero
        fout<<strand<<"\t";							// strand
        fout<<".\t";								            // frame
        fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
        fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
        fout<<"exon \""<<++cnt<<"\"; "<<endl;
        //fout<<"coverage \""<<coverage<<"\";"<<endl;
    }

    return 0;
}

int circular_transcript::print(int id)
{
    //printf("circRNA %d - ", id);
    /*printf("circRNA_id:%s, chrm:%s, bundle_size:%lu, transcript_id:%s, start:%d, end:%d, coverage:%d, fake_supple:%d, supple_len:%d, path_score: %lf, path_type: %d, path: ",circRNA_id.c_str(), seqname.c_str(), bundle_size, transcript_id.c_str(), start, end, coverage, (fake_supple == true)?1:0, supple_len, path_score, path_type);
    for(int i = 0; i < circ_path.size() - 1; i++)
	{
		printf("%d, ", circ_path[i]);
	}
	printf("%d\n", circ_path[circ_path.size() - 1]);
    printf("path coordinates= ");
    for(int i=0;i<circ_path_regions.size();i++)
    {
        printf("[%d, %d) ",circ_path_regions[i].lpos,circ_path_regions[i].rpos); 
    }
    printf("\n");
    printf("merged coordinates= ");
        for(int i=0;i<merged_regions.size();i++)
    {
        printf("[%d, %d) ",merged_regions[i].lpos,merged_regions[i].rpos); 
    }
    printf("\n");*/
    printf("circRNA_id:%s, chrm:%s, bundle_size:%lu, coverage:%d, fake_count:%d, supple_len:%d, candidate_path_count:%d, path_score:%lf, path_count_1:%d, path_count_2:%d, path_count_3:%d, path_count_4:%d, exon_count:%d, total_exon_len:%d, max_exon_len:%d, min_exon_len:%d, avg_exon_len:%lf\n",circRNA_id.c_str(), seqname.c_str(), bundle_size, coverage, fake_count, supple_len, candidate_path_count, path_score, path_count_1, path_count_2, path_count_3, path_count_4, exon_count, total_exon_length, max_exon_length, min_exon_length, avg_exon_length);
    return 0;

}

circular_transcript::~circular_transcript()
{
}

