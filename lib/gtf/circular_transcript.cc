#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "circular_transcript.h"
#include "util.h"

circular_transcript::circular_transcript()
{
    //chrm_id = "";
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

    junc_reads = 0;
    non_junc_reads = 0;
    circ_path.clear();
    circ_path_regions.clear();
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path)
{
    transcript_id = circRNA_ID;
    seqname = chrm_id;
    start = start;
    end = end;
    circ_path.clear();
    this->circ_path.insert(this->circ_path.begin(), circ_path.begin(),circ_path.end());
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path, int32_t junc_reads, int32_t non_junc_reads)
{
    transcript_id = circRNA_ID;
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
    fout<<"seqname="<<seqname.c_str()<<"\t";
    fout<<"source="<<source.c_str()<<"\t";
    fout<<"feature="<<feature.c_str()<<"\t";
    fout<<"start="<<start<<"\t";
    fout<<"end="<<end<<"\t";
	fout<<"score="<<score<<"\t";							// score, for now as zero
	fout<<"strand="<<strand<<"\t";							// strand
	fout<<".\t";								            // frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
    fout<<"coverage \""<<coverage<<"\";"<<endl;
    
    fout<<"path vertices= ( ";
    for(int i=0;i<circ_path.size();i++)
    {
        fout<<circ_path[i]<<" "; 
    }
    fout<<")";
    fout<<endl;

    join_interval_map jmap;
	for(int k = 0; k < circ_path_regions.size() ; k++)
	{
		int32_t p1 = circ_path_regions[k].lpos;
		int32_t p2 = circ_path_regions[k].rpos;
		jmap += make_pair(ROI(p1, p2), 1);
	}

	int cnt = 0;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
        fout<<"seqname="<<seqname.c_str()<<"\t";
        fout<<"source="<<source.c_str()<<"\t";
        fout<<"feature=exon\t";
        fout<<"start="<<lower(it->first)<<"\t"; //is +1 needed here?
        fout<<"end="<<upper(it->first)<<"\t";
        fout<<"score="<<score<<"\t";							// score, for now as zero
        fout<<"strand="<<strand<<"\t";							// strand
        fout<<".\t";								            // frame
        fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
        fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon_number \""<<++cnt<<"\"; ";
		fout<<"coverage \""<<coverage<<"\";"<<endl;
	}

    fout<<"path coordinates= ";
    for(int i=0;i<circ_path_regions.size();i++)
    {
        fout<<"["<<circ_path_regions[i].lpos<<","<<circ_path_regions[i].rpos<<") "; 
        /*if(circ_path_regions.size() == 1)
        {
            fout<<"["<<start<<","<<end<<") ";
        }
        else
        {
            if(i == 0)
            {
                fout<<"["<<start<<","<<circ_path_regions[i].rpos<<") ";
            }
            if(i == circ_path_regions.size()-1)
            {
                fout<<"["<<circ_path_regions[i].lpos<<","<<end<<") ";
            }
            if(i > 0 && i < circ_path_regions.size()-1)
            {
                fout<<"["<<circ_path_regions[i].lpos<<","<<circ_path_regions[i].rpos<<") "; 
            }
        }*/
    }
    fout<<endl<<endl;

    return 0;
}

int circular_transcript::print(int id)
{
    printf("circRNA %d - ", id);
    printf("seqname: %s, transcript_id: %s, start: %d, end: %d, path: ",seqname.c_str(), transcript_id.c_str(), start, end);
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

    return 0;

}

circular_transcript::~circular_transcript()
{
}

