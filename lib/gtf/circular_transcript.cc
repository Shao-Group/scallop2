#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "circular_transcript.h"
#include "util.h"

circular_transcript::circular_transcript()
{
    circRNA_ID = "";
    chrm_id = "";
	start = 0;
	end = 0;

    junc_reads = 0;
    non_junc_reads = 0;
    circ_path.clear();
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path)
{
    circRNA_ID = circRNA_ID;
    chrm_id = chrm_id;
    start = start;
    end = end;
    circ_path.clear();
    this->circ_path.insert(this->circ_path.begin(), circ_path.begin(),circ_path.end());
}

circular_transcript::circular_transcript(string circRNA_ID, string chrm_id, int32_t start, int32_t end, vector<int> circ_path, int32_t junc_reads, int32_t non_junc_reads)
{
    circRNA_ID = circRNA_ID;
    chrm_id = chrm_id;
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

    fout<<"ID="<<circRNA_ID.c_str()<<"\t";
    fout<<"chrm="<<chrm_id.c_str()<<"\t";
    fout<<"start="<<start<<"\t";
    fout<<"end="<<end<<"\t";
    fout<<"path vertices= ( ";
    for(int i=0;i<circ_path.size();i++)
    {
        fout<<circ_path[i]<<" "; 
    }
    fout<<")";
    fout<<endl;
    //fout<< "end";

    return 0;
}

int circular_transcript::print(int id)
{
    printf("circRNA %d - ", id);
    printf("ID: %s, chrm_ID: %s, start: %d, end: %d, path: ",circRNA_ID.c_str(), chrm_id.c_str(), start, end);
    for(int i = 0; i < circ_path.size() - 1; i++)
	{
		printf("%d, ", circ_path[i]);
	}
	printf("%d\n", circ_path[circ_path.size() - 1]);

    return 0;

}

circular_transcript::~circular_transcript()
{
}

