/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "previewer.h"
#include "assembler.h"

using namespace std;

int init_files()
{
	string tss_file_name = berth_folder + string("tss_merged_features.tsv");
	string tes_file_name = berth_folder + string("tes_merged_features.tsv");
	// string read_st_end_file_name = berth_folder + string("read_st_end.tsv");
	
	// cout << tss_file_name << tes_file_name << endl;
    ofstream tss_file(tss_file_name);
	ofstream tes_file(tes_file_name);
	ofstream anchor_file(anchor_file_name);
	cout << "hit_id\t" << "left_clipped_seq\t" << "left_adapter_pos\t" << "right_clipped_seq\t" << "right_adapter_pos" << endl;
	// ofstream read_st_end_file(read_st_end_file_name);

	tss_file.close();
	tes_file.close();
	anchor_file.close();
	// read_st_end_file.close();
	return 0;
}

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		print_copyright();
		print_help();
		printf("\n");
		return 0;
	}

	parse_arguments(argc, argv);

	if(verbose >= 1)
	{
		print_copyright();
		printf("\n");
		print_command_line(argc, argv);
		printf("\n");
	}

	// previewer pv;
	// pv.preview();

	if(preview_only == true) return 0;

	init_files();
	assembler asmb;
	asmb.assemble();

	return 0;
}
