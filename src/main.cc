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
	print_parameters();

	if(verbose >= 1)
	{
		print_copyright();
		printf("\n");
		print_command_line(argc, argv);
		printf("\n");
	}

	previewer pv; 
	pv.preview(); //resolve strandness and estimate fragment length distirbution

	if(preview_only == true) return 0;

	assembler asmb;
	asmb.assemble();

	printf("Outward count = %d\n", asmb.outward_count);
	printf("Non outward count = %d\n", asmb.non_outward_count);

	return 0;
}
