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

int nt;		//total fragments #
int nnp;	//total paired end fr
int nnu;	//total umi-lined only fr
int nnb;	//total both fr
int np;		//bridged paired-end fr
int nu;		//bridged UMI-linked only fr
int nb;		//bridged both fr

int main(int argc, const char **argv)
{
	srand(time(0));

	// TODO: bridged fragment type count
	nt = nnp = nnu = nnb = np = nu = nb = 0;

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

	previewer pv;
	pv.preview();

	if(preview_only == true) return 0;

	assembler asmb;
	asmb.assemble();

	printf("Total fragments = %d, paired-end = %d, UMI-linked only = %d, both = %d\nBridged fragments: total = %d, paired-end = %d, UMI-linked only = %d, both = %d\n", nt, nnp, nnu, nnb, np+nu, np, nu, nb);

	return 0;
}
