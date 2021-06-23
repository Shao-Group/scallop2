/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.

Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
(c) 2019 by Mingfu Shao and The Pennsylvania State University.

See LICENSE for licensing.
*/

#ifndef __PREVIEWER_H__
#define __PREVIEWER_H__

#include "hit.h"
#include "bundle_base.h"

#include <fstream>
#include <string>

using namespace std;

class previewer
{
private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

public:
	int preview();

private:
	int open_file();
	int close_file();
	int solve_strandness();
	int solve_insertsize();
	int process_bundle(bundle_base& bb, map<int32_t, int>& m);
};

#endif
