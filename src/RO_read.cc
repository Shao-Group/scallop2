/*
Part of Coral
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "RO_read.h"

RO_read::RO_read()
{
    read_name = "";
    chrm = "";
    BSJ_pos = 0;
    BSJ_rpos = 0;
}

RO_read::~RO_read()
{
}

int RO_read::print()
{
    printf("read:%s, chrm:%s, BSJ:%d-%d\n",read_name.c_str(),chrm.c_str(),BSJ_pos,BSJ_rpos);
    return 0;
}
