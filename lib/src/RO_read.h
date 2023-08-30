#ifndef __RO_READ_H__
#define __RO_READ_H__

#include <string>

using namespace std;

class RO_read
{
    public:
        string read_name;
        string chrm;
        int32_t BSJ_pos;
        int32_t BSJ_rpos;

    public:
        RO_read();
        ~RO_read();
        int print();
        
};

#endif