#ifndef __TSS_TES_H__
#define __TSS_TES_H__

#include <vector>
#include <cstdint>
#include "hit.h"
#include "junction.h"

using namespace std;

class tss_tes
{
public:
    tss_tes(int type = 0);
    tss_tes(int type, int32_t pos, int weight_berth, int weight_sg);
    void calculate_junction_cnt(vector<junction> &sorted_junctions_start, vector<junction> &sorted_junctions_end);
    void calculate_clip_length(vector<hit> &hits, char bb_strand);
    void calculate_anchor_features(vector<hit> &hits);
    void soft_clip_entropy(vector<hit> &hits);

public:
    int type;                          // 0: TSS, 1: TES
    int32_t pos;                       // position
    int weight_berth;                  // weight from berth
    int weight_sg;                     // weight from splice graph
    int read_density;                  // count of read starting/ending in the neighborhood
    int spanning_reads_cnt;            // number of reads passing through the neighborhood
    float leading_clip_length;         // average leading soft clip length of reads starting/ending in the neighborhood
    float trailing_clip_length;        // average trailing soft clip length of reads starting/ending in the neighborhood
    int junction_start_cnt;            // count of junction starting in the neighborhood
    int junction_end_cnt;              // count of junction ending in the neighborhood
    int junction_cross_cnt;            // count of junction crossing the neighborhood
    int left_anchor_cnt;               // count of anchors in leftclip of reads in the neighborhood
    int right_anchor_cnt;              // count of anchors in rightclip of reads in the neighborhood
    int left_anchor_padding_mean;      // average anchor padding length of leftclip of reads starting/ending in the neighborhood
    int right_anchor_padding_mean;     // average anchor padding length of rightclip reads starting/ending in the neighborhood
    double left_soft_clip_entropy;     // entropy of left soft clip length distribution
    double right_soft_clip_entropy;    // entropy of right soft clip length distribution
    int coverage_before;               // coverage before the site
    int coverage_after;                // coverage after the site
    int delta_coverage;                // coverage before the site minus coverage after the site within the neighborhood
};

#endif 