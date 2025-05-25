#ifndef __TSS_TES_H__
#define __TSS_TES_H__

#include <vector>
#include <cstdint>
#include "hit.h"
#include "junction.h"
#include "bundle_base.h"

using namespace std;

class tss_tes
{
public:
    tss_tes(int type = 0);
    tss_tes(int type, int32_t pos, int weight_berth, int weight_sg);
    void calculate_junction_cnt(vector<junction> &sorted_junctions_start, vector<junction> &sorted_junctions_end);
    void calculate_clip_length(const vector<hit> &hits);
    void calculate_anchor_features(const vector<hit> &hits);
    void calculate_coverage_features(bundle_base &bb, int32_t pos);
    void build(bundle_base &bb, const vector<hit> &hits, vector<junction> &sorted_junctions_start, vector<junction> &sorted_junctions_end);
    // void soft_clip_entropy(vector<hit> &hits);

public:
    int type;                          // 0: TSS, 1: TES
    int32_t pos;                       // position
    int weight_berth;                  // weight from berth
    int weight_sg;                     // weight from splice graph
    int read_density;                  // count of read starting/ending in the neighborhood
    int spanning_reads_cnt;            // number of reads passing through the neighborhood
    double mean_clip_length;         // average leading soft clip length of reads starting/ending in the neighborhood
    double std_clip_length;        // average trailing soft clip length of reads starting/ending in the neighborhood
    int junction_start_cnt;            // count of junction starting in the neighborhood
    int junction_end_cnt;              // count of junction ending in the neighborhood
    int junction_cross_cnt;            // count of junction crossing the neighborhood
    int anchor_cnt;               // count of anchors in leftclip of reads in the neighborhood
    double mean_anchor_padding;      // average anchor padding length of leftclip of reads starting/ending in the neighborhood
    double stddev_anchor_padding;     // average anchor padding length of rightclip reads starting/ending in the neighborhood
    double soft_clip_entropy;     // entropy of left soft clip length distribution
    int coverage_before;               // coverage before the site
    int coverage_after;                // coverage after the site
    int delta_coverage;                // coverage before the site minus coverage after the site within the neighborhood

protected:
    // Helper function to determine which padding to use based on strand and type
    bool should_use_left_padding(char strand) const {
        return (strand == '+' && type == 0) || (strand == '-' && type == 1);
    }

    // Helper function to get appropriate anchor index based on strand and type
    int get_anchor_index(char strand) const {
        return should_use_left_padding(strand) ? 0 : 1;
    }

    // Helper function to get appropriate padding based on strand and type
    int get_appropriate_padding(const hit& h, int left_padding, int right_padding) const {
        return should_use_left_padding(h.strand) ? left_padding : right_padding;
    }

    int calculate_window_coverage(bundle_base &bb, int32_t window_start, int32_t window_end) {
        if(window_start >= window_end) return -1;
        if(bb.hits.size() == 0) return 0;

        // Use the pre-calculated mmap from bundle_base
        PSIMI pei = locate_boundary_iterators(bb.mmap, window_start, window_end);
        SIMI lit = pei.first, rit = pei.second;

        if(lit == bb.mmap.end()) return 0;

        int32_t window_size = window_end - window_start;
        int coverage_sum = compute_sum_overlap(bb.mmap, lit, rit);
        return coverage_sum;
    }
};

#endif 