#include "tss_tes.h"
#include <unordered_map>
#include "config.h"
#include <cmath>
#include <map>

// Global constant for berth neighborhood size
// const int berth_neighborhood = 50;

tss_tes::tss_tes(int type)
{
    this->type = type;
}

tss_tes::tss_tes(int type, int32_t pos, int weight_sg, int weight_berth)
{
    this->type = type;
    this->pos = pos;
    this->weight_sg = weight_sg;
    this->weight_berth = weight_berth;
}

void tss_tes::calculate_clip_length(vector<hit> &sorted_hits_compatible, char bb_strand)
{
    double mean_clip_length = 0, stddev_clip_length = 0, soft_clip_entropy = 0;
    vector<int> clip_lengths;
    clip_lengths.reserve(sorted_hits_compatible.size());
    for(int i = 0; i < sorted_hits_compatible.size(); i++)
    {
        hit &h = sorted_hits_compatible[i];
        float lc = abs(h.itvc1.second - h.itvc1.first) ;
        float rc = abs(h.itvc2.second - h.itvc2.first) ;
        if(h.strand == '+' && this->type == 0) clip_lengths.push_back(lc);
        else if (h.strand == '+' && this->type == 1) clip_lengths.push_back(rc);
        else if (h.strand == '-' && this->type == 0) clip_lengths.push_back(rc);
        else if (h.strand == '-' && this->type == 1) clip_lengths.push_back(lc);
    }
    
    unordered_map<int, double> clip_length_freq;
    clip_length_freq.reserve(clip_lengths.size());
    for (auto &l : clip_lengths)
    {
        mean_clip_length += l / clip_lengths.size();
        clip_length_freq[l] += 1 / clip_lengths.size();
    }
    for (auto &l : clip_lengths)
    {
        stddev_clip_length += (l - mean_clip_length) * (l - mean_clip_length);
    }
    stddev_clip_length = sqrt(stddev_clip_length / clip_lengths.size());

    // calculate soft clip entropy
    for (auto &p : clip_length_freq)
    {
        soft_clip_entropy -= p.second * log2(p.second);
    }
    
    this->mean_clip_length = mean_clip_length;
    this->std_clip_length = stddev_clip_length;
    this->soft_clip_entropy = soft_clip_entropy;

}

void tss_tes::calculate_junction_cnt(vector<junction> &sorted_junctions_start, vector<junction> &sorted_junctions_end)
{
    auto junction_start_low = lower_bound(sorted_junctions_start.begin(), sorted_junctions_start.end(), 
        this->pos - berth_neighborhood, [](const junction &a, const int32_t &b) { return a.lpos < b; });
    auto junction_start_high = upper_bound(sorted_junctions_start.begin(), sorted_junctions_start.end(), 
        this->pos + berth_neighborhood, [](const int32_t &a, const junction &b) { return a < b.lpos; });
    this->junction_start_cnt = std::distance(junction_start_low, junction_start_high);
    
    auto junction_end_low = lower_bound(sorted_junctions_end.begin(), sorted_junctions_end.end(), 
        this->pos - berth_neighborhood, [](const junction &a, const int32_t &b) { return a.rpos < b; });
    auto junction_end_high = upper_bound(sorted_junctions_end.begin(), sorted_junctions_end.end(), 
        this->pos + berth_neighborhood, [](const int32_t &a, const junction &b) { return a < b.rpos; });
    this->junction_end_cnt = std::distance(junction_end_low, junction_end_high);

    // Count junctions that cross the neighborhood
    this->junction_cross_cnt = 0;
    for(auto &j : sorted_junctions_start)
    {
        if(j.lpos >= this->pos - berth_neighborhood) break;
        if(j.rpos <= this->pos + berth_neighborhood) continue;
        if(j.lpos < this->pos - berth_neighborhood && j.rpos > this->pos + berth_neighborhood)
        {
            this->junction_cross_cnt++;
        }
    }
}

void tss_tes::calculate_anchor_features(vector<hit> &hits)
{
    int anchor_count = 0;
    double mean_anchor_padding = 0, stddev_anchor_padding = 0;
    vector<int> anchor_padding;
    int left_padding, right_padding;
    for(auto &h : hits)
    {
        left_padding = h.is_anchor_satisfactory(0, 1000) ? h.left_anchor_padding : -1;
        right_padding = h.is_anchor_satisfactory(1, 1000) ? h.right_anchor_padding : 0;
        if (h.strand == '+' && this->type == 0) anchor_padding.push_back(left_padding);
        else if (h.strand == '+' && this->type == 1) anchor_padding.push_back(right_padding);
        else if (h.strand == '-' && this->type == 0) anchor_padding.push_back(right_padding);
        else if (h.strand == '-' && this->type == 1) anchor_padding.push_back(left_padding);
    }

    for (auto &p : anchor_padding)
    {
        mean_anchor_padding += p / anchor_padding.size();
    }
    for (auto &p : anchor_padding)
    {
        stddev_anchor_padding += (p - mean_anchor_padding) * (p - mean_anchor_padding);
    }
    stddev_anchor_padding = sqrt(stddev_anchor_padding / anchor_padding.size());
}
