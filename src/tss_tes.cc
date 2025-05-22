#include "tss_tes.h"
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
    float avg_leading_clip_length = 0, avg_trailing_clip_length = 0;
    for(int i = 0; i < sorted_hits_compatible.size(); i++)
    {
        hit &h = sorted_hits_compatible[i];
        float lc = abs(h.itvc1.second - h.itvc1.first) / sorted_hits_compatible.size();
        float rc = abs(h.itvc2.second - h.itvc2.first) / sorted_hits_compatible.size();
        if(h.strand == '+' || (h.strand == '.' && bb_strand == '+'))
        {
            avg_leading_clip_length += lc;
            avg_trailing_clip_length += rc;
        }
        else
        {
            avg_leading_clip_length += rc;
            avg_trailing_clip_length += lc;
        }
    }
    
    this->leading_clip_length = avg_leading_clip_length;
    this->trailing_clip_length = avg_trailing_clip_length;
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
    int left_anchor_count = 0, right_anchor_count = 0;
    float left_anchor_padding_sum = 0, right_anchor_padding_sum = 0;
    
    for(auto &h : hits)
    {
        if(h.is_anchor_satisfactory(0, 1000))
        {
            left_anchor_count++;
            left_anchor_padding_sum += h.left_anchor_padding;
        }
        if(h.is_anchor_satisfactory(1, 1000))
        {
            right_anchor_count++;
            right_anchor_padding_sum += h.right_anchor_padding;
        }
    }

    this->left_anchor_cnt = left_anchor_count;
    this->right_anchor_cnt = right_anchor_count;
    this->left_anchor_padding_mean = left_anchor_count > 0 ? left_anchor_padding_sum / left_anchor_count : 0;
    this->right_anchor_padding_mean = right_anchor_count > 0 ? right_anchor_padding_sum / right_anchor_count : 0;
}

void tss_tes::soft_clip_entropy(vector<hit> &hits)
{
    vector<int> left_soft_clip_lengths, right_soft_clip_lengths;
    map<int, int> left_freq, right_freq;
    
    // Collect lengths and frequencies
    for(auto &h : hits)
    {
        int left_len = h.left_s_seq.length();
        int right_len = h.right_s_seq.length();
        
        left_freq[left_len]++;
        right_freq[right_len]++;
        
        left_soft_clip_lengths.push_back(left_len);
        right_soft_clip_lengths.push_back(right_len);
    }
    
    // Calculate entropy for left clips
    double left_entropy = 0;
    for(auto &pair : left_freq)
    {
        double p = (double)pair.second / left_soft_clip_lengths.size();
        left_entropy -= p * log2(p);
    }
    
    // Calculate entropy for right clips
    double right_entropy = 0;
    for(auto &pair : right_freq)
    {
        double p = (double)pair.second / right_soft_clip_lengths.size();
        right_entropy -= p * log2(p);
    }
    
    this->left_soft_clip_entropy = left_entropy;
    this->right_soft_clip_entropy = right_entropy;
} 