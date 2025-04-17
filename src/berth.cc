/*
Part of Berth
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and Pennsylvania State University.
See LICENSE for licensing.
*/

#include "berth.h"
#include <cassert>
#include <cstdint>
#include <fstream> 

berth::berth(bundle_base &b)
	: bb(b), hits(b.hits)
{
    status = 0;
    berths.clear();
}

bool berth::is_empty() const
{
    return berths.empty();
}


const vector<PI32>& berth::get_berths() const
{
    assert (status != 0);
    return berths;
}

const vector<int>& berth::get_weights() const
{
    assert (status != 0);
    assert (berths.size() == weights.size());
    return weights;
}

/**
 * @param side - 0: left side of berth
 *             - 1: right side of berth
 * @return A map of left/right side pos of berth & weights
 */
map<int32_t, int>  berth::get_berth_side(int side) const
{
    assert (side == 0 || side == 1);
    assert (berths.size() == weights.size());
    map<int32_t, int> side_berth;
    for (int i = 0; i < berths.size(); i++)
    {

        int p;
        if (bb.strand == '+')
        {
            p = (side == 0) ? berths[i].first : berths[i].second;
        }
        else
        {   
            assert(bb.strand == '-');
            p = (side == 0) ? berths[i].second : berths[i].first;
        }
    
        if (side_berth.find(p) == side_berth.end())
        {
            side_berth[p] = weights[i];
        }
        else
        {
            side_berth[p] += weights[i];
        }        
    }
    return side_berth;
}


int berth::refine_berths(vector<int>& ss)
{
    assert (status != 0);
    if (status == 1 || is_empty()) return 0;
    refine_berths_by_splice_sites(ss);
    status = 3;
    return 0;
}

int berth::build_berths()
{
    assert (status == 0);
    
    init_sites();
    pick_peaks();
    assign_hit_berth();

    if (is_empty()) 
    {
        status = 1;
    }
    else 
    {
        status = 2;
    }
    
    this->print();
    string berth_file = berth_folder + string("berth.bed");
    string berth_file_one_end = berth_folder + string("b_one_end.bed");
    write_bed(berth_file.c_str());
    write_bed_one_end(berth_file_one_end.c_str());
    
    return 0;
}

int berth::init_sites()
{
    vector<int32_t> ss;  // TSS; may contain duplicates at the same position
	vector<int32_t> sl;  // left range of TSS
	vector<int32_t> tt;  // TTS
	vector<int32_t> tr;  // right range of TTS
    
    for (const auto &hit : hits) 
    {
        sl.push_back(hit.itvc1.first);
        ss.push_back(hit.itvc1.second);
        tt.push_back(hit.itvc2.first);
        tr.push_back(hit.itvc2.second);
    }
    assert (sl.size() == ss.size());
    assert (tt.size() == tr.size());

    vector<vector<int32_t>*> vecs = {&ss, &sl, &tt, &tr};
    int i = 0;
    for (auto v : vecs)
    {
        //CLEAN:
        i ++;
        // printf("un sorted [%d] ss sl tt tr: ", i);  
        // printv(*v);
        // printf("\n");  

        sort(v->begin(), v->end());
        
        // printf("is sorted [%d] ss sl tt tr: ", i);  
        // printv(*v);
        // printf("\n");
    }

    vector<map<int32_t, int>*> maps = {&ssc, &slc, &ttc, &trc};
    for (int i = 0; i < vecs.size(); i++)
    {
        const vector<int32_t> &v = *(vecs[i]);
        map<int32_t, int>     &m = *(maps[i]);
        for (int32_t val : v)
        {
            if (m.find(val) == m.end())
            {
                m[val] = 1;
            }
            else
            {
                m[val] ++;
            }
        }
    }

    //CLEAN:
    // for (int i = 0; i < maps.size(); i++)
    // {
    //     printf("map ssc etc [%d]....: \n", i);
    //     map<int32_t, int>     &m = *(maps[i]);
    //     for (const auto& [key, value] : m) 
    //         std::cout << key << ": " << value << '\n';
    // }
    // cout << "map ssc printed....: " << i << endl;

    return 0;
}

int berth::pick_peaks()
{
    // OR filter_gaussian();
    cout << "window start CLEAN: " << endl;
    filter_window(&ssc, false);
    cout << "window1 done CLEAN: " << endl;
    filter_window(&slc, false);
    cout << "window2 done  CLEAN: " << endl;
    filter_window(&ttc, true);
    cout << "window3 done  CLEAN: " << endl;
    filter_window(&trc, true);
    cout << "window4 done CLEAN: " << endl;
    return 0;
}

// dir: 0 for ss, 1 for tt
// find local maxima of a sliding window
// default SLIDINGD_WINDOW_SIZE = 10
int berth::filter_window(map<int32_t, int>* site_counts, bool dir)  
{
    if (site_counts->size() <= 0) return 0;
    
    map<int32_t, int> peakc;
    
    for (auto it = site_counts->begin(); it != site_counts->end(); ++it) 
    {
        int window_left, window_right;
        if (!dir)   // ss tail on right side   
        {    
            window_left = it->first;
            window_right = window_left + SLIDINGD_WINDOW_SIZE - 1;
        }
        else        // tt tail on left side
        {   
            window_right = it->first;
            window_left = window_right - SLIDINGD_WINDOW_SIZE + 1;
        }
        assert (window_left <= window_right);

        int sum = 0;
        for (auto jt = site_counts->lower_bound(window_left); jt != site_counts->upper_bound(window_right); ++jt) {
            sum += jt->second;
        }
        peakc[it->first] = sum;
    }

    assert(peakc.size() >= 1);

    {   //CLEAN:
        // std::cout << "peakc: " << '\n'; 
        // for (const auto& [key, value] : peakc) cout << key << ": " << value << '\n';
    }

    auto _tmp = find_local_or_isolated_max(peakc);
    {   //CLEAN:
        // std::cout << "peakc after local_or_iso_max: " << '\n'; 
        // for (const auto& [key, value] : _tmp) cout << key << ": " << value << '\n';
    }
    *site_counts = _tmp;
    // *site_counts = find_local_or_isolated_max(peakc);
    return 0;
}

// Find local maxima and isolated maxima 
// default: 500 bp away considered isolated
// default: count 3 or more for ss/tt
map<int32_t, int> berth::find_local_or_isolated_max(const map<int32_t, int>& peakc)
{
    map<int32_t, int> local_max;
    if (peakc.size() == 0) 
    {
        return local_max;
    }

    cout << "compute local_max ." << endl;

    if (peakc.size() == 1) 
    {
        if(peakc.begin()->second >= MIN_BOUND_COUNT) {
            local_max[peakc.begin()->first] = peakc.begin()->second;
        }
        return local_max;
    }
    
    cout << "compute local_max first pos" << endl;

    // first point
    if (peakc.size() >= 2) 
    {
        auto it1st = peakc.begin();
        auto it2nd = std::next(it1st);
        if ((it1st->second > it2nd->second) || (abs(it2nd->first - it1st->first) > ISOLATED_DISTANCE)) 
        {
            if (it1st->second >= MIN_BOUND_COUNT)
            {
                local_max[it1st->first] = it1st->second;
            }
        }
    }

    cout << "local_max size 1:" << local_max.size() << endl;
    // middle points
    for (auto it = std::next(peakc.begin()); it != prev(peakc.end()); ++it) 
    {
        if (it->second < MIN_BOUND_COUNT) continue;

        auto pre = std::prev(it);
        auto nxt = std::next(it);
        bool greater_than_or_isolated_from_pre = (it->second > pre->second) || (it->first - pre->first > ISOLATED_DISTANCE);
        bool greater_than_or_isolated_from_nxt = (it->second > nxt->second) || (nxt->first - it->first > ISOLATED_DISTANCE);
        if (! greater_than_or_isolated_from_pre) continue;
        if (! greater_than_or_isolated_from_nxt) continue;
        local_max[it->first] = it->second;
    }
    cout << "local_max size 2:" << local_max.size() << endl;
    
    // last point
    if (peakc.size() >= 3)  // otherwise added in the beginning
    {
        auto last1st = (--peakc.end());  
        auto last2nd = (--peakc.end());  
        if ((last1st->second > last2nd->second) || (abs(last2nd->first - last1st->first) > ISOLATED_DISTANCE)) 
        {
            if (last1st->second >= MIN_BOUND_COUNT) {
                local_max[last1st->first] = last1st->second;
            }
        }
    }
    cout << "local_max size 3: " << local_max.size() << endl;
    return local_max;
}

int berth::assign_hit_berth()
{
    for (int i = 0; i < hits.size(); ++i)
    {
        const auto & h = hits[i];

        //todo sl and tr
        int32_t htss = h.itvc1.second;
        int32_t htts = h.itvc2.first;
        auto it1 = closest_site(&ssc, htss);
        auto it2 = closest_site(&ttc, htts);
        
        int32_t tss = (it1 == ssc.end()) ? -1 : it1->first;
        int32_t tts = (it2 == ttc.end()) ? INT32_MAX : it2->first;
        
        auto b = make_pair(tss, tts);
        if (berth2hit.find(b) == berth2hit.end())
        {
            vector<int> v;
            v.push_back(i);
            berth2hit[b] = v;
        }
        else
        {
            berth2hit[b].push_back(i);
        }
    }

    for (const auto& [b, v]  : berth2hit)
    {
        int weight = v.size();
        if (weight < MIN_BERTH_HIT_COUNT) continue;
        if (b.first <= 1 || b.second >= INT32_MAX - 1) 
        {
            one_end_berths.push_back(b);
            one_end_weights.push_back(weight);
        }
        else
        {
            berths.push_back(b);
            weights.push_back(weight);
        }
    }
    return 0;
}

// get closest site->first in ssmap to val; if gt HITBERTH DISTANCE, return end()
// return the closest site iterator
map<int32_t, int>::const_iterator berth::closest_site(const map<int32_t, int>* ssmap, int val)
{
    if (ssmap->size() == 0) return ssmap->end();
    
    auto closest = ssmap->end();

    auto it = ssmap->lower_bound(val);
    if (it == ssmap->begin()) 
    {
        closest = it;
    }
    else if (it == ssmap->end()) 
    {
        closest = std::prev(it);
    }
    else    // it is in the middle
    {
        auto pre = std::prev(it);
        if (abs(val - pre->first) < abs(it->first - val)) 
        {
            closest = pre;
        }
        else 
        {
            closest =  it;
        }
    }

    if(abs(val - closest->first) < CLOSESITE_DISTANCE) 
    {
        return closest;
    }
    return ssmap->end();
}

// refine. edit berth in-place 
// At most one berth between each pair of ss
int berth::refine_berths_by_splice_sites(vector<int>& ss)
{
    //TODO:
    return 0;
}

/* Write bed format
1. chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
2. chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
3. chromEnd - End position of the feature in standard chromosomal coordinates
4. name - Name of the line in the BED file
5. score - A score between 0 and 1000. 
6. strand - defined as + (forward) or - (reverse) or '.'
7. thickStart - coordinate at which to start drawing the feature as a solid rectangle
8. thickEnd - coordinate at which to stop drawing the feature as a solid rectangle
9. itemRgb - an RGB colour value (e.g. 0,0,255). Only used if there is a track line with the value of itemRgb set to "on" (case-insensitive).
10. blockCount - the number of sub-elements (e.g. exons) within the feature
11. blockSizes - List of values separated by commas corresponding to the size of the blocks (the number of values must correspond to that of the "blockCount")
12. blockStarts - List of values separated by commas corresponding to the starting coordinates of the blocks, coordinates calculated relative to those present in the chromStart column (the number of values must correspond to that of the "blockCount")
*/
int berth::write_bed(string filename = "berth.bed") const
{
    ofstream bed_file(filename, ios::app);
    if (!bed_file.is_open())
    {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return 0;
    }

    assert (berths.size() == weights.size());
    for (int i = 0; i < berths.size(); i++)
    {
        const auto &berth = berths[i];
        const auto &weight = weights[i];
        int st = berth.first;
        int ed = berth.second;
        string name = "berth" + to_string(i) + "_" + to_string(st) + "-" + to_string(ed) + "_w" + to_string(weight);
        bed_file << bb.chrm << "\t";
        bed_file << st << "\t";
        bed_file << ed << "\t";
        bed_file << name << "\t";       // 4. name
        bed_file << weight << "\t";     // 5. score
        bed_file << bb.strand << "\t";  // 6. strand
        // bed_file << st << "\t";         // 7. thickStart
        // bed_file << ed << "\t";         // 8. thickEnd
        // bed_file << "0,0,255" << "\t";  // 9. itemRgb
        // bed_file << "0" << "\t";        // 10. blockCount 
        // bed_file << "\t";               // 11. blockSizes
        // bed_file << "\t";               // 12. blockStarts 
        bed_file << "\n";
    }
    bed_file.close();
    return 0;
}

int berth::write_bed_one_end(string filename = "berth_one_end.bed") const
{
    ofstream bed_file(filename, ios::app);
    if (!bed_file.is_open())
    {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return 0;
    }

    for (const auto &berth : one_end_berths)
    {
        int st = berth.first;
        int ed = berth.second;
        if (st <= 1 && ed >= INT32_MAX-1) continue;
        if (st <= 1)                      st = ed -1;
        if (ed >= INT32_MAX-1)            ed = st + -1;
        bed_file << bb.chrm << "\t" << st << "\t" << ed << "\n";
    }
    bed_file.close();
    return 0;
}

int berth::print(int index) const
{
    printf("Berth %d (size %d, status) :  \n", index, berths.size(), status);

    assert (berths.size() == weights.size());
    for (int i = 0; i < berths.size(); i++)
    {
        printf("\tberth [%d, %d], w = %d \n", berths[i].first, berths[i].second, weights[i]);
    }
    
    printf("berththit: \n");
    for (const auto& bh: berth2hit)
    {
        printf("\tberth [%d, %d]: ", bh.first.first, bh.first.second);
        printv(bh.second);
        printf("\n");
    }
    
    if (is_empty()) 
    {
        printf("Empty berth printing ssc etc....\n");
        printf("ssc: \n");
        for (const auto& [key, value] : ssc) std::cout << "\t" << key << ": " << value << '\n';
        
        printf("slc: \n");
        for (const auto& [key, value] : slc) std::cout << "\t" << key << ": " << value << '\n';
        
        printf("ttc: \n");
        for (const auto& [key, value] : ttc) std::cout << "\t" << key << ": " << value << '\n';
        
        printf("trc: \n");
        for (const auto& [key, value] : trc) std::cout << "\t" << key << ": " << value << '\n';
    }
    return 0;
}

