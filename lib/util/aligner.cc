/*
Part of Berth
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and Pennsylvania State University.
See LICENSE for licensing.
*/

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <unordered_map>
#include "aligner.h"
#include <climits>


/**
 * Assumes seq1 is shorter than seq2
 * @param nm number of maximal allowed mismatches
 * @param indel_penalty 
 * @param mis_penalty
 * @return position of seq2's first subsequence that can be transformed into seq1
 * @return -1 if unfeasible (or seq1 is longer than seq2)
 */
pair<int, int> subseq_pos(const std::string& seq1, const std::string& seq2, int nm, int indel_penalty, int mis_penalty) 
{
    if (seq1.length() > seq2.length())  return make_pair(-1, -1);

    if (! validate_dna_seq(seq1, seq2)) 
    {
        cerr << "WARNING:\t sequenes have non-ATCG bases. Proceed anyway." << endl;
    }

    int end_pos = -1;
    int min_penalty = INT_MAX;

    // Initialize
    vector<vector<int>> mx(seq1.length() + 1, vector<int>(seq2.length() + 1, 0));
    mx[0][0] = 0;
    for (size_t i = 0; i < seq2.length() + 1; ++i) mx[0][i] = 0;
    for (size_t i = 1; i < seq1.length() + 1; ++i) mx[i][0] = i * indel_penalty;

    // DP body
    for (size_t j = 1; j <= seq2.length(); ++j) 
    {
        for (size_t i = 1; i <= seq1.length(); ++i) 
        {
            int examine_match = (seq1[i - 1] == seq2[j - 1]) ? 0 : mis_penalty;
            int consume_both = mx[i-1][j-1] + examine_match;
            int consume_seq1 = mx[i][j-1]   + indel_penalty;
            int consume_seq2 = mx[i-1][j]   + indel_penalty;
            
            mx[i][j] = min({consume_both, consume_seq1, consume_seq2});
        }
        
        // Check if we found a match within allowed edit distance
        // if (mx[seq1.length()][j] <= nm) return j; //FIXME: return highest score position
        if (mx[seq1.length()][j] <= nm && mx[seq1.length()][j] < min_penalty)  //FIXED: saving min dist position
        {
            end_pos = j;
            min_penalty = mx[seq1.length()][j];
        }
    }

    if(end_pos == -1) return make_pair(-1,-1);
    int start_pos = end_pos;
    int i = seq1.length();
    while( i > 1)
    {
        if(mx[i-1][start_pos-1] < mx[i-1][start_pos] && mx[i-1][start_pos-1] < mx[i][start_pos-1]) 
        {
            start_pos = start_pos - 1;
            i--;
        }
        else if (mx[i-1][start_pos] < mx[i][start_pos-1]) i--;
        else start_pos--;
        if (start_pos < 1) cout << seq1 << " " << seq2 << endl;
        assert(start_pos >= 1);

    }

    return make_pair(start_pos-1, end_pos-1);

}

// returns the length of polyT with a local search
// truncate when the ratio of T's in the sequence is less than polyT_ratio
int polyT(const std::string& seq1, int start, double polyT_ratio, int min_polyT_len) 
{
    if (seq1.length() < min_polyT_len) return 0;

    int count = 0;
    int lastT = start;
    for (size_t i = start; i < seq1.length(); ++i) 
    {
        cout << seq1[i];
        if (seq1[i] == 'T' || seq1[i] == 't') 
        {
            count++;
            lastT = i;
        }
        else if (count <= min_polyT_len || (count / (double)(i - start + 1) >= polyT_ratio)) continue;
        else break;
    }
    cout << endl;

    if (count >= min_polyT_len) return lastT - start + 1; // not counts of T but length of T tail
    else return 0;
}


// Reverse complement, convert to uppercase
char revcomp_char(const char& c)
{
    static const unordered_map<char, char> REV_COMP_TABLE = { 
        {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'},
        {'U', 'A'}, {'R', 'Y'}, {'Y', 'R'}, {'S', 'W'},
        {'W', 'S'}, {'K', 'M'}, {'M', 'K'}, {'B', 'V'},
        {'D', 'H'}, {'H', 'D'}, {'V', 'B'}, {'N', 'N'},
        {'*', '*'} 
    };
    char rc = '*';
    auto it = REV_COMP_TABLE.find(toupper(c));
    if (it == REV_COMP_TABLE.end()) throw runtime_error("Invalid character in sequence: " + string(1, c));
    else rc = it->second;
    return rc;
}

/**
 * Returns the reverse complement of a DNA/RNA sequence
 * @param seq Input sequence
 * @return Reverse complement sequence
 */
std::string revcomp(const string& seq) 
{
    string rc_seq = seq;
    
    for (char& c : rc_seq) c = revcomp_char(c);
    
    reverse(rc_seq.begin(), rc_seq.end());
    
    return rc_seq;
}