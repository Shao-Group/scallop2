/*
Part of Berth
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BERTH_H__
#define __BERTH_H__

#include <vector>
#include <map>
#include <random>
#include "hit.h"
#include "util.h"
#include "bundle_base.h"

using namespace std;

// Gaussian component for mixture model
struct GaussianComponent {
    double mean;
    double variance;
    double weight;
    double log_likelihood;
};

// Directional coverage data structure
struct DirectionalCoverage {
    map<int32_t, int> forward_coverage;  // coverage in transcription direction
    map<int32_t, int> reverse_coverage;  // coverage opposite to transcription
    map<int32_t, double> directional_score; // forward - reverse normalized
};

class berth
{
public:
	berth(bundle_base &b);
	int  clear();
	bool is_empty() const;
	int build_berths();
	int refine_berths(vector<int>& ss);

public:
	const vector<PI32>& get_berths() const;
	const vector<int>&  get_weights() const;
	map<int32_t, int>   get_berth_side(int side) const;
	int write_bed(string filename) const;
	int write_bed_one_end(string filename) const;
	int print(int index = 0) const;

private:	
	int status;          // 0: unbuilt, 1: built empty, 2: built populated, 3: refined populated
	const bundle_base &bb;
	const vector<hit> &hits;
	map<int32_t, int> ssc;  // <pos, count> of TSS 
	map<int32_t, int> slc;  // <pos, count> of left range of TSS
	map<int32_t, int> ttc;  // <pos, count> of TTS
	map<int32_t, int> trc;  // <pos, count> of right range of TTS
	map<PI32, vector<int>> berth2hit;  // <TSS, TES> -> hit index. may incld illegal berths

	vector<PI32> berths;          // a vector of legal <TSS, TES> pairs
	vector<int>  weights;         // a vector of weights of the berths
	vector<PI32> one_end_berths;  // a vector of one-end berth <TSS, TES>, the other end marked -1 or INT32_MAX
	vector<int>  one_end_weights; // a vector of weights of one-end berths

	// Enhanced data structures for directional analysis and GMM
	DirectionalCoverage tss_directional_cov;
	DirectionalCoverage tes_directional_cov;
	vector<GaussianComponent> tss_gmm_components;
	vector<GaussianComponent> tes_gmm_components;

	// default parameters
	const int SLIDINGD_WINDOW_SIZE = 100;
	const int ISOLATED_DISTANCE    = 2000;
	const int CLOSESITE_DISTANCE   = 250;
	const int MIN_BOUND_COUNT      = 5;
	const int MIN_BERTH_HIT_COUNT  = 3;
	
	// New parameters for enhanced analysis
	const int DIRECTIONAL_WINDOW_SIZE = 50;
	const double MIN_DIRECTIONAL_SCORE = 0.3;
	const int GMM_MAX_COMPONENTS = 10;
	const int GMM_MAX_ITERATIONS = 100;
	const double GMM_CONVERGENCE_THRESHOLD = 1e-6;

private:
	int init_sites();
	int pick_peaks();
	
	// Enhanced methods
	int init_directional_coverage();
	int pick_peaks_gmm();
	int filter_window(map<int32_t, int> * v, bool dir);	// dir: 0 for ss, 1 for tt
	map<int32_t, int> find_local_or_isolated_max(const map<int32_t, int>& peakc);
	
	// Directional coverage analysis methods
	DirectionalCoverage compute_directional_coverage(const map<int32_t, int>& site_counts, bool is_tss);
	double calculate_directional_score(int forward_count, int reverse_count);
	map<int32_t, int> filter_by_directionality(const map<int32_t, int>& sites, 
	                                           const DirectionalCoverage& dir_cov);
	
	// GMM peak calling methods
	vector<GaussianComponent> fit_gaussian_mixture(const map<int32_t, int>& coverage_data, int max_components);
	vector<int32_t> extract_peaks_from_gmm(const vector<GaussianComponent>& components, 
	                                       const map<int32_t, int>& original_data);
	double gaussian_pdf(double x, double mean, double variance);
	double calculate_bic(const vector<GaussianComponent>& components, const map<int32_t, int>& data);
	void expectation_maximization(vector<GaussianComponent>& components, 
	                             const vector<pair<int32_t, int>>& data_points);
	
	int assign_hit_berth();
	map<int32_t, int>::const_iterator closest_site(const map<int32_t, int>* ssmap, int val);
	int refine_berths_by_splice_sites(vector<int>& ss); 	// refine. At most one berth between each pair of ss
};

#endif
