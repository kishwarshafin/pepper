//
// Created by Kishwar Shafin on 10/24/18.
//

#ifndef PEPPER_VARIANT_CANDIDATE_FINDER_HP_H
#define PEPPER_VARIANT_CANDIDATE_FINDER_HP_H

#include <cmath>
#include <numeric>
#include <utility>
#include <vector>
#include "bam_handler.h"
#include "candidate_finder.h"
using namespace std;

namespace CandidateFinderHP_options {
    static constexpr int min_mapping_quality = 1;
    static constexpr int min_base_quality = 1;
    static constexpr int freq_threshold = 5;
    static constexpr int min_count_threshold = 2;
};

namespace ONTLinearRegression {
    static constexpr double SNP_ALT_FREQ_COEF = 0;
    static constexpr double SNP_NON_REF_PROB_COEF = -0.002397;
    static constexpr double SNP_ALLELE_WEIGHT_COEF = 1.008378;
    static constexpr double SNP_BIAS_TERM = 0.001291;
    static constexpr double SNP_THRESHOLD = 0.01;
    static constexpr double SNP_LOWER_FREQ_THRESHOLD = 0.10;
    static constexpr double SNP_UPPER_FREQ = 0.4;

    static constexpr double INSERT_ALT_FREQ_COEF = 0;
    static constexpr double INSERT_NON_REF_PROB_COEF = 0.239488;
    static constexpr double INSERT_ALLELE_WEIGHT_COEF = 0.822283;
    static constexpr double INSERT_BIAS_TERM = 0.000298;
    static constexpr double INSERT_THRESHOLD = 0.2;
    static constexpr double IN_LOWER_FREQ_THRESHOLD = 0.10;
    static constexpr double IN_UPPER_FREQ = 0.5;

    static constexpr double DELETE_ALT_FREQ_COEF = 0;
    static constexpr double DELETE_NON_REF_PROB_COEF = 0.039434;
    static constexpr double DELETE_ALLELE_WEIGHT_COEF = 0.765909;
    static constexpr double DELETE_BIAS_TERM = -0.003304;
    static constexpr double DELETE_THRESHOLD = 0.15;
    static constexpr double DEL_LOWER_FREQ_THRESHOLD = 0.10;
    static constexpr double DEL_UPPER_FREQ_THRESHOLD = 0.5;
}

class CandidateFinderHP {
    long long region_start;
    long long region_end;
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<Candidate, int> AlleleFrequencyMap;
    vector< set<Candidate> > AlleleMap;
public:
    CandidateFinderHP(string reference_sequence,
                      string chromosome_name,
                      long long region_start,
                      long long region_end,
                      long long ref_start,
                      long long ref_end);
    void add_read_alleles(type_read &read, vector<int> &coverage);
    vector<PositionalCandidateRecord> find_candidates(vector <type_read>& reads,
                                                      vector<long long> positions,
                                                      vector<int>indices,
                                                      const vector< vector<int> >& base_predictions_h1,
                                                      const vector< vector<int> >& base_predictions_h2,
                                                      bool freq_based,
                                                      double freq);
    static bool filter_candidate(const Candidate& candidate, bool freq_based, double freq);
};



#endif //PEPPER_VARIANT_CANDIDATE_FINDER_H
