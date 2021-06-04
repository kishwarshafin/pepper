//
// Created by Kishwar Shafin on 10/24/18.
//

#ifndef PEPPER_VARIANT_CANDIDATE_FINDER_H
#define PEPPER_VARIANT_CANDIDATE_FINDER_H

#include <cmath>
#include <numeric>
#include <utility>
#include <vector>
#include "bam_handler.h"
using namespace std;

namespace CandidateFinder_options {
    static constexpr int min_mapping_quality = 1;
    static constexpr int min_base_quality = 0;
    static constexpr int freq_threshold = 4;
    static constexpr int min_count_threshold = 2;
    static constexpr bool report_indels = true;
};

namespace AlleleType {
    static constexpr int SNP_ALLELE = 1;
    static constexpr int INSERT_ALLELE = 2;
    static constexpr int DELETE_ALLELE = 3;
};

namespace Genotype {
    static constexpr int HOM = 0;
    static constexpr int HET = 1;
    static constexpr int HOM_ALT = 2;
};

namespace LinearRegression {
    static constexpr double SNP_ALT_FREQ_COEF = 0;
    static constexpr double SNP_NON_REF_PROB_COEF = 0.145468;
    static constexpr double SNP_ALLELE_WEIGHT_COEF = 1.049319;
    static constexpr double SNP_ALT_PROB1_COEF = 0.01494;
    static constexpr double SNP_ALT_PROB2_COEF = 0.014063;
    static constexpr double SNP_BIAS_TERM = -0.230597;
    static constexpr double SNP_THRESHOLD = 0.1;
    static constexpr double SNP_LOWER_FREQ_THRESHOLD = 0.10;
    static constexpr double SNP_UPPER_FREQ = 0.4;

    static constexpr double INSERT_ALT_FREQ_COEF = 0;
    static constexpr double INSERT_NON_REF_PROB_COEF = 0;
    static constexpr double INSERT_ALLELE_WEIGHT_COEF = 0;
    static constexpr double INSERT_BIAS_TERM = 0;
    static constexpr double INSERT_THRESHOLD = 0;
    static constexpr double IN_LOWER_FREQ_THRESHOLD = 0;
    static constexpr double IN_UPPER_FREQ = 0;

    static constexpr double DELETE_ALT_FREQ_COEF = 0;
    static constexpr double DELETE_NON_REF_PROB_COEF = 0;
    static constexpr double DELETE_ALLELE_WEIGHT_COEF = 0;
    static constexpr double DELETE_BIAS_TERM = 0;
    static constexpr double DELETE_THRESHOLD = 0;
    static constexpr double DEL_LOWER_FREQ_THRESHOLD = 0;
    static constexpr double DEL_UPPER_FREQ_THRESHOLD = 0;
}

struct CandidateAllele{
    string ref;
    string alt;
    int alt_type;

    CandidateAllele(string ref, string alt, int alt_type) {
        this->ref = std::move(ref);
        this->alt = std::move(alt);
        this->alt_type = alt_type;
    }

    CandidateAllele() {}
};

struct Candidate {
    long long pos;
    long long pos_end;
    CandidateAllele allele;
    int genotype;
    int depth;
    int read_support;
    double allele_probability;
    double genotype_probability;

    // will be deprecated
    double alt_prob;
    double alt_prob_h1;
    double alt_prob_h2;
    double non_ref_prob;

    Candidate() {
        this->genotype = 0;
    }

    void set_genotype(int genotype_) {
        this->genotype = genotype_;
    }

    void set_depth_values(int depth_, int read_support_) {
        this->depth = depth_;
        this->read_support = read_support_;
    }

    Candidate(long long pos_start, long long pos_end, string ref_, string alt_, int alt_type_) {
        this->pos = pos_start;
        this->pos_end = pos_end;
        this->allele.alt = std::move(alt_);
        this->allele.ref = std::move(ref_);
        this->allele.alt_type = alt_type_;
        this->genotype = 0;
    }
    bool operator< (const Candidate& that ) const {
        if(this->pos != that.pos) return this->pos < that.pos;
        if(this->pos_end != that.pos_end) return this->pos_end < that.pos_end;
        if(this->allele.alt != that.allele.alt) return this->allele.alt < that.allele.alt;
        if(this->allele.ref != that.allele.ref) return this->allele.ref < that.allele.ref;
        if(this->allele.alt_type != that.allele.alt_type) return this->allele.alt_type < that.allele.alt_type;
        return this->pos < that.pos;
    }

    bool operator==(const Candidate& that ) const {
        if(this->pos == that.pos &&
           this->pos_end == that.pos_end &&
           this->allele.ref == that.allele.ref &&
           this->allele.alt == that.allele.alt &&
           this->allele.alt_type == that.allele.alt_type)
            return true;
        return false;
    }

    void print() {
        cout<<this->pos<<" "<<this->pos_end<<" "<<this->allele.ref<<" "<<this->allele.alt<<" "<<this->allele.alt_type<<" "<<this->genotype<<endl;
    }

};

struct PositionalCandidateRecord{
    string chromosome_name;
    long long pos_start;
    long long pos_end;
    int depth;
    vector<Candidate> candidates;

    PositionalCandidateRecord(string chromosome_name, long long pos_start, long long pos_end, int depth) {
        this->chromosome_name = std::move(chromosome_name);
        this->pos_start = pos_start;
        this->pos_end = pos_end;
        this->depth = depth;
    }
    PositionalCandidateRecord() {
    }
    bool operator< (const PositionalCandidateRecord& that ) const {
        if(this->chromosome_name != that.chromosome_name) return this->chromosome_name < that.chromosome_name;
        if(this->pos_start != that.pos_start) return this->pos_start < that.pos_start;
        if(this->pos_end != that.pos_end) return this->pos_end < that.pos_end;
        return this->depth > that.depth;
    }
};

class CandidateFinder {
    long long region_start;
    long long region_end;
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<Candidate, int> AlleleFrequencyMap;
    vector< set<Candidate> > AlleleMap;
public:
    CandidateFinder(string reference_sequence,
                    string chromosome_name,
                    long long region_start,
                    long long region_end,
                    long long ref_start,
                    long long ref_end);
    void add_read_alleles(type_read &read, vector<int> &coverage);
    void add_read_alleles_consensus(type_read &read, vector<int> &coverage, vector<int> &insert_count, vector<int> &delete_count, vector<int> &snp_count);
//    vector<PositionalCandidateRecord> find_candidates(vector<type_read>& reads,
//                                                      vector<long long> positions,
//                                                      vector< vector<float> > predictions,
//                                                      vector< vector<float> > type_predictions,
//                                                      vector<int> base_labels,
//                                                      vector<int> type_labels,
//                                                      bool freq_based,
//                                                      double freq);
    vector<PositionalCandidateRecord> find_candidates(vector<type_read>& reads,
                                                      vector<long long> positions,
                                                      vector< vector<float> > predictions,
                                                      vector<int> base_labels,
                                                      bool freq_based,
                                                      double freq);
    vector<long long> find_candidates_consensus(vector <type_read>& reads, double snp_freq_threshold, double insert_freq_threshold, double delete_freq_threshold);
    static bool filter_candidate(const Candidate& candidate, bool freq_based, double freq);
};



#endif //PEPPER_VARIANT_CANDIDATE_FINDER_H
