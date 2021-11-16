//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_REGION_SUMMARY_H
#define PEPPER_PRIVATE_REGION_SUMMARY_H

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>

namespace ImageOptionsRegion {
    static constexpr int MAX_COLOR_VALUE = 125;
    static constexpr int MIN_COLOR_VALUE = -125;
    static constexpr int MISMATCH_COLOR_START = 128;
    static constexpr int REFERENCE_INDEX_START = 0;
    static constexpr int REFERENCE_INDEX_SIZE = 0;
    static constexpr int SUPPORT_INDEX_SIZE = 0;
    static constexpr int BASE_INDEX_START = 11;
    static constexpr int BASE_INDEX_SIZE = 14;
    vector<string> column_values{"REFB:",
                                 "SNPS:",
                                 "INSS:",
                                 "DELS:",
                                 "REFF:",
                                 "SNPF:",
                                 "INSF:",
                                 "DELF:",
                                 "AFRW:",
                                 "CFRW:",
                                 "GFRW:",
                                 "TFRW:",
                                 "IFRW:",
                                 "DFRW:",
                                 "*FRW:",
                                 "REFR:",
                                 "SNPR:",
                                 "INSR:",
                                 "DELR:",
                                 "AREV:",
                                 "CREV:",
                                 "GREV:",
                                 "TREV:",
                                 "IREV:",
                                 "DREV:",
                                 "*REV:"};

    static constexpr bool GENERATE_INDELS = false;
};

struct type_truth_record{
    string contig;
    long long pos_start;
    long long pos_end;
    string ref;
    string alt;

    type_truth_record(string contig, long long pos, long long pos_end, string ref, string alt) {
        this->contig = std::move(contig);
        this->pos_start = pos;
        this->pos_end = pos_end;
        this->ref = std::move(ref);
        this->alt = std::move(alt);
    }
};


namespace VariantTypes {

    static constexpr int HOM_REF = 0;
    static constexpr int SNP = 1;
    static constexpr int INSERT = 2;
    static constexpr int DELETE = 3;
};

struct RegionalImageSummary {
    vector< vector< vector<uint8_t> > > chunked_image_matrix;
    vector< vector<int64_t> > chunked_positions;
    vector< vector<int32_t> > chunked_index;
    vector< vector<uint8_t> > chunked_labels;
    vector< vector<uint8_t> > chunked_type_labels;
    vector<int> chunked_ids;
};


struct CandidateImageSummary {
    string contig;
    long long position;
    vector< vector<int> > image_matrix;
    vector <string> candidates;
    vector <int> candidate_frequency;
    int depth;
    uint8_t base_label;
    uint8_t type_label;
    CandidateImageSummary() {
    }

    // string &, long long &, int &, vector<string> &, vector<int> &, vector<vector<int> > &, int &, int &
    CandidateImageSummary(string contig, long long position, int depth, vector<string> candidates, vector<int> candidate_frequency, vector< vector<int> > image, int base_label, int type_label) {
        this->contig = std::move(contig);
        this->position = position;
        this->image_matrix = std::move(image);
        this->candidates = std::move(candidates);
        this->candidate_frequency = std::move(candidate_frequency);
        this->depth = depth;
        this->base_label = base_label;
        this->type_label = type_label;
    }
};


struct CandidateImagePrediction {
    string contig;
    long long position;
    int depth;
    vector <string> candidates;
    vector <int> candidate_frequency;
    vector<float> prediction_base;
    vector <float> prediction_type;

    CandidateImagePrediction() {
    }

    // string &, long long &, int &, vector<string> &, vector<int> &, vector<vector<int> > &, int &, int &
    CandidateImagePrediction(string contig, long long position, int depth, vector<string> candidates, vector<int> candidate_frequency, const vector<float>& prediction_base, const vector <float>& prediction_type) {
        this->contig = std::move(contig);
        this->position = position;
        this->depth = depth;
        this->candidates = std::move(candidates);
        this->candidate_frequency = std::move(candidate_frequency);
        this->prediction_base = prediction_base;
        this->prediction_type = prediction_type;
    }
};

class RegionalSummaryGenerator {
    string contig;
    long long ref_start;
    long long ref_end;
    string reference_sequence;
    vector< vector<type_truth_record> > hp1_truth_alleles;
    vector< vector<type_truth_record> > hp2_truth_alleles;
    vector<char> labels_hp1;
    vector<char> labels_hp2;
    vector<int> variant_type_labels_hp1;
    vector<int> variant_type_labels_hp2;
    vector<char> ref_at_labels;
public:
    vector<uint16_t> labels;
    vector<uint16_t> labels_variant_type;
    vector<uint64_t> max_observed_insert;
    vector<uint64_t> positions;
    vector<uint32_t> index;
    vector<uint64_t> cumulative_observed_insert;
    uint64_t total_observered_insert_bases;

    RegionalSummaryGenerator(string contig, long long region_start, long long region_end, string reference_sequence);

    void generate_max_insert_observed(const type_read& read);

    void generate_max_insert_summary(vector <type_read> &reads);

    static int get_reference_feature_index(char base);

    void encode_reference_bases(vector< vector<int> >& image_matrix);

    void generate_labels(const vector<type_truth_record>& hap1_records, const vector<type_truth_record>& hap2_records);

    void populate_summary_matrix(vector< vector<int> >& image_matrix,
                                 int *coverage_vector,
                                 int *snp_count,
                                 int *insert_count,
                                 int *delete_count,
                                 vector< map<string, int> > &AlleleFrequencyMap,
                                 vector< map<string, int> > &AlleleFrequencyMapFwdStrand,
                                 vector< map<string, int> > &AlleleFrequencyMapRevStrand,
                                 vector< set<string> > &AlleleMap,
                                 type_read read,
                                 double min_snp_baseq,
                                 double min_indel_baseq);


    static int get_feature_index(char ref_base, char base, bool is_reverse);

    void debug_print_matrix(vector<vector<int> > image_matrix, bool train_mode);

    void debug_candidate_summary(CandidateImageSummary candidate, int small_chunk_size, bool train_mode);

    vector<CandidateImageSummary> generate_summary(vector <type_read> &reads,
                                                   double min_snp_baseq,
                                                   double min_indel_baseq,
                                                   double snp_freq_threshold,
                                                   double insert_freq_threshold,
                                                   double delete_freq_threshold,
                                                   double min_coverage_threshold,
                                                   double snp_candidate_freq_threshold,
                                                   double indel_candidate_freq_threshold,
                                                   double candidate_support_threshold,
                                                   bool skip_indels,
                                                   long long candidate_region_start,
                                                   long long candidate_region_end,
                                                   int candidate_window_size,
                                                   int feature_size,
                                                   bool train_mode);

    int get_reference_feature_value(char base);

};

#endif //PEPPER_PRIVATE_REGION_SUMMARY_H
