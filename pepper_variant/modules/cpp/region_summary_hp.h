//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_REGION_SUMMARY_HP_H
#define PEPPER_PRIVATE_REGION_SUMMARY_HP_H

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>

namespace ImageOptionsRegionHP {
    static constexpr int MAX_COLOR_VALUE = 125;
    static constexpr int MIN_COLOR_VALUE = -125;
    static constexpr int MISMATCH_COLOR_START = 128;
    static constexpr int REFERENCE_INDEX_START = 0;
    static constexpr int REFERENCE_INDEX_SIZE = 0;
    static constexpr int SUPPORT_INDEX_SIZE = 0;
    static constexpr int BASE_INDEX_START = 11;
    static constexpr int BASE_INDEX_SIZE = 14;
    vector<string> column_values{"REFB1:", //0
                                 "SNPS1:",
                                 "INSS1:",
                                 "DELS1:",
                                 "REFF1:",
                                 "SNPF1:",
                                 "INSF1:",
                                 "DELF1:", //7
                                 "AFRW1:",
                                 "CFRW1:",
                                 "GFRW1:",
                                 "TFRW1:",
                                 "IFRW1:",
                                 "DFRW1:",
                                 "*FRW1:",
                                 "REFR1:",
                                 "SNPR1:",
                                 "INSR1:",
                                 "DELR1:", //18
                                 "AREV1:",
                                 "CREV1:",
                                 "GREV1:",
                                 "TREV1:",
                                 "IREV1:",
                                 "DREV1:",
                                 "*REV1:",
                                 "REFF2:",
                                 "SNPF2:",
                                 "INSF2:",
                                 "DELF2:", //29
                                 "AFRW2:",
                                 "CFRW2:",
                                 "GFRW2:",
                                 "TFRW2:",
                                 "IFRW2:",
                                 "DFRW2:",
                                 "*FRW2:",
                                 "REFR2:",
                                 "SNPR2:",
                                 "INSR2:",
                                 "DELR2:", //40
                                 "AREV2:",
                                 "CREV2:",
                                 "GREV2:",
                                 "TREV2:",
                                 "IREV2:",
                                 "DREV2:",
                                 "*REV2:"};

    static constexpr bool GENERATE_INDELS = false;
};

struct type_truth_recordHP {
    string contig;
    long long pos_start;
    long long pos_end;
    string ref;
    string alt;

    type_truth_recordHP(string contig, long long pos, long long pos_end, string ref, string alt) {
        this->contig = std::move(contig);
        this->pos_start = pos;
        this->pos_end = pos_end;
        this->ref = std::move(ref);
        this->alt = std::move(alt);
    }
};


namespace VariantTypesHP {

    static constexpr int HOM_REF = 0;
    static constexpr int SNP = 1;
    static constexpr int INSERT = 2;
    static constexpr int DELETE = 3;
};

struct CandidateImageSummaryHP {
    string contig;
    long long position;
    vector< vector<int> > image_matrix;
    vector <string> candidates;
    vector <int> candidate_frequency;
    int depth;
    uint8_t base_label;
    uint8_t type_label;
    CandidateImageSummaryHP() {
    }

    // string &, long long &, int &, vector<string> &, vector<int> &, vector<vector<int> > &, int &, int &
    CandidateImageSummaryHP(string contig, long long position, int depth, vector<string> candidates, vector<int> candidate_frequency, vector< vector<int> > image, int base_label, int type_label) {
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

class RegionalSummaryGeneratorHP {
    string contig;
    long long ref_start;
    long long ref_end;
    string reference_sequence;
    vector< vector<type_truth_recordHP> > hp1_truth_alleles;
    vector< vector<type_truth_recordHP> > hp2_truth_alleles;
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

    RegionalSummaryGeneratorHP(string contig, long long region_start, long long region_end, string reference_sequence);

    void generate_max_insert_observed(const type_read& read);

    void generate_max_insert_summary(vector <type_read> &reads);

    static int get_reference_feature_index(char base);

    void encode_reference_bases(vector< vector<int> >& image_matrix);

    void generate_labels(const vector<type_truth_recordHP>& hap1_records, const vector<type_truth_recordHP>& hap2_records);

    void populate_summary_matrix(vector< vector<int> >& image_matrix,
                                 int *coverage_vector,
                                 int *snp_count,
                                 int *insert_count,
                                 int *delete_count,
                                 vector< map<string, int> > &AlleleFrequencyMap,
                                 vector< map<string, int> > &AlleleFrequencyMapFwdStrandHP1,
                                 vector< map<string, int> > &AlleleFrequencyMapFwdStrandHP2,
                                 vector< map<string, int> > &AlleleFrequencyMapRevStrandHP1,
                                 vector< map<string, int> > &AlleleFrequencyMapRevStrandHP2,
                                 vector< set<string> > &AlleleMap,
                                 type_read read,
                                 double min_snp_baseq,
                                 double min_indel_baseq);


    static int get_feature_index(char ref_base, char base, bool is_reverse, int hp_tag);

    void debug_print_matrix_hp(vector<vector<int> > image_matrix, bool train_mode);

    void debug_candidate_summary_hp(CandidateImageSummaryHP candidate, int small_chunk_size, bool train_mode);

    vector<CandidateImageSummaryHP> generate_summary(vector <type_read> &reads,
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
