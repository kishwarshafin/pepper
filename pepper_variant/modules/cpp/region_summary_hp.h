//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_REGION_SUMMARY_HP_H
#define PEPPER_PRIVATE_REGION_SUMMARY_HP_H

#include <cstring>
#include <iostream>
#include <algorithm>

namespace ImageOptionsRegionHP {
    static constexpr int MAX_COLOR_VALUE = 254;
    static constexpr bool GENERATE_INDELS = true;
};


struct RegionalImageSummaryHP {
    vector< vector< vector<uint8_t> > > chunked_image_matrix_hp1;
    vector< vector< vector<uint8_t> > > chunked_image_matrix_hp2;
    vector< vector<int64_t> > chunked_positions;
    vector< vector<int32_t> > chunked_index;
    vector< vector<uint8_t> > chunked_labels_hp1;
    vector< vector<uint8_t> > chunked_labels_hp2;
    vector<int> chunked_ids;
};

class RegionalSummaryGeneratorHP {
    long long ref_start;
    long long ref_end;
    string reference_sequence;
    vector<char> labels_hp1;
    vector<char> labels_hp2;
public:
    vector<uint16_t> labels_hp1_report;
    vector<uint16_t> labels_hp2_report;
    vector<uint64_t> max_observed_insert;
    vector<uint64_t> positions;
    vector<uint32_t> index;
    vector<uint64_t> cumulative_observed_insert;
    uint64_t total_observered_insert_bases;

    RegionalSummaryGeneratorHP(long long region_start, long long region_end, string reference_sequence);

    void generate_max_insert_observed(const type_read& read);

    void generate_max_insert_summary(vector <type_read> &reads);

    void generate_labels_from_truth_read(type_read read, int hp_tag);

    void generate_labels(const type_read& truth_read_hp1, const type_read& truth_read_hp2);

    void populate_summary_matrix(int **image_matrix_hp1,
                                 int **image_matrix_hp2,
                                 int *coverage_hp1,
                                 int *coverage_hp2,
                                 int *coverage,
                                 type_read read);

    static int get_feature_index(char base, bool is_reverse);

    void debug_print_matrix(int** image_matrix_hp1, int** image_matrix_hp2, bool train_mode);

    RegionalImageSummaryHP generate_summary(vector <type_read> &reads,
                                            int chunk_overlap,
                                            int smaller_chunk_size,
                                            int feature_size,
                                            int chunk_id_start,
                                            bool train_mode);
};

#endif //PEPPER_PRIVATE_REGION_SUMMARY_H
