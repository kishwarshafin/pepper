//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_REGION_SUMMARY_H
#define PEPPER_PRIVATE_REGION_SUMMARY_H

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ImageOptionsRegion {

    static constexpr int MAX_COLOR_VALUE = 256;
    static constexpr int MAX_REF_COLOR_VALUE = 100;
    static constexpr int MISMATCH_COLOR_START = 128;
    static constexpr int REFERENCE_INDEX_START = 0;
    static constexpr int REFERENCE_INDEX_SIZE = 5;
    static constexpr int BASE_INDEX_START = 5;
    static constexpr int BASE_INDEX_SIZE = 14;
    vector<string> column_values{"AREF:",
                                 "CREF:",
                                 "GREF:",
                                 "TREF:",
                                 "*REF:",
                                 "AFRW:",
                                 "CFRW:",
                                 "GFRW:",
                                 "TFRW:",
                                 "*FRW:",
                                 "IFRW:",
                                 "DFRW:",
                                 "AREV:",
                                 "CREV:",
                                 "GREV:",
                                 "TREV:",
                                 "*REV:",
                                 "IREV:",
                                 "DREV:"};

    static constexpr bool GENERATE_INDELS = true;
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

class RegionalSummaryGenerator {
    long long ref_start;
    long long ref_end;
    string reference_sequence;
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

    RegionalSummaryGenerator(long long region_start, long long region_end, string reference_sequence);

    void generate_max_insert_observed(const type_read& read);

    void generate_max_insert_summary(vector <type_read> &reads);

    void generate_labels_from_truth_read(type_read read, int hp_tag);

    static int get_reference_feature_index(char base);

    void encode_reference_bases(int **image_matrix);

    void generate_labels(const type_read& truth_read_hp1, const type_read& truth_read_hp2);

    void populate_summary_matrix(int **image_matrix,
                                 int *coverage_vector,
                                 type_read read);

    static int get_feature_index(char base, bool is_reverse);

    void debug_print_matrix(int** image_matrix, bool train_mode);

    RegionalImageSummary generate_summary(vector <type_read> &reads,
                                          int chunk_overlap,
                                          int smaller_chunk_size,
                                          int feature_size,
                                          int chunk_id_start,
                                          bool train_mode);
};

#endif //PEPPER_PRIVATE_REGION_SUMMARY_H
