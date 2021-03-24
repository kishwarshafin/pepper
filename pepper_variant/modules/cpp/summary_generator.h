//
// Created by Kishwar Shafin on 25/2/19.
//

#ifndef PEPPER_VARIANT_SUMMARY_GENERATOR_H
#define PEPPER_VARIANT_SUMMARY_GENERATOR_H

#include <math.h>
#include <algorithm>
#include <iomanip>
#include <assert.h>
using namespace std;
#include "bam_handler.h"


namespace ImageOptions {
    static constexpr int MAX_COLOR_VALUE = 254;
};

struct ImageSummary {
    vector< vector< vector<uint8_t> > > images;
    vector< vector<pair<long long, int> > > positions;
    vector< vector<uint8_t> > refs;
    vector< vector<uint8_t> > labels;
    vector<int> chunk_ids;
};

class SummaryGenerator {
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;

    map< pair< pair<long long, int>, int>, double> insert_summaries;
    map< pair<long long, int>, double> base_summaries;

    map<long long, double> coverage;

    map< pair<long long, int>, char> insert_labels_hp1;
    map< pair<long long, int>, char> insert_labels_hp2;
    map< long long, char> base_labels_hp1;
    map< long long, char> base_labels_hp2;

public:
    map<long long, long long> longest_insert_count;

    vector< vector<uint8_t> > image;
    vector<uint8_t> labels;

    vector<pair<long long, int> > genomic_pos;
    vector<int> bad_label_positions;
    vector<uint8_t> ref_image;

    SummaryGenerator(string reference_sequence,
                     string chromosome_name,
                     long long ref_start,
                     long long ref_end);
    void generate_summary(vector <type_read> &reads,
                          long long start_pos,
                          long long end_pos);

    void generate_train_summary(vector <type_read> &reads,
                                long long start_pos,
                                long long end_pos,
                                type_read truth_read_hp1,
                                type_read truth_read_hp2);

    void iterate_over_read(type_read read, long long region_start, long long region_end);
    int get_sequence_length(long long start_pos, long long end_pos);
    void generate_labels(type_read read, long long region_start, long long region_end, int hp_tag);
    void generate_ref_features();
    void debug_print(long long start_pos, long long end_pos);
    void generate_image(long long start_pos, long long end_pos);
    ImageSummary chunk_image(int chunk_size, int chunk_overlap, int image_height);
    ImageSummary chunk_image_train(int chunk_size, int chunk_overlap, int image_height, int chunk_id_start);
};


#endif //PEPPER_VARIANT_SUMMARY_GENERATOR_H
