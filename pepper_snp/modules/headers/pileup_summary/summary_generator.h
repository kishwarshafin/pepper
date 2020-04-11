//
// Created by Kishwar Shafin on 25/2/19.
//

#ifndef HELEN_SUMMARY_GENERATOR_H
#define HELEN_SUMMARY_GENERATOR_H

#include <math.h>
#include <algorithm>
#include <iomanip>
#include <assert.h>
using namespace std;
#include "../dataio/bam_handler.h"


namespace ImageOptions {
    static constexpr int MAX_COLOR_VALUE = 254;
    static constexpr bool GENERATE_INSERTS = false;
};

class SummaryGenerator {
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;

    map< pair< pair<long long, int>, int>, double> insert_summaries;
    map< pair<long long, int>, double> base_summaries;
    map<long long, long long> longest_insert_count;
    map<long long, double> coverage;

    map< pair<long long, int>, char> insert_labels_h1;
    map< pair<long long, int>, char> insert_labels_h2;
    map< long long, char> base_labels_h1;
    map< long long, char> base_labels_h2;
public:
    vector< vector<uint8_t> > image;
    vector<uint8_t> labels;
    vector<int> ref_image;
    vector<pair<long long, int> > genomic_pos;
    vector<int> bad_label_positions;

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
                                vector<type_read> &truth_reads_h1,
                                vector<type_read> &truth_reads_h2);

    void iterate_over_read(type_read read, long long region_start, long long region_end);
    int get_sequence_length(long long start_pos, long long end_pos);
    void generate_labels(type_read truth_read, long long region_start, long long region_end, int hp_tag);
    void generate_ref_features();
    void debug_print(long long start_pos, long long end_pos);
    void generate_image(long long start_pos, long long end_pos);
};


#endif //HELEN_SUMMARY_GENERATOR_H
