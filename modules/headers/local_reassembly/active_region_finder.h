//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef HELEN_ACTIVE_REGION_FINDER_H
#define HELEN_ACTIVE_REGION_FINDER_H

#include "../../headers/dataio/bam_handler.h"
#include "../../headers/dataio/fasta_handler.h"

#include <vector>
#include <algorithm>

using namespace std;

namespace ActiveRegionFinder_options {
    static constexpr float bias = -0.7;
    static constexpr float w_match = -0.06;
    static constexpr float w_mismatch = -0.08;
    static constexpr float w_insert = 2.5;
    static constexpr float w_delete = 1.8;
    static constexpr float w_soft_clip = 3.0;
    static constexpr float threshold = 3.0;

    static constexpr int min_region_size = 80;
    static constexpr int min_mapping_quality = 20;
    static constexpr int min_base_quality = 20;
};

class ActiveRegionFinder {
    long long region_start;
    long long region_end;
    string chromosome_name;
    string reference_sequence;
public:
    ActiveRegionFinder(string reference_sequence,
                       string chromosome_name,
                       long long region_start,
                       long long region_end);
    void update_weights(type_read &read, vector<float> &weight_vector);
    vector< pair<long long, long long> > find_active_region(vector<type_read> reads);
    ~ActiveRegionFinder();
};
#endif //FRIDAY_ACTIVE_REGION_FINDER_H
