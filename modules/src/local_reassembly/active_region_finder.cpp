//
// Created by Kishwar Shafin on 10/18/18.
//

#include "../../headers/local_reassembly/active_region_finder.h"

ActiveRegionFinder::ActiveRegionFinder(string reference_sequence,
                                       string chromosome_name,
                                       long long region_start,
                                       long long region_end) {
    this->reference_sequence = reference_sequence;
    this->region_start = region_start;
    this->region_end = region_end;
    this->chromosome_name = chromosome_name;
}

void ActiveRegionFinder::update_weights(type_read &read, vector<float> &weight_vector) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long w_pos_start = 0;
    long long w_pos_end = 0;
    int base_quality = 0;

    for(auto &cigar: read.cigar_tuples) {
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if(ref_position < region_start) {
                    cigar_index = min(region_start - ref_position, (long long)cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for(int i=cigar_index; i < cigar.length ; i++) {
                    if(ref_position <= region_end &&
                       reference_sequence[ref_position - region_start] == read.sequence[read_index] &&
                       read.base_qualities[ref_position - region_start] >= ActiveRegionFinder_options::min_base_quality)

                        weight_vector[ref_position-region_start] += ActiveRegionFinder_options::w_match;
                    else if(ref_position <= region_end &&
                            read.base_qualities[ref_position - region_start] >= ActiveRegionFinder_options::min_base_quality)
                        weight_vector[ref_position-region_start] += ActiveRegionFinder_options::w_mismatch;

                    read_index += 1;
                    ref_position += 1;
                }
                break;

            case CIGAR_OPERATIONS::IN:
                w_pos_start = max(region_start, ref_position - cigar.length + 1);
                w_pos_end = max(region_start, min(region_end, ref_position + cigar.length));
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                if(base_quality >= ActiveRegionFinder_options::min_base_quality) {
                    for (int pos = w_pos_start; pos <= w_pos_end; pos++) {
                        weight_vector[pos - region_start] += ActiveRegionFinder_options::w_insert;
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::DEL:
                w_pos_start = max(region_start, ref_position + 1);
                w_pos_end = max(region_start, min(region_end, ref_position + cigar.length));

                base_quality = read.base_qualities[max(0, read_index - 1)];
                if(base_quality >= ActiveRegionFinder_options::min_base_quality) {
                    for (int pos = w_pos_start; pos <= w_pos_end; pos++) {
                        weight_vector[pos - region_start] += ActiveRegionFinder_options::w_delete;
                    }
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                w_pos_start = max(region_start, ref_position - cigar.length + 1);
                w_pos_end = max(region_start, min(region_end, ref_position + cigar.length));
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                if(base_quality >= ActiveRegionFinder_options::min_base_quality) {
                    for (int pos = w_pos_start; pos <= w_pos_end; pos++) {
                        weight_vector[pos - region_start] += ActiveRegionFinder_options::w_soft_clip;
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;

        }
    }
}

vector< pair<long long, long long> > ActiveRegionFinder::find_active_region(vector<type_read> reads) {
    vector<float> weight_vector(region_end - region_start + 1, ActiveRegionFinder_options::bias);

    for(auto &read:reads) {
        if(read.mapping_quality < ActiveRegionFinder_options::min_mapping_quality) continue;
        update_weights(read, weight_vector);
    }

    vector<long long> positions;
    // get all the positions that pass the threshold
    for(int i=0; i < weight_vector.size(); i++) {
        if(weight_vector[i] > ActiveRegionFinder_options::threshold) {
            positions.push_back(region_start+i);
        }
    }

    // create windows
    vector< pair<long long, long long> > windows;

    long long start_pos = -1;
    long long end_pos = -1;

    for(auto &pos: positions) {
        if(start_pos == -1) {
            start_pos = pos;
            end_pos = pos;
        } else if(pos > end_pos + ActiveRegionFinder_options::min_region_size) {
            windows.push_back(make_pair(start_pos - ActiveRegionFinder_options::min_region_size,
                                        end_pos + ActiveRegionFinder_options::min_region_size));
            start_pos = pos;
            end_pos = pos;
        }
        else {
            end_pos = pos;
        }
    }
    if(start_pos != -1) {
        windows.push_back(make_pair(start_pos - ActiveRegionFinder_options::min_region_size,
                                    end_pos + ActiveRegionFinder_options::min_region_size));
    }
    sort(windows.begin(), windows.end());

    return windows;
}

ActiveRegionFinder::~ActiveRegionFinder(){
}