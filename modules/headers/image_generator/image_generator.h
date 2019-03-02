//
// Created by Kishwar Shafin on 11/1/18.
//

#ifndef FRIDAY_IMAGE_GENERATOR_H
#define FRIDAY_IMAGE_GENERATOR_H

#include <algorithm>
using namespace std;
#include "../dataio/bam_handler.h"
#include "../candidate_finding/candidate_finder.h"

namespace PileupPixels {
    static constexpr int MAX_COLOR_VALUE = 254;
    static constexpr int BASE_QUALITY_CAP = 40;
    static constexpr int MAP_QUALITY_CAP = 60;
    static constexpr int REF_ROW_BAND = 1;
    static constexpr int IMAGE_HEIGHT = 100;
    static constexpr int CONTEXT_SIZE = 10;
};


namespace Genotypes {
    static constexpr int HOM = 0;
    static constexpr int HET = 1;
    static constexpr int HOM_ALT = 2;
};

struct PileupImage {
    string chromosome_name;
    long long start_pos;
    long long end_pos;
    vector<vector<vector<int> > >image;
    vector<int> label;
    void set_values(string chromosome_name, long long start_pos, long long end_pos) {
        this->chromosome_name = chromosome_name;
        this->start_pos = start_pos;
        this->end_pos = end_pos;
    }
};

class ImageGenerator {
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<long long, PositionalCandidateRecord> all_positional_candidates;
    map<char, int> global_base_color;
    map<long long, vector<type_positional_vcf_record> > pos_vcf;
public:
    ImageGenerator(string reference_sequence,
                   string chromosome_name,
                   long long ref_start,
                   long long ref_end,
                   map<long long, PositionalCandidateRecord> all_positional_candidates);
    void set_positional_vcf(map<long long, vector<type_positional_vcf_record> > pos_vcf);

    string get_reference_sequence(long long st_pos, long long end_pos);
    vector<vector<int> > read_to_image_row(type_read read, long long &read_start, long long &read_end);
    vector<vector<int> > get_reference_row(string ref_seq);
    int get_image_label(int gt1, int gt2);
    vector<int> get_window_labels(pair<long long, long long> window);
    vector<PileupImage> create_window_pileups(vector<pair<long long, long long> > windows,
                                              vector<type_read> reads,
                                              bool train_mode);
    int get_which_allele(long long pos, string ref, string alt, int alt_type);
    long long overlap_length_between_ranges(pair<long long, long long> range_a,
                                            pair<long long, long long> range_b);
    void assign_read_to_window(PileupImage& pileup_image,
                               vector<vector<int> >& image_row,
                               long long read_start,
                               long long read_end);
};


#endif //FRIDAY_IMAGE_GENERATOR_H
