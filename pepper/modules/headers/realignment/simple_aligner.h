//
// Created by Kishwar Shafin on 10/21/19.
//

#ifndef PEPPER_SIMPLE_ALIGNER_H
#define PEPPER_SIMPLE_ALIGNER_H

#include "ssw_cpp.h"
#include "ssw.h"
#include "../dataio/bam_handler.h"
#include <unordered_map>
#include <list>
#include <regex>
#include <string>
#include <sstream>
#include <memory>
using namespace StripedSmithWaterman;

namespace Aligner_options {
    static constexpr int match = 4;
    static constexpr int mismatch = 6;
    static constexpr int gap_open_penalty = 8;
    static constexpr int gap_extend_penalty = 2;
    static constexpr int max_number_of_mismatches = 2;
};

class LibSSWPairwiseAligner {
    Filter filter;
    Aligner* ssw_aligner;
public:
    LibSSWPairwiseAligner();
    void set_reference(string reference);
    Alignment align(string query);
};

class ReadAligner {
    string reference_sequence;
    int region_start;
    int region_end;
    LibSSWPairwiseAligner SSWAligner;
    list<CigarOp> CigarStringToVector(const string& cigar);
public:
    ReadAligner(int ref_start, int ref_end, string ref_seq);
    vector<type_read> align_reads_to_reference(vector<type_read> reads);
};

#endif //PEPPER_SIMPLE_ALIGNER_H
