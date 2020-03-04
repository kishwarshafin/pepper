//
// Created by Kishwar Shafin on 10/21/18.
//

// REUSE DECLARATION: THIS SCRIPT HEAVILY REUSES THE METHODS USED IN DEEPVARIANT https://github.com/google/deepvariant/
// LICENSE INCLUDED IN: third_party/deepVariant.LICENSE

#ifndef HELEN_SSW_ALIGNER_H
#define HELEN_SSW_ALIGNER_H

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

struct LinearAlignment {
    long long pos;
    int mapping_quality;
    vector<CigarOp> cigar_tuples;

    void copy_from_read(const type_read& that){
        this->pos = that.pos;
        this->cigar_tuples = that.cigar_tuples;
        this->mapping_quality = that.mapping_quality;
    }

    void clear_cigar() {
        cigar_tuples.clear();
    }

    void set_position(long long pos){
        this->pos = pos;
    }

    void add_cigar(CigarOp cigar_unit) {
        this->cigar_tuples.push_back(cigar_unit);
    }
};


struct ReadAlignment {
    ReadAlignment() : position(0), cigar(""), score(0) {}

    ReadAlignment(uint16_t position_param, const string& cigar_param, int score_param)
            : position(position_param), cigar(cigar_param), score(score_param) {}

    bool operator==(const ReadAlignment& that) const {
        return score == that.score && position == that.position &&
               cigar == that.cigar;
    }

    void reset() {
        score = 0;
        position = 0;
        cigar = "";
    }

    long long position;
    string cigar;
    int score;
};


struct HaplotypeReadsAlignment {
    HaplotypeReadsAlignment() : haplotype_index(0), haplotype_score(0) {}
    HaplotypeReadsAlignment(size_t haplotype_index, int score, const vector<ReadAlignment>& read_alignment_scores)
            : haplotype_index(haplotype_index), haplotype_score(score) {
        this->read_alignment_scores.assign(read_alignment_scores.begin(),
                                           read_alignment_scores.end());
    }

    bool operator==(const HaplotypeReadsAlignment& that) const {
        return haplotype_index == that.haplotype_index &&
               haplotype_score == that.haplotype_score &&
               read_alignment_scores == that.read_alignment_scores &&
               cigar == that.cigar && cigar_ops == that.cigar_ops &&
               is_reference == that.is_reference &&
               hap_to_ref_positions_map == that.hap_to_ref_positions_map;
    }

    bool operator<(const HaplotypeReadsAlignment& that) const {
        return haplotype_score < that.haplotype_score;
    }
    // halpotype index in haplotypes member of FastPassAligner class
    size_t haplotype_index;

    // sum of all aligned read scores. Each read's
    // haplotype_score = number of matched bases.
    int haplotype_score;

    // Simple alignment haplotype_score for each read.
    std::vector<ReadAlignment> read_alignment_scores;

    // Cigar string for haplotype against reference alignment.
    string cigar;

    // Cigar as a list of CigarOp
    std::list<CigarOp> cigar_ops;

    // Hypolotype to reference position.
    long long ref_pos;

    // Map of shifts that allow to easily calculate read to reference position
    // from read to haplotype position. read pos = read_to_haplotype_pos
    //    + haplotype_to_ref_position
    //    + hap_to_ref_positions_map[read_to_haplotype_pos]
    std::vector<int> hap_to_ref_positions_map;

    // If true the haplotype is a reference.
    bool is_reference;
};




struct KmerMap {
    KmerMap() {}
    KmerMap(long long  read_id, long long pos) : read_id(read_id), read_pos(pos) {}

    bool operator==(const KmerMap& b) const {
        return read_id == b.read_id && read_pos == b.read_pos;
    }

    long long read_id;
    long long read_pos;
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
    vector<string> read_sequences;
    vector<string> haplotypes;
    int kmer_size = 32;
    float similarity_threshold = 0.17;
    double alignment_score_threshold;
    unordered_map<string, vector<KmerMap> > KmerIndexMap;
    vector<HaplotypeReadsAlignment> read_to_haplotype_alignments;
    LibSSWPairwiseAligner SSWAligner;

    void CreateKmerMap();
    void AlignReadsToHaplotypes();
    void AlignReadsToHaplotype( const string& haplotype,
                                int* haplotype_score,
                                vector<ReadAlignment>* haplotype_read_alignment_scores);
    int AlignStrings(string s1, string s2, int max_mismatches, int* num_of_mismatches) const;
    void AlignHaplotypesToReference();
    list<CigarOp> CigarStringToVector(const string& cigar);
    void CalculatePositionMaps();
    void SswAlignReadsToHaplotypes(int score_threshold);
    void RealignReadsToReference(const vector<type_read>& reads,
                                 vector<type_read>& realigned_reads);
    bool GetBestReadAlignment(size_t readId, int* best_hap_index) const;
    void CalculateReadToRefAlignment(size_t read_index, const ReadAlignment& read_to_haplotype_alignment,
                                     const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
                                     list<CigarOp>* read_to_ref_cigar_ops);
    list<CigarOp> LeftTrimHaplotypeToRefAlignment(const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
                                                  int read_to_haplotype_pos);
    void MergeCigarOp(const CigarOp& op, int read_len, list<CigarOp>* cigar);
    int AlignedLength(const list<CigarOp>& cigar);
public:
    ReadAligner(int ref_start, int ref_end, string ref_seq);
    vector<type_read> align_reads(vector<string> haplotypes, vector<type_read> reads);
    vector<type_read> align_reads_to_reference(vector<type_read> reads);
    vector<type_read> align_haplotypes_to_reference(vector<string> haplotypes);
};

inline bool BothOpsAreMatch(const CigarOp& op1, const CigarOp& op2);
inline bool OneOfOpsIsSoftClip(const CigarOp& op1, const CigarOp& op2);
inline bool DelAndMatch(const CigarOp& op1, const CigarOp& op2);
inline bool BothOpsAreDel(const CigarOp& op1, const CigarOp& op2);
inline bool InsAndMatch(const CigarOp& op1, const CigarOp& op2);
inline bool BothOpsAreIns(const CigarOp& op1, const CigarOp& op2);
inline void PushFrontIfNotEmpty(const CigarOp& op, list<CigarOp>* cigar);
#endif //FRIDAY_SSW_ALIGNER_H
