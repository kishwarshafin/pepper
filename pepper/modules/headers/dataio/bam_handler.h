//
// Created by Kishwar Shafin on 6/12/18.
//

#ifndef FRIDAY_CPP_BAM_HANDLER_H
#define FRIDAY_CPP_BAM_HANDLER_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "sam.h"
#include "hts.h"
#include "cram.h"
#include "hts_endian.h"

using namespace std;

#define SNP_TYPE 1
#define INSERT_TYPE 2
#define DELETE_TYPE 3

//typedef struct type_read type_read;
//typedef struct LinearAlignment LinearAlignment;

struct CigarOp {
    CigarOp() : operation(-1), length(0) {}
    CigarOp(int op, int len) : operation(op), length(len) {}

    bool operator==(const CigarOp& that) const {
        return operation == that.operation && length == that.length;
    }

    void operator=(const CigarOp& that) {
        this->operation = that.operation;
        this->length = that.length;
    }

    void set_operation(int op) {
        operation = op;
    }

    void set_length(int len) {
        length = len;
    }

    int operation;
    int length;
};

struct type_read_flags{
    bool is_paired;
    bool is_proper_pair;
    bool is_unmapped;
    bool is_mate_unmapped;
    bool is_reverse;
    bool is_mate_is_reverse;
    bool is_read1;
    bool is_read2;
    bool is_secondary;
    bool is_qc_failed;
    bool is_duplicate;
    bool is_supplementary;
    type_read_flags(){
        is_paired= 0;
        is_proper_pair= 0;
        is_unmapped= 0;
        is_mate_unmapped= 0;
        is_reverse= 0;
        is_mate_is_reverse= 0;
        is_read1= 0;
        is_read2= 0;
        is_secondary= 0;
        is_qc_failed= 0;
        is_duplicate= 0;
        is_supplementary= 0;
    }
    void operator=(const type_read_flags& that) {
        this->is_paired = that.is_paired;
        this->is_proper_pair = that.is_proper_pair;
        this->is_unmapped = that.is_unmapped;
        this->is_mate_unmapped = that.is_mate_unmapped;
        this->is_reverse = that.is_reverse;
        this->is_mate_is_reverse = that.is_mate_is_reverse;
        this->is_read1 = that.is_read1;
        this->is_read2 = that.is_read2;
        this->is_secondary = that.is_secondary;
        this->is_qc_failed = that.is_qc_failed;
        this->is_duplicate = that.is_duplicate;
        this->is_supplementary = that.is_supplementary;

    }
};

struct type_read{
    long long pos;
    long long pos_end;
    string query_name;
    type_read_flags flags;
    string sequence;
    vector <CigarOp> cigar_tuples;
    vector <int> bad_indicies;
    int mapping_quality;
    vector <int> base_qualities;
    int read_id;
    int hp_tag;

    void set_read_id(int id) {
        this->read_id = id;
    }

    void set_position(long long pos){
        this->pos = pos;
    }
    void set_end_position(long long pos){
        this->pos_end = pos;
    }

    void set_cigar_tuples(vector <CigarOp> cigar_tuples) {
        this->cigar_tuples = cigar_tuples;
    }

    void operator=(const type_read& that) {
        this->pos = that.pos;
        this->pos_end = that.pos_end;
        this->query_name = that.query_name;
        this->flags = that.flags;
        this->sequence = that.sequence;
        this->cigar_tuples = that.cigar_tuples;
        this->bad_indicies = that.bad_indicies;
        this->mapping_quality = that.mapping_quality;
        this->base_qualities = that.base_qualities;
        this->hp_tag = that.hp_tag;
    }

    bool operator<(const type_read& that) {
        if(this->pos == that.pos) {
            return this->pos_end < that.pos_end;
        } else {
            return this->pos < that.pos;
        }
    }
};

typedef struct{
    string sequence_name;
    int sequence_length;
} type_sequence;


class CIGAR_OPERATIONS {
public:
    static constexpr int MATCH = 0;
    static constexpr int IN = 1;
    static constexpr int DEL = 2;
    static constexpr int REF_SKIP = 3;
    static constexpr int SOFT_CLIP = 4;
    static constexpr int HARD_CLIP = 5;
    static constexpr int PAD = 6;
    static constexpr int EQUAL = 7;
    static constexpr int DIFF = 8;
    static constexpr int BACK = 9;
    static constexpr int UNSPECIFIED = -1;
};

class BAM_handler {
    public:
        htsFile* hts_file;
        hts_idx_t* idx;
        bam_hdr_t* header;

        BAM_handler(string path);

        // this will divide reads in haplotype bins and then return
        vector<type_read> get_reads(string region,
                                    long long start,
                                    long long stop,
                                    bool include_supplementary,
                                    int min_mapq,
                                    int min_baseq);
        vector<string> get_chromosome_sequence_names();
        vector<type_sequence> get_chromosome_sequence_names_with_length();
        set<string> get_sample_names();
        type_read_flags get_read_flags(int flag);

    	~BAM_handler();
};

#endif // BAM_HANDLER_H
