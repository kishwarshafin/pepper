//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_READ_H
#define PEPPER_PRIVATE_READ_H

#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <algorithm>
#include "cigar.h"
using namespace std;

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



#endif //PEPPER_PRIVATE_READ_H
