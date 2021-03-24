//
// Created by Kishwar Shafin on 3/21/21.
//

#ifndef PEPPER_PRIVATE_CIGAR_H
#define PEPPER_PRIVATE_CIGAR_H

#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <algorithm>
using namespace std;

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

#endif //PEPPER_PRIVATE_CIGAR_H
