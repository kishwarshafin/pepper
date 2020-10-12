//
// Created by Kishwar Shafin on 6/14/18.
//

#ifndef PEPPER_VARIANT_FASTA_HANDLER_H
#define PEPPER_VARIANT_FASTA_HANDLER_H

#include "faidx.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class FASTA_handler {
    public:
        FASTA_handler(string path);
        // this start and stop is zero-based. Should we change everything to one-based?
        string get_reference_sequence(string region, long long start, long long stop);
        int get_chromosome_sequence_length(string chromosome_name);
        vector<string> get_chromosome_names();
        ~FASTA_handler();
    private:
        faidx_t* fasta;
};


#endif //FRIDAY_CPP_FASTA_HANDLER_H
