//
// Created by Kishwar Shafin on 6/12/18.
//

#ifndef PEPPER_VARIANT_BAM_HANDLER_H
#define PEPPER_VARIANT_BAM_HANDLER_H

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
#include "read.h"
#include "cigar.h"
#include "sequence.h"
using namespace std;

#define SNP_TYPE 1
#define INSERT_TYPE 2
#define DELETE_TYPE 3

class BAM_handler {
    public:
        htsFile* hts_file;
        hts_idx_t* idx;
        bam_hdr_t* header;

        // initialize a bam file
        BAM_handler(string path);

        // get reads from a bam file given a region
        vector<type_read> get_reads(string region,
                                    long long start,
                                    long long stop,
                                    bool include_supplementary,
                                    int min_mapq,
                                    int min_baseq);

        // get sequence names from a bam file
        vector<string> get_chromosome_sequence_names();

        // get contig names and lengths from the bam header
        vector<type_sequence> get_chromosome_sequence_names_with_length();

        // get sample name from the bam header
        set<string> get_sample_names();

        // decipher the read flag
        type_read_flags get_read_flags(int flag);

    	~BAM_handler();
};

#endif // BAM_HANDLER_H
