//
// Created by Kishwar Shafin on 6/14/18.
//

#include "../../headers/dataio/fasta_handler.h"

FASTA_handler::FASTA_handler(string path) {
    this->fasta = fai_load(path.c_str());
    if (fasta == NULL) {
        cerr<<"INVALID FASTA FILE. PLEASE CHECK IF PATH IS CORRECT AND FILE IS INDEXED: "<<path<<endl;
        exit (EXIT_FAILURE);
    }
}

FASTA_handler::~FASTA_handler() {
    fai_destroy(this->fasta);
}

vector<string> FASTA_handler::get_chromosome_names() {
    vector<string> chromosome_names;
    int number_of_chromosome_sequences = faidx_nseq(this->fasta);

    for(int i=0; i < number_of_chromosome_sequences; i++) {
        string chromosome_name = faidx_iseq(this->fasta, i);
        chromosome_names.push_back(chromosome_name);
    }

    return chromosome_names;
}

string FASTA_handler::get_reference_sequence(string region, long long start, long long stop) {
    // the fetch is zero-based, in one-based system the previous base would be the start.
    // this is done to make the reference sequence compatible with the BAM file's read sequences.
    // start += 1;
    // stop += 1;

    int len = 0;
    string sequence;

    sequence = faidx_fetch_seq(fasta, region.c_str(), start, stop - 1, &len);
    //-2 if c_name not present, -1 general error
    if(len == -2){
        cerr<<"CHROMOSOME NAME NOT PRESENT IN REFERENCE FASTA FILE: "<<region<<" "<<start<<" "<<stop<<endl;
        return NULL;
    } else if(len == -1){
        cerr<<"ENCOUNTERED ERROR IN FETCHING REFERENCE FASTA FILE: "<<region<<" "<<start<<" "<<stop<<endl;
        return NULL;
    }
    return sequence;
}

int FASTA_handler::get_chromosome_sequence_length(string chromosome_name) {
    int len = faidx_seq_len(this->fasta, chromosome_name.c_str());
    return len;
}