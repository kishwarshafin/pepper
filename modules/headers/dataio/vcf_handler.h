//
// Created by Kishwar Shafin on 6/14/18.
//

#ifndef FRIDAY_CPP_VCF_HANDLER_H
#define FRIDAY_CPP_VCF_HANDLER_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "vcf.h"
#include "tbx.h"

using namespace std;

#define SNP_TYPE 1
#define INSERT_TYPE 2
#define DELETE_TYPE 3

typedef struct{
    string ref;
    string alt_allele;
    int alt_type;

    string get_ref() {
        return ref;
    }
    string get_alt_allele() {
        return alt_allele;
    }
    int get_alt_type(){
        return alt_type;
    }
} type_alt_allele;

typedef struct{
    string chromosome_name;
    long long start_pos;
    long long end_pos;
    string id;
    double qual;
    bool is_phased;
    bool is_filter_pass;
    string sample_name;
    vector<int> genotype;
    vector<string> filters;
    vector<string> alleles;
//    vector<type_alt_allele> alt_allele;
} type_vcf_record;


typedef struct{
    string chromosome_name;
    uint64_t start_pos;
    uint64_t end_pos;
    string id;
    double qual;
    bool is_phased;
    bool is_filter_pass;
    string sample_name;
    vector<int> genotype;
    vector<string> filters;
    vector<type_alt_allele> alt_allele;
} type_positional_vcf_record;

class VCF_handler {
    public:
        VCF_handler(string file_path);
        map<long long, vector<type_positional_vcf_record> > get_positional_vcf_records(string chromosome_name,
                                                                                       long long start,
                                                                                       long long stop);
        vector<type_vcf_record> get_vcf_records(string chromosome_name,
                                                long long start,
                                                long long stop) ;
        int process_insert_allele(string ref, string alt, type_alt_allele &alt_allele);
        int process_delete_allele(string ref, string alt, type_alt_allele &alt_allele);
        int process_snp_allele(string ref, string alt, type_alt_allele &alt_allele);
        ~VCF_handler();
    private:
        htsFile * vcf_file;
        bcf_hdr_t * vcf_header;
        tbx_t * vcf_index;
};

#endif //FRIDAY_CPP_VCF_HANDLER_H
