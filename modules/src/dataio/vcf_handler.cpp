//
// Created by Kishwar Shafin on 6/14/18.
//
#include "../../headers/dataio/vcf_handler.h"

VCF_handler::VCF_handler(string file_path)
{
    // Opening the VCF file
    this->vcf_file = bcf_open(file_path.c_str(), "r");
    if(this->vcf_file == NULL) {
        cerr<<"INVALID VCF FILE. PLEASE CHECK IF PATH IS CORRECT AND FILE IS INDEXED: "<<file_path<<endl;
        exit (EXIT_FAILURE);
    }
    this->vcf_header = bcf_hdr_read(this->vcf_file);
    if(this->vcf_header == NULL) {
        cerr<<"INVALID VCF HEADER. PLEASE CHECK IF PATH IS CORRECT AND FILE IS INDEXED: "<<file_path<<endl;
        exit (EXIT_FAILURE);
    }

    this->vcf_index = tbx_index_load(file_path.c_str());
    if(this->vcf_index == NULL){
        cerr<<"INVALID VCF INDEX FILE. PLEASE CHECK IF FILE IS INDEXED: "<<file_path<<endl;
        exit (EXIT_FAILURE);
    }
}

int VCF_handler::process_insert_allele(string ref, string alt, type_alt_allele &alt_allele) {
    // make sure it's an insert
    if(ref.length() < alt.length()) {
        // it's an insert
        string processed_alt;
        int insert_length = alt.length() - ref.length() + 1;
        for(int i=0; i<insert_length; i++) processed_alt += alt[i];
        alt_allele.ref = ref[0];
        alt_allele.alt_allele = processed_alt;
        alt_allele.alt_type = INSERT_TYPE;
        return 0;
    }
    return -1;
}

int VCF_handler::process_delete_allele(string ref, string alt, type_alt_allele &alt_allele) {
    // make sure it's a delete
    if(ref.length() > alt.length()) {
        // it's a delete
        string processed_alt;
        int insert_length = ref.length() - alt.length() + 1;
        for(int i=0; i<insert_length; i++) processed_alt += ref[i];
        alt_allele.ref = ref[0];
        alt_allele.alt_allele = processed_alt;
        alt_allele.alt_type = DELETE_TYPE;
        return 0;
    }
    return -1;
}

int VCF_handler::process_snp_allele(string ref, string alt, type_alt_allele &alt_allele) {
    // see if there's a snp
    if(ref[0] != alt[0]){
        alt_allele.ref = ref[0];
        alt_allele.alt_allele = alt[0];
        alt_allele.alt_type = SNP_TYPE;
        return 0;
    }
    return -1;
}

map<long long, vector<type_positional_vcf_record> > VCF_handler::get_positional_vcf_records(string chromosome_name,
                                                                                 long long start,
                                                                                 long long stop) {
    map<long long, vector<type_positional_vcf_record> > positional_vcf_records;

    // get the id of the chromosome or sequence name
    int tid = tbx_name2id(this->vcf_index, chromosome_name.c_str());

    // get the iterator
    hts_itr_t *iter  = tbx_itr_queryi(this->vcf_index, tid, start, stop);

    // initialize an alignment
    kstring_t str_ = {0, 0, nullptr};
    // cout<<"Chr\t"<<"start\t"<<"end\t"<<"id\t\t"<<"Nalt\t"<<"Ref\tAlts\t"<<"Qual\t"<<"Filter\t"<<"Sample\t"<<"GT\t"<<"GT\t"<<"PD\t"<<endl;
    while(tbx_itr_next(this->vcf_file, this->vcf_index, iter, &str_) > -1) {
        type_positional_vcf_record vcf_rec;
        bcf1_t *bcf1_ = bcf_init();

        if(vcf_parse(&str_, this->vcf_header, bcf1_) > -1) {
            bcf_unpack(bcf1_, BCF_UN_ALL);
            vcf_rec.chromosome_name = bcf_hdr_id2name(this->vcf_header, bcf1_->rid);
            vcf_rec.start_pos = bcf1_->pos;
            vcf_rec.end_pos = (bcf1_->pos + bcf1_->rlen);

            if (bcf1_->d.id && strcmp(bcf1_->d.id, ".") != 0) {
                vcf_rec.id = bcf1_->d.id;
            } else {
                vcf_rec.id = '.';
            }

            // allele specific processing
            if (bcf1_->n_allele > 0) {
                string ref_allele = bcf1_->d.allele[0];
                for (int i = 1; i < bcf1_->n_allele; i++) {
                    string alt_allele = bcf1_->d.allele[i];
                    if(ref_allele.length() < alt_allele.length()){
                        type_alt_allele in_allele;
                        if(process_insert_allele(ref_allele, alt_allele, in_allele) >= 0){
                            vcf_rec.alt_allele.push_back(in_allele);
                        }
                    }
                    else if(ref_allele.length() > alt_allele.length()){
                        type_alt_allele del_allele;
                        if(process_delete_allele(ref_allele, alt_allele, del_allele) >= 0){
                            vcf_rec.alt_allele.push_back(del_allele);
                        }
                    }
                    type_alt_allele snp_allele;
                    if(process_snp_allele(ref_allele, alt_allele, snp_allele) >= 0){
                        vcf_rec.alt_allele.push_back(snp_allele);
                    }
                }
            }

            if (bcf_float_is_missing(bcf1_->qual)) {
                vcf_rec.qual = -1;
            } else {
                vcf_rec.qual = bcf1_->qual;
            }
            // get filters
            vcf_rec.is_filter_pass = false;
            for (int i = 0; i < bcf1_->d.n_flt; ++i) {
                string filter = this->vcf_header->id[BCF_DT_ID][bcf1_->d.flt[i]].key;
                vcf_rec.filters.push_back(filter);
                if(filter.compare("PASS") == 0) {
                    vcf_rec.is_filter_pass = true;
                }
            }

            // Parse the calls of the variant.
            if (bcf1_->n_sample > 0) {
                int *gt_array = nullptr;
                int ploidy, n_gts = 0;
                if (bcf_get_genotypes(this->vcf_header, bcf1_, &gt_array, &n_gts) < 0) {
                    free(gt_array);
                    cerr<<"VCF ERROR: COULD NOT PARSE GT"<<endl;
                    continue;
                }
                ploidy = n_gts / bcf1_->n_sample;

                for (int i = 0; i < bcf1_->n_sample; i++) {
                    vcf_rec.sample_name = this->vcf_header->samples[i];
                    // Get the GT calls, if requested and available.
                    bool gt_is_phased = false;
                    for (int j = 0; j < ploidy; j++) {
                        int gt_idx = gt_array[i * ploidy + j];
                        int gt = bcf_gt_allele(gt_idx);
                        gt_is_phased = gt_is_phased || bcf_gt_is_phased(gt_idx);
                        vcf_rec.genotype.push_back(gt);
                    }
                    vcf_rec.is_phased = gt_is_phased;
                }
            }
            positional_vcf_records[vcf_rec.start_pos].push_back(vcf_rec);

            //TODO: parse the info and format field. We don't need it now.
        }
    }
    free(str_.s);
    tbx_itr_destroy(iter);
    return positional_vcf_records;
}


vector<type_vcf_record> VCF_handler::get_vcf_records(string chromosome_name, long long start, long long stop) {
    vector<type_vcf_record> vcf_records;

    // get the id of the chromosome or sequence name
    int tid = tbx_name2id(this->vcf_index, chromosome_name.c_str());

    // get the iterator
    hts_itr_t *iter  = tbx_itr_queryi(this->vcf_index, tid, start, stop);

    // initialize an alignment
    kstring_t str_ = {0, 0, nullptr};
    // cout<<"Chr\t"<<"start\t"<<"end\t"<<"id\t\t"<<"Nalt\t"<<"Ref\tAlts\t"<<"Qual\t"<<"Filter\t"<<"Sample\t"<<"GT\t"<<"GT\t"<<"PD\t"<<endl;
    while(tbx_itr_next(this->vcf_file, this->vcf_index, iter, &str_) > -1) {
        type_vcf_record vcf_rec;
        bcf1_t *bcf1_ = bcf_init();

        if(vcf_parse(&str_, this->vcf_header, bcf1_) > -1) {
            bcf_unpack(bcf1_, BCF_UN_ALL);
            vcf_rec.chromosome_name = bcf_hdr_id2name(this->vcf_header, bcf1_->rid);
            vcf_rec.start_pos = bcf1_->pos;
            vcf_rec.end_pos = (bcf1_->pos + bcf1_->rlen);

            if (bcf1_->d.id && strcmp(bcf1_->d.id, ".") != 0) {
                vcf_rec.id = bcf1_->d.id;
            } else {
                vcf_rec.id = '.';
            }

            // allele specific processing
            if (bcf1_->n_allele > 0) {
                string ref_allele = bcf1_->d.allele[0];
                vcf_rec.alleles.push_back(ref_allele);
                for (int i = 1; i < bcf1_->n_allele; i++) {
                    string alt_allele = bcf1_->d.allele[i];
                    vcf_rec.alleles.push_back(alt_allele);
                }
            }

            if (bcf_float_is_missing(bcf1_->qual)) {
                vcf_rec.qual = -1;
            } else {
                vcf_rec.qual = bcf1_->qual;
            }
            // get filters
            vcf_rec.is_filter_pass = false;
            for (int i = 0; i < bcf1_->d.n_flt; ++i) {
                string filter = this->vcf_header->id[BCF_DT_ID][bcf1_->d.flt[i]].key;
                vcf_rec.filters.push_back(filter);
                if(filter.compare("PASS") == 0) {
                    vcf_rec.is_filter_pass = true;
                }
            }

            // Parse the calls of the variant.
            if (bcf1_->n_sample > 0) {
                int *gt_array = nullptr;
                int ploidy, n_gts = 0;
                if (bcf_get_genotypes(this->vcf_header, bcf1_, &gt_array, &n_gts) < 0) {
                    free(gt_array);
                    cerr<<"VCF ERROR: COULD NOT PARSE GT"<<endl;
                    continue;
                }
                ploidy = n_gts / bcf1_->n_sample;

                for (int i = 0; i < bcf1_->n_sample; i++) {
                    vcf_rec.sample_name = this->vcf_header->samples[i];
                    // Get the GT calls, if requested and available.
                    bool gt_is_phased = false;
                    for (int j = 0; j < ploidy; j++) {
                        int gt_idx = gt_array[i * ploidy + j];
                        int gt = bcf_gt_allele(gt_idx);
                        gt_is_phased = gt_is_phased || bcf_gt_is_phased(gt_idx);
                        vcf_rec.genotype.push_back(gt);
                    }
                    vcf_rec.is_phased = gt_is_phased;
                }
            }
            vcf_records.push_back(vcf_rec);

            //TODO: parse the info and format field. We don't need it now.
        }
    }
    free(str_.s);
    tbx_itr_destroy(iter);
    return vcf_records;
}

VCF_handler::~VCF_handler() {
    hts_close(this->vcf_file);
    bcf_hdr_destroy(this->vcf_header);
    tbx_destroy(this->vcf_index);
}
