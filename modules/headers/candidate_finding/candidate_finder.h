//
// Created by Kishwar Shafin on 10/24/18.
//

#ifndef FRIDAY_CANDIDATE_FINDER_H
#define FRIDAY_CANDIDATE_FINDER_H

#include <cmath>
#include "../dataio/bam_handler.h"
using namespace std;

namespace CandidateFinder_options {
    static constexpr int min_mapping_quality = 10;
    static constexpr int min_base_quality = 10;
    static constexpr int freq_threshold = 12;
    static constexpr int min_count_threshold = 2;
};

namespace AlleleType {
    static constexpr int SNP_ALLELE = 1;
    static constexpr int INSERT_ALLELE = 2;
    static constexpr int DELETE_ALLELE = 3;
};

namespace Genotype {
    static constexpr int HOM = 0;
    static constexpr int HET = 1;
    static constexpr int HOM_ALT = 2;
};


struct PositionalCandidateRecord{
    string chromosome_name;
    long long pos;
    long long pos_end;
    string ref;
    string alt1;
    string alt2;
    int alt1_type;
    int alt2_type;


    bool labeled;
    int alt1_gt;
    int alt2_gt;
    vector<int> genotype {0, 0};

    PositionalCandidateRecord(string chr_name, long long pos, long long pos_end, string ref, string alt1, string alt2,
    int alt1_type, int alt2_type) {
        this->chromosome_name = chr_name;
        this->pos = pos;
        this->pos_end = pos_end;
        this->ref = ref;
        this->alt1 = alt1;
        this->alt2 = alt2;
        this->alt1_type = alt1_type;
        this->alt2_type = alt2_type;
    }

    PositionalCandidateRecord() {
        this->alt1_type = 0;
        this->alt2_type = 0;
        this->labeled = 0;
        this->alt1_gt = 0;
        this->alt2_gt = 0;
        this->alt1 = '.';
        this->alt2 = '.';
    }

    void set_positions(string chromosme_name, long long pos, long long pos_end) {
        this->chromosome_name = chromosme_name;
        this->pos = pos;
        this->pos_end = pos_end;
    }

    void set_reference(string ref) {
        this->ref = ref;
    }

    void set_alt1(string alt1, int alt1_type) {
        this->alt1 = alt1;
        this->alt1_type = alt1_type;
    }

    void set_alt2(string alt2, int alt2_type, long long pos_end) {
        this->alt2 = alt2;
        this->alt2_type = alt2_type;
        this->pos_end = max(this->pos_end, pos_end);
    }

    void set_genotype(vector<int> genotype) {
        this->genotype = genotype;
        this->labeled = true;
        int rec_genotype = Genotype::HOM;
        if(genotype[0] != 0 or genotype[1] != 0) {
            if(genotype[0] == genotype[1]) {
                rec_genotype = Genotype::HOM_ALT;
            } else {
                rec_genotype = Genotype::HET;
            }
        }

        if(genotype[0] == 1 or genotype[1] == 1) {
            this->alt1_gt = rec_genotype;
        }
        if(genotype[0] == 2 or genotype[1] == 2) {
            this->alt2_gt = rec_genotype;
        }
    }

    string get_type(int type) {
        if(type == AlleleType::SNP_ALLELE) return "SNP";
        if(type == AlleleType::INSERT_ALLELE) return "IN";
        if(type == AlleleType::DELETE_ALLELE) return "DEL";
        return "";
    }

    string get_gt(int gt){
        if(gt == Genotype::HOM) return "HOM";
        if(gt == Genotype::HET) return "HET";
        if(gt == Genotype::HOM_ALT) return "HOM_ALT";
        return "";
    }

    void print() {
        cout<<this->chromosome_name<<" "<<this->pos<<" "<<this->pos_end<<" "<<this->ref;
        cout<<" "<<this->alt1<<" "<<this->alt2<<" "<<get_type(this->alt1_type)<<" "<<get_type(this->alt2_type);
        if(this->labeled) cout<<" ("<<genotype[0]<<", "<<genotype[1]<<")";
        if(this->labeled) cout<<" "<<get_gt(this->alt1_gt)<<" "<<get_gt(this->alt2_gt);
        cout<<endl;
    }

    vector<string> get_candidate_record() {
        vector<string> records;

        string rec_alt1 = this->chromosome_name + "\t" + to_string(this->pos) + "\t" + to_string(this->pos_end);
        rec_alt1 = rec_alt1 + "\t" + this->ref + "\t" + this->alt1 + "\t" + '.' ;
        rec_alt1 = rec_alt1 + "\t" + to_string(this->alt1_type) + "\t" + '0' + "\t" + to_string(this->alt1_gt);
        records.push_back(rec_alt1);

        if(this->alt2.compare(".") == 0) {
            return records;
        }

        string rec_alt2 = this->chromosome_name + "\t" + to_string(this->pos) + "\t" + to_string(this->pos_end);
        rec_alt2 = rec_alt2 + "\t" + this->ref + "\t" + this->alt2 + "\t" + '.' ;
        rec_alt2 = rec_alt2 + "\t" + to_string(this->alt2_type) + "\t" + '0' + "\t" + to_string(this->alt2_gt);
        records.push_back(rec_alt2);

        int rec_genotype = Genotype::HOM;
        if(this->genotype[0] != 0 or this->genotype[1] != 0) {
            if(this->genotype[0] == this->genotype[1]) {
                rec_genotype = Genotype::HOM_ALT;
            } else if(this->genotype[0] != 0 and this->genotype[1] != 0) {
                rec_genotype = Genotype::HOM_ALT; // 2/1 or 1/2
            } else {
                rec_genotype = Genotype::HET;
            }
        }
        string rec_com = this->chromosome_name + "\t" + to_string(this->pos) + "\t" + to_string(this->pos_end);
        rec_com = rec_com + "\t" + this->ref + "\t" + this->alt1 + "\t" + this->alt2 ;
        rec_com = rec_com + "\t" + to_string(this->alt1_type) + "\t" + to_string(this->alt2_type) + "\t" + to_string(rec_genotype);
        records.push_back(rec_com);

        return records;
    }
};

struct CandidateAllele{
    string ref;
    string alt;
    int alt_type;

    CandidateAllele(string ref, string alt, int alt_type) {
        this->ref = ref;
        this->alt = alt;
        this->alt_type = alt_type;
    }

    CandidateAllele() {}
};

struct Candidate{
    long long pos;
    long long pos_end;
    CandidateAllele allele;

    Candidate(long long pos_start, long long pos_end, string ref, string alt, int alt_type) {
        this->pos = pos_start;
        this->pos_end = pos_end;
        this->allele.alt = alt;
        this->allele.ref = ref;
        this->allele.alt_type = alt_type;
    }
    bool operator< (const Candidate& that ) const {
        if(this->pos != that.pos) return this->pos < that.pos;
        if(this->pos_end != that.pos_end) return this->pos_end < that.pos_end;
        if(this->allele.alt != that.allele.alt) return this->allele.alt < that.allele.alt;
        if(this->allele.ref != that.allele.ref) return this->allele.ref < that.allele.ref;
        if(this->allele.alt_type != that.allele.alt_type) return this->allele.alt_type < that.allele.alt_type;
        return this->pos < that.pos;
    }

    bool operator==(const Candidate& that ) const {
        if(this->pos == that.pos &&
           this->pos_end == that.pos_end &&
           this->allele.ref == that.allele.ref &&
           this->allele.alt == that.allele.alt &&
           this->allele.alt_type == that.allele.alt_type)
            return true;
        return false;
    }

};

class CandidateFinder {
    long long region_start;
    long long region_end;
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<Candidate, int> AlleleFrequencyMap;
    vector< set<Candidate> > AlleleMap;
public:
    CandidateFinder(string reference_sequence,
                    string chromosome_name,
                    long long region_start,
                    long long region_end,
                    long long ref_start,
                    long long ref_end);
    void add_read_alleles(type_read &read, vector<int> &coverage);
    pair<set<long long>, map<long long, PositionalCandidateRecord>  >  find_candidates(vector<type_read> reads);
};



#endif //FRIDAY_CANDIDATE_FINDER_H
