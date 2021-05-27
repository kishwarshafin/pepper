//
// Created by Kishwar Shafin
//

#include "candidate_finder.h"

#include <utility>


CandidateFinder::CandidateFinder(string reference_sequence,
                                 string chromosome_name,
                                 long long region_start,
                                 long long region_end,
                                 long long ref_start,
                                 long long ref_end) {
    this->reference_sequence = std::move(reference_sequence);
    this->region_start = region_start;
    this->region_end = region_end;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    this->chromosome_name = std::move(chromosome_name);
    AlleleMap.resize(region_end - region_start + 1);
}

void CandidateFinder::add_read_alleles(type_read &read, vector<int> &coverage) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index, region_index;

    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < region_start) {
                    cigar_index = min(region_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
                    region_index = ref_position - region_start;

                    if (ref_position >= region_start && ref_position <= region_end &&
                        reference_sequence[reference_index] != read.sequence[read_index] &&
                        read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        //look forward and make sure this is not an anchor base
                        bool check_this_base = true;
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size() && CandidateFinder_options::report_indels) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage[region_index] += 1;
                                check_this_base = false;
                            }
                        }
                        if(check_this_base) {
                            // process the SNP allele here
                            string ref(1, reference_sequence[reference_index]);
                            string alt(1, read.sequence[read_index]);
                            Candidate candidate_alt(ref_position, ref_position + 1, ref, alt, AlleleType::SNP_ALLELE);

                            if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                                AlleleFrequencyMap[candidate_alt] += 1;
                            } else {
                                AlleleFrequencyMap[candidate_alt] = 1;
                            }

                            if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                                AlleleMap[region_index].insert(candidate_alt);

                            coverage[region_index] += 1;
                        }
                    } else if (ref_position <= region_end &&
                               read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        coverage[region_index] += 1;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;


                if (ref_position - 1 >= region_start &&
                    ref_position - 1 <= region_end && CandidateFinder_options::report_indels) {
                    // process insert allele here
                    string ref = reference_sequence.substr(reference_index, 1);
                    string alt;
                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref + read.sequence.substr(read_index, cigar.length);

                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }

                    if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                        AlleleMap[region_index].insert(candidate_alt);

//                    cout<<"INSERT: "<<ref_position-1<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<base_quality<<endl;
                }
                read_index += cigar.length;
                break;

            case CIGAR_OPERATIONS::DEL:
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;

                if (ref_position - 1 >= region_start && ref_position - 1 <= region_end &&
                    ref_position + cigar.length < ref_end && CandidateFinder_options::report_indels) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, cigar.length + 1);
                    string alt;

                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, 1);
                    else alt = reference_sequence.substr(ref_position - ref_start - 1, 1);

                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length + 1, ref, alt,
                                            AlleleType::DELETE_ALLELE);

                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }

//                    cout<<"DEL: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;

                    if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                        AlleleMap[region_index].insert(candidate_alt);
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;

        }
    }
}

int get_index_from_base(char base) {
    if(base == '*')
        return 0;
    if(base == 'A')
        return 1;
    if(base == 'C')
        return 2;
    if(base == 'G')
        return 3;
    if(base == 'T')
        return 4;
    return  -1;
}

int get_index_from_type(char base) {
    if(base == 'R')
        return 0;
    if(base == 'S')
        return 1;
    if(base == 'I')
        return 2;
    if(base == 'D')
        return 3;
    return  -1;
}

int get_genotype(string type_predicted) {
    if(type_predicted[0] == 'R' || type_predicted[1] == 'R') {
        if(type_predicted[0] == type_predicted[1]) return 0; // hom-var
        if(type_predicted[0] != type_predicted[1]) return 1; // het
    } else {
        if(type_predicted[0] == type_predicted[1]) return 2; // hom-alt
        if(type_predicted[0] != type_predicted[1]) return 1; // het
    }
    return 0;
}

int get_genotype_from_base(char ref_base, char base_prediction1, char base_prediction2) {
    if(ref_base == base_prediction1 || ref_base == base_prediction2) {
        if(base_prediction1 == base_prediction2) return 0; // hom-var
        else return 1;
    } else if(base_prediction1 == base_prediction2) {
        return 2; // hom-alt
    } else {
        return 1; // het
    }
}

bool CandidateFinder::filter_candidate(const Candidate& candidate, bool freq_based, double freq) {
    double allele_frequency = candidate.read_support / max(1.0, double(candidate.depth));
    return false;
}


void CandidateFinder::add_read_alleles_consensus(type_read &read, vector<int> &coverage, vector<int> &insert_count, vector<int> &delete_count, vector<int> &snp_count) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index, region_index;

    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < region_start) {
                    cigar_index = min(region_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
                    region_index = ref_position - region_start;

                    if (ref_position >= region_start && ref_position <= region_end &&
                        reference_sequence[reference_index] != read.sequence[read_index] &&
                        read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        //look forward and make sure this is not an anchor base
                        bool check_this_base = true;
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size()) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage[region_index] += 1;
                                check_this_base = false;
                            }
                        }
                        if(check_this_base) {
                            // process the SNP allele here
                            snp_count[region_index] += 1;
                            coverage[region_index] += 1;
                        }
                    } else if (ref_position <= region_end && read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        coverage[region_index] += 1;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
                                                 read.base_qualities.begin() + (read_index + cigar.length));
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;
                insert_count[region_index] += 1;

                if (ref_position - 1 >= region_start &&
                    ref_position - 1 <= region_end) {
                    // process insert allele here
                    string ref = reference_sequence.substr(reference_index, 1);
                    string alt;
                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref + read.sequence.substr(read_index, cigar.length);


//                    cout<<"INSERT: "<<ref_position-1<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<base_quality<<endl;
                }
                read_index += cigar.length;
                break;

            case CIGAR_OPERATIONS::DEL:
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;
                if(ref_position - 1 >= region_start && ref_position - 1 <= region_end) insert_count[region_index] += 1;

                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;

        }
    }
}

vector<long long> CandidateFinder::find_candidates_consensus(vector <type_read>& reads, double freq) {

    // populate all the prediction maps
    vector<int> coverage(region_end - region_start + 1, 0);
    vector<int> allele_ends(region_end - region_start + 1, 0);
    vector<int> insert_count(region_end - region_start + 1, 0);
    vector<int> snp_count(region_end - region_start + 1, 0);
    vector<int> delete_count(region_end - region_start + 1, 0);
    vector<PositionalCandidateRecord> all_records;
    int read_index = 0;

    for (auto &read:reads) {
        add_read_alleles_consensus(read, coverage, insert_count, delete_count, snp_count);
        read_index += 1;
    }

    vector<long long> positions;
    int local_region_size = (int)(region_end - region_start);

    for(int pos_index = 0; pos_index < local_region_size; pos_index++) {
        if(coverage[pos_index] == 0)continue;
//        cout<<region_start+pos_index<<" SNP: "<<snp_count[pos_index]<<" IN: "<<insert_count[pos_index]<<" DEL: "<<delete_count[pos_index]<<" "<<coverage[pos_index]<<endl;
        double insert_frequency = (double) insert_count[pos_index] / (double) coverage[pos_index];
        double delete_frequency = (double) delete_count[pos_index] / (double) coverage[pos_index];
        double snp_frequency = (double) snp_count[pos_index] / (double) coverage[pos_index];
        if(insert_frequency >= freq || delete_frequency >= freq || snp_frequency >= freq) {
            positions.push_back(region_start+pos_index);
        }
    }
    return positions;
}

vector<PositionalCandidateRecord> CandidateFinder::find_candidates(vector <type_read>& reads,
                                                                   vector<long long> positions,
                                                                   vector< vector<float> > predictions,
                                                                   vector< vector<float> > type_predictions,
                                                                   vector<int> base_labels,
                                                                   vector<int> type_labels,
                                                                   bool freq_based,
                                                                   double freq) {

    // populate all the prediction maps
    vector<string> decoded_base_lables {"##", "AA", "AC", "AT", "AG", "A*", "A#", "CC", "CT", "CG", "C*", "C#", "TT", "TG", "T*", "T#", "GG", "G*", "G#", "**", "*#"};
    vector<string> decoded_type_lables {"RR", "RS", "RI", "RD", "SS", "SI", "SD", "II", "ID", "DD" };

    // first go over and see how many base positions are there.
    // all of the positions will be relative to the positions vector so we need to count where it starts and ends
    long long local_region_start = region_start;
    long long local_region_end = region_end;

    int local_region_size = (int) (local_region_end - local_region_start + 1);

    vector<int> max_observed_insert;
    vector<uint64_t> cumulative_observed_insert;
    uint64_t total_observered_insert_bases = 0;

    max_observed_insert.resize(local_region_size + 1, 0);
    cumulative_observed_insert.resize(local_region_size + 1, 0);

//    for(int i=0; i < (int) positions.size(); i++) {
//        if(positions[i] < 0) continue;
//        long long pos = positions[i];
////        max_observed_insert[pos-local_region_start] = max(max_observed_insert[pos-local_region_start]);
//    }

//    total_observered_insert_bases = max_observed_insert[0];
//    for(int i=1;i < max_observed_insert.size(); i++) {
//        cumulative_observed_insert[i] = cumulative_observed_insert[i-1] + max_observed_insert[i-1];
//        total_observered_insert_bases += max_observed_insert[i];
//    }

    // create a prediction map, this is for all positions and 5 base predictions
    vector<int> prediction_base_map(local_region_size + total_observered_insert_bases + 1);
    vector<int> prediction_type_map(local_region_size + total_observered_insert_bases + 1);

    vector< vector<float> > prediction_base_values_map;
    vector< vector<float> > prediction_type_values_map;

    prediction_base_values_map.resize(local_region_size + total_observered_insert_bases + 1, vector<float>((int)decoded_base_lables.size(), 0));
    prediction_type_values_map.resize(local_region_size + total_observered_insert_bases + 1, vector<float>((int)decoded_type_lables.size(), 0));


    for(int i=0; i< positions.size(); i++) {
        long long position = positions[i];
        if(position < 0) continue;
        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start]);
        int predicted_base_label = base_labels[i];
        int predicted_type_label = type_labels[i];

        prediction_base_map[position_index] = predicted_base_label;
        prediction_type_map[position_index] = predicted_type_label;

        prediction_base_values_map[position_index] = predictions[i];
        prediction_type_values_map[position_index] = type_predictions[i];
//        cout<<position<<" "<<predicted_type_label<<" "<<decoded_type_lables[predicted_type_label]<<" "<<position_index<<" "<<predicted_base_label<<" "<<predicted_type_label<<endl;
    }

    // all prediction maps populated now do candidate finding
    set<long long> filtered_candidate_positions;

    vector<int> coverage(region_end - region_start + 1, 0);
    vector<int> allele_ends(region_end - region_start + 1, 0);
    vector<PositionalCandidateRecord> all_records;
    int read_index = 0;
    for (auto &read:reads) {
        add_read_alleles(read, coverage);
        read_index += 1;
    }

    // allele maps are ready now filter through candidates
    int ref_buffer = (int) (region_start - ref_start);
    for (long long i = 0; i < coverage.size(); i++) {
        allele_ends[i] = 1;
        int max_del_length = 0;

        // first figure out the longest delete and figure out the end position of the allele
        for (auto &candidate: AlleleMap[i]) {
            double freq_can = 0.0;
            if (coverage[i] > 0)
                freq_can = 100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]);

            if(candidate.allele.alt_type == DELETE_TYPE) {
                allele_ends[i] = max(allele_ends[i], (int) candidate.allele.ref.length());
                max_del_length = max(max_del_length, (int) candidate.allele.ref.length());
            }
        }

        PositionalCandidateRecord positional_record;
        positional_record.chromosome_name = this->chromosome_name;
        positional_record.pos_start = this->region_start + i;
        positional_record.pos_end = positional_record.pos_start + allele_ends[i];

        bool candidate_found = false;

        for (auto &can: AlleleMap[i]) {
            Candidate candidate = can;

            if(candidate.pos > local_region_end || candidate.pos < local_region_start) continue;

            double alt_freq =(int) (((double) AlleleFrequencyMap[candidate] / max(1.0, (double) coverage[i])) * 100.0);
            int supported_reads = AlleleFrequencyMap[candidate];


//            if(alt_freq < CandidateFinder_options::freq_threshold || supported_reads < CandidateFinder_options::min_count_threshold) continue;
//            cout<<alt_freq<<"\t"<<supported_reads<<"\t"<<candidate.pos<<" "<<candidate.pos_end<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<candidate.allele.alt_type<<" "<<candidate.genotype<<endl;

            candidate_found = true;
            filtered_candidate_positions.insert(i + this->region_start);
            // allele, depth and frequency
            candidate.set_depth_values(coverage[i], AlleleFrequencyMap[candidate]);

            double non_ref_prob = 0.0;
            double alt_prob = 0.0;
            double alt_prob_h1 = 0.0;
            double alt_prob_h2 = 0.0;
            // all set, now calculate allele_prob and non_ref_prob
            if(candidate.allele.alt_type == SNP_TYPE) {
                long long position = candidate.pos;
                long long index = 0;

                int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);
                int alt_allele_index = get_index_from_base((char)candidate.allele.alt[0]);

                int base_prediction_index = prediction_base_map[position_index];
                int type_prediction_index = prediction_type_map[position_index];

                string bases_predicted = decoded_base_lables[base_prediction_index];
                string types_predicted = decoded_type_lables[type_prediction_index];

                float base_prediction_value = prediction_base_values_map[position_index][base_prediction_index];
                float type_prediction_value = prediction_type_values_map[position_index][type_prediction_index];

//                int predicted_genotype = get_genotype(types_predicted);
                int predicted_genotype = get_genotype_from_base(candidate.allele.ref[0], bases_predicted[0], bases_predicted[1]);

                if(bases_predicted[0] == candidate.allele.alt[0] || bases_predicted[1] == candidate.allele.alt[0]) {
                    candidate.allele_probability = base_prediction_value;
                    candidate.genotype_probability = type_prediction_value;
                    candidate.genotype = predicted_genotype;
//                    cout<<"SNP: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<bases_predicted<<" "<<base_prediction_value<<" "<<types_predicted<<" "<<type_prediction_value<<" "<<predicted_genotype<<endl;
                    positional_record.candidates.push_back(candidate);
                }
            }
            else if(candidate.allele.alt_type == INSERT_TYPE) {
//                string alt_allele = candidate.allele.alt;
//                long long pos = candidate.pos;
//                int length = 0;
//                alt_prob = 0.0;
//                non_ref_prob = 0.0;
//
//
////                cout<<"###################"<<endl;
//                for(int index=1; index <= max_observed_insert[candidate.pos - local_region_start]; index++) {
//
//                    if(index < candidate.allele.alt.length()) {
//                        int alt_allele_index = 0;
//                        long long position = candidate.pos;
//                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);
//
////                        cout << "INDEX: " << index << " " << candidate.allele.alt[index] << endl;
//                        alt_allele_index = get_index_from_base(candidate.allele.alt[index]);
//                        int alt_type_index = get_index_from_type('I');
//
//                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0.0)));
//                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0.0)));
//
//                        double prob_alt_h1 = prediction_map_h1[position_index][alt_allele_index] / max(1.0, sum_h1_probs);
//                        double prob_alt_h2 = prediction_map_h2[position_index][alt_allele_index] / max(1.0, sum_h2_probs);
//
//                        double sum_h1_type_probs = max(1.0, double(accumulate(prediction_map_type_h1[position_index].begin(), prediction_map_type_h1[position_index].end(), 0.0)));
//                        double sum_h2_type_probs = max(1.0, double(accumulate(prediction_map_type_h2[position_index].begin(), prediction_map_type_h2[position_index].end(), 0.0)));
//
//                        double prob_alt_type_h1 = prediction_map_type_h1[position_index][alt_type_index] / max(1.0, sum_h1_type_probs);
//                        double prob_alt_type_h2 = prediction_map_type_h2[position_index][alt_type_index] / max(1.0, sum_h2_type_probs);
//
//
////                        cout << "BASE PREDICTION VECTORS: " << endl;
////                        for (int j = 0; j < prediction_map_h1[position_index].size(); j++)cout << prediction_map_h1[position_index][j] << " ";
////                        cout << endl;
////                        for (int j = 0; j < prediction_map_h2[position_index].size(); j++)cout << prediction_map_h2[position_index][j] << " ";
////                        cout << endl;
////
////
////                        cout << "TYPE PREDICTION VECTORS: " << endl;
////                        for (int j = 0; j < prediction_map_type_h1[position_index].size(); j++)cout << prediction_map_type_h1[position_index][j] << " ";
////                        cout << endl;
////                        for (int j = 0; j < prediction_map_type_h2[position_index].size(); j++)cout << prediction_map_type_h2[position_index][j] << " ";
////                        cout << endl;
////                        cout<<sum_h1_type_probs<<" "<<sum_h2_type_probs<<endl;
////                        cout<<sum_h1_probs<<" "<<sum_h2_probs<<endl;
//                        alt_prob += max(prob_alt_h1, prob_alt_h2);
//                        non_ref_prob += max(prob_alt_type_h1, prob_alt_type_h2);
//                        length += 1;
////                        cout<<"ALT PROB CALCULATION (inside): "<<position<<" "<<alt_prob<<" ("<<prob_alt_h1<<","<<prob_alt_h2<<")"<<"\t\t";
////                        cout<<"NRF PROB CALCULATION (inside): "<<position<<" "<<non_ref_prob<<" ("<<prob_alt_type_h1<<", "<<prob_alt_type_h2<<")"<<endl;
//                    }
////                    else {
////                        int alt_allele_index = 0;
////                        long long position = candidate.pos;
////                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);
////
//////                        cout << "INDEX: " << index << endl;
////                        alt_allele_index = get_index_from_base('*');
////                        int alt_type_index = get_index_from_type('I');
////
////                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0.0)));
////                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0.0)));
////
////                        double prob_alt_h1 = prediction_map_h1[position_index][alt_allele_index] / max(1.0, sum_h1_probs);
////                        double prob_alt_h2 = prediction_map_h2[position_index][alt_allele_index] / max(1.0, sum_h2_probs);
////
////                        double sum_h1_type_probs = max(1.0, double(accumulate(prediction_map_type_h1[position_index].begin(), prediction_map_type_h1[position_index].end(), 0.0)));
////                        double sum_h2_type_probs = max(1.0, double(accumulate(prediction_map_type_h2[position_index].begin(), prediction_map_type_h2[position_index].end(), 0.0)));
////
////                        double prob_alt_type_h1 = prediction_map_type_h1[position_index][alt_type_index] / max(1.0, sum_h1_type_probs);
////                        double prob_alt_type_h2 = prediction_map_type_h2[position_index][alt_type_index] / max(1.0, sum_h2_type_probs);
////
////
//////                        cout << "BASE PREDICTION VECTORS: " << endl;
//////                        for (int j = 0; j < prediction_map_h1[position_index].size(); j++)cout << prediction_map_h1[position_index][j] << " ";
//////                        cout << endl;
//////                        for (int j = 0; j < prediction_map_h2[position_index].size(); j++)cout << prediction_map_h2[position_index][j] << " ";
//////                        cout << endl;
//////
//////
//////                        cout << "TYPE PREDICTION VECTORS: " << endl;
//////                        for (int j = 0; j < prediction_map_type_h1[position_index].size(); j++)cout << prediction_map_type_h1[position_index][j] << " ";
//////                        cout << endl;
//////                        for (int j = 0; j < prediction_map_type_h2[position_index].size(); j++)cout << prediction_map_type_h2[position_index][j] << " ";
//////                        cout << endl;
//////                        cout<<sum_h1_type_probs<<" "<<sum_h2_type_probs<<endl;
//////                        cout<<sum_h1_probs<<" "<<sum_h2_probs<<endl;
////                        alt_prob += max(prob_alt_h1, prob_alt_h2);
////                        // penalize the insertion based on the likelihood of observing an insertion outside the contained index
////                        non_ref_prob -= max(prob_alt_type_h1, prob_alt_type_h2);
////                        length += 1;
//////                        cout<<"ALT PROB CALCULATION (outside): "<<position<<" "<<alt_prob<<" ("<<prob_alt_h1<<","<<prob_alt_h2<<")"<<"\t\t";
//////                        cout<<"NRF PROB CALCULATION (outside): "<<position<<" "<<non_ref_prob<<" ("<<prob_alt_type_h1<<", "<<prob_alt_type_h2<<")"<<endl;
////
////                    }
//
//                }
////                cout<<"###################"<<endl;
//                alt_prob = alt_prob / max(1.0, (double)length);
//                non_ref_prob = non_ref_prob / max(1.0, (double)length);
//                candidate.alt_prob = alt_prob;
//                candidate.non_ref_prob = non_ref_prob;
////                cout<<"IN: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<alt_prob<<" "<<non_ref_prob<<endl;
////                cout<<"-----------------------"<<endl;
            }
            else if(candidate.allele.alt_type == DELETE_TYPE) {
                // The inference mode we are using for delete works like this:
                // we have two prediction vectors : base_prediction, type_prediction
                // for the positions that fall within the delete allele we take the prediction of seeing a deletion
                // for the positions outside, we take the maximum likelihood of the other bases and use that to quantify the likelihood of the current allele
                // we do that until we hit the maximum length of the likelihood.
                // Finally, normalize by length.

//                int length = 0;
//                alt_prob = 0.0;
//                alt_prob_h1 = 0.0;
//                alt_prob_h2 = 0.0;
//                non_ref_prob = 0.0;

//                cout<<"##################"<<endl;
//                for(long long pos=candidate.pos; pos < candidate.pos + max_del_length; pos++) {
//                    if(candidate.pos < pos && pos < candidate.pos_end) {
//                        int del_allele_index = get_index_from_base('*');
//                        int del_type_index = get_index_from_type('D');
//
//                        long long position = pos;
//                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start]);
//
//                        double sum_h1_base_probs = 0.0;
//                        double sum_h2_base_probs = 0.0;
//
//                        for(int j=0; j<prediction_map_h1[position_index].size(); j++) {
//                            sum_h1_base_probs += prediction_map_h1[position_index][j];
//                            sum_h2_base_probs += prediction_map_h2[position_index][j];
//                        }
//
//                        double sum_h1_type_probs = 0.0;
//                        double sum_h2_type_probs = 0.0;
//
//                        for(int j=0; j<prediction_map_type_h1[position_index].size(); j++) {
//                            sum_h1_type_probs += prediction_map_type_h1[position_index][j];
//                            sum_h2_type_probs += prediction_map_type_h2[position_index][j];
//                        }
//
//                        double prob_del_h1 = (prediction_map_h1[position_index][del_allele_index]) / max(1.0, sum_h1_base_probs);
//                        double prob_del_h2 = (prediction_map_h2[position_index][del_allele_index]) / max(1.0, sum_h2_base_probs);
//
//                        double prob_del_type_h1 = (prediction_map_type_h1[position_index][del_type_index]) / max(1.0, sum_h1_type_probs);
//                        double prob_del_type_h2 = (prediction_map_type_h2[position_index][del_type_index]) / max(1.0, sum_h2_type_probs);
//                        cout<<del_type_index<<endl;
//                        cout<<prediction_map_type_h1[position_index][del_type_index]<<" "<<max(1.0, sum_h1_type_probs)<<endl;
//                        cout<<prediction_map_type_h2[position_index][del_type_index]<<" "<<max(1.0, sum_h2_type_probs)<<endl;


//                        alt_prob += max(prob_del_h1, prob_del_h2);
//                        non_ref_prob += max(prob_del_type_h1, prob_del_type_h2);
////                        cout<<"ALT PROB CALCULATION ( inside): "<<position<<" "<<alt_prob<<" ("<<prob_del_h1<<","<<prob_del_h2<<")"<<"\t\t";
////                        cout<<"NRF PROB CALCULATION ( inside): "<<position<<" "<<non_ref_prob<<" ("<<prob_del_type_h1<<", "<<prob_del_type_h2<<")"<<endl;
//
//                        length += 1;
//                    }
//                    else if(pos>=candidate.pos_end) {
//                        // the the index for the Deleteion predictions from the base and the type prediction vectors
//                        int del_allele_index = get_index_from_base('*');
//                        int del_type_index = get_index_from_type('D');
//
//                        // get the position
//                        long long position = pos;
//                        // calculate the index
//                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start]);
//
//                        // now calculate the sum of the entire prediction vectors
//                        double sum_h1_base_probs = 0.0;
//                        double sum_h2_base_probs = 0.0;
//
//                        for(int j=0; j<prediction_map_h1[position_index].size(); j++) {
//                            sum_h1_base_probs += prediction_map_h1[position_index][j];
//                            sum_h2_base_probs += prediction_map_h2[position_index][j];
//                        }
//                        double sum_h1_type_probs = 0.0;
//                        double sum_h2_type_probs = 0.0;
//
//                        for(int j=0; j<prediction_map_type_h1[position_index].size(); j++) {
//                            sum_h1_type_probs = sum_h1_type_probs + prediction_map_type_h1[position_index][j];
//                            sum_h2_type_probs = sum_h2_type_probs + prediction_map_type_h2[position_index][j];
//                        }
//
//                        // now calcucate what is the maximum value you can have instead of a deleteion at those positions
//                        double max_alt_base_prob_h1 = 0.0;
//                        double max_alt_base_prob_h2 = 0.0;
//                        for(int j=0; j<prediction_map_h1[position_index].size();j++) {
//                            if(j==del_allele_index)continue;
//                            max_alt_base_prob_h1 = max(max_alt_base_prob_h1, (double) prediction_map_h1[position_index][j]);
//                            max_alt_base_prob_h2 = max(max_alt_base_prob_h2, (double) prediction_map_h2[position_index][j]);
//                        }
//
//                        double prob_del_h1 = max_alt_base_prob_h1 / max(1.0, sum_h1_base_probs);
//                        double prob_del_h2 = max_alt_base_prob_h2 / max(1.0, sum_h2_base_probs);
//
//                        // And also calculate the maximum value if not a delete type
//                        double max_alt_type_prob_h1 = 0.0;
//                        double max_alt_type_prob_h2 = 0.0;
//                        for(int j=0; j<prediction_map_type_h1[position_index].size();j++) {
//                            if(j==del_type_index)continue;
//                            max_alt_type_prob_h1 = max(max_alt_type_prob_h1, (double) prediction_map_type_h1[position_index][j]);
//                            max_alt_type_prob_h2 = max(max_alt_type_prob_h2, (double) prediction_map_type_h2[position_index][j]);
//                        }
//
//                        double prob_del_type_h1 = max_alt_type_prob_h1 / max(1.0, sum_h1_type_probs);
//                        double prob_del_type_h2 = max_alt_type_prob_h2 / max(1.0, sum_h2_type_probs);
////                        cout<<del_type_index<<endl;
////                        cout<<prediction_map_type_h1[position_index][del_type_index]<<" "<<max(1.0, sum_h1_type_probs)<<endl;
////                        cout<<prediction_map_type_h2[position_index][del_type_index]<<" "<<max(1.0, sum_h2_type_probs)<<endl;
//                        double non_del_alt_prob = max(prob_del_h1, prob_del_h2);
//                        double non_del_type_prob = max(prob_del_type_h1, prob_del_type_h2);
//
//                        alt_prob += non_del_alt_prob;
//                        non_ref_prob += non_del_type_prob;
////                        cout<<"ALT PROB CALCULATION (outside): "<<position<<" "<<alt_prob<<" ("<<prob_del_h1<<","<<prob_del_h2<<")"<<"\t\t";
////                        cout<<"NRF PROB CALCULATION (outside): "<<position<<" "<<non_ref_prob<<" ("<<prob_del_type_h1<<", "<<prob_del_type_h2<<")"<<endl;
//
//                        length += 1;
//                    }
//                }
//                cout<<"LENGTH: "<<length<<endl;
//                cout<<"##################"<<endl;

//                alt_prob = alt_prob / max(1.0, (double)length);
//                non_ref_prob = non_ref_prob / max(1.0, (double)length);
//                non_ref_prob_h1 = non_ref_prob_h1 / max(1.0, non_ref_length);
//                non_ref_prob_h2 = non_ref_prob_h2 / max(1.0, non_ref_length);
//                non_ref_prob = max(non_ref_prob_h1, non_ref_prob_h2);

//                cout<<"DEL: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" Alt-prob: "<<alt_prob<<" non-ref-prob: "<<non_ref_prob<<" Read support: "<<AlleleFrequencyMap[candidate]<<" Allele freq: "<<alt_freq<<endl;
//                cout<<"-----------------------"<<endl;
//                cout<<"DEL: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<alt_prob_h1<<" "<<alt_prob_h2<<" "<<non_ref_prob<<endl;
//                candidate.alt_prob = alt_prob;
//                candidate.alt_prob_h1 = alt_prob_h1;
//                candidate.alt_prob_h2 = alt_prob_h2;
//                candidate.non_ref_prob = non_ref_prob;
            }

//            if(filter_candidate(candidate, freq_based, freq)) positional_record.candidates.push_back(candidate);
        }

        if (!candidate_found) continue;

        all_records.push_back(positional_record);
    }

    return all_records;
}