//
// Created by Kishwar Shafin
//

#include "../../headers/candidate_finding/candidate_finder.h"


CandidateFinder::CandidateFinder(string reference_sequence,
                                 string chromosome_name,
                                 long long region_start,
                                 long long region_end,
                                 long long ref_start,
                                 long long ref_end) {
    this->reference_sequence = reference_sequence;
    this->region_start = region_start;
    this->region_end = region_end;
    this->chromosome_name = chromosome_name;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
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
                        bool check_this_base = 1;
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size()) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage[region_index] += 1;
                                check_this_base = 0;
                            }
                        }
                        if(check_this_base == 1) {
                            // process the SNP allele here
                            string ref(1, reference_sequence[reference_index]);
                            string alt(1, read.sequence[read_index]);
                            Candidate candidate_alt(ref_position, ref_position + 1, ref, alt, AlleleType::SNP_ALLELE);

                            if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                                AlleleFrequencyMap[candidate_alt] += 1;

                                if(read.hp_tag == 0) {
                                    AlleleFrequencyMapH0[candidate_alt] += 1;
                                } else if(read.hp_tag == 1) {
                                    AlleleFrequencyMapH1[candidate_alt] += 1;
                                } else if(read.hp_tag == 2) {
                                    AlleleFrequencyMapH2[candidate_alt] += 1;
                                }
                            } else {
                                AlleleFrequencyMap[candidate_alt] = 1;
                                AlleleFrequencyMapH0[candidate_alt] = 0;
                                AlleleFrequencyMapH1[candidate_alt] = 0;
                                AlleleFrequencyMapH2[candidate_alt] = 0;

                                if(read.hp_tag == 0) {
                                    AlleleFrequencyMapH0[candidate_alt] += 1;
                                } else if(read.hp_tag == 1) {
                                    AlleleFrequencyMapH1[candidate_alt] += 1;
                                } else if(read.hp_tag == 2) {
                                    AlleleFrequencyMapH2[candidate_alt] += 1;
                                }
                            }

                            if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                                AlleleMap[region_index].insert(candidate_alt);

                            coverage[region_index] += 1;
//                            cout<<"SNP: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<read.base_qualities[read_index]<<" "<<read.mapping_quality<<" "<<read.flags.is_supplementary<<endl;
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
                    ref_position - 1 <= region_end) {
                    // process insert allele here
                    string ref = reference_sequence.substr(reference_index, 1);
                    string alt;
                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref + read.sequence.substr(read_index, cigar.length);

                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;

                        if(read.hp_tag == 0) {
                            AlleleFrequencyMapH0[candidate_alt] += 1;
                        } else if(read.hp_tag == 1) {
                            AlleleFrequencyMapH1[candidate_alt] += 1;
                        } else if(read.hp_tag == 2) {
                            AlleleFrequencyMapH2[candidate_alt] += 1;
                        }
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;

                        AlleleFrequencyMapH0[candidate_alt] = 0;
                        AlleleFrequencyMapH1[candidate_alt] = 0;
                        AlleleFrequencyMapH2[candidate_alt] = 0;

                        if(read.hp_tag == 0) {
                            AlleleFrequencyMapH0[candidate_alt] += 1;
                        } else if(read.hp_tag == 1) {
                            AlleleFrequencyMapH1[candidate_alt] += 1;
                        } else if(read.hp_tag == 2) {
                            AlleleFrequencyMapH2[candidate_alt] += 1;
                        }
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
                    ref_position + cigar.length < ref_end) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, cigar.length + 1);
                    string alt;

                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, 1);
                    else alt = reference_sequence.substr(ref_position - ref_start - 1, 1);

                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length + 1, ref, alt,
                                            AlleleType::DELETE_ALLELE);

                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                        if(read.hp_tag == 0) {
                            AlleleFrequencyMapH0[candidate_alt] += 1;
                        } else if(read.hp_tag == 1) {
                            AlleleFrequencyMapH1[candidate_alt] += 1;
                        } else if(read.hp_tag == 2) {
                            AlleleFrequencyMapH2[candidate_alt] += 1;
                        }
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;

                        AlleleFrequencyMapH0[candidate_alt] = 0;
                        AlleleFrequencyMapH1[candidate_alt] = 0;
                        AlleleFrequencyMapH2[candidate_alt] = 0;

                        if(read.hp_tag == 0) {
                            AlleleFrequencyMapH0[candidate_alt] += 1;
                        } else if(read.hp_tag == 1) {
                            AlleleFrequencyMapH1[candidate_alt] += 1;
                        } else if(read.hp_tag == 2) {
                            AlleleFrequencyMapH2[candidate_alt] += 1;
                        }
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


bool CandidateFinder::filter_candidate(Candidate candidate) {
    double allele_frequency = candidate.read_support / max(1.0, double(candidate.depth));

    // CONDITIONS FOR INSERT
    if(candidate.allele.alt_type == SNP_TYPE) {
//        if(allele_frequency >= 0.10) return true;
//        else return false;

        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < LinearRegression::SNP_LOWER_FREQ_THRESHOLD) {
            if(allele_frequency >= 0.05) {
                if(allele_weight >= 0.5) return true;
                else return false;
            }
            return false;
        }

        double predicted_val = allele_weight * LinearRegression::SNP_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * LinearRegression::SNP_NON_REF_PROB_COEF + LinearRegression::SNP_BIAS_TERM;

        if(predicted_val >= LinearRegression::SNP_THRESHOLD) return true;
        if(allele_frequency >= LinearRegression::SNP_UPPER_FREQ && allele_weight >= 0.4) return true;
        return false;
    }
    // CONDITIONS FOR INSERT
    else if (candidate.allele.alt_type == INSERT_TYPE) {
//        if(allele_frequency >= 0.10) return true;
//        else return false;

        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < LinearRegression::IN_LOWER_FREQ_THRESHOLD) {
            if(allele_frequency >= 0.05 && allele_weight >= 0.8) return true;
            return false;
        }
        double predicted_val = allele_frequency * LinearRegression::INSERT_ALT_FREQ_COEF + allele_weight * LinearRegression::INSERT_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * LinearRegression::INSERT_NON_REF_PROB_COEF + LinearRegression::INSERT_BIAS_TERM;

        if(predicted_val >= LinearRegression::INSERT_THRESHOLD) return true;
        if(allele_frequency >= LinearRegression::IN_UPPER_FREQ && allele_weight >= 0.6) return true;
        return false;
    }
    // CONDITIONS FOR DELETE
    else if (candidate.allele.alt_type == DELETE_TYPE) {
//        if(allele_frequency >= 0.20) return true;
//        else return false;

        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < LinearRegression::DEL_LOWER_FREQ_THRESHOLD) {
            if(allele_frequency >= 0.10 && allele_weight >= 0.5) return true;
            if(allele_frequency >= 0.05 && allele_weight >= 0.8) return true;

            return false;
        }

        double predicted_val = allele_weight * LinearRegression::DELETE_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * LinearRegression::DELETE_NON_REF_PROB_COEF + LinearRegression::DELETE_BIAS_TERM;

        if(predicted_val >= LinearRegression::DELETE_THRESHOLD) return true;
        if(allele_frequency >= LinearRegression::DEL_UPPER_FREQ_THRESHOLD && allele_weight >= 0.6) return true;

        return false;
    }
    return false;
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

vector<PositionalCandidateRecord> CandidateFinder::find_candidates(vector <type_read>& reads, vector<long long> positions, vector<int>indices, vector< vector<int> > base_predictions_h1, vector< vector<int> > base_predictions_h2) {

    // populate all the prediction maps

    map < pair<long long, int>, vector<int> > prediction_map_h1;
    map < pair<long long, int>, vector<int> > prediction_map_h2;
    map < long long, int > max_index_map;

    for(int i=0; i<positions.size(); i++) {
        pair<long long, int> pos_pair (positions[i], indices[i]);
        prediction_map_h1[pos_pair] = base_predictions_h1[i];
        prediction_map_h2[pos_pair] = base_predictions_h2[i];

        if(max_index_map.find(positions[i]) != max_index_map.end()) {
            max_index_map[positions[i]] = max(indices[i], max_index_map[positions[i]]);
        } else{
            max_index_map[positions[i]] = indices[i];
        }
    }

//    // all prediction maps populated now do candidate finding
    map<long long, vector <Candidate> > all_positional_candidates;
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
    int ref_buffer = region_start - ref_start;
    for (long long i = 0; i < coverage.size(); i++) {
        allele_ends[i] = 1;
        int max_del_length = 0;
        // first figure out the longest delete and figure out the end position of the allele
        for (auto &candidate: AlleleMap[i]) {
            double freq_can = 0.0;
            if (coverage[i] > 0)
                freq_can = 100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]);

            if (freq_can >= CandidateFinder_options::freq_threshold &&
                AlleleFrequencyMap[candidate] >= CandidateFinder_options::min_count_threshold) {
                if(candidate.allele.alt_type == DELETE_TYPE) {
                    allele_ends[i] = max(allele_ends[i], (int) candidate.allele.ref.length());
                    max_del_length = max(max_del_length, (int) candidate.allele.ref.length());
                }
            }
        }

        PositionalCandidateRecord positional_record;
        positional_record.chromosome_name = this->chromosome_name;
        positional_record.pos_start = this->region_start + i;
        positional_record.pos_end = positional_record.pos_start + allele_ends[i];

        bool candidate_found = false;

        for (auto &can: AlleleMap[i]) {
            Candidate candidate = can;
            double alt_freq =(int) (((double) AlleleFrequencyMap[candidate] / max(1.0, (double) coverage[i])) * 100.0);
            int supported_reads = AlleleFrequencyMap[candidate];
            if(alt_freq < CandidateFinder_options::freq_threshold || supported_reads < CandidateFinder_options::min_count_threshold) continue;

            candidate_found = true;
            filtered_candidate_positions.insert(i + this->region_start);
            // allele, depth and frequency
            candidate.set_depth_values(coverage[i], AlleleFrequencyMap[candidate], AlleleFrequencyMapH0[candidate], AlleleFrequencyMapH1[candidate], AlleleFrequencyMapH2[candidate]);

            double non_ref_prob;
            double alt_prob_h1;
            double alt_prob_h2;
            // all set, now calculate allele_prob and non_ref_prob
            if(candidate.allele.alt_type == SNP_TYPE) {

                pair<long long, int> pos_pair (candidate.pos, 0); // index is 0 for SNPs
                int alt_allele_index = get_index_from_base((char)candidate.allele.alt[0]);
                non_ref_prob = 0.0;
                if(prediction_map_h1.find(pos_pair)==prediction_map_h1.end() || prediction_map_h2.find(pos_pair)==prediction_map_h2.end()) {
                    non_ref_prob = 0.0;
                    alt_prob_h1 = 0.0;
                    alt_prob_h2 = 0.0;
                } else {
                    non_ref_prob = 0.0;
                    alt_prob_h1 = 0.0;
                    alt_prob_h2 = 0.0;
                    double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[pos_pair].begin(), prediction_map_h1[pos_pair].end(), 0)));
                    double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[pos_pair].begin(), prediction_map_h2[pos_pair].end(), 0)));
                    double prob_alt_h1 = prediction_map_h1[pos_pair][alt_allele_index] / sum_h1_probs;
                    double prob_alt_h2 = prediction_map_h2[pos_pair][alt_allele_index] / sum_h2_probs;

                    for(int index=0; index <= max_index_map[candidate.pos]; index++) {

                        int ref_allele_index = 0;

                        if(index == 0) ref_allele_index = get_index_from_base(candidate.allele.ref[0]);
                        else ref_allele_index = get_index_from_base('*');

                        pair<long long, int> ref_pos_pair (candidate.pos, index);
                        sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ref_pos_pair].begin(), prediction_map_h1[ref_pos_pair].end(), 0)));
                        sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ref_pos_pair].begin(), prediction_map_h2[ref_pos_pair].end(), 0)));

                        double non_ref_prob_h1 = (sum_h1_probs - prediction_map_h1[ref_pos_pair][ref_allele_index]) / sum_h1_probs;
                        double non_ref_prob_h2 = (sum_h2_probs - prediction_map_h2[ref_pos_pair][ref_allele_index]) / sum_h2_probs;
                        non_ref_prob = max(non_ref_prob, max(non_ref_prob_h1, non_ref_prob_h2));
                    }
                    alt_prob_h1 = max(0.0001, prob_alt_h1);
                    alt_prob_h2 = max(0.0001, prob_alt_h2);
                }

                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
//                cout<<"SNP: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" Alt-prob-1: "<<alt_prob_h1<<" alt-prob-2: "<<alt_prob_h2<<" non-ref-prob: "<<non_ref_prob<<" Read support: "<<AlleleFrequencyMap[candidate]<<" Allele freq: "<<alt_freq<<endl;
//                cout<<"------------------"<<endl;
            }
            else if(candidate.allele.alt_type == INSERT_TYPE) {
                string alt_allele = candidate.allele.alt;
                long long pos = candidate.pos;
                int length = 0;
                alt_prob_h1 = 1.0;
                alt_prob_h2 = 1.0;
//                cout<<"IN: "<<candidate.pos<<" "<<candidate.allele.alt<<endl;
//                cout<<"max index: "<<max_index_map[candidate.pos]<<endl;

                for(int index=1; index <= max_index_map[candidate.pos]; index++) {
                    int alt_allele_index = 0;
                    pair<long long, int> ins_pos_pair (candidate.pos, index);

                    if(index < candidate.allele.alt.length()) alt_allele_index = get_index_from_base(candidate.allele.alt[index]);
                    else alt_allele_index = get_index_from_base('*');

                    double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ins_pos_pair].begin(), prediction_map_h1[ins_pos_pair].end(), 0)));
                    double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ins_pos_pair].begin(), prediction_map_h2[ins_pos_pair].end(), 0)));

                    double prob_alt_h1 = (prediction_map_h1[ins_pos_pair][alt_allele_index]) / max(1.0, sum_h1_probs);
                    double prob_alt_h2 = (prediction_map_h2[ins_pos_pair][alt_allele_index]) / max(1.0, sum_h2_probs);

                    alt_prob_h1 *= max(0.0001, prob_alt_h1);
                    alt_prob_h2 *= max(0.0001, prob_alt_h2);
                    length += 1;
                }
                alt_prob_h1 = max(0.0001, alt_prob_h1);
                alt_prob_h2 = max(0.0001, alt_prob_h2);
                // now calculate non-ref-prob
                length = 0;
                double non_ref_prob_h1 = 0.0;
                double non_ref_prob_h2 = 0.0;
                for(int index=0; index <= min(max_index_map[candidate.pos], int(candidate.allele.alt.size() - 1)); index++) {
                    int ref_allele_index = 0;
                    if(index==0) ref_allele_index = get_index_from_base(candidate.allele.ref[0]);
                    else ref_allele_index = get_index_from_base('*');

                    pair<long long, int> ref_pos_pair (candidate.pos, index);

                    double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ref_pos_pair].begin(), prediction_map_h1[ref_pos_pair].end(), 0)));
                    double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ref_pos_pair].begin(), prediction_map_h2[ref_pos_pair].end(), 0)));

                    non_ref_prob_h1 = non_ref_prob_h1 + ((sum_h1_probs - prediction_map_h1[ref_pos_pair][ref_allele_index]) / max(1.0, sum_h1_probs));
                    non_ref_prob_h2 = non_ref_prob_h2 + ((sum_h2_probs - prediction_map_h2[ref_pos_pair][ref_allele_index]) / max(1.0, sum_h2_probs));
                    length += 1;
                }
                non_ref_prob_h1 = non_ref_prob_h1 / max(1, length);
                non_ref_prob_h2 = non_ref_prob_h2 / max(1, length);
                non_ref_prob = max(non_ref_prob_h1, non_ref_prob_h2);

                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
//                cout<<"IN: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" Alt-prob-1: "<<alt_prob_h1<<" alt-prob-2: "<<alt_prob_h2<<" non-ref-prob: "<<non_ref_prob<<" Read support: "<<AlleleFrequencyMap[candidate]<<" Allele freq: "<<alt_freq<<endl;
//                cout<<"------------------"<<endl;
            }
            else if(candidate.allele.alt_type == DELETE_TYPE) {
                int length = 0;
                double non_ref_length = 0.0;
                double non_ref_prob_h1 = 0.0;
                double non_ref_prob_h2 = 0.0;
                double alt_prob_h1 = 1.0;
                double alt_prob_h2 = 1.0;
                for(int pos=candidate.pos; pos < candidate.pos + max_del_length; pos++) {
                    if(candidate.pos < pos && pos < candidate.pos_end) {
                        int ref_allele_index = get_index_from_base(candidate.allele.ref[pos - candidate.pos]);

                        pair<long long, int> ref_pos_pair(pos, 0);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ref_pos_pair].begin(), prediction_map_h1[ref_pos_pair].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ref_pos_pair].begin(), prediction_map_h2[ref_pos_pair].end(), 0)));

                        non_ref_prob_h1 = non_ref_prob_h1 + ((sum_h1_probs - prediction_map_h1[ref_pos_pair][ref_allele_index]) / max(1.0, sum_h1_probs));
                        non_ref_prob_h2 = non_ref_prob_h2 + ((sum_h2_probs - prediction_map_h2[ref_pos_pair][ref_allele_index]) / max(1.0, sum_h2_probs));
                        non_ref_length += 1.0;
                    }

                    if(candidate.pos < pos && pos < candidate.pos_end) {
                        int del_allele_index = get_index_from_base('*');

                        pair<long long, int> ref_pos_pair(pos, 0);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ref_pos_pair].begin(), prediction_map_h1[ref_pos_pair].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ref_pos_pair].begin(), prediction_map_h2[ref_pos_pair].end(), 0)));

                        double prob_del_h1 = (prediction_map_h1[ref_pos_pair][del_allele_index] + 0.1) / max(1.0, sum_h1_probs);
                        double prob_del_h2 = (prediction_map_h2[ref_pos_pair][del_allele_index] + 0.1) / max(1.0, sum_h2_probs);
                        alt_prob_h1 *= max(0.0001, max(prob_del_h1, prob_del_h2));
                        alt_prob_h2 *= max(0.0001, max(prob_del_h1, prob_del_h2));
                        length += 1;
                    } else if(pos>=candidate.pos_end) {
                        int del_allele_index = get_index_from_base('*');

                        pair<long long, int> ref_pos_pair(pos, 0);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[ref_pos_pair].begin(), prediction_map_h1[ref_pos_pair].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[ref_pos_pair].begin(), prediction_map_h2[ref_pos_pair].end(), 0)));

                        double prob_non_del_h1 = (sum_h1_probs - prediction_map_h1[ref_pos_pair][del_allele_index]) / max(1.0, sum_h1_probs);
                        double prob_non_del_h2 = (sum_h2_probs - prediction_map_h2[ref_pos_pair][del_allele_index]) / max(1.0, sum_h2_probs);
                        alt_prob_h1 *= max(0.0001, prob_non_del_h1);
                        alt_prob_h2 *= max(0.0001, prob_non_del_h2);
                        length += 1;
                    }
                }

                alt_prob_h1 = alt_prob_h1;
                alt_prob_h2 = alt_prob_h2;

                non_ref_prob_h1 = non_ref_prob_h1 / max(1.0, non_ref_length);
                non_ref_prob_h2 = non_ref_prob_h2 / max(1.0, non_ref_length);
                non_ref_prob = max(non_ref_prob_h1, non_ref_prob_h2);

//                cout<<"DEL: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" Alt-prob-1: "<<alt_prob_h1<<" alt-prob-2: "<<alt_prob_h2<<" non-ref-prob: "<<non_ref_prob<<" Read support: "<<AlleleFrequencyMap[candidate]<<" Allele freq: "<<alt_freq<<endl;
//                cout<<"------------------"<<endl;

//                cout<<"DEL: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<alt_prob_h1<<" "<<alt_prob_h2<<" "<<non_ref_prob<<endl;
                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
            }

            if(filter_candidate(candidate)) positional_record.candidates.push_back(candidate);
        }

        if (!candidate_found) continue;

        all_records.push_back(positional_record);
    }

    return all_records;
}