//
// Created by Kishwar Shafin
//

#include "candidate_finder_hp.h"

#include <utility>


CandidateFinderHP::CandidateFinderHP(string reference_sequence,
                                     string chromosome_name,
                                     long long region_start,
                                     long long region_end,
                                     long long ref_start,
                                     long long ref_end) {
    this->reference_sequence = std::move(reference_sequence);
    this->region_start = region_start;
    this->region_end = region_end;
    this->chromosome_name = std::move(chromosome_name);
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    AlleleMap.resize(region_end - region_start + 1);
}

void CandidateFinderHP::add_read_alleles(type_read &read, vector<int> &coverage) {
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
                        if(check_this_base == true) {
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
                    ref_position - 1 <= region_end) {
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

bool CandidateFinderHP::filter_candidate(const Candidate& candidate, bool freq_based, double freq) {
    double allele_frequency = candidate.read_support / max(1.0, double(candidate.depth));
    if(freq_based) {
        if(allele_frequency >= freq) return true;
        return false;
    }

    // CONDITIONS FOR INSERT
    if(candidate.allele.alt_type == SNP_TYPE) {
        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < ONTLinearRegression::SNP_LOWER_FREQ_THRESHOLD) {
            //if(allele_frequency >= 0.05) {
            //    if(allele_weight >= 0.5) return true;
            //    else return false;
            //}
            return false;
        }

        double predicted_val = allele_weight * ONTLinearRegression::SNP_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * ONTLinearRegression::SNP_NON_REF_PROB_COEF + ONTLinearRegression::SNP_BIAS_TERM;

        if(predicted_val >= ONTLinearRegression::SNP_THRESHOLD) return true;
        // if(allele_frequency >= LinearRegression::SNP_UPPER_FREQ && allele_weight >= 0.4) return true;
        return false;
    }
        // CONDITIONS FOR INSERT
    else if (candidate.allele.alt_type == INSERT_TYPE) {
        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < ONTLinearRegression::IN_LOWER_FREQ_THRESHOLD) {
            // if(allele_frequency >= 0.05 && allele_weight >= 0.8) return true;
            return false;
        }
        double predicted_val = allele_weight * ONTLinearRegression::INSERT_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * ONTLinearRegression::INSERT_NON_REF_PROB_COEF + ONTLinearRegression::INSERT_BIAS_TERM;

        if(predicted_val >= ONTLinearRegression::INSERT_THRESHOLD) return true;
        // if(allele_frequency >= LinearRegression::IN_UPPER_FREQ && allele_weight >= 0.6) return true;
        return false;
    }
        // CONDITIONS FOR DELETE
    else if (candidate.allele.alt_type == DELETE_TYPE) {
        double allele_weight = max(candidate.alt_prob_h1, candidate.alt_prob_h2);

        if(allele_frequency < ONTLinearRegression::DEL_LOWER_FREQ_THRESHOLD) {
            // if(allele_frequency >= 0.10 && allele_weight >= 0.5) return true;
            // if(allele_frequency >= 0.05 && allele_weight >= 0.8) return true;
            return false;
        }

        double predicted_val = allele_weight * ONTLinearRegression::DELETE_ALLELE_WEIGHT_COEF + candidate.non_ref_prob * ONTLinearRegression::DELETE_NON_REF_PROB_COEF + ONTLinearRegression::DELETE_BIAS_TERM;

        if(predicted_val >= ONTLinearRegression::DELETE_THRESHOLD) return true;
        // if(allele_frequency >= LinearRegression::DEL_UPPER_FREQ_THRESHOLD && allele_weight >= 0.6) return true;

        return false;
    }
    return false;
}


vector<PositionalCandidateRecord> CandidateFinderHP::find_candidates(vector <type_read>& reads,
                                                                     vector<long long> positions,
                                                                     vector<int>indices,
                                                                     const vector< vector<int> >& base_predictions_h1,
                                                                     const vector< vector<int> >& base_predictions_h2,
                                                                     bool freq_based,
                                                                     double freq) {
    // first go over and see how many base positions are there.
    // all of the positions will be relative to the positions vector so we need to count where it starts and ends
    long long local_region_start = positions[0];
    long long local_region_end = positions[0];

    for(long long & position : positions) {
        if(position < 0) continue;
        local_region_start = min(local_region_start, position);
        local_region_end = max(local_region_end, position);
    }

    int local_region_size = (int) (local_region_end - local_region_start + 1);

    vector<int> max_observed_insert;
    vector<uint64_t> cumulative_observed_insert;
    uint64_t total_observered_insert_bases = 0;

    max_observed_insert.resize(local_region_size + 1, 0);
    cumulative_observed_insert.resize(local_region_size + 1, 0);

    for(int i=0; i < (int) positions.size(); i++) {
        if(positions[i] < 0) continue;
        long long pos = positions[i];
        int index = indices[i];
        max_observed_insert[pos-local_region_start] = max(max_observed_insert[pos-local_region_start], index);
    }

    total_observered_insert_bases = max_observed_insert[0];
    for(int i=1;i < max_observed_insert.size(); i++) {
        cumulative_observed_insert[i] = cumulative_observed_insert[i-1] + max_observed_insert[i-1];
        total_observered_insert_bases += max_observed_insert[i];
    }

    // create a prediction map, this is for all positions and 5 base predictions
    vector< vector<int> > prediction_map_h1;
    vector< vector<int> > prediction_map_h2;

    prediction_map_h1.resize(local_region_size + total_observered_insert_bases + 1, vector<int>(5, 0));
    prediction_map_h2.resize(local_region_size + total_observered_insert_bases + 1, vector<int>(5, 0));


    for(int i=0; i< positions.size(); i++) {
        long long position = positions[i];
        int index = indices[i];
        if(position < 0) continue;

        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

        prediction_map_h1[position_index] = base_predictions_h1[i];
        prediction_map_h2[position_index] = base_predictions_h2[i];
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

            if(candidate.pos > local_region_end || candidate.pos < local_region_start) continue;

            double alt_freq =(int) (((double) AlleleFrequencyMap[candidate] / max(1.0, (double) coverage[i])) * 100.0);
            int supported_reads = AlleleFrequencyMap[candidate];
            if(alt_freq < CandidateFinder_options::freq_threshold || supported_reads < CandidateFinder_options::min_count_threshold) continue;

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

                non_ref_prob = 0.0;
                alt_prob_h1 = 0.0;
                alt_prob_h2 = 0.0;
                double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));
                double prob_alt_h1 = prediction_map_h1[position_index][alt_allele_index] / sum_h1_probs;
                double prob_alt_h2 = prediction_map_h2[position_index][alt_allele_index] / sum_h2_probs;

                for(index=0; index <= max_observed_insert[candidate.pos - local_region_start]; index++) {

                    int ref_allele_index = 0;

                    if(index == 0) ref_allele_index = get_index_from_base(candidate.allele.ref[0]);
                    else ref_allele_index = get_index_from_base('*');

                    pair<long long, int> ref_pos_pair (candidate.pos, index);
                    position_index = (int) (candidate.pos - local_region_start + cumulative_observed_insert[candidate.pos - local_region_start] + index);

                    sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                    sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                    double non_ref_prob_h1 = (sum_h1_probs - prediction_map_h1[position_index][ref_allele_index]) / sum_h1_probs;
                    double non_ref_prob_h2 = (sum_h2_probs - prediction_map_h2[position_index][ref_allele_index]) / sum_h2_probs;
                    non_ref_prob = max(non_ref_prob, max(non_ref_prob_h1, non_ref_prob_h2));
                }
                alt_prob_h1 = max(0.0001, prob_alt_h1);
                alt_prob_h2 = max(0.0001, prob_alt_h2);

                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
            }
            else if(candidate.allele.alt_type == INSERT_TYPE) {
                string alt_allele = candidate.allele.alt;
                long long pos = candidate.pos;
                int length = 0;
                alt_prob = 1.0;
                alt_prob_h1 = 1.0;
                alt_prob_h2 = 1.0;

                for(int index=1; index <= max_observed_insert[candidate.pos - local_region_start]; index++) {

                    int alt_allele_index = 0;
                    long long position = candidate.pos;
                    int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

                    if(index < candidate.allele.alt.length()) alt_allele_index = get_index_from_base(candidate.allele.alt[index]);
                    else alt_allele_index = get_index_from_base('*');

                    double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                    double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                    double prob_alt_h1 = (prediction_map_h1[position_index][alt_allele_index] + 0.1) / max(1.0, sum_h1_probs);
                    double prob_alt_h2 = (prediction_map_h2[position_index][alt_allele_index] + 0.1) / max(1.0, sum_h2_probs);

                    alt_prob *= max(0.0001, max(prob_alt_h1, prob_alt_h2));
                    alt_prob_h1 *= max(0.0001, prob_alt_h1);
                    alt_prob_h2 *= max(0.0001, prob_alt_h2);
                    length += 1;

                }
                alt_prob = max(0.0001, alt_prob);
                alt_prob_h1 = max(0.0001, alt_prob_h1);
                alt_prob_h2 = max(0.0001, alt_prob_h2);
                // now calculate non-ref-prob
                length = 0;
                double non_ref_prob_h1 = 0.0;
                double non_ref_prob_h2 = 0.0;
                for(int index=0; index <= min(max_observed_insert[candidate.pos - local_region_start], int(candidate.allele.alt.size() - 1)); index++) {
                    int ref_allele_index = 0;
                    if(index==0) ref_allele_index = get_index_from_base(candidate.allele.ref[0]);
                    else ref_allele_index = get_index_from_base('*');

                    long long position = candidate.pos;
                    int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

                    double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                    double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                    non_ref_prob_h1 = non_ref_prob_h1 + ((sum_h1_probs - prediction_map_h1[position_index][ref_allele_index]) / max(1.0, sum_h1_probs));
                    non_ref_prob_h2 = non_ref_prob_h2 + ((sum_h2_probs - prediction_map_h2[position_index][ref_allele_index]) / max(1.0, sum_h2_probs));
                    length += 1;
                }
                non_ref_prob_h1 = non_ref_prob_h1 / max(1, length);
                non_ref_prob_h2 = non_ref_prob_h2 / max(1, length);
                non_ref_prob = max(non_ref_prob_h1, non_ref_prob_h2);

                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
//                cout<<"IN: "<<candidate.pos<<" "<<candidate.allele.ref<<" "<<candidate.allele.alt<<" "<<alt_prob_h1<<" "<<alt_prob_h2<<" "<<non_ref_prob<<endl;
//                cout<<"-----------------------"<<endl;
            }
            else if(candidate.allele.alt_type == DELETE_TYPE) {
                int length = 0;
                double non_ref_length = 0.0;
                double non_ref_prob_h1 = 0.0;
                double non_ref_prob_h2 = 0.0;
                alt_prob = 1.0;
                alt_prob_h1 = 1.0;
                alt_prob_h2 = 1.0;

                for(long long pos=candidate.pos; pos < candidate.pos + max_del_length; pos++) {
                    if(candidate.pos < pos && pos < candidate.pos_end) {
                        int ref_allele_index = get_index_from_base(candidate.allele.ref[pos - candidate.pos]);

                        long long position = pos;
                        int index = 0;
                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                        non_ref_prob_h1 = non_ref_prob_h1 + ((sum_h1_probs - prediction_map_h1[position_index][ref_allele_index]) / max(1.0, sum_h1_probs));
                        non_ref_prob_h2 = non_ref_prob_h2 + ((sum_h2_probs - prediction_map_h2[position_index][ref_allele_index]) / max(1.0, sum_h2_probs));
                        non_ref_length += 1.0;
                    }

                    if(candidate.pos < pos && pos < candidate.pos_end) {
                        int del_allele_index = get_index_from_base('*');

                        long long position = pos;
                        int index = 0;
                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                        double prob_del_h1 = (prediction_map_h1[position_index][del_allele_index] + 0.1) / max(1.0, sum_h1_probs);
                        double prob_del_h2 = (prediction_map_h2[position_index][del_allele_index] + 0.1) / max(1.0, sum_h2_probs);
                        alt_prob *= max(0.0001, max(prob_del_h1, prob_del_h2));
                        alt_prob_h1 *= max(0.0001, max(prob_del_h1, prob_del_h2));
                        alt_prob_h2 *= max(0.0001, max(prob_del_h1, prob_del_h2));
                        length += 1;
                    } else if(pos>=candidate.pos_end) {
                        int del_allele_index = get_index_from_base('*');

                        long long position = pos;
                        int index = 0;
                        int position_index = (int) (position - local_region_start + cumulative_observed_insert[position - local_region_start] + index);

                        double sum_h1_probs = max(1.0, double(accumulate(prediction_map_h1[position_index].begin(), prediction_map_h1[position_index].end(), 0)));
                        double sum_h2_probs = max(1.0, double(accumulate(prediction_map_h2[position_index].begin(), prediction_map_h2[position_index].end(), 0)));

                        double prob_non_del_h1 = (sum_h1_probs - prediction_map_h1[position_index][del_allele_index]) / max(1.0, sum_h1_probs);
                        double prob_non_del_h2 = (sum_h2_probs - prediction_map_h2[position_index][del_allele_index]) / max(1.0, sum_h2_probs);
                        alt_prob *= max(0.0001, max(prob_non_del_h1, prob_non_del_h2));
                        alt_prob_h1 *= max(0.0001, prob_non_del_h1);
                        alt_prob_h2 *= max(0.0001, prob_non_del_h2);
                        length += 1;
                    }
                }

                alt_prob = max(0.0000001, alt_prob);
                alt_prob_h1 = max(0.0000001, alt_prob_h1);
                alt_prob_h2 = max(0.0000001, alt_prob_h2);

                non_ref_prob_h1 = non_ref_prob_h1 / max(1.0, non_ref_length);
                non_ref_prob_h2 = non_ref_prob_h2 / max(1.0, non_ref_length);
                non_ref_prob = max(non_ref_prob_h1, non_ref_prob_h2);

//                cout << "DEL: " << candidate.pos << " " << candidate.allele.ref << " " << candidate.allele.alt << " Alt-prob-1: " << alt_prob_h1 << " alt-prob-2: " << alt_prob_h2 << " non-ref-prob: " << non_ref_prob << " Read support: " << AlleleFrequencyMap[candidate] << " Allele freq: " << alt_freq << endl;
//                cout << "DEL: " << candidate.pos << " " << candidate.allele.ref << " " << candidate.allele.alt << " " << alt_prob_h1 << " " << alt_prob_h2 << " " << non_ref_prob << endl;

                candidate.alt_prob_h1 = alt_prob_h1;
                candidate.alt_prob_h2 = alt_prob_h2;
                candidate.non_ref_prob = non_ref_prob;
            }

            if(filter_candidate(candidate, freq_based, freq)) positional_record.candidates.push_back(candidate);
        }

        if (!candidate_found) continue;

        all_records.push_back(positional_record);
    }

    return all_records;
}