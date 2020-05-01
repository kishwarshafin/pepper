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
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size()) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage[region_index] += 1;
                            }
                        } else {
                            // process the SNP allele here
                            string ref(1, reference_sequence[reference_index]);
                            string alt(1, read.sequence[read_index]);
                            Candidate candidate_alt(ref_position, ref_position + 1, ref, alt, AlleleType::SNP_ALLELE);

//                            cout<<"SNP "<<ref_position<<" "<<ref_position + 1<<" "<<ref<<" "<<alt<<endl;

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


vector<PositionalCandidateRecord> CandidateFinder::find_candidates(
        vector <type_read>& reads) {

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

    int ref_buffer = region_start - ref_start;
    // get all the positions that pass the threshold
    for (long long i = 0; i < coverage.size(); i++) {
        allele_ends[i] = 1;
        // first figure out the longest delete
        for (auto &candidate: AlleleMap[i]) {
            double freq_can = 0.0;
            if (coverage[i] > 0)
                freq_can = 100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]);

            if (freq_can >= CandidateFinder_options::freq_threshold &&
                AlleleFrequencyMap[candidate] >= CandidateFinder_options::min_count_threshold) {
                if(candidate.allele.alt_type == DELETE_TYPE) {
                    allele_ends[i] = max(allele_ends[i], (int) candidate.allele.ref.length());
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
            candidate.set_depth_values(coverage[i], AlleleFrequencyMap[candidate], AlleleFrequencyMapH0[candidate], AlleleFrequencyMapH1[candidate], AlleleFrequencyMapH2[candidate]);

            // allele, depth and frequency
            positional_record.candidates.push_back(candidate);
        }

        if (!candidate_found) continue;

        all_records.push_back(positional_record);
    }

    return all_records;
}