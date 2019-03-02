//
// Created by Kishwar Shafin on 10/24/18.
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
//    cout<<CIGAR_OPERATIONS::SOFT_CLIP<<" "<<CIGAR_OPERATIONS::IN<<endl;
//    cout<<read.pos<<" "<<read.pos_end<<endl;
//    for (auto &cigar: read.cigar_tuples) {
//
//        cout<<"("<<cigar.operation<<", "<<cigar.length<<"), ";
//    }
//    cout<<endl;
    for (auto &cigar: read.cigar_tuples) {
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
//                        cout<<read.pos<<" "<<ref_position<<" "<<region_start<<" "<<i<<" ";
//                        cout<<"SNP: "<<read.query_name<<" "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<" "<<AlleleMap[region_index].size()<<endl;
                        coverage[region_index] += 1;
                    } else if (ref_position <= region_end &&
                               read.base_qualities[read_index] >= CandidateFinder_options::min_base_quality) {
                        coverage[region_index] += 1;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
//            case CIGAR_OPERATIONS::SOFT_CLIP:
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
//                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    alt = ref + read.sequence.substr(read_index, cigar.length);

                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
                    if (AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
                        AlleleFrequencyMap[candidate_alt] += 1;
                    } else {
                        AlleleFrequencyMap[candidate_alt] = 1;
                    }

                    if (AlleleMap[region_index].find(candidate_alt) == AlleleMap[region_index].end())
                        AlleleMap[region_index].insert(candidate_alt);

//                    cout<<"INSERT: "<<ref_position-1<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
                }
                read_index += cigar.length;
                break;

            case CIGAR_OPERATIONS::DEL:
//                base_quality = read.base_qualities[max(0, read_index)];
                reference_index = ref_position - ref_start - 1;
                region_index = ref_position - region_start - 1;

                if (ref_position - 1 >= region_start && ref_position - 1 <= region_end &&
                    ref_position + cigar.length < ref_end) {
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, 1);
                    string alt = ref + reference_sequence.substr(ref_position - ref_start, cigar.length);
                    Candidate candidate_alt(ref_position - 1, ref_position - 1 + cigar.length, ref, alt,
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

//                if (ref_position >= region_start && ref_position <= region_end) {
//                    for (long long pos = ref_position; pos <= min(region_end, ref_position + cigar.length); pos++)
//                        coverage[pos - region_start] += 1;
//                }

                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
//                base_quality = *std::min_element(read.base_qualities.begin() + read_index,
//                                                 read.base_qualities.begin() + (read_index + cigar.length));
//                reference_index = ref_position - ref_start - 1;
//                region_index = ref_position - region_start - 1;
//
//                if(ref_position - 1 >= region_start &&
//                   ref_position - 1 <= region_end &&
//                   base_quality >= ActiveRegionFinder_options::min_base_quality) {
//                    // process insert allele here
//                    string ref = reference_sequence.substr(reference_index, 1);
//                    string alt;
//                    if(read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
//                    else alt = ref + read.sequence.substr(read_index, cigar.length);
//
//                    Candidate candidate_alt(ref_position - 1, ref_position, ref, alt, AlleleType::INSERT_ALLELE);
//                    if(AlleleFrequencyMap.find(candidate_alt) != AlleleFrequencyMap.end()) {
//                        AlleleFrequencyMap[candidate_alt] += 1;
//                    } else {
//                        AlleleFrequencyMap[candidate_alt] = 1;
//                    }
//
//                    if(AlleleMap[region_index].find(candidate_alt) ==  AlleleMap[region_index].end())
//                        AlleleMap[region_index].insert(candidate_alt);
//
//                    cout<<"SOFT CLIP: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
//                }
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

pair <set<long long>, map<long long, PositionalCandidateRecord>> CandidateFinder::find_candidates(
        vector <type_read> reads) {

    map<long long, PositionalCandidateRecord> all_positional_candidates;
    set<long long> filtered_candidate_positions;

    vector<int> coverage(region_end - region_start + 1, 0);
    for (auto &read:reads) {
        add_read_alleles(read, coverage);
    }

    // get all the positions that pass the threshold
    for (long long i = 0; i < coverage.size(); i++) {
        vector <pair<double, Candidate>> positional_candidates;

        for (auto &candidate: AlleleMap[i]) {
            int freq_can = 0;
            if (coverage[i] > 0)
                freq_can = (int) ceil(100.0 * ((double) AlleleFrequencyMap[candidate] / (double) coverage[i]));

            if (freq_can >= CandidateFinder_options::freq_threshold &&
                AlleleFrequencyMap[candidate] >= CandidateFinder_options::min_count_threshold) {
                filtered_candidate_positions.insert(i + this->region_start);
                positional_candidates.push_back(make_pair(freq_can, candidate));
//                cout << "CANDIDATE: " << i + this->region_start << " " << candidate.allele.ref << " "
//                     << candidate.allele.alt << " " << candidate.allele.alt_type << " " << AlleleFrequencyMap[candidate]
//                     << " " << coverage[i] << " "<< freq_can << endl;
            } else {
//                cout << "SKIPPED CANDIDATE: " << i + this->region_start << " " << candidate.allele.ref << " "
//                     << candidate.allele.alt << " " << candidate.allele.alt_type << " " << AlleleFrequencyMap[candidate]
//                     << " " << coverage[i] << " "<< freq_can << endl;
            }
        }
        if (positional_candidates.empty()) continue;

        sort(positional_candidates.begin(), positional_candidates.end());

        PositionalCandidateRecord pos_candidate;
        pair<double, Candidate> candidates = positional_candidates.back();
        positional_candidates.pop_back();
        pos_candidate.set_positions(this->chromosome_name, candidates.second.pos, candidates.second.pos_end);
        pos_candidate.set_reference(candidates.second.allele.ref);
        pos_candidate.set_alt1(candidates.second.allele.alt, candidates.second.allele.alt_type);

        if (!positional_candidates.empty()) {
            candidates = positional_candidates.back();
            positional_candidates.pop_back();
            pos_candidate.set_alt2(candidates.second.allele.alt,
                                   candidates.second.allele.alt_type,
                                   candidates.second.pos_end);
        }

        all_positional_candidates[i + this->region_start] = pos_candidate;
    }

    return make_pair(filtered_candidate_positions, all_positional_candidates);
}