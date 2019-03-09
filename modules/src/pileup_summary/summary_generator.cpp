//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/pileup_summary/summary_generator.h"

SummaryGenerator::SummaryGenerator(string reference_sequence, string chromosome_name, long long ref_start,
                                   long long ref_end) {
    this->reference_sequence = reference_sequence;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
}


int get_feature_index(char base, bool is_reverse, bool is_ref_base) {
    base = toupper(base);
    if(is_ref_base) {
        if (base == 'A') return 0;
        if (base == 'C') return 1;
        if (base == 'G') return 2;
        if (base == 'T') return 3;
        return 4;
    } else if(is_reverse) {
        if (base == 'A') return 5;
        if (base == 'C') return 6;
        if (base == 'G') return 7;
        if (base == 'T') return 8;
        return 13;
    } else {
        if (base == 'A') return 9;
        if (base == 'C') return 10;
        if (base == 'G') return 11;
        if (base == 'T') return 12;
        return 14;
    }
}


int get_labels(char base) {
    base = toupper(base);
    if (base == 'A') return 1;
    if (base == 'C') return 2;
    if (base == 'G') return 3;
    if (base == 'T') return 4;
    if (base == '*') return 5;
    if (base == '#') return 0;
    return 0;
}


void SummaryGenerator::iterate_over_read(type_read read, long long region_start, long long region_end) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > region_end) break;
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
                    //read.base_qualities[read_index] base quality
                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];

                        // update the summary of base
                        base_summaries[make_pair(ref_position,
                                                 get_feature_index(base, read.flags.is_reverse, false))] += 1.0;
                        coverage[ref_position] += 1.0;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
//                base_qualities = read.base_qualities.begin() + read_index, read.base_qualities.begin() + (read_index + cigar.length);
                reference_index = ref_position - ref_start - 1;

                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    alt = read.sequence.substr(read_index, cigar.length);
                    for (int i = 0; i < cigar.length; i++) {
                        pair<long long, int> position_pair = make_pair(ref_position - 1, i);
                        insert_summaries[make_pair(position_pair,
                                                   get_feature_index(alt[i], read.flags.is_reverse, false))] += 1.0;
                    }
                    longest_insert_count[ref_position - 1] = std::max(longest_insert_count[ref_position - 1],
                                                                      (long long) alt.length());
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
//                base_quality = read.base_qualities[max(0, read_index)];
                reference_index = ref_position - ref_start - 1;
                // process delete allele here
                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        // update the summary of base
                        base_summaries[make_pair(ref_position + i,
                                                 get_feature_index('*', read.flags.is_reverse, false))] += 1.0;
                        coverage[ref_position + i] += 1.0;
                    }
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }
}

int SummaryGenerator::get_sequence_length(long long start_pos, long long end_pos) {
    int length = 0;
    for (int i = start_pos; i <= end_pos; i++) {
        length += 1;
        if (longest_insert_count[i] > 0) {
            length += longest_insert_count[i];
        }
    }
    return length;
}

//bool check_base(char base) {
//    if(base=='A' || base=='a' || base=='C' || base=='c' || base=='T' || base=='t' || base =='G' || base==;)
//}

void SummaryGenerator::generate_labels(type_read read, long long region_start, long long region_end) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > region_end) break;
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;

                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        base_labels[ref_position] = base;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                reference_index = ref_position - ref_start - 1;
                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    alt = read.sequence.substr(read_index, cigar.length);

                    for (int i = 0; i < longest_insert_count[ref_position - 1]; i++) {
                        char base = '#';
                        if (i < alt.length()) {
                            base = alt[i];
                        }
                        pair<long long, int> position_pair = make_pair(ref_position - 1, i);
                        insert_labels[position_pair] = base;
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
                reference_index = ref_position - ref_start - 1;

                if (ref_position >= ref_start && ref_position <= ref_end) {
                    // process delete allele here
                    for (int i = 0; i < cigar.length; i++) {
                        if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                            // DELETE
                            char base = '*';
                            base_labels[ref_position + i] = '*';
                        }
                    }
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }
}


void SummaryGenerator::debug_print(long long start_pos, long long end_pos) {
    cout << setprecision(2);
    for (int i = start_pos; i <= end_pos; i++) {
        cout << reference_sequence[i - start_pos] << "\t";
        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++)cout << "*" << "\t";
        }
    }
    cout << endl;
    for (int i = start_pos; i <= end_pos; i++) {
        cout << base_labels[i] << "\t";
        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                if (insert_labels[make_pair(i, ii)]) cout << insert_labels[make_pair(i, ii)] << "\t";
                else cout << "*" << "\t";
            }
        }
    }
    cout << endl;
    for (int i = 0; i < labels.size(); i++) {
        cout << labels[i] << "\t";
    }
    cout << endl;
    for (int i = 0; i < labels.size(); i++) {
        cout << "(" << genomic_pos[i].first << "," << genomic_pos[i].second << ")\t";
    }
    cout << endl;

    cout << "-------------" << endl;
    for (int i = 0; i < image[0].size(); i++) {
        for (int j = 0; j < image.size(); j++) {
            printf("%.3lf\t", image[j][i]);
        }
        cout << endl;
    }

}

void SummaryGenerator::generate_ref_features() {
    for(int i=ref_start; i<=ref_end; i++) {
        char base = reference_sequence[i-ref_start];
        base_summaries[make_pair(i, get_feature_index(base, false, true))] = 1.0;
    }
}


void SummaryGenerator::generate_train_summary(vector <type_read> &reads,
                                              long long start_pos,
                                              long long end_pos,
                                              type_read truth_read,
                                              bool train_mode) {
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            iterate_over_read(read, start_pos, end_pos);
        }
    }

    generate_ref_features();

    // this populates base_labels and insert_labels dictionaries
    generate_labels(truth_read, start_pos, end_pos);

    // after all the dictionaries are populated, we can simply walk through the region and generate a sequence
    for (long long pos = start_pos; pos <= end_pos; pos++) {
        labels.push_back(get_labels(base_labels[pos]));
        genomic_pos.push_back(make_pair(pos, 0));
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++) {
                genomic_pos.push_back(make_pair(pos, ii + 1));
                if (insert_labels[make_pair(pos, ii)])
                    labels.push_back(get_labels(insert_labels[make_pair(pos, ii)]));
                else labels.push_back(get_labels('#'));
            }
        }
    }
    assert(labels.size() == genomic_pos.size());

    // at this point labels and positions are generated, now generate the pileup
    for (long long i = start_pos; i <= end_pos; i++) {
        vector<double> row;
        // iterate through the summaries
        for(int j = 5; j <= 14; j++) {
            if (j > 4) row.push_back(base_summaries[make_pair(i, j)] / max(1.0, coverage[i]));
            else row.push_back(base_summaries[make_pair(i, j)]);
        }
        assert(row.size() == 10);
        image.push_back(row);

        if (longest_insert_count[i] > 0) {

            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                vector<double> ins_row;

                // iterate through the summaries
                for(int j = 5; j <= 14; j++)
                    ins_row.push_back(insert_summaries[make_pair(make_pair(i, ii), j)] / max(1.0, coverage[i]));

                assert(ins_row.size() == 10);
                image.push_back(ins_row);
            }

        }

    }
    assert(image.size() == genomic_pos.size());
//     at this point everything should be generated
//    debug_print(start_pos, end_pos);
}


void SummaryGenerator::generate_summary(vector <type_read> &reads,
                                        long long start_pos,
                                        long long end_pos) {
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            iterate_over_read(read, start_pos, end_pos);
        }
    }

    // after all the dictionaries are populated, we can simply walk through the region and generate a sequence
    // otherwise only generate the positions
    for (long long pos = start_pos; pos <= end_pos; pos++) {
        int indx = 0;
        genomic_pos.push_back(make_pair(pos, 0));
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++)
                genomic_pos.push_back(make_pair(pos, ii + 1));
        }
    }

    generate_ref_features();

    // at this point labels and positions are generated, now generate the pileup
    for (long long i = start_pos; i <= end_pos; i++) {
        vector<double> row;
        // iterate through the summaries
        for(int j = 0; j <= 14; j++) {
            if (j > 4) row.push_back(base_summaries[make_pair(i, j)] / max(1.0, coverage[i]));
            else row.push_back(base_summaries[make_pair(i, j)]);
        }
        assert(row.size() == 15);
        image.push_back(row);

        if (longest_insert_count[i] > 0) {

            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                vector<double> ins_row {0.0, 0.0, 0.0, 0.0, 1.0};

                // iterate through the summaries
                for(int j = 5; j <= 14; j++)
                    ins_row.push_back(insert_summaries[make_pair(make_pair(i, ii), j)] / max(1.0, coverage[i]));

                assert(ins_row.size() == 15);
                image.push_back(ins_row);
            }

        }

    }
    assert(image.size() == genomic_pos.size());
//     at this point everything should be generated
//    debug_print(start_pos, end_pos);
}