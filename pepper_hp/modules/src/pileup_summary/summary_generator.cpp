//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/pileup_summary/summary_generator.h"

SummaryGenerator::SummaryGenerator(string reference_sequence, string chromosome_name, long long ref_start,
                                   long long ref_end) {
    this->reference_sequence = reference_sequence;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    this->chromosome_name = chromosome_name;
}


int get_feature_index(char base, bool is_reverse) {
    base = toupper(base);
    if (is_reverse) {
        if (base == 'A') return 0;
        if (base == 'C') return 1;
        if (base == 'G') return 2;
        if (base == 'T') return 3;
        return 8;
    } else {
        // tagged and forward
        if (base == 'A') return 4;
        if (base == 'C') return 5;
        if (base == 'G') return 6;
        if (base == 'T') return 7;
        return 9;
    }
}

int get_reference_feature_index(char base) {
    base = toupper(base);
    if (base == 'A') return 1;
    if (base == 'C') return 2;
    if (base == 'G') return 3;
    if (base == 'T') return 4;
    return 0;
}


uint8_t get_labels(char base) {
    base = toupper(base);
    if (base == 'A') return 1;
    if (base == 'C') return 2;
    if (base == 'G') return 3;
    if (base == 'T') return 4;
    if (base == '*') return 0; // this is for deleted bases, but the number is so small that it creates confusion
    if (base == '#') return 0;
    return 0;
}


void SummaryGenerator::iterate_over_read(type_read read, long long region_start, long long region_end) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;
    // scaled mapping quality
    float mapping_quality = min(60, read.mapping_quality) / 60.0;
    float base_quality = 0;
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
                        base_quality = min(100, read.base_qualities[read_index]) / 100.0;

                        // update the summary of base
                        if(read.hp_tag == 0) {
                            base_summaries_hp1[make_pair(ref_position, get_feature_index(base, read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            base_summaries_hp2[make_pair(ref_position, get_feature_index(base, read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp1[ref_position] += 1;
                            coverage_hp2[ref_position] += 1;
                        }
                        else if(read.hp_tag == 1) {
                            base_summaries_hp1[make_pair(ref_position, get_feature_index(base, read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp1[ref_position] += 1;
                        }
                        else if(read.hp_tag == 2) {
                            base_summaries_hp2[make_pair(ref_position, get_feature_index(base, read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp2[ref_position] += 1;
                        }

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

                        base_quality = min(100, read.base_qualities[read_index + i]) / 100.0;

                        if(read.hp_tag == 0) {
                            insert_summaries_hp1[make_pair(position_pair, get_feature_index(alt[i], read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            insert_summaries_hp2[make_pair(position_pair, get_feature_index(alt[i], read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                        }
                        else if(read.hp_tag == 1) {
                            insert_summaries_hp1[make_pair(position_pair, get_feature_index(alt[i], read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                        }
                        else if(read.hp_tag == 2)    {
                            insert_summaries_hp2[make_pair(position_pair, get_feature_index(alt[i], read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                        }
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
                base_quality = min(100, read.base_qualities[read_index]) / 100.0;
                // process delete allele here
                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        // update the summary of base
                        if(read.hp_tag == 0) {
                            base_summaries_hp1[make_pair(ref_position + i, get_feature_index('*', read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            base_summaries_hp2[make_pair(ref_position + i, get_feature_index('*', read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp1[ref_position] += 1;
                            coverage_hp2[ref_position] += 1;
                        } else if(read.hp_tag == 1) {
                            base_summaries_hp1[make_pair(ref_position + i, get_feature_index('*', read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp1[ref_position] += 1;
                        } else if(read.hp_tag == 2) {
                            base_summaries_hp2[make_pair(ref_position + i, get_feature_index('*', read.flags.is_reverse))] += (base_quality + mapping_quality) / 2.0;
                            coverage_hp2[ref_position] += 1;
                        }
                        coverage[ref_position] += 1.0;
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

bool check_base(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base =='G' || base=='g' || base == '*' || base == '#') return true;
    return false;
}

void SummaryGenerator::generate_labels(type_read read, long long region_start, long long region_end, int hp_tag) {
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
//                    cout<<ref_position<<" "<<ref_end<<" "<<region_end<<endl;
                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        if(hp_tag == 1)
                            base_labels_hp1[ref_position] = base;
                        else
                            base_labels_hp2[ref_position] = base;
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
                        if(hp_tag == 1)
                            insert_labels_hp1[position_pair] = base;
                        else
                            insert_labels_hp2[position_pair] = base;
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
                            if(hp_tag == 1)
                                base_labels_hp1[ref_position + i] = '*';
                            else
                                base_labels_hp2[ref_position + i] = '*';
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
    cout << setprecision(3);
    for (int i = start_pos; i <= end_pos; i++) {
        if(i==start_pos) cout<<"REF:\t";
        cout << "  " << reference_sequence[i - start_pos] << "\t";
        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++)cout << "  *" << "\t";
        }
    }
    cout << endl;
    for (int i = start_pos; i <= end_pos; i++) {
        if(i==start_pos) cout<<"TRH1:\t";
        cout << "  " <<base_labels_hp1[i] << "\t";
        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                if (insert_labels_hp1[make_pair(i, ii)]) cout << "  "<< insert_labels_hp1[make_pair(i, ii)] << "\t";
                else cout << "  *" << "\t";
            }
        }
    }
    cout << endl;
    for (int i = 0; i < labels_hp1.size(); i++) {
        if(i==0) cout<<"LBL1:\t";
        printf("%3d\t", labels_hp1[i]);
    }
    cout << endl;

    for (int i = start_pos; i <= end_pos; i++) {
        if(i==start_pos) cout<<"TRH2:\t";
        cout << "  " <<base_labels_hp2[i] << "\t";
        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                if (insert_labels_hp2[make_pair(i, ii)]) cout << "  "<< insert_labels_hp2[make_pair(i, ii)] << "\t";
                else cout << "  *" << "\t";
            }
        }
    }
    cout << endl;
    for (int i = 0; i < labels_hp1.size(); i++) {
        if(i==0) cout<<"LBL2:\t";
        printf("%3d\t", labels_hp2[i]);
    }
    cout << endl;

    cout<<ref_image.size()<<endl;
    for (int i = 0; i < ref_image.size(); i++) {
        if (i == 0) cout << "REF2:\t";
        printf("%3d\t", ref_image[i]);
    }
    cout << endl;

    cout<<"POS:\t";
    for (int i = 0; i < genomic_pos.size(); i++) {

//        cout << "(" << genomic_pos[i].first << "," << genomic_pos[i].second << ")\t";
        printf("%3lld\t", (genomic_pos[i].first) % 100);
    }
    cout << endl;

    cout << "-------------" << endl;
    for (int i = 0; i < image_hp1[0].size(); i++) {
        if(i==0)cout<<"AFW1:\t";
        if(i==1)cout<<"CFW1:\t";
        if(i==2)cout<<"GFW1:\t";
        if(i==3)cout<<"TFW1:\t";
        if(i==4)cout<<"ARV1:\t";
        if(i==5)cout<<"CRV1:\t";
        if(i==6)cout<<"GRV1:\t";
        if(i==7)cout<<"TRV1:\t";
        if(i==8)cout<<"GFW1:\t";
        if(i==9)cout<<"GRV1:\t";

        for (int j = 0; j < image_hp1.size(); j++) {
            printf("%3d\t", image_hp1[j][i]);
        }
        cout << endl;
    }
    cout << "-------------" << endl;
    for (int i = 0; i < image_hp2[0].size(); i++) {
        if(i==0)cout<<"AFW2:\t";
        if(i==1)cout<<"CFW2:\t";
        if(i==2)cout<<"GFW2:\t";
        if(i==3)cout<<"TFW2:\t";
        if(i==4)cout<<"ARV2:\t";
        if(i==5)cout<<"CRV2:\t";
        if(i==6)cout<<"GRV2:\t";
        if(i==7)cout<<"TRV2:\t";
        if(i==8)cout<<"GFW2:\t";
        if(i==9)cout<<"GRV2:\t";

        for (int j = 0; j < image_hp2.size(); j++) {
            printf("%3d\t", image_hp2[j][i]);
        }
        cout << endl;
    }

}

void SummaryGenerator::generate_image(long long start_pos, long long end_pos) {
    // at this point labels and positions are generated, now generate the pileup
    for (long long i = start_pos; i <= end_pos; i++) {

        vector<uint8_t> row_h1;
        vector<uint8_t> row_h2;
        uint8_t pixel_value_h1 = 0;
        uint8_t pixel_value_h2 = 0;
        // iterate through the summaries
        for(int j = 0; j <= 9; j++) {
            pixel_value_h1 = (base_summaries_hp1[make_pair(i, j)] / max(1.0, coverage_hp1[i])) * ImageOptions::MAX_COLOR_VALUE;
            pixel_value_h2 = (base_summaries_hp2[make_pair(i, j)] / max(1.0, coverage_hp2[i])) * ImageOptions::MAX_COLOR_VALUE;
            row_h1.push_back(pixel_value_h1);
            row_h2.push_back(pixel_value_h2);
        }

        assert(row_h1.size() == 10);
        assert(row_h2.size() == 10);
        image_hp1.push_back(row_h1);
        image_hp2.push_back(row_h2);

        if (longest_insert_count[i] > 0) {

            for (int ii = 0; ii < longest_insert_count[i]; ii++) {
                vector<uint8_t> ins_row_hp1;
                vector<uint8_t> ins_row_hp2;

                // iterate through the summaries
                for(int j = 0; j <= 9; j++) {
                    pixel_value_h1 = (insert_summaries_hp1[make_pair(make_pair(i, ii), j)] /
                            max(1.0, coverage_hp1[i])) * ImageOptions::MAX_COLOR_VALUE ;
                    pixel_value_h2 = (insert_summaries_hp2[make_pair(make_pair(i, ii), j)] /
                                   max(1.0, coverage_hp2[i])) * ImageOptions::MAX_COLOR_VALUE ;
                    ins_row_hp1.push_back(pixel_value_h1);
                    ins_row_hp2.push_back(pixel_value_h2);
                }
                assert(ins_row_hp1.size() == 10);
                assert(ins_row_hp2.size() == 10);
                image_hp1.push_back(ins_row_hp1);
                image_hp2.push_back(ins_row_hp2);
            }
        }
    }

    assert(image_hp1.size() == genomic_pos.size());
    assert(image_hp1.size() == image_hp2.size());
}



void SummaryGenerator::generate_train_summary(vector <type_read> &reads,
                                              long long start_pos,
                                              long long end_pos,
                                              type_read truth_read_hp1,
                                              type_read truth_read_hp2) {


    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        iterate_over_read(read, start_pos, end_pos);
    }

    // this populates base_labels and insert_labels dictionaries
    generate_labels(truth_read_hp1, start_pos, end_pos + 1, 1);
    generate_labels(truth_read_hp2, start_pos, end_pos + 1, 2);

    // after all the dictionaries are populated, we can simply walk through the region and generate a sequence
    for (long long pos = start_pos; pos <= end_pos; pos++) {

        if(coverage[pos] > 0) {
            labels_hp1.push_back(get_labels(base_labels_hp1[pos]));
            labels_hp2.push_back(get_labels(base_labels_hp2[pos]));
        } else {
            labels_hp1.push_back(get_labels('*'));
            labels_hp2.push_back(get_labels('*'));
        }

        // if the label contains anything but ACTG
        if(!check_base(base_labels_hp1[pos]) || !check_base(base_labels_hp2[pos])) {
//            cerr<<"INFO: INVALID REFERENCE BASE INDEX FOUND: ["<<chromosome_name<<":"<<start_pos<<"-"<<end_pos<<"] " <<
//                pos<<" "<<" "<<base_labels[pos]<<endl;
            bad_label_positions.push_back(labels_hp1.size());
        }

        genomic_pos.push_back(make_pair(pos, 0));
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++) {
                genomic_pos.push_back(make_pair(pos, ii + 1));
                if (insert_labels_hp1[make_pair(pos, ii)]) {
                    labels_hp1.push_back(get_labels(insert_labels_hp1[make_pair(pos, ii)]));

                    // if the label contains anything but ACTG
                    if(!check_base(insert_labels_hp1[make_pair(pos, ii)])) {
//                        cerr<<"INFO: INVALID REFERENCE INSERT BASE INDEX FOUND: "<<chromosome_name<<" "<<
//                            pos<<" "<<insert_labels[make_pair(pos, ii)]<<endl;
                        bad_label_positions.push_back(labels_hp1.size());
                    }
                }
                else labels_hp1.push_back(get_labels('#'));
                if (insert_labels_hp2[make_pair(pos, ii)]) {
                    labels_hp2.push_back(get_labels(insert_labels_hp2[make_pair(pos, ii)]));

                    // if the label contains anything but ACTG
                    if(!check_base(insert_labels_hp2[make_pair(pos, ii)])) {
//                        cerr<<"INFO: INVALID REFERENCE INSERT BASE INDEX FOUND: "<<chromosome_name<<" "<<
//                            pos<<" "<<insert_labels[make_pair(pos, ii)]<<endl;
                        bad_label_positions.push_back(labels_hp2.size());
                    }
                }
                else labels_hp2.push_back(get_labels('#'));
            }
        }
    }
    bad_label_positions.push_back(labels_hp1.size());
    assert(labels_hp1.size() == genomic_pos.size());
    assert(labels_hp1.size() == labels_hp2.size());

    // generate reference sequence
    for (int i = start_pos; i <= end_pos; i++) {
        ref_image.push_back(get_reference_feature_index(reference_sequence[i - start_pos]));

        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++)
                ref_image.push_back(get_reference_feature_index('*'));
        }
    }

    generate_image(start_pos, end_pos);
//     at this point everything should be generated
//    debug_print(start_pos, end_pos);
//    exit(1);
}


void SummaryGenerator::generate_summary(vector <type_read> &reads,
                                        long long start_pos,
                                        long long end_pos) {
    int read_count = 0;
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        read_count += 1;
        iterate_over_read(read, start_pos, end_pos);
    }

    // generate reference sequence
    for (int i = start_pos; i <= end_pos; i++) {
        ref_image.push_back(get_reference_feature_index(reference_sequence[i - start_pos]));

        if (longest_insert_count[i] > 0) {
            for (int ii = 0; ii < longest_insert_count[i]; ii++)
                ref_image.push_back(get_reference_feature_index('*'));
        }
    }

    // after all the dictionaries are populated, we can simply walk through the region and generate a sequence
    for (long long pos = start_pos; pos <= end_pos; pos++) {
        genomic_pos.push_back(make_pair(pos, 0));
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++) {
                genomic_pos.push_back(make_pair(pos, ii + 1));
            }
        }
    }

    generate_image(start_pos, end_pos);
//     at this point everything should be generated
//    debug_print(start_pos, end_pos);
}
