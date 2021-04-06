//
// Created by Kishwar Shafin on 3/21/21.
//

#include "region_summary_hp.h"

#include <utility>

RegionalSummaryGeneratorHP::RegionalSummaryGeneratorHP(long long region_start, long long region_end, string reference_sequence) {
    this->ref_start = region_start;
    this->ref_end = region_end;
    this->reference_sequence = std::move(reference_sequence);
    this->total_observered_insert_bases = 0;
    this->max_observed_insert.resize(region_end-region_start+1, 0);
    this->cumulative_observed_insert.resize(region_end-region_start+1, 0);
}

void RegionalSummaryGeneratorHP::generate_max_insert_observed(const type_read& read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > ref_end) break;
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
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                reference_index = ref_position - ref_start - 1;

                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {
                    max_observed_insert[reference_index] = std::max(max_observed_insert[reference_index], (uint64_t) cigar.length);
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
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


void RegionalSummaryGeneratorHP::generate_max_insert_summary(vector <type_read> &reads) {
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        generate_max_insert_observed(read);
    }

    cumulative_observed_insert[0] = 0;
    total_observered_insert_bases += max_observed_insert[0];

    positions.push_back(ref_start);
    index.push_back(0);

    for(int j=1; j <= max_observed_insert[0]; j++) {
        positions.push_back(ref_start);
        index.push_back(j);
    }

    for(int i=1;i < max_observed_insert.size(); i++) {
        cumulative_observed_insert[i] = cumulative_observed_insert[i-1] + max_observed_insert[i-1];
        total_observered_insert_bases += max_observed_insert[i];
        positions.push_back(ref_start + i);
        index.push_back(0);
        for(int j=1; j <= max_observed_insert[i]; j++) {
            positions.push_back(ref_start + i);
            index.push_back(j);
        }
    }
}


char check_truth_base_hp(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base=='G' || base=='g' ||
       base=='*' || base=='#') return base;
    return '*';
}


void RegionalSummaryGeneratorHP::generate_labels_from_truth_read(type_read read, int hp_tag) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > ref_end) break;
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

                        int base_index = (int)(ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
                        if(hp_tag == 1)
                            labels_hp1[base_index] = check_truth_base_hp(base);
                        else
                            labels_hp2[base_index] = check_truth_base_hp(base);
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

                    for (int i = 0; i < max_observed_insert[(ref_position -1) - ref_start]; i++) {
                        char base = '*';
                        if (i < alt.length()) {
                            base = alt[i];
                        }
                        int base_index = (int)((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start] + (i + 1));

                        if(hp_tag == 1)
                            labels_hp1[base_index] = check_truth_base(base);
                        else
                            labels_hp2[base_index] = check_truth_base(base);
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
                            int base_index = (int)(ref_position - ref_start + i + cumulative_observed_insert[ref_position - ref_start + i]);

                            if(hp_tag == 1)
                                labels_hp1[base_index] = base;
                            else
                                labels_hp2[base_index] = base;
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


void RegionalSummaryGeneratorHP::generate_labels(const type_read& truth_read_hp1,
                                               const type_read& truth_read_hp2) {
    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);
    labels_hp1.resize(region_size + 1, '*');
    labels_hp2.resize(region_size + 1, '*');

    generate_labels_from_truth_read(truth_read_hp1, 1);
    generate_labels_from_truth_read(truth_read_hp2, 2);
}


int RegionalSummaryGeneratorHP::get_feature_index(char base, bool is_reverse) {
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

void RegionalSummaryGeneratorHP::populate_summary_matrix(int **image_matrix_hp1,
                                                         int **image_matrix_hp2,
                                                         int *coverage_hp1,
                                                         int *coverage_hp2,
                                                         int *coverage_vector,
                                                         type_read read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > ref_end) break;
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
                    //read.base_qualities[read_index] base quality
                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];

                        int base_index = (int)(ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
                        int feature_index = get_feature_index(base, read.flags.is_reverse);

                        // update the summary of base
                        if(ref_position >= ref_start && ref_position <= ref_end) {

                            if(read.hp_tag == 0) {
                                image_matrix_hp1[base_index][feature_index] += 1;
                                image_matrix_hp2[base_index][feature_index] += 1;
                                coverage_hp1[ref_position - ref_start] += 1;
                                coverage_hp2[ref_position - ref_start] += 1;
                            }
                            else if(read.hp_tag == 1) {
                                image_matrix_hp1[base_index][feature_index] += 1;
                                coverage_hp1[ref_position - ref_start] += 1;
                            }
                            else if(read.hp_tag == 2) {
                                image_matrix_hp2[base_index][feature_index] += 1;
                                coverage_hp2[ref_position - ref_start] += 1;
                            }

                            coverage_vector[ref_position - ref_start] += 1;
                        }

                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    alt = read.sequence.substr(read_index, cigar.length);
                    for (int i = 0; i < cigar.length; i++) {
                        int base_index = (int)((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start] + (i + 1));
                        int feature_index = get_feature_index(alt[i], read.flags.is_reverse);

                        if(read.hp_tag == 0) {
                            image_matrix_hp1[base_index][feature_index] += 1;
                            image_matrix_hp2[base_index][feature_index] += 1;
                        }
                        else if(read.hp_tag == 1) {
                            image_matrix_hp1[base_index][feature_index] += 1;
                        }
                        else if(read.hp_tag == 2) {
                            image_matrix_hp2[base_index][feature_index] += 1;
                        }
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
                // process delete allele here
                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        // update the summary of base
                        int base_index = (int)(ref_position - ref_start + i + cumulative_observed_insert[ref_position - ref_start + i]);
                        int feature_index = get_feature_index('*', read.flags.is_reverse);

                        if(read.hp_tag == 0) {
                            image_matrix_hp1[base_index][feature_index] += 1;
                            image_matrix_hp2[base_index][feature_index] += 1;
                            coverage_hp1[ref_position - ref_start] += 1;
                            coverage_hp2[ref_position - ref_start] += 1;
                        }
                        else if(read.hp_tag == 1) {
                            image_matrix_hp1[base_index][feature_index] += 1;
                            coverage_hp1[ref_position - ref_start] += 1;
                        }
                        else if(read.hp_tag == 2) {
                            image_matrix_hp2[base_index][feature_index] += 1;
                            coverage_hp2[ref_position - ref_start] += 1;
                        }

                        coverage_vector[ref_position - ref_start] += 1;
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

RegionalImageSummaryHP RegionalSummaryGeneratorHP::generate_summary(vector <type_read> &reads,
                                                                    int chunk_overlap,
                                                                    int smaller_chunk_size,
                                                                    int feature_size,
                                                                    int chunk_id_start,
                                                                    bool train_mode) {
    RegionalImageSummaryHP summary{};

    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

    // Generate a cover vector of chunk size. Chunk size = 10kb defining the region
    int coverage_vector[ref_end - ref_start + 1];
    int coverage_hp1[ref_end - ref_start + 1];
    int coverage_hp2[ref_end - ref_start + 1];

    // generate the image matrix of chunk_size (10kb) * feature_size (10)
    int** image_matrix_hp1 = new int*[region_size + 1];
    int** image_matrix_hp2 = new int*[region_size + 1];

    for (int i = 0; i < region_size + 1; i++)
    {
        image_matrix_hp1[i] = new int[feature_size];
        image_matrix_hp2[i] = new int[feature_size];
        for (int j = 0; j < feature_size; j++) {
            image_matrix_hp1[i][j] = 0;
            image_matrix_hp2[i][j] = 0;
        }
    }

    memset(coverage_vector, 0, sizeof(coverage_vector));
    memset(coverage_hp1, 0, sizeof(coverage_hp1));
    memset(coverage_hp2, 0, sizeof(coverage_hp2));

    // now iterate over all of the reads and populate the image matrix and coverage matrix
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            populate_summary_matrix(image_matrix_hp1, image_matrix_hp2, coverage_hp1, coverage_hp2, coverage_vector, read);
        }
    }
//    cout<<"SUMMARY MATRIX GENERATED"<<endl;

    // once the image matrix is generated, scale the counted values.
    for(int i=0;i<region_size;i++){
        for(int j=0;j<feature_size;j++){
            image_matrix_hp1[i][j] = (int) (((double)image_matrix_hp1[i][j] / max(1.0, (double) coverage_hp1[positions[i]-ref_start])) * ImageOptionsRegion::MAX_COLOR_VALUE);
            image_matrix_hp2[i][j] = (int) (((double)image_matrix_hp2[i][j] / max(1.0, (double) coverage_hp2[positions[i]-ref_start])) * ImageOptionsRegion::MAX_COLOR_VALUE);
        }
    }

//    debug_print_matrix(image_matrix_hp1, image_matrix_hp2, train_mode);

    labels_hp1_report.resize(region_size + 1, 0);
    labels_hp2_report.resize(region_size + 1, 0);
    // check if train mode, if yes, then generate labels
    if(train_mode) {
        for (int i = 0; i < labels_hp1.size(); i++) {
            labels_hp1_report[i] = get_labels(labels_hp1[i]);
        }
        for (int i = 0; i < labels_hp2.size(); i++) {
            labels_hp2_report[i] = get_labels(labels_hp2[i]);
        }
    }

    // Now chunk the image into 1kb smaller matrices. Making them 10 * 1000 * feature_size
    int64_t chunk_start = 0;
    int64_t chunk_end = smaller_chunk_size;
    int total_chunks = 0;

    // first count how many chunks will there be
    while(true) {
        total_chunks += 1;
        if (chunk_end == (int) region_size) {
            break;
        }
        chunk_start = chunk_end - chunk_overlap;
        chunk_end = min((int64_t) region_size, chunk_start + smaller_chunk_size);

    }
    // once we know how many chunks there will be, we can generate a pre-defined vector and do the chunking
    summary.chunked_image_matrix_hp1.resize(total_chunks, vector<vector<uint8_t> >(smaller_chunk_size,vector<uint8_t>(feature_size)));
    summary.chunked_image_matrix_hp2.resize(total_chunks, vector<vector<uint8_t> >(smaller_chunk_size,vector<uint8_t>(feature_size)));
    summary.chunked_labels_hp1.resize(total_chunks, vector<uint8_t>(smaller_chunk_size));
    summary.chunked_labels_hp2.resize(total_chunks, vector<uint8_t>(smaller_chunk_size));
    summary.chunked_positions.resize(total_chunks, vector<int64_t>(smaller_chunk_size));
    summary.chunked_index.resize(total_chunks, vector<int32_t>(smaller_chunk_size));
    summary.chunked_ids.resize(total_chunks);

    chunk_start = 0;
    chunk_end = smaller_chunk_size;
    int current_chunk = 0;

    while(true) {
        summary.chunked_ids[current_chunk] = chunk_id_start + current_chunk;
        for(int64_t i = chunk_start; i < chunk_end; i++) {

            if(i<region_size) {
                summary.chunked_positions[current_chunk][i - chunk_start] = positions[i];
                summary.chunked_index[current_chunk][i - chunk_start] = index[i];
                summary.chunked_labels_hp1[current_chunk][i - chunk_start] = labels_hp1_report[i];
                summary.chunked_labels_hp2[current_chunk][i - chunk_start] = labels_hp1_report[i];
            } else {
                summary.chunked_positions[current_chunk][i - chunk_start] = -1;
                summary.chunked_index[current_chunk][i - chunk_start] = -1;
                summary.chunked_labels_hp1[current_chunk][i - chunk_start] = 0;
                summary.chunked_labels_hp2[current_chunk][i - chunk_start] = 0;
            }
            for(int j=0; j < feature_size; j++) {
                if(i<=region_size) {
                    summary.chunked_image_matrix_hp1[current_chunk][i - chunk_start][j] = image_matrix_hp1[i][j];
                    summary.chunked_image_matrix_hp2[current_chunk][i - chunk_start][j] = image_matrix_hp2[i][j];
                } else {
                    summary.chunked_image_matrix_hp1[current_chunk][i - chunk_start][j] = 0;
                    summary.chunked_image_matrix_hp2[current_chunk][i - chunk_start][j] = 0;
                }
            }
        }
        if (chunk_end >= region_size) {
            break;
        }

        chunk_start = chunk_end - chunk_overlap;
        chunk_end = chunk_start + smaller_chunk_size;

        // in train mode, we don't go outside the window
        if(train_mode) {
            if(chunk_end > region_size) {
                int padding_required = (int)(chunk_end - region_size);
                chunk_start -= padding_required;
                chunk_end = chunk_start + smaller_chunk_size;
                if(chunk_start < 0) break;
            }
        }
        current_chunk += 1;
    }

    // debug_print_matrix(image_matrix_hp1, image_matrix_hp2, train_mode);

    for (int i = 0; i < region_size + 1; i++)
    {
        free(image_matrix_hp1[i]);
        free(image_matrix_hp2[i]);
    }
    free(image_matrix_hp1);
    free(image_matrix_hp2);


    return summary;
}


void RegionalSummaryGeneratorHP::debug_print_matrix(int** image_matrix_hp1, int** image_matrix_hp2, bool train_mode) {
    cout << "------------- IMAGE MATRIX" << endl;

    cout << setprecision(3);
    for (long long i = ref_start; i <= ref_end; i++) {
        if(i==ref_start) cout<<"REF:\t";
        cout << "  " << reference_sequence[i - ref_start] << "\t";
        if (max_observed_insert[i - ref_start] > 0) {
            for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
        }
    }
    cout << endl;

    if(train_mode) {
        for (int i = 0; i <= ref_end - ref_start + total_observered_insert_bases; i++) {
            if (i == 0) cout << "TRH1:\t";
            cout << "  " << labels_hp1[i] << "\t";
        }
        cout << endl;

        for (int i = 0; i <= ref_end - ref_start + total_observered_insert_bases; i++) {
            if (i == 0) cout << "TRH2:\t";
            cout << "  " << labels_hp2[i] << "\t";
        }
        cout << endl;
    }

//    for (int i = 0; i < labels.size(); i++) {
//        if(i==0) cout<<"LBL:\t";
//        printf("%3d\t", labels[i]);
//    }
//    cout << endl;

    cout<<"POS:\t";
    for(int i=0; i < positions.size(); i++ ) {
        printf("%3lld\t", positions[i] % 100);
    }
    cout << endl;
    for (int i = 0; i < 10; i++) {
        if(i==0)cout<<"AFW:\t";
        if(i==1)cout<<"CFW:\t";
        if(i==2)cout<<"GFW:\t";
        if(i==3)cout<<"TFW:\t";
        if(i==4)cout<<"ARV:\t";
        if(i==5)cout<<"CRV:\t";
        if(i==6)cout<<"GRV:\t";
        if(i==7)cout<<"TRV:\t";
        if(i==8)cout<<"*FW:\t";
        if(i==9)cout<<"*RV:\t";

        int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", image_matrix_hp1[j][i]);
        }
        cout << endl;
    }

    for (int i = 0; i < 10; i++) {
        if(i==0)cout<<"AFW:\t";
        if(i==1)cout<<"CFW:\t";
        if(i==2)cout<<"GFW:\t";
        if(i==3)cout<<"TFW:\t";
        if(i==4)cout<<"ARV:\t";
        if(i==5)cout<<"CRV:\t";
        if(i==6)cout<<"GRV:\t";
        if(i==7)cout<<"TRV:\t";
        if(i==8)cout<<"*FW:\t";
        if(i==9)cout<<"*RV:\t";

        int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", image_matrix_hp2[j][i]);
        }
        cout << endl;
    }

}