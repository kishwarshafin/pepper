//
// Created by Kishwar Shafin on 3/21/21.
//

#include "region_summary.h"

#include <utility>

RegionalSummaryGenerator::RegionalSummaryGenerator(string contig, long long region_start, long long region_end, string reference_sequence) {
    this->contig = contig;
    this->ref_start = region_start;
    this->ref_end = region_end;
    this->reference_sequence = std::move(reference_sequence);
    this->total_observered_insert_bases = 0;
    this->max_observed_insert.resize(region_end-region_start+1, 0);
    this->cumulative_observed_insert.resize(region_end-region_start+1, 0);
}

void RegionalSummaryGenerator::generate_max_insert_observed(const type_read& read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;
    if(ImageOptionsRegion::GENERATE_INDELS == true) {
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
                    reference_index = ref_position - ref_start - 1;
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
}


void RegionalSummaryGenerator::generate_max_insert_summary(vector <type_read> &reads) {
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


char check_truth_base(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base=='G' || base=='g' ||
       base=='*' || base=='#') return base;
    return '*';
}

uint8_t get_label_index(char base_h1, char base_h2) {
    base_h1 = toupper(base_h1);
    base_h2 = toupper(base_h2);
    // AA
    if (base_h1 == 'A' && base_h2 == 'A') return 1;
    // AC
    if (base_h1 == 'A' && base_h2 == 'C') return 2;
    if (base_h1 == 'C' && base_h2 == 'A') return 2;
    // AT
    if (base_h1 == 'A' && base_h2 == 'T') return 3;
    if (base_h1 == 'T' && base_h2 == 'A') return 3;
    // AG
    if (base_h1 == 'A' && base_h2 == 'G') return 4;
    if (base_h1 == 'G' && base_h2 == 'A') return 4;
    // A*
    if (base_h1 == 'A' && base_h2 == '*') return 5;
    if (base_h1 == '*' && base_h2 == 'A') return 5;
    // A#
    if (base_h1 == 'A' && base_h2 == '#') return 6;
    if (base_h1 == '#' && base_h2 == 'A') return 6;
    // CC
    if (base_h1 == 'C' && base_h2 == 'C') return 7;
    // CT
    if (base_h1 == 'C' && base_h2 == 'T') return 8;
    if (base_h1 == 'T' && base_h2 == 'C') return 8;
    // CG
    if (base_h1 == 'C' && base_h2 == 'G') return 9;
    if (base_h1 == 'G' && base_h2 == 'C') return 9;
    // C*
    if (base_h1 == 'C' && base_h2 == '*') return 10;
    if (base_h1 == '*' && base_h2 == 'C') return 10;
    // C#
    if (base_h1 == 'C' && base_h2 == '#') return 11;
    if (base_h1 == '#' && base_h2 == 'C') return 11;
    // TT
    if (base_h1 == 'T' && base_h2 == 'T') return 12;
    // TG
    if (base_h1 == 'T' && base_h2 == 'G') return 13;
    if (base_h1 == 'G' && base_h2 == 'T') return 13;
    // T*
    if (base_h1 == 'T' && base_h2 == '*') return 14;
    if (base_h1 == '*' && base_h2 == 'T') return 14;
    // T#
    if (base_h1 == 'T' && base_h2 == '#') return 15;
    if (base_h1 == '#' && base_h2 == 'T') return 15;
    // GG
    if (base_h1 == 'G' && base_h2 == 'G') return 16;
    // G*
    if (base_h1 == 'G' && base_h2 == '*') return 17;
    if (base_h1 == '*' && base_h2 == 'G') return 17;
    // G#
    if (base_h1 == 'G' && base_h2 == '#') return 18;
    if (base_h1 == '#' && base_h2 == 'G') return 18;
    // **
    if (base_h1 == '*' && base_h2 == '*') return 19;
    // *#
    if (base_h1 == '*' && base_h2 == '#') return 20;
    if (base_h1 == '#' && base_h2 == '*') return 20;
    // ##
    if (base_h1 == '#' && base_h2 == '#') return 0;
    return 0;
}


uint8_t get_variant_type_label_index(int type_h1, int type_h2) {

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::HOM_REF) return 0;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::SNP) return 1;
    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::HOM_REF) return 1;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::INSERT) return 2;
    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::HOM_REF) return 2;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::DELETE) return 3;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::HOM_REF) return 3;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::SNP) return 4;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::INSERT) return 5;
    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::SNP) return 5;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::DELETE) return 6;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::SNP) return 6;

    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::INSERT) return 7;

    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::DELETE) return 8;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::INSERT) return 8;

    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::DELETE) return 9;

    cout<<"ERROR: VARIANT LABEL NOT DEFINED: "<<type_h1<<" "<<type_h2<<endl;
    exit(1);
}


int RegionalSummaryGenerator::get_reference_feature_index(char base) {
    base = toupper(base);
    if (base == 'A') return ImageOptionsRegion::REFERENCE_INDEX_START;
    if (base == 'C') return ImageOptionsRegion::REFERENCE_INDEX_START + 1;
    if (base == 'G') return ImageOptionsRegion::REFERENCE_INDEX_START + 2;
    if (base == 'T') return ImageOptionsRegion::REFERENCE_INDEX_START + 3;
    return ImageOptionsRegion::REFERENCE_INDEX_START + 4;
}

void RegionalSummaryGenerator::encode_reference_bases(vector< vector<int> >& image_matrix) {
    for (long long ref_position = ref_start; ref_position <= ref_end; ref_position++) {
        // encode the C base
        int base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
        int feature_index = get_reference_feature_index(reference_sequence[ref_position - ref_start]);
        image_matrix[base_index][feature_index] = 255;

        for(int i = 1; i <= max_observed_insert[ref_position - ref_start]; i++) {
            base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]) + i;
            feature_index = get_reference_feature_index('*');
            image_matrix[base_index][feature_index] = 255;
        }
    }
}

int RegionalSummaryGenerator::get_feature_index(char base, bool is_reverse) {
    base = toupper(base);
    if (is_reverse) {
        if (base == 'A') return ImageOptionsRegion::BASE_INDEX_START;
        if (base == 'C') return ImageOptionsRegion::BASE_INDEX_START + 1;
        if (base == 'G') return ImageOptionsRegion::BASE_INDEX_START + 2;
        if (base == 'T') return ImageOptionsRegion::BASE_INDEX_START + 3;
        if (base == 'I') return ImageOptionsRegion::BASE_INDEX_START + 5;
        if (base == 'D') return ImageOptionsRegion::BASE_INDEX_START + 6;
        return ImageOptionsRegion::BASE_INDEX_START + 4;
    } else {
        // tagged and forward
        if (base == 'A') return ImageOptionsRegion::BASE_INDEX_START + 7;
        if (base == 'C') return ImageOptionsRegion::BASE_INDEX_START + 8;
        if (base == 'G') return ImageOptionsRegion::BASE_INDEX_START + 9;
        if (base == 'T') return ImageOptionsRegion::BASE_INDEX_START + 10;
        if (base == 'I') return ImageOptionsRegion::BASE_INDEX_START + 12;
        if (base == 'D') return ImageOptionsRegion::BASE_INDEX_START + 13;
        return ImageOptionsRegion::BASE_INDEX_START + 11;
    }
}



void RegionalSummaryGenerator::generate_labels(const vector<type_truth_record>& hap1_records, const vector<type_truth_record>& hap2_records) {
    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);
    labels_hp1.resize(region_size + 1, '*');
    labels_hp2.resize(region_size + 1, '*');
    variant_type_labels_hp1.resize(region_size + 1, VariantTypes::HOM_REF);
    variant_type_labels_hp2.resize(region_size + 1, VariantTypes::HOM_REF);

    for(long long pos = ref_start; pos <= ref_end; pos++) {
        int base_index = (int)(pos - ref_start + cumulative_observed_insert[pos - ref_start]);
        labels_hp1[base_index] = reference_sequence[base_index];
        labels_hp2[base_index] = reference_sequence[base_index];
    }

    for(const auto& truth_record : hap1_records) {
        if(truth_record.ref.length() > truth_record.alt.length()) {
            //it's a delete
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::DELETE;
                labels_hp1[base_index] = '#';
            }
            // dont expand to the rest of the bases

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    if (pos - truth_record.pos_start < truth_record.alt.length()) {
//                        labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                    } else {
//                        labels_hp1[base_index] = '*';
//                    }
//                }
//            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::INSERT;
                labels_hp1[base_index] = '*';
            }

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                }
//            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::SNP;
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
                    labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
                }
            }
        }
    }

    for(const auto& truth_record : hap2_records) {
        if(truth_record.ref.length() > truth_record.alt.length()) {
            //it's a delete
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::DELETE;
                labels_hp2[base_index] = '#';
            }
            // dont expand to the rest of the bases

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    if (pos - truth_record.pos_start < truth_record.alt.length()) {
//                        labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                    } else {
//                        labels_hp2[base_index] = '*';
//                    }
//                }
//            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::INSERT;
                labels_hp2[base_index] = '*';
            }

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                }
//            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::SNP;
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
                    labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
                }
            }
        }
    }
}


void RegionalSummaryGenerator::populate_summary_matrix(vector< vector<int> >& image_matrix,
                                                       int *coverage_vector,
                                                       int *snp_count,
                                                       int *insert_count,
                                                       int *delete_count,
                                                       type_read read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    double base_quality = 1.0 - pow(10.0 , (-1.0 * read.base_qualities[read_index])/10);
    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];
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
//                    double base_quality = 1.0 - pow(10.0 , (-1.0 * read.base_qualities[read_index])/10);

                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        char ref_base = reference_sequence[ref_position - ref_start];

                        int base_index = (int)(ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
                        int feature_index = get_feature_index(base, read.flags.is_reverse);

                        // update the summary of base
                        if (ref_position >= ref_start && ref_position <= ref_end) {
                            image_matrix[base_index][feature_index] += 1;
                            coverage_vector[ref_position - ref_start] += 1;

                            if(ref_base != base) {
                                snp_count[ref_position - ref_start] += 1;
                            }
                        }



                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                if (ref_position - 1 >= ref_start && ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start]);
                    int insert_count_index =  get_feature_index('I', read.flags.is_reverse);
                    image_matrix[base_index][insert_count_index] += 1;
                    insert_count[ref_position - 1 - ref_start] += 1;

//                    if(ImageOptionsRegion::GENERATE_INDELS == true) {
//                        alt = read.sequence.substr(read_index, cigar.length);
//                        for (int i = 0; i < cigar.length; i++) {
//                            base_index = (int) ((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start] + (i + 1));
//                            int feature_index = get_feature_index(alt[i], read.flags.is_reverse);
//                            image_matrix[base_index][feature_index] += 1;
//                        }
//                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
                // process delete allele here

                if (ref_position -1 >= ref_start && ref_position - 1 <= ref_end) {
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)(ref_position - 1 - ref_start + cumulative_observed_insert[ref_position - 1 - ref_start]);
                    int delete_count_index =  get_feature_index('D', read.flags.is_reverse);
                    image_matrix[base_index][delete_count_index] += 1.0;
                    delete_count[ref_position - 1 - ref_start] += 1;

//                    int feature_index = get_feature_index('*', read.flags.is_reverse);
//                    image_matrix[base_index][feature_index] += 1.0;
                }
                // dont' expand to the full delete length, rather just mount everything to the anchor

                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        // update the summary of base
                        int base_index = (int) (ref_position - ref_start + i + cumulative_observed_insert[ref_position - ref_start + i]);
                        char ref_base = reference_sequence[ref_position - ref_start + i];
                        int feature_index = get_feature_index('*', read.flags.is_reverse);

                        image_matrix[base_index][feature_index] += 1.0;
//                        coverage_vector[ref_position - ref_start + i] += 1.0;
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

vector<CandidateImageSummary> RegionalSummaryGenerator::generate_summary(vector <type_read> &reads,
                                                                         double snp_freq_threshold,
                                                                         double insert_freq_threshold,
                                                                         double delete_freq_threshold,
                                                                         double min_coverage_threshold,
                                                                         long long candidate_region_start,
                                                                         long long candidate_region_end,
                                                                         int candidate_window_size,
                                                                         int feature_size,
                                                                         bool train_mode) {


    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

    // Generate a cover vector of chunk size. Chunk size = 10kb defining the region
    int coverage_vector[ref_end - ref_start + 1];
    int snp_count[ref_end - ref_start + 1];
    int insert_count[ref_end - ref_start + 1];
    int delete_count[ref_end - ref_start + 1];

    // generate the image matrix of chunk_size (10kb) * feature_size (10)
    vector< vector<int> > image_matrix;

    image_matrix.resize(region_size + 1, vector<int>(feature_size));

    for (int i = 0; i < region_size + 1; i++) {
        for (int j = 0; j < feature_size; j++)
            image_matrix[i][j] = 0;
    }

    memset(coverage_vector, 0, sizeof(coverage_vector));
    memset(snp_count, 0, sizeof(snp_count));
    memset(insert_count, 0, sizeof(insert_count));
    memset(delete_count, 0, sizeof(delete_count));

    encode_reference_bases(image_matrix);

    // now iterate over all of the reads and populate the image matrix and coverage matrix
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            populate_summary_matrix(image_matrix, coverage_vector, snp_count, insert_count, delete_count, read);
        }
    }

    vector<long long> filtered_candidate_positions;
    // once the image matrix is generated, scale the counted values.
    for(int i=0;i<region_size;i++){
        double snp_fraction = snp_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);
        double insert_fraction = insert_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);
        double delete_fraction = delete_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);

        if(snp_fraction >= snp_freq_threshold || insert_fraction >= insert_freq_threshold || delete_fraction >= delete_freq_threshold) {
            if(positions[i] >= candidate_region_start && positions[i] <= candidate_region_end && coverage_vector[positions[i]-ref_start] >= min_coverage_threshold) {
                filtered_candidate_positions.push_back(positions[i]);
            }
        }
        // normalize things
//        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++){
//            image_matrix[i][j] = (int) (((double)image_matrix[i][j] / max(1.0, (double) coverage_vector[positions[i]-ref_start])) * ImageOptionsRegion::MAX_COLOR_VALUE);
//        }

        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++){
            image_matrix[i][j] = (int) min(image_matrix[i][j], ImageOptionsRegion::MAX_COLOR_VALUE);
        }
//        int fwd_feature_index = get_feature_index(ref_at_labels[i], false);
//        int rev_feature_index = get_feature_index(ref_at_labels[i], true);
//
//        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++) {
//            if(image_matrix[i][j] == 0) continue;
//
//            if (j == fwd_feature_index || j == rev_feature_index){
//                image_matrix[i][j] = min(ImageOptionsRegion::MAX_REF_COLOR_VALUE, image_matrix[i][j]);
//            } else {
//                image_matrix[i][j] = min(ImageOptionsRegion::MAX_COLOR_VALUE, ImageOptionsRegion::MISMATCH_COLOR_START + image_matrix[i][j]);
//            }
//        }
    }

    labels.resize(region_size + 1, 0);
    labels_variant_type.resize(region_size + 1, 0);
    // check if train mode, if yes, then generate labels
    if(train_mode) {
        for (int i = 0; i < labels_hp1.size(); i++) {
            labels[i] = get_label_index(labels_hp1[i], labels_hp2[i]);
            labels_variant_type[i] = get_variant_type_label_index(variant_type_labels_hp1[i], variant_type_labels_hp2[i]);
        }
    }

    vector<CandidateImageSummary> all_candidate_images;
    // at this point all of the images are generated. So we can create the images for each candidate position.
    for(long long candidate_position : filtered_candidate_positions) {
//        cout<<candidate_position<<endl;
        int base_index = (int)(candidate_position - ref_start + cumulative_observed_insert[candidate_position - ref_start]);

        CandidateImageSummary candidate_summary;
        candidate_summary.contig = contig;
        candidate_summary.position = candidate_position;
        if(train_mode) {
            candidate_summary.base_label = labels[base_index];
            candidate_summary.type_label = labels_variant_type[base_index];
        } else {
            candidate_summary.base_label = 0;
            candidate_summary.type_label = 0;
        }
        candidate_summary.region_start = ref_start;
        candidate_summary.region_stop = ref_end;

        int base_left = base_index - candidate_window_size / 2;
        int base_right = base_index + candidate_window_size / 2;

//        cout<<base_left<<" "<<base_right<<endl;
        candidate_summary.image_matrix.resize(candidate_window_size + 1, vector<uint8_t>(feature_size));
        // now copy the entire feature matrix

        for(int i=base_left; i<=base_right;i++) {
            for (int j = 0; j < feature_size; j++) {
                candidate_summary.image_matrix[i-base_left][j] = image_matrix[i][j];
            }
        }
        all_candidate_images.push_back(candidate_summary);
//        debug_candidate_summary(candidate_summary, candidate_window_size, train_mode);
    }



//    debug_print_matrix(image_matrix, train_mode);
//    for (int i = 0; i < region_size + 1; i++)
//        free(image_matrix[i]);
//    free(image_matrix);
//    exit(0);


    return all_candidate_images;
}


void RegionalSummaryGenerator::debug_print_matrix(vector<vector<int> > image_matrix, bool train_mode) {
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

    for (int i = 0; i < labels.size(); i++) {
        if(i==0) cout<<"LBL:\t";
        printf("%3d\t", labels[i]);
    }
    cout << endl;

    for (int i = 0; i < labels_variant_type.size(); i++) {
        if(i==0) cout<<"TYP:\t";
        printf("%3d\t", labels_variant_type[i]);
    }
    cout << endl;

    cout<<"POS:\t";
    for(int i=0; i < positions.size(); i++ ) {
        printf("%3lld\t", positions[i] % 100);
    }
    cout << endl;
    int image_size = ImageOptionsRegion::REFERENCE_INDEX_SIZE + ImageOptionsRegion::BASE_INDEX_SIZE;
    for (int i = 0; i < image_size; i++) {
        cout<< ImageOptionsRegion::column_values[i] <<"\t";
        int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", image_matrix[j][i]);
        }
        cout << endl;
    }

}


void RegionalSummaryGenerator::debug_candidate_summary(CandidateImageSummary candidate, int small_chunk_size, bool train_mode) {
    vector<string> decoded_base_lables {"##", "AA", "AC", "AT", "AG", "A*", "A#", "CC", "CT", "CG", "C*", "C#", "TT", "TG", "T*", "T#", "GG", "G*", "G#", "**", "*#"};
    vector<string> decoded_type_lables {"RR", "RS", "RI", "RD", "SS", "SI", "SD", "II", "ID", "DD" };
    cout << "------------- CANDIDATE PILEUP" << endl;
    cout<<"Contig: "<<candidate.contig<<endl;
    cout<<"Position: "<<candidate.position<<endl;
    cout<<"Type label: "<<(int)candidate.type_label<<" "<<decoded_type_lables[candidate.type_label]<<endl;
    cout<<"Base label: :"<<(int)candidate.base_label<<" "<<decoded_base_lables[candidate.base_label]<<endl;

    long long candidate_ref_start = candidate.position - small_chunk_size / 2;
    long long candidate_ref_end = candidate.position + small_chunk_size / 2;
    cout << setprecision(3);
    for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
        if (i == candidate_ref_start) cout << "POS:\t";
        printf("%3lld\t", (i - candidate_ref_start) % 100);
    }
    cout << endl;

    for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
        if(i==candidate_ref_start) cout<<"REF:\t";
        cout << "  " << reference_sequence[i - ref_start] << "\t";
        if (max_observed_insert[i - ref_start] > 0) {
            for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
        }
    }
    cout << endl;


    if(train_mode) {

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRH1:\t";
            cout << "  " << labels_hp1[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRH2:\t";
            cout << "  " << labels_hp2[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRT1:\t";
            cout << "  " << variant_type_labels_hp1[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRT2:\t";
            cout << "  " << variant_type_labels_hp2[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;
    }

    int image_size = candidate.image_matrix[0].size();
    for (int i = 0; i < image_size; i++) {
        cout<< ImageOptionsRegion::column_values[i] <<"\t";
        int region_size = candidate.image_matrix.size();

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", candidate.image_matrix[j][i]);
        }
        cout << endl;
    }

}