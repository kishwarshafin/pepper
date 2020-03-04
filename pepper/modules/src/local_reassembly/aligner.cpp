//
// Created by Kishwar Shafin on 10/21/18.
//


// REUSE DECLARATION: THIS SCRIPT HEAVILY REUSES THE METHODS USED IN DEEPVARIANT https://github.com/google/deepvariant/
// LICENSE INCLUDED IN: third_party/deepVariant.LICENSE
#include "../../headers/realignment/aligner.h"
#include "ssw_cpp.cpp"
#include "ssw.c"

LibSSWPairwiseAligner::LibSSWPairwiseAligner() {
    Filter filter;
    ssw_aligner = new Aligner(Aligner_options::match,
                              Aligner_options::mismatch,
                              Aligner_options::gap_open_penalty,
                              Aligner_options::gap_extend_penalty);
}

Alignment LibSSWPairwiseAligner::align(string query){
    Alignment alignment;
    bool status = ssw_aligner->Align_cpp(query.c_str(), filter, &alignment, 0);
    return alignment;
}

void LibSSWPairwiseAligner::set_reference(string reference) {
    ssw_aligner->SetReferenceSequence(reference.c_str(), reference.length());
}

ReadAligner::ReadAligner(int ref_start, int ref_end, string ref_seq) {
    reference_sequence = ref_seq;
    region_start = ref_start;
    region_end = ref_end;
}

void ReadAligner::CreateKmerMap() {
    int read_hash_id = 0;
    for (const auto& read : read_sequences) {
        if (read.length() <= kmer_size) {
            continue;
        }
        auto last_pos = read.length() - kmer_size;
        for (int pos = 0; pos <= last_pos; pos++) {
            KmerIndexMap[read.substr(pos, kmer_size)].push_back(KmerMap(read_hash_id, pos));
        }
        read_hash_id += 1;
    }
}

int ReadAligner::AlignStrings(string s1, string s2, int max_mismatches, int* num_of_mismatches) const {
    int num_of_matches = 0;
    *num_of_mismatches = 0;

    for (int i = 0; i < s1.size(); i++) {
        const auto& c1 = s1[i];
        const auto& c2 = s2[i];
        if (c1 != c2 && (c1 != 'N' && c2 != 'N')) {
            if (c1 != c2) {
                (*num_of_mismatches)++;
            }
            if (*num_of_mismatches == max_mismatches) {
                return 0;
            }
        } else {
            num_of_matches++;
        }
    }
    return num_of_matches * Aligner_options::match - *num_of_mismatches * Aligner_options::mismatch;
}

void ReadAligner::AlignReadsToHaplotype( const string& haplotype,
                                         int* haplotype_score,
                                         vector<ReadAlignment>* haplotype_read_alignment_scores) {
    // In the loop we try to align reads for each position in haplotype up to
    // lastPos.
    const auto& lastPos = haplotype.length() - kmer_size;
    for (int i = 0; i <= lastPos; i++) {
        // get all reads that are aligned against i-th position
        auto index_it = KmerIndexMap.find(haplotype.substr(i, kmer_size));
        if (index_it == KmerIndexMap.end()) {
            continue;
        }
        // Iterate through all the reads that are found in the index for the current
        // kmer.
        for (const auto& it : index_it->second) {
            uint64_t read_id_index = static_cast<uint64_t>(it.read_id);
            if(read_id_index >= read_sequences.size()) continue;

            size_t target_start_pos = std::max(
                    static_cast<int64_t>(0),
                    static_cast<int64_t>(i) - static_cast<int64_t>(it.read_pos));
            size_t cur_read_size = read_sequences[read_id_index].size();
            size_t span = cur_read_size;
            if (target_start_pos + cur_read_size > haplotype.length()) {
                continue;
            }
            auto& read_alignment =
                    (*haplotype_read_alignment_scores)[read_id_index];

            // This read is already aligned, skip it.
            if (read_alignment.position != 0
                && read_alignment.position == target_start_pos) {
                continue;
            }
            int num_of_mismatches = 0;
            int new_read_alignment_score = AlignStrings(haplotype.substr(target_start_pos, span),
                                                        read_sequences[read_id_index],
                                                        Aligner_options::max_number_of_mismatches + 1,
                                                        &num_of_mismatches);

            // For reads that cannot be aligned with fast alignment we want to avoid
            // tying them over and over. In order to do that we set position for the
            // read even if the read could not be aligned. This way we know that the
            // read was already tried at this position and we can skip it. Doing so
            // reduces a number of checks per read 10 times.
            // If score is not zero we cannot change position without fist checking
            // the score.
            if (read_alignment.score == 0) {
                read_alignment.position = target_start_pos;
            }

            if (num_of_mismatches <= Aligner_options::max_number_of_mismatches) {
                int oldScore = read_alignment.score;
                if (oldScore < new_read_alignment_score) {
                    read_alignment.score = new_read_alignment_score;
                    *haplotype_score -= oldScore;
                    *haplotype_score += read_alignment.score;
                    read_alignment.position = target_start_pos;
                    read_alignment.cigar = std::to_string(cur_read_size) + "=";
                }
            }
        }  // for (matching reads)
    }    // for (all k-mer positions)
}


void ReadAligner::AlignReadsToHaplotypes() {
    vector<ReadAlignment> read_alignment_scores(read_sequences.size());
    for (int i = 0; i < haplotypes.size(); i++) {
        const auto& haplotype = haplotypes[i];
        int haplotype_score = 0;
        for (auto& readAlignment : read_alignment_scores) {
            readAlignment.reset();
        }
        AlignReadsToHaplotype(haplotype,
                              &haplotype_score,
                              &read_alignment_scores);
        read_to_haplotype_alignments.push_back(HaplotypeReadsAlignment(i, haplotype_score, read_alignment_scores));
    }
}


int CigarOperationFromChar(char op) {
    switch (op) {
        case '=':
        case 'X':
            return CIGAR_OPERATIONS::MATCH;
        case 'S':
            return CIGAR_OPERATIONS::SOFT_CLIP;
        case 'D':
            return CIGAR_OPERATIONS::DEL;
        case 'I':
            return CIGAR_OPERATIONS::IN;
        default:
            return CIGAR_OPERATIONS::UNSPECIFIED;
    }
}


list<CigarOp> ReadAligner::CigarStringToVector(const string& cigar) {
    list<CigarOp> cigarOps;
    istringstream parser(cigar);
    char cigar_opch;
    int cigar_len;

    while(parser >> cigar_len >> cigar_opch) {
        int cigar_op = CigarOperationFromChar(cigar_opch);
        cigarOps.push_back(CigarOp(cigar_op, cigar_len));
    }
    return cigarOps;
}


void ReadAligner::AlignHaplotypesToReference() {
    this->SSWAligner.set_reference(reference_sequence);

    // Initialize read_to_haplotype_alignments if it is not initialized yet.
    if (read_to_haplotype_alignments.empty()) {
        for (int i = 0; i < haplotypes.size(); i++) {
            read_to_haplotype_alignments.push_back(HaplotypeReadsAlignment(i, -1,
                                                                           vector<ReadAlignment>(read_sequences.size())));
        }
    }

    for (auto& haplotype_alignment : read_to_haplotype_alignments) {
        Alignment alignment = SSWAligner.align(haplotypes[haplotype_alignment.haplotype_index]);
        auto hap_len = haplotypes[haplotype_alignment.haplotype_index].size();
        if (alignment.sw_score > 0) {
            string all_equal_cigar = to_string(hap_len) + "=";
            haplotype_alignment.is_reference = alignment.cigar_string == all_equal_cigar;
            haplotype_alignment.cigar = alignment.cigar_string;
            haplotype_alignment.cigar_ops = CigarStringToVector(haplotype_alignment.cigar);
            haplotype_alignment.ref_pos = alignment.ref_begin;
        }
    }
}


void ReadAligner::CalculatePositionMaps() {
    for (auto& hyplotype_alignment : read_to_haplotype_alignments) {
        hyplotype_alignment.hap_to_ref_positions_map.resize(haplotypes[hyplotype_alignment.haplotype_index].size());
        int cur_shift = 0;
        int haplotype_pos = 0;
        int last_pos = 0;


        istringstream parser(hyplotype_alignment.cigar);
        char op;
        int operation_len;

        while (parser >> operation_len >> op) {
            switch (op) {
                case '=':
                case 'X':
                    last_pos = haplotype_pos + operation_len;
                    while (haplotype_pos != last_pos) {
                        hyplotype_alignment.hap_to_ref_positions_map[haplotype_pos] = cur_shift;
                        haplotype_pos++;
                    }
                    break;
                case 'S':
                    last_pos = haplotype_pos + operation_len;
                    cur_shift -= operation_len;
                    while (haplotype_pos != last_pos) {
                        hyplotype_alignment.hap_to_ref_positions_map[haplotype_pos] = cur_shift;
                        haplotype_pos++;
                    }
                    break;
                case 'D':
                    cur_shift += operation_len;
                    break;
                case 'I':
                    last_pos = haplotype_pos + operation_len;
                    while (haplotype_pos != last_pos) {
                        hyplotype_alignment.hap_to_ref_positions_map[haplotype_pos] = cur_shift;
                        cur_shift--;
                        haplotype_pos++;
                    }
                    break;
            }
        }
    }
}


void ReadAligner::SswAlignReadsToHaplotypes(int score_threshold) {
    // For each read
    for (int i = 0; i < read_sequences.size(); i++) {
        bool has_at_least_one_alignment = false;
        // Check if this read is aligned to at least one haplotype
        for (const auto& hap_alignment : read_to_haplotype_alignments) {
            if (hap_alignment.read_alignment_scores[i].score > 0) {
                has_at_least_one_alignment = true;
                break;
            }
        }
//        // If this read is not aligned to any of the haplotypes we try SSW.
        if (!has_at_least_one_alignment) {
            for (auto& hap_alignment : read_to_haplotype_alignments) {
                SSWAligner.set_reference(haplotypes[hap_alignment.haplotype_index]);
                Alignment alignment = SSWAligner.align(read_sequences[i]);
                if (alignment.sw_score > 0) {
                    if (alignment.sw_score >= score_threshold) {
                        if (hap_alignment.read_alignment_scores[i].score < alignment.sw_score) {
                            hap_alignment.read_alignment_scores[i].score = alignment.sw_score;
                            hap_alignment.read_alignment_scores[i].cigar = alignment.cigar_string;
                            hap_alignment.read_alignment_scores[i].position = alignment.ref_begin;
                        }
                    }
                }
            }
        }
    }  // for all reads
}


bool ReadAligner::GetBestReadAlignment(size_t readId, int* best_hap_index) const {
    int best_score = 0;
    bool best_haplotype_found = false;
    for (int hap_index = 0; hap_index < haplotypes.size(); hap_index++) {
        if (read_to_haplotype_alignments[hap_index].read_alignment_scores[readId].score > best_score
            // If compared scores are equal preference is given to a read alignment
            // to a non-reference haplotype.
            || (best_score > 0 &&
                read_to_haplotype_alignments[hap_index].read_alignment_scores[readId].score == best_score &&
                !read_to_haplotype_alignments[hap_index].is_reference)) {
            best_score = read_to_haplotype_alignments[hap_index].read_alignment_scores[readId].score;
            *best_hap_index = hap_index;
            best_haplotype_found = true;
        }
    }
    return best_haplotype_found;
}


void ReadAligner::RealignReadsToReference(const vector<type_read>& reads, vector<type_read>& realigned_reads) {
    // Loop through all reads
    for (size_t read_index = 0; read_index < reads.size(); read_index++) {
        const type_read& read = reads[read_index];
        type_read realigned_read;
        realigned_read = read;
        int best_hap_index = -1;
        // See if we have a better alignment
        if (GetBestReadAlignment(read_index, &best_hap_index)) {
            const HaplotypeReadsAlignment& bestHaplotypeAlignments = read_to_haplotype_alignments[best_hap_index];
            unique_ptr<LinearAlignment> new_alignment(new LinearAlignment());
            new_alignment->copy_from_read(read);
            new_alignment->clear_cigar();
            // Calculate new alignment position.
            long long new_position = read.pos;
            auto read_to_hap_pos = bestHaplotypeAlignments.read_alignment_scores[read_index].position;
            // We only change position of original read alignment and don't change
            // chromosome, it shouldn't change anyway!
            new_position = region_start + bestHaplotypeAlignments.ref_pos + read_to_hap_pos
                    + bestHaplotypeAlignments.hap_to_ref_positions_map[read_to_hap_pos];
            new_alignment->set_position(new_position);
            list<CigarOp> readToRefCigarOps;
            // Calculate new cigar by merging read to haplotype and haplotype to ref
            // alignments.
            CalculateReadToRefAlignment(read_index, bestHaplotypeAlignments.read_alignment_scores[read_index],
                                        bestHaplotypeAlignments.cigar_ops, &readToRefCigarOps);

            for (auto& op : readToRefCigarOps) {
                new_alignment->add_cigar(op);
            }
            if (!readToRefCigarOps.empty()) {
                realigned_read.set_cigar_tuples(new_alignment->cigar_tuples);
                realigned_read.set_position(new_alignment->pos);
//                realigned_read.set_alignment_from_linear_alignment(new_alignment);
            }
            realigned_reads.push_back(realigned_read);
        } else {  // keep original alignment
            realigned_reads.push_back(realigned_read);
        }
    }  // for
}


list<CigarOp> ReadAligner::LeftTrimHaplotypeToRefAlignment(const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
                                              int read_to_haplotype_pos) {
    int cur_pos = 0;
    list<CigarOp> haplotype_to_ref_cigar_ops(haplotype_to_ref_cigar_ops_input);

    while (cur_pos != read_to_haplotype_pos) {
        CigarOp cur_hap_op = haplotype_to_ref_cigar_ops.front();
        haplotype_to_ref_cigar_ops.pop_front();
        if (cur_hap_op.operation == CIGAR_OPERATIONS::MATCH ||
            cur_hap_op.operation == CIGAR_OPERATIONS::HARD_CLIP ||
            cur_hap_op.operation == CIGAR_OPERATIONS::SOFT_CLIP ||
            cur_hap_op.operation == CIGAR_OPERATIONS::IN) {
                if (cur_hap_op.length + cur_pos > read_to_haplotype_pos) {
                    haplotype_to_ref_cigar_ops.push_front(CigarOp(cur_hap_op.operation,
                                                                  cur_hap_op.length - (read_to_haplotype_pos - cur_pos)));
            }
            cur_pos = min(cur_hap_op.length + cur_pos, read_to_haplotype_pos);
        }
    }

    // If after trimming the first operation is DEL we need to remove it,
    // because read alignment cannot start with DEL.
    if (haplotype_to_ref_cigar_ops.front().operation == CIGAR_OPERATIONS::DEL) {
        haplotype_to_ref_cigar_ops.pop_front();
    }

    return haplotype_to_ref_cigar_ops;
}


int ReadAligner::AlignedLength(const list<CigarOp>& cigar) {
    int len = 0;
    for (auto& op : cigar) {
        if (op.operation != CIGAR_OPERATIONS::DEL) {
            len += op.length;
        }
    }
    return len;
}

void ReadAligner::MergeCigarOp(const CigarOp& op, int read_len, std::list<CigarOp>* cigar) {
    int last_cigar_op = cigar->empty() ? CIGAR_OPERATIONS::UNSPECIFIED : cigar->back().operation;
    int aligned_length_before_merge = AlignedLength(*cigar);
    int new_op_length = 0;
    if (op.operation != CIGAR_OPERATIONS::DEL) {
        new_op_length = min(op.length, read_len - aligned_length_before_merge);
    } else {
        new_op_length = op.length;
    }

    // Nothing is merged if we already aligned all positions of the read.
    if (new_op_length <= 0 || aligned_length_before_merge == read_len) {
        return;
    }

    // If the op we are adding is the same as the last one on the list, directly
    // add the length to it.
    if (op.operation == last_cigar_op) {
        cigar->back().length += new_op_length;

        // If the op we are adding is not the same as the last, set the proper
        // length and add it to the list.
    }  else {
        cigar->push_back(CigarOp(op.operation, new_op_length));
    }
}

void ReadAligner::CalculateReadToRefAlignment(size_t read_index, const ReadAlignment& read_to_haplotype_alignment,
                                              const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
                                              list<CigarOp>* read_to_ref_cigar_ops) {
    int read_len = read_sequences[read_index].length();
    int read_to_haplotype_pos = read_to_haplotype_alignment.position;
    std::list<CigarOp> read_to_haplotype_cigar_ops =
            CigarStringToVector(read_to_haplotype_alignment.cigar);

    // Left trim haplotype to reference cigar to match read to haplotype
    // alignment position.
    list<CigarOp> haplotype_to_ref_cigar_ops = LeftTrimHaplotypeToRefAlignment(haplotype_to_ref_cigar_ops_input,
                                                                               read_to_haplotype_pos);

    // Skip heading soft clips.
    if (!read_to_haplotype_cigar_ops.empty() &&
        read_to_haplotype_cigar_ops.front().operation == CIGAR_OPERATIONS::SOFT_CLIP) {

        MergeCigarOp(CigarOp(CIGAR_OPERATIONS::SOFT_CLIP, read_to_haplotype_cigar_ops.front().length),
                     read_len, read_to_ref_cigar_ops);
        read_to_haplotype_cigar_ops.pop_front();
    }


    // Build read to reference cigar by iterating CigarOp overlaps.
    while ((!read_to_haplotype_cigar_ops.empty() || !haplotype_to_ref_cigar_ops.empty()) &&
           AlignedLength(*read_to_ref_cigar_ops) < read_len) {
        // redacted
        // This can happen if read was aligned to hyplotype partially. In this case
        // The tail (or head) of read to haplotype alignment would be soft-clipped.
        if (!read_to_haplotype_cigar_ops.empty() && haplotype_to_ref_cigar_ops.empty()) {
            MergeCigarOp(read_to_haplotype_cigar_ops.front(), read_len, read_to_ref_cigar_ops);
            read_to_haplotype_cigar_ops.pop_front();
            continue;
        }

        // Read is aligned completely, we are done.
        if (read_to_haplotype_cigar_ops.empty() && !haplotype_to_ref_cigar_ops.empty()) {
            break;
        }

        // Assign current Cigar Ops for each alignment.
        CigarOp cur_read_to_hap_op = read_to_haplotype_cigar_ops.front();
        read_to_haplotype_cigar_ops.pop_front();
        CigarOp cur_hap_to_ref_op = haplotype_to_ref_cigar_ops.front();
        haplotype_to_ref_cigar_ops.pop_front();

        // We look at cur_read_to_hap_op, cur_hap_to_ref_op.
        // For each of the op, they can be either MATCH(M), INS(I), or DEL(D).
        // In addition first or last read operation can be a SOFT CLIP (S)
        // As a result, we need to consider 9 combinations (soft clips are treated
        // the same way as match).
        // Out of 9 combinations we don not consider INS/DEL and DEL/INS. Those
        // cases are skipped due to an ambiguity in general case as well as due to
        // a very low impact. Reads that contain INS/DEL at the same position are
        // not realined.
        // cur_read_to_hap_op, cur_hap_to_ref_op = M|S, M|S
        if (BothOpsAreMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
            int new_op_len = min(cur_read_to_hap_op.length, cur_hap_to_ref_op.length);

            if (OneOfOpsIsSoftClip(cur_read_to_hap_op, cur_hap_to_ref_op))
                MergeCigarOp(CigarOp(CIGAR_OPERATIONS::SOFT_CLIP, new_op_len), read_len, read_to_ref_cigar_ops);
            else
                MergeCigarOp(CigarOp(CIGAR_OPERATIONS::MATCH, new_op_len), read_len, read_to_ref_cigar_ops);

            cur_read_to_hap_op.length -= new_op_len;
            PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);
            cur_hap_to_ref_op.length -= new_op_len;
            PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);

            // cur_read_to_hap_op, cur_hap_to_ref_op = D, M
        } else if (DelAndMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::DEL, cur_read_to_hap_op.length), read_len, read_to_ref_cigar_ops);
            cur_hap_to_ref_op.length -= cur_read_to_hap_op.length;
            PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);
            // cur_read_to_hap_op, cur_hap_to_ref_op = M, D
        } else if (DelAndMatch(cur_hap_to_ref_op, cur_read_to_hap_op)) {
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::DEL, cur_hap_to_ref_op.length), read_len, read_to_ref_cigar_ops);
            PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);

            // cur_read_to_hap_op, cur_hap_to_ref_op = D, D
        } else if (BothOpsAreDel(cur_read_to_hap_op, cur_hap_to_ref_op)) {
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::DEL, cur_hap_to_ref_op.length + cur_read_to_hap_op.length),
                         read_len, read_to_ref_cigar_ops);
            // cur_read_to_hap_op, cur_hap_to_ref_op = I, M
        } else if (InsAndMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
            cur_read_to_hap_op.length = min(read_len - AlignedLength(*read_to_ref_cigar_ops), cur_read_to_hap_op.length);
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::IN, cur_read_to_hap_op.length), read_len, read_to_ref_cigar_ops);
            PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);

            // cur_read_to_hap_op, cur_hap_to_ref_op = M, I
        } else if (InsAndMatch(cur_hap_to_ref_op, cur_read_to_hap_op)) {
            cur_hap_to_ref_op.length = min(read_len - AlignedLength(*read_to_ref_cigar_ops), cur_hap_to_ref_op.length);
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::IN, cur_hap_to_ref_op.length), read_len, read_to_ref_cigar_ops);
            // We need to decrease the length of cur_read_to_hap_op by INS length
            cur_read_to_hap_op.length = max(0, cur_read_to_hap_op.length - cur_hap_to_ref_op.length);
            PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);

            // cur_read_to_hap_op, cur_hap_to_ref_op = I, I
        } else if (BothOpsAreIns(cur_hap_to_ref_op, cur_read_to_hap_op)) {
            cur_hap_to_ref_op.length = cur_hap_to_ref_op.length + cur_read_to_hap_op.length;
            MergeCigarOp(CigarOp(CIGAR_OPERATIONS::IN, cur_hap_to_ref_op.length), read_len, read_to_ref_cigar_ops);
            // In all other cases read realignment is discarded.
            // redacted
        } else {
//            cerr << "read " << static_cast<int>(read_index)
//                 << ", could not be aligned, alignedLength=" << AlignedLength(*read_to_ref_cigar_ops);
            read_to_ref_cigar_ops->clear();
            return;
        }
    }  // while
}


vector<type_read> ReadAligner::align_reads_to_reference(vector<type_read> reads) {
    // set the reference sequence
    SSWAligner.set_reference(reference_sequence);
    vector<type_read> realigned_reads;
    for(auto &read: reads) {
        // go through each of the alignments and try to align it
        Alignment alignment = SSWAligner.align(read.sequence);

        if(alignment.sw_score > 1) {
            // create a new read
            type_read realigned_read;
            realigned_read = read;

            // get all the cigar operations
            list <CigarOp> cigarOps = CigarStringToVector(alignment.cigar_string);
            vector<CigarOp> cigar_tuples;
            for (auto& op : cigarOps) {
                cigar_tuples.push_back(op);
            }
            // if not empty then add this to realigned reads
            if (!cigar_tuples.empty()) {
                realigned_read.set_cigar_tuples(cigar_tuples);
                realigned_read.set_position(region_start + alignment.ref_begin);
                realigned_read.set_end_position(region_start + alignment.ref_end);
            }
            realigned_reads.push_back(realigned_read);
        } else {
            realigned_reads.push_back(read);
        }
    }
    return realigned_reads;
}


vector<type_read> ReadAligner::align_haplotypes_to_reference(vector<string> haplotypes) {
    // set the reference sequence
    SSWAligner.set_reference(reference_sequence);
    vector<type_read> realigned_reads;
    for(auto &haplotype_sequence: haplotypes) {
        // go through each of the alignments and try to align it
        Alignment alignment = SSWAligner.align(haplotype_sequence);

        if(alignment.sw_score > 1) {
            // create a new read
            type_read realigned_read;
            realigned_read.sequence = haplotype_sequence;

            // get all the cigar operations
            list <CigarOp> cigarOps = CigarStringToVector(alignment.cigar_string);
            vector<CigarOp> cigar_tuples;
            for (auto& op : cigarOps) {
                cigar_tuples.push_back(op);
            }
            // if not empty then add this to realigned reads
            if (!cigar_tuples.empty()) {
                realigned_read.set_cigar_tuples(cigar_tuples);
                realigned_read.set_position(region_start + alignment.ref_begin);
                realigned_read.set_end_position(region_start + alignment.ref_end);
            }
            realigned_reads.push_back(realigned_read);
        } else {
            cerr<<"HAPLOTYPE ALIGNMENT FAILED: "<<region_start<<" "<<region_end<<endl;
        }
    }
    return realigned_reads;
}


vector<type_read> ReadAligner::align_reads(vector<string> haplotypes, vector<type_read> reads) {
    this->haplotypes = haplotypes;

    // copy reads
    for(auto &read: reads) {
        read_sequences.push_back(read.sequence);
    }
    // infer alignment score threshold
    int read_size = read_sequences[0].length();
    alignment_score_threshold = Aligner_options::match * read_size * similarity_threshold
                                - Aligner_options::mismatch * read_size * (1 - similarity_threshold);
    if (alignment_score_threshold < 0) {
        alignment_score_threshold = 1;
    }
    // create kmer map
    CreateKmerMap();

    AlignReadsToHaplotypes();

    AlignHaplotypesToReference();

    CalculatePositionMaps();

    SswAlignReadsToHaplotypes(alignment_score_threshold);

    sort(read_to_haplotype_alignments.begin(), read_to_haplotype_alignments.end());

    vector<type_read> realigned_reads;

    RealignReadsToReference(reads, realigned_reads);

    return realigned_reads;
}


inline bool BothOpsAreMatch(const CigarOp& op1, const CigarOp& op2) {
    return (op1.operation == CIGAR_OPERATIONS::MATCH || op1.operation == CIGAR_OPERATIONS::SOFT_CLIP) &&
           (op2.operation == CIGAR_OPERATIONS::MATCH || op2.operation == CIGAR_OPERATIONS::SOFT_CLIP);
}

inline bool OneOfOpsIsSoftClip(const CigarOp& op1, const CigarOp& op2) {
    return op1.operation == CIGAR_OPERATIONS::SOFT_CLIP || op2.operation == CIGAR_OPERATIONS::SOFT_CLIP;
}

inline bool DelAndMatch(const CigarOp& op1, const CigarOp& op2) {
    return op1.operation == CIGAR_OPERATIONS::DEL &&
           (op2.operation == CIGAR_OPERATIONS::MATCH || op2.operation == CIGAR_OPERATIONS::SOFT_CLIP);
}

inline bool BothOpsAreDel(const CigarOp& op1, const CigarOp& op2) {
    return op1.operation == CIGAR_OPERATIONS::DEL && op2.operation == CIGAR_OPERATIONS::DEL;
}

inline bool InsAndMatch(const CigarOp& op1, const CigarOp& op2) {
    return op1.operation == CIGAR_OPERATIONS::IN &&
           (op2.operation == CIGAR_OPERATIONS::MATCH || op2.operation == CIGAR_OPERATIONS::SOFT_CLIP);
}

inline bool BothOpsAreIns(const CigarOp& op1, const CigarOp& op2) {
    return (op1.operation == CIGAR_OPERATIONS::IN && op2.operation == CIGAR_OPERATIONS::IN);
}

inline void PushFrontIfNotEmpty(const CigarOp& op, list<CigarOp>* cigar) {
    if (cigar == nullptr) {
        return;
    }
    if (op.length > 0) {
        cigar->push_front(op);
    }
}