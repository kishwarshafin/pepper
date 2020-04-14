//
// Created by Kishwar Shafin on 10/21/18.
//


// REUSE DECLARATION: THIS SCRIPT HEAVILY REUSES THE METHODS USED IN DEEPVARIANT https://github.com/google/deepvariant/
// LICENSE INCLUDED IN: third_party/deepVariant.LICENSE
#include "../../headers/realignment/simple_aligner.h"
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

ReadAligner::ReadAligner(int ref_start, int ref_end, string ref_seq) {
    reference_sequence = ref_seq;
    region_start = ref_start;
    region_end = ref_end;
}

vector<type_read> ReadAligner::align_reads_to_reference(vector<type_read> reads) {
    // set the reference sequence
//    SSWAligner.set_reference(reference_sequence);
    vector<type_read> realigned_reads;
    for(auto &read: reads) {
//        cout<<read.pos<<" "<<read.pos_end<<" "<<region_start<<" "<<region_end<<" "<<endl;
        if(read.pos < region_start) {
            // read should be between the selected region.
            cerr<<"READ STARTS BEFORE REGION: "<<read.pos<<" "<<region_start<<" "<<region_end<<endl;
            continue;
        }
        long long reference_start_index = read.pos - region_start;
        // do a left trim of the sequence and set as a reference sequence
        SSWAligner.set_reference(reference_sequence.substr(reference_start_index));

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
                realigned_read.set_position(read.pos + alignment.ref_begin);
                realigned_read.set_end_position(read.pos + alignment.ref_end);
            }
            realigned_reads.push_back(realigned_read);
        } else {
            realigned_reads.push_back(read);
        }
    }
    return realigned_reads;
}