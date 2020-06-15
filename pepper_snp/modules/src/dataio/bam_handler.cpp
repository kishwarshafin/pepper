//
// Created by Kishwar Shafin on 6/14/18.
//
#include "../../headers/dataio/bam_handler.h"

BAM_handler::BAM_handler(string path) {
    this->hts_file = sam_open(path.c_str(), "r");
    // Used to check if valid bam file
    if(this->hts_file == NULL){
        cerr<<"INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT: "<<path<<endl;
        exit (EXIT_FAILURE);
    }

    // Setting up the indexing
    this->idx = sam_index_load(this->hts_file, path.c_str());
    if(this->idx == NULL){
        cerr<<"INVALID BAM INDEX FILE. PLEASE CHECK IF FILE IS INDEXED: "<<path<<endl;
        exit (EXIT_FAILURE);
    }

    // setting up the header
    this->header = sam_hdr_read(this->hts_file);
    if(this->header == 0){
        cerr<<"HEADER ERROR: INVALID BAM FILE. PLEASE CHECK IF PATH IS CORRECT."<<endl;
        exit (EXIT_FAILURE);
    }

}

set<string> BAM_handler::get_sample_names() {
    int l_text = header->l_text;
    char *text = header->text;
    set<string> samples;
    stringstream full_header_tokenizer(text);
    string line;
    while(getline(full_header_tokenizer, line, '\n')){
        istringstream header_line_tokenizer(line);
        string token;
        getline(header_line_tokenizer, token, '\t');
        if(token == "@RG") {
            while(getline(header_line_tokenizer, token, '\t')){
                istringstream tag_tokenizer(token);
                string tag;
                getline(tag_tokenizer, tag, ':');
                if(tag == "SM"){
                    string sample_name;
                    getline(tag_tokenizer, sample_name, ':');
                    samples.insert(sample_name);
                }
            }
        }
    }
    bam_hdr_destroy(header);
    return samples;
}

static inline int HtslibAuxSize(uint8_t type) {
    switch (type) {
        case 'A': case 'c': case 'C':
            return 1;
        case 's': case 'S':
            return 2;
        case 'f': case 'i': case 'I':
            return 4;
        default:
            return -1;
    }
}

type_read_flags BAM_handler::get_read_flags(int flag) {
    type_read_flags flags;
    if ( flag&BAM_FPAIRED ) flags.is_paired = 1;
    if ( flag&BAM_FPROPER_PAIR ) flags.is_proper_pair = 1;
    if ( flag&BAM_FUNMAP ) flags.is_unmapped = 1;
    if ( flag&BAM_FMUNMAP ) flags.is_mate_unmapped = 1;
    if ( flag&BAM_FREVERSE ) flags.is_reverse = 1;
    if ( flag&BAM_FMREVERSE ) flags.is_mate_is_reverse = 1;
    if ( flag&BAM_FREAD1 ) flags.is_read1 = 1;
    if ( flag&BAM_FREAD2 ) flags.is_read2 = 1;
    if ( flag&BAM_FSECONDARY ) flags.is_secondary = 1;
    if ( flag&BAM_FQCFAIL ) flags.is_qc_failed = 1;
    if ( flag&BAM_FDUP ) flags.is_duplicate = 1;
    if ( flag&BAM_FSUPPLEMENTARY ) flags.is_supplementary = 1;
    return flags;
}

vector<type_sequence> BAM_handler::get_chromosome_sequence_names_with_length() {
    // Get all the sequence names. These are the chromosome names from the bed header file.
    vector<type_sequence> sequence_names;
    int total_targets = header->n_targets;
    for (int i=0; i < total_targets; i++){
        type_sequence sq_info;
        sq_info.sequence_name = header->target_name[i];
        sq_info.sequence_length = header->target_len[i];
        sequence_names.push_back(sq_info);
    }

    bam_hdr_destroy(header);

    return sequence_names;
}

vector<string> BAM_handler::get_chromosome_sequence_names() {
    // Get all the sequence names. These are the chromosome names from the bed header file.
    vector<string> sequence_names;
    int total_targets = header->n_targets;
    for (int i=0; i < total_targets; i++){
        string seq_name = header->target_name[i];
        sequence_names.push_back(seq_name);
    }

    return sequence_names;
}

vector<type_read> BAM_handler::get_reads(string chromosome,
                                                   long long start,
                                                   long long stop,
                                                   bool include_supplementary,
                                                   int min_mapq=0,
                                                   int min_baseq = 0) {
    // safe bases
//    stop += 0;

    vector <type_read> all_reads;

    // get the id of the chromosome
    const int tid = bam_name2id(this->header, chromosome.c_str());

    // get the iterator
    hts_itr_t *iter  = sam_itr_queryi(this->idx, tid, start, stop);

    // initialize an alignment
    bam1_t* alignment = bam_init1();

    while(sam_itr_next(this->hts_file, iter, alignment) >= 0) {
        //get read flags
        type_read_flags read_flags = get_read_flags(alignment->core.flag);

        if (read_flags.is_qc_failed || read_flags.is_duplicate || read_flags.is_secondary
            || read_flags.is_unmapped) {
            continue;
        }
        if (!include_supplementary && read_flags.is_supplementary) {
            continue;
        }

        // mapping quality
        if (alignment->core.qual < min_mapq) {
            continue;
        }

        // get query name
        string query_name = bam_get_qname(alignment);

        // get the position where it gets aligned
        int32_t pos = alignment->core.pos;
        uint32_t len = alignment->core.l_qseq;

        // get the base qualities and sequence bases
        uint8_t *seqi = bam_get_seq(alignment);
        uint8_t *qual = bam_get_qual(alignment);

        vector<int> base_qualities;
        vector<int> bad_bases;
        string read_seq;

        // get the cigar operations of the alignment
        uint32_t *cigar = bam_get_cigar(alignment);
        string str_cigar;
        vector <CigarOp> cigar_tuples;
        long long pos_start = -1;
        long long pos_end = -1;

        long long current_read_pos = pos;
        int current_read_index = 0;
        int running_sequence_index = 0;

        // this is a bit ambitious, we are cutting all the reads to desired regions so we don't have to deal with ultra-long reads
        // I am not sure if there are any downside to it, but I would really love the speed-up
        for (int k = 0; k < alignment->core.n_cigar; k++) {
            // we are going on all cigar operations to cut the reads short
            int cigar_op = bam_cigar_op(cigar[k]);
            int cigar_len = bam_cigar_oplen(cigar[k]);
            int modified_cigar_length;
            int cigar_index;
            if (current_read_pos > stop) {
                break;
            }

            switch (cigar_op) {
                case BAM_CMATCH:
                case BAM_CDIFF:
                case BAM_CEQUAL:
                    cigar_index = 0;
                    // if the current read position is to the left then we jump forward
                    if (current_read_pos < start) {
                        // jump as much as we can but not more than the boundary of start
                        cigar_index = min(start - current_read_pos, (long long) cigar_len);
                        current_read_index += cigar_index;
                        current_read_pos += cigar_index;
                    }
                    // once we've made the jump now add each of the elements
                    modified_cigar_length = 0;
                    for (int i = cigar_index; i < cigar_len; i++) {
                        if (current_read_pos <= stop) {
                            if (pos_start == -1) {
                                pos_start = current_read_pos;
                                pos_end = pos_start;
                            }
                            // we are adding the base and quality
                            int base_quality = (int) qual[current_read_index];
                            base_qualities.push_back(base_quality);
                            char base = ::toupper(seq_nt16_str[bam_seqi(seqi, current_read_index)]);
                            read_seq += base;

                            if (base_quality < min_baseq or
                                (base != 'A' &&
                                 base != 'C' &&
                                 base != 'G' &&
                                 base != 'T')) {
                                bad_bases.push_back(running_sequence_index);
                            }
                            running_sequence_index += 1;
                            modified_cigar_length += 1;
                            pos_end += 1;
                        } else break;

                        current_read_index += 1;
                        current_read_pos += 1;
                    }
                    if (modified_cigar_length > 0) {
                        // save the cigar tuple now
                        CigarOp cigar_instance;

                        cigar_instance.operation = cigar_op;
                        cigar_instance.length = modified_cigar_length;
                        cigar_tuples.push_back(cigar_instance);
                    }
                    break;
                case BAM_CSOFT_CLIP:
                case BAM_CINS:
                    modified_cigar_length = 0;
                    // this only happens after the first position, I am also forcing an anchor
                    if (current_read_pos >= start && current_read_pos <= stop && pos_start != -1) {
                        for (int i = 0; i < cigar_len; i++) {
                            // we are adding the base and quality
                            int base_quality = (int) qual[current_read_index];
                            base_qualities.push_back(base_quality);
                            char base = ::toupper(seq_nt16_str[bam_seqi(seqi, current_read_index)]);
                            read_seq += base;

                            if (base_quality < min_baseq or
                                (base != 'A' &&
                                 base != 'C' &&
                                 base != 'G' &&
                                 base != 'T')) {
                                bad_bases.push_back(running_sequence_index);
                            }
                            running_sequence_index += 1;
                            modified_cigar_length += 1;
                            current_read_index += 1;
                        }

                    } else {
                        current_read_index += cigar_len;
                    }
                    if (modified_cigar_length > 0) {
                        // save the cigar tuple now
                        CigarOp cigar_instance;

                        cigar_instance.operation = cigar_op;
                        cigar_instance.length = modified_cigar_length;
                        cigar_tuples.push_back(cigar_instance);
                    }
                    break;
                case BAM_CREF_SKIP:
                case BAM_CDEL:
                    modified_cigar_length = 0;
                    if (current_read_pos >= start && current_read_pos <= stop && pos_start != -1) {
                        modified_cigar_length = 0;
                        for (int i = 0; i < cigar_len; i++) {
                            if (current_read_pos <= stop) {
                                modified_cigar_length += 1;
                                pos_end += 1;
                            } else break;

                            current_read_pos += 1;
                        }

                    } else {
                        current_read_pos += cigar_len;
                    }
                    if (modified_cigar_length > 0) {
                        // save the cigar tuple now
                        CigarOp cigar_instance;

                        cigar_instance.operation = cigar_op;
                        cigar_instance.length = modified_cigar_length;
                        cigar_tuples.push_back(cigar_instance);
                    }
                    break;
                case BAM_CHARD_CLIP:
                    break;

            }
        }
        bad_bases.push_back(read_seq.length() + 1);

        // mapping quality
        int map_quality = alignment->core.qual;

        // handle auxiliary data
        uint8_t *s = bam_get_aux(alignment);
        const uint8_t *aux_end = alignment->data + alignment->l_data;
        int HP_tag = 0;
        // FORMAT OF TAG IS: TAG:TYPE:VALUE
        // WE ARE ONLY INTERESTED IN HP TAG WHICH SHOULD ALWAYS HAVE INTEGER VALUE

        bool tag_state_ok = true;
        while (aux_end - s >= 4 && tag_state_ok) {
            // Each block is encoded like (each element is a byte):
            // [tag char 1, tag char 2, type byte, ...]
            // where the ... contents depends on the 2-character tag and type.
            const string tag = string(reinterpret_cast<char *>(s), 2);
            s += 2;
            const uint8_t tag_type = *s++;

            switch (tag_type) {
                // An 'A' is just a single character string.
                case 'A': {
                    // Safe since we know s is at least 4 bytes from the end.
                    const string value = string(reinterpret_cast<char *>(s), 1);
                    s += 1;
                    if (tag == "HP")
                        cerr << "HP TAG HAS INVALID TYPE: " << tag_type << " " << query_name << endl;
                }
                    break;
                    // These are all different byte-sized integers.
                case 'C':
                case 'c':
                case 'S':
                case 's':
                case 'I':
                case 'i': {
                    // MOST OF THE TIMES HP WILL BE OF THIS TYPE
                    const int size = HtslibAuxSize(tag_type);
                    if (size < 0 || aux_end - s < size) {
                        tag_state_ok = false;
                        cerr << "INVALID TAG: " << tag << endl;
                        break;
                    }
                    errno = 0;
                    const int value = bam_aux2i(s - 1);
                    if (value == 0 && errno == EINVAL) {
                        tag_state_ok = false;
                        cerr << "INVALID TAG: " << tag << endl;
                        break;
                    }
                    // VALID TAG

                    if (tag == "HP") {
                        HP_tag = value;
                    }

                    s += size;
                }
                    break;
                    // A 4-byte floating point.
                case 'f': {
                    if (aux_end - s < 4) {
                        tag_state_ok = false;
                        cerr << "INVALID TAG: " << tag << endl;
                        break;
                    }
                    const float value = le_to_float(s);
                    // VALID TAG
                    if (tag == "HP")
                        cerr << "HP TAG HAS INVALID TYPE: " << tag_type << " " << query_name << endl;
                    s += 4;
                }
                    break;
                    // Z and H are null-terminated strings.
                case 'Z':
                case 'H': {
                    char *value = reinterpret_cast<char *>(s);
                    for (; s < aux_end && *s; ++s) {}  // Loop to the end.
                    if (s >= aux_end) {
                        tag_state_ok = false;
                        cerr << "INVALID TAG: " << tag << endl;
                    }
                    s++;
                    // VALID TAG
                    if (tag == "HP")
                        cerr << "HP TAG HAS INVALID TYPE: " << tag_type << " " << query_name << endl;
//                    if (type == 'Z') {
//                        SetInfoField(tag, value, read_message);
//                    }
                }
                    break;
                    // B is an array of atomic types (strings, ints, floats).
                case 'B': {
                    if (tag == "HP")
                        cerr << "HP TAG HAS INVALID TYPE: " << tag_type << " " << query_name << endl;

                    const uint8_t sub_type = *s++;
                    const int element_size = HtslibAuxSize(sub_type);
                    if (element_size < 0) {
                        tag_state_ok = false;
                        cerr << "SIZE == 0 for TAG: " << tag << endl;
                        break;
                    }
                    // Prevents us from reading off the end of our buffer with le_to_u32.
                    if (aux_end - s < 4) {
                        tag_state_ok = false;
                        break;
                    }
                    const int n_elements = le_to_u32(s);
                    if (n_elements == 0) cerr << "READ TAG: n_elements is zero" << endl;
                    s += 4 + n_elements * element_size;
                }
                    break;
                default: {
                    tag_state_ok = false;
                    cerr << "UNKNOWN TAG: " << tag << endl;
                    break;
                }
            }
        }

        // set all fetched attributes
        type_read read;
        if (read_seq.length() > 0) {
            read.query_name = query_name;
            read.pos = pos_start;
            read.pos_end = pos_end;
            read.sequence = read_seq;
            read.flags = read_flags;
            read.mapping_quality = map_quality;
            read.base_qualities = base_qualities;
            read.cigar_tuples = cigar_tuples;
            read.bad_indicies = bad_bases;
            read.hp_tag = HP_tag;

            all_reads.push_back(read);
        }
    }
    bam_destroy1(alignment);
    hts_itr_destroy(iter);

    return all_reads;
}

BAM_handler::~BAM_handler() {
    // destroy everything
    hts_idx_destroy( this->idx );
    sam_close( this->hts_file );
    bam_hdr_destroy( this->header );
}