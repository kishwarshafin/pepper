from pepper.build import PEPPER
import argparse
from collections import defaultdict

"""
    static constexpr int MATCH = 0;
    static constexpr int EQUAL = 7;
    static constexpr int DIFF = 8;
    
    static constexpr int IN = 1;
    
    static constexpr int DEL = 2;
    static constexpr int REF_SKIP = 3;
    static constexpr int PAD = 6;
    
    static constexpr int SOFT_CLIP = 4;
    static constexpr int HARD_CLIP = 5;
    
    static constexpr int BACK = 9;
    static constexpr int UNSPECIFIED = -1;
"""
MATCH_CIGARS = [0, 7, 8]
IN_CIGARS = [1]
DEL_CIGARS = [2, 3, 6]
CLIP_CIGARS = [4, 5]

def pileup_from_reads(reference, start_pos, end_pos, reads):

    longest_insert_count = defaultdict(int)
    insert_dict = defaultdict(lambda: defaultdict(int))
    base_dict = defaultdict(lambda: defaultdict(int))
    read_start_pos = defaultdict(int)
    read_end_pos = defaultdict(int)
    query_names = list()

    for read in reads:
        read_start_pos[read.query_name] = read.pos
        read_end_pos[read.query_name] = read.pos_end
        cigar_tuples = read.cigar_tuples
        read_index = 0
        reference_position = read.pos
        query_names.append(read.query_name)

        for cigar_tup in cigar_tuples:
            # match
            if cigar_tup.cigar_op in MATCH_CIGARS:
                for i in range(cigar_tup.cigar_len):
                    base_dict[read.query_name][reference_position] = read.sequence[read_index]
                    read_index += 1
                    reference_position += 1
            # insert
            if cigar_tup.cigar_op in IN_CIGARS:
                longest_insert_count[reference_position] = max(longest_insert_count[reference_position], cigar_tup.cigar_len)
                in_allele = ""
                for i in range(cigar_tup.cigar_len):
                    in_allele += read.sequence[read_index]
                    read_index += 1

                insert_dict[read.query_name][reference_position] = in_allele
            # delete
            if cigar_tup.cigar_op in DEL_CIGARS:
                for i in range(cigar_tup.cigar_len):
                    base_dict[read.query_name][reference_position] = '-'
                    reference_position += 1
            # soft-clip
            if cigar_tup.cigar_op in CLIP_CIGARS:
                read_index += cigar_tup.cigar_len

    for i in range(start_pos, end_pos):
        ref_base = reference[i - start_pos]
        print(ref_base, end='')
        if longest_insert_count[i]:
            print("*" * longest_insert_count[i], end='')
    print()

    for read_name in query_names:
        for i in range(start_pos, end_pos):
            # if read_start_pos[read_name] < i:
            #     print(' ', end='')
            #     continue
            # if read_end_pos[read_name] < i:
            #     break

            if base_dict[read_name][i]:
                read_base = base_dict[read_name][i]
            else:
                read_base = ' '
            print(read_base, end='')

            if insert_dict[read_name][i] and i < read_end_pos[read_name]:
                print(insert_dict[read_name][i], end='')
                star_needed = longest_insert_count[i] - int(len(insert_dict[read_name][i]))
                if star_needed > 0:
                    print("*" * star_needed, end='')
            elif longest_insert_count[i] and i < read_end_pos[read_name]:
                print("*" * longest_insert_count[i], end='')
        print(read_name)


def msa_align(sequences):
    import pyabpoa as pa
    a = pa.msa_aligner(aln_mode='g')
    res = a.msa(sequences, out_cons=True, out_msa=True)
    res.print_msa()

    return res.cons_seq

def generate_pileup(bam_file_path, fasta_file_path, region):
    chr_name, start_end = region.rstrip().split(':')
    start_position, end_position = start_end.rstrip().split('-')
    start_position = int(start_position)
    end_position = int(end_position)
    bam_handler = PEPPER.BAM_handler(bam_file_path)
    fasta_handler = PEPPER.FASTA_handler(fasta_file_path)
    all_reads = bam_handler.get_reads(chr_name,
                                      start_position,
                                      end_position,
                                      False,
                                      60,
                                      0)

    min_start = start_position
    max_end = end_position
    for read in all_reads:
        min_start = min(read.pos, min_start)
        max_end = max(read.pos_end, max_end)


    print("ORIGINAL ALIGNMENT:")
    # MINIMAP2 ALIGNMENT
    # reference_sequence = fasta_handler.get_reference_sequence(chr_name, start_position, end_position)
    reference_sequence = fasta_handler.get_reference_sequence(chr_name, min_start, max_end)
    pileup_from_reads(reference_sequence, start_position, end_position, all_reads)

    # MSA ALIGNMENT
    all_sequences = [reference_sequence]
    for read in all_reads:
        all_sequences.append(read.sequence)

    print("MSA ALIGNMENT:")
    consensus_sequence = msa_align(all_sequences)

    print("REFERENCE SEQUENCE")
    print(reference_sequence)

    print("CONSENSUS SEQUENCE")
    print(consensus_sequence[0])

    reference_sequence_file = open("ref_seq.fa", 'w')
    reference_sequence_file.write(">ref\n")
    reference_sequence_file.write(reference_sequence + "\n")

    consensus_sequence_file = open("consensus_seq.fa", 'w')
    consensus_sequence_file.write(">consensus\n")
    consensus_sequence_file.write(consensus_sequence[0] + "\n")

    # import pyabpoa as pa
    # all_sequences.append(consensus_sequence)
    # a = pa.msa_aligner(aln_mode='g')
    # res = a.msa(all_sequences, out_cons=True, out_msa=True)
    # res.print_msa()

    # print("AFTER PAIRWISE ALIGNMENT")
    # # PAIRWISE ALIGNMENT

    aligner = PEPPER.Aligner()
    filter = PEPPER.Filter()
    alignment = PEPPER.Alignment()
    aligner.SetReferenceSequence(consensus_sequence[0], len(consensus_sequence[0]))
    for read in all_reads:
        print(read.pos)
        aligner.Align_cpp(read.sequence, filter, alignment, 0)
        print(alignment.cigar_string)
        print(alignment.reference_begin)


    # realigned_reads = aligner.align_reads_to_reference(all_reads)
    # pileup_from_reads(reference_sequence, start_position, end_position, realigned_reads)








if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="realigner\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--bam",
        required=True,
        type=str,
        help="path to bam."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=str,
        help="path to fasta."
    )
    parser.add_argument(
        "--region",
        required=True,
        type=str,
        help="region in the format chrX:XXX-XXX"
    )

    FLAGS, unparsed = parser.parse_known_args()
    generate_pileup(FLAGS.bam, FLAGS.fasta, FLAGS.region)
