from build import PEPPER
from collections import defaultdict


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
            if cigar_tup.cigar_op == 0:
                for i in range(cigar_tup.cigar_len):
                    base_dict[read.query_name][reference_position] = read.sequence[read_index]
                    read_index += 1
                    reference_position += 1
            # insert
            if cigar_tup.cigar_op == 1:
                longest_insert_count[reference_position] = max(longest_insert_count[reference_position], cigar_tup.cigar_len)
                in_allele = ""
                for i in range(cigar_tup.cigar_len):
                    in_allele += read.sequence[read_index]
                    read_index += 1

                insert_dict[read.query_name][reference_position] = in_allele
            # delete
            if cigar_tup.cigar_op == 2:
                for i in range(cigar_tup.cigar_len):
                    base_dict[read.query_name][reference_position] = '*'
                    reference_position += 1
            # soft-clip
            if cigar_tup.cigar_op == 4:
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

