import h5py
import argparse
import sys
from modules.python.TextColor import TextColor
from modules.python.CandidateFinder import find_candidates
from modules.python.VcfWriter import VCFWriter
from modules.python.ImageGenerationUI import UserInterfaceSupport


def candidates_to_variants(candidates, contig):
    alleles = []
    ref_start = candidates[0][0]
    ref_end = candidates[0][1]
    ref_seq = candidates[0][2]

    if len(candidates) > 1:
        if candidates[1][1] > candidates[0][1]:
            ref_end = candidates[1][1]
            ref_seq = candidates[1][2]

    for candidate in candidates:
        changed_candidate = list(candidate)
        if candidate[1] < ref_end:
            bases_needed = ref_end - candidate[1]
            ref_suffix = ref_seq[-bases_needed:]
            changed_candidate[1] = ref_end
            changed_candidate[2] = ref_seq
            changed_candidate[3] = candidate[3] + ref_suffix
        alleles.append(changed_candidate[3])

    genotype = [0, 0]
    if len(alleles) == 1:
        genotype = [0, 1]
    elif alleles[0] == alleles[1]:
        alleles = list(set(alleles))
        genotype = [1, 1]
    else:
        genotype = [1, 2]

    return contig, ref_start, ref_end, ref_seq, alleles, genotype


def write_vcf(contig, all_candidate_positions, candidate_positional_map, vcf_file):
    # print(candidate_map)
    # candidate_map = {2931716: {(2931716, 2931719, 'CTT', 'C', 1, 'DEL'), (2931716, 2931718, 'CT', 'C', 2, 'DEL')}}
    for pos in sorted(all_candidate_positions):
        candidates = list()
        if pos in candidate_positional_map:
            candidates.append(candidate_positional_map[pos])

        if len(candidates) == 0:
            continue

        if len(candidates) > 2:
            print("ERROR: OBSERVED MORE THAN 2 CANDIDATES AT SITE: ", contig, pos, candidates)
            exit(1)

        variant = candidates_to_variants(list(candidates), contig)
        vcf_file.write_vcf_records(variant)


def candidate_finder(hdf_file_path, reference_file, sample_name, output_path, threads):
    with h5py.File(hdf_file_path, 'r') as hdf5_file:
        contigs = list(hdf5_file['predictions'].keys())

    # candidates_file = open(output_path+'candidates.tsv', 'w')
    vcf_file = VCFWriter(reference_file, contigs, sample_name, output_path, "candidates_as_variants")

    for contig in contigs:
        sys.stderr.write(TextColor.YELLOW + "INFO: PROCESSING CONTIG: " + contig + "\n" + TextColor.END)

        with h5py.File(hdf_file_path, 'r') as hdf5_file:
            chunk_keys = sorted(hdf5_file['predictions'][contig].keys())

        all_candidate_positions, candidate_positional_map = \
            find_candidates(hdf_file_path,  reference_file, contig, chunk_keys, threads)

        sys.stderr.write(TextColor.BLUE + "INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                         + str(len(all_candidate_positions)) + " " + ".\n" + TextColor.END)

        # output TSV
        # header = ["CONTIG", "START", "END_POS", "REF", "ALT", "HP", "TYPE"]
        # header_out = '\t'.join([str(item) for item in list(header)])
        # candidates_file.write(header_out + "\n")

        # if all_candidate_positions is not None and len(all_candidate_positions) > 0:
        #     for candidate_pos in sorted(all_candidate_positions):
        #         if candidate_pos in candidate_positional_map_h1:
        #             candidate = candidate_positional_map_h1[candidate_pos]
        #             candidate_out = '\t'.join([contig] + [str(item) for item in list(candidate)])
        #             candidates_file.write(candidate_out + "\n")
        #
        #         if candidate_pos in candidate_positional_map_h2:
        #             candidate = candidate_positional_map_h2[candidate_pos]
        #             candidate_out = '\t'.join([contig] + [str(item) for item in list(candidate)])
        #             candidates_file.write(candidate_out + "\n")

        write_vcf(contig, all_candidate_positions, candidate_positional_map, vcf_file)

    hdf5_file.close()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="3_pepper_stitch.py performs the final stitching to generate  "
                                                 "the polished sequences.")
    parser.add_argument(
        "-i",
        "--input_hdf",
        type=str,
        required=True,
        help="Input hdf prediction file."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Input reference fasta file."
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=True,
        help="Name of sample."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads."
    )

    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = UserInterfaceSupport.handle_output_directory(FLAGS.output_dir)
    candidate_finder(FLAGS.input_hdf, FLAGS.reference, FLAGS.sample_name, FLAGS.output_dir, FLAGS.threads)
