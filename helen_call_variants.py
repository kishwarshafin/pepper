import h5py
import argparse
import sys
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator
from modules.python.Haplotype2Variant import Haplotype2Variant
from modules.python.VcfWriter import VCFWriter
from modules.python.TextColor import TextColor

label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
BASE_ERROR_RATE = 0.1


def chunks(file_names, threads):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(file_names), threads):
        chunks.append(file_names[i:i + threads])
    return chunks


def find_variants(file_name, ref_file_path, contig, small_chunk_keys):
    # for chunk_key in small_chunk_keys:
    all_variants = list()

    for chunk_name in small_chunk_keys:
        with h5py.File(file_name, 'r') as hdf5_file:
            other_keys = ['contig', 'region_start', 'region_end']
            haplotypes = set(hdf5_file['predictions'][contig][chunk_name].keys()) - set(other_keys)
            contig = hdf5_file['predictions'][contig][chunk_name]['contig'][()]
            region_start = hdf5_file['predictions'][contig][chunk_name]['region_start'][()]
            region_end = hdf5_file['predictions'][contig][chunk_name]['region_end'][()]

        haplotype_sequences = list()
        for haplotype in haplotypes:
            with h5py.File(file_name, 'r') as hdf5_file:
                chunk_ids = list(hdf5_file['predictions'][contig][chunk_name][haplotype].keys())

            positions = set()
            prediction_dict = defaultdict()
            for chunk in chunk_ids:
                with h5py.File(file_name, 'r') as hdf5_file:
                    predictions = hdf5_file['predictions'][contig][chunk_name][haplotype][chunk]['predictions'][()]
                    position = hdf5_file['predictions'][contig][chunk_name][haplotype][chunk]['position'][()]
                    index = hdf5_file['predictions'][contig][chunk_name][haplotype][chunk]['index'][()]

                position = np.array(position, dtype=np.int64)
                predictions = np.array(predictions, dtype=np.int)
                index = np.array(index, dtype=np.int)

                for pos, indx, pred in zip(position, index, predictions):
                    if (pos, indx) not in prediction_dict:
                        prediction_dict[(pos, indx)] = pred
                        positions.add((pos, indx))

            pos_list = sorted(list(positions), key=lambda element: (element[0], element[1]))
            dict_fetch = operator.itemgetter(*pos_list)
            predicted_labels = list(dict_fetch(prediction_dict))
            sequence = ''.join([label_decoder[x] for x in predicted_labels])
            haplotype_sequences.append(sequence)

        hap_2_variant = Haplotype2Variant(ref_file_path, contig, region_start, region_end)
        variant_set = hap_2_variant.generate_records_from_haplotypes(haplotype_sequences)
        all_variants.extend(variant_set)

    return all_variants


def get_variants_from_contig(hdf5_file_path, ref_file_path, contig, sequence_chunk_keys, threads):
    contig_variants = list()
    record_set = set()

    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_keys, 1)
        futures = [executor.submit(find_variants, hdf5_file_path, ref_file_path, contig, file_chunk)
                   for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                variant_set = fut.result()

                for variant in variant_set:
                    chromosome_name, pos, pos_end, ref, alternate_alleles, genotype = variant
                    if (chromosome_name, pos) in record_set:
                        continue

                    record_set.add((chromosome_name, pos))
                    contig_variants.append(variant)
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL " + str(len(contig_variants)) + " VARIANTS FOUND IN " + str(contig)
                     + "\n" + TextColor.END)

    return contig_variants


def process_marginpolish_h5py(bam_file_path, sample_name, hdf_file_path, ref_path_path, output_path, threads):
    with h5py.File(hdf_file_path, 'r') as hdf5_file:
        contigs = sorted(list(hdf5_file['predictions'].keys()))

    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL " + str(len(contigs)) + " CONTING FOUND." + "\n" + TextColor.END)

    vcf_file = VCFWriter(bam_file_path, sample_name, output_path, contigs)

    for contig in contigs:
        sys.stderr.write(TextColor.GREEN + "INFO: CALLING VARIANT ON " + str(contig) + "." + "\n" + TextColor.END)
        with h5py.File(hdf_file_path, 'r') as hdf5_file:
            chunk_keys = sorted(hdf5_file['predictions'][contig].keys())

        contig_variants = get_variants_from_contig(hdf_file_path, ref_path_path, contig, chunk_keys, threads)
        sys.stderr.write(TextColor.GREEN + "INFO: WRITING VARIANTS OF " + str(contig) + "." + "\n" + TextColor.END)
        vcf_file.add_variants(contig_variants)
        sys.stderr.write(TextColor.GREEN + "INFO: FINISHED CALLING VARIANTS ON " + str(contig) + "." + "\n"
                         + TextColor.END)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequence_hdf",
        type=str,
        required=True,
        help="H5PY file generated by HELEN."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference FASTA file."
    )
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file path."
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        required=True,
        help="Sample name."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="CONSENSUS output directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )

    FLAGS, unparsed = parser.parse_known_args()
    process_marginpolish_h5py(FLAGS.bam, FLAGS.sample_name, FLAGS.sequence_hdf,
                              FLAGS.ref, FLAGS.output_dir, FLAGS.threads)
    # read_marginpolish_h5py(FLAGS.marginpolish_h5py_dir, FLAGS.output_h5py_dir, FLAGS.train_mode, FLAGS.threads)
