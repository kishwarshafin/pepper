import h5py
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict


def decode_ref_base(ref_code):
    if ref_code == 0:
        return 'N'
    elif ref_code == 1:
        return 'A'
    elif ref_code == 2:
        return 'C'
    elif ref_code == 3:
        return 'G'
    elif ref_code == 4:
        return 'T'


def decode_bases(pred_code):
    if pred_code == 0:
        return ['*', '*']
    elif pred_code == 1:
        return ['A', 'A']
    elif pred_code == 2:
        return ['A', 'C']
    elif pred_code == 3:
        return ['A', 'T']
    elif pred_code == 4:
        return ['A', 'G']
    elif pred_code == 5:
        return ['A', '*']
    elif pred_code == 6:
        return ['C', 'C']
    elif pred_code == 7:
        return ['C', 'T']
    elif pred_code == 8:
        return ['C', 'G']
    elif pred_code == 9:
        return ['C', '*']
    elif pred_code == 10:
        return ['T', 'T']
    elif pred_code == 11:
        return ['T', 'G']
    elif pred_code == 12:
        return ['T', '*']
    elif pred_code == 13:
        return ['G', 'G']
    elif pred_code == 14:
        return ['G', '*']


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path)
                  if isfile(join(directory_path, file)) and file[-2:] == 'h5']
    return file_paths


def chunks(file_names, threads):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(file_names), threads):
        chunks.append(file_names[i:i + threads])
    return chunks


def chunks_alignment_sequence(alignment_sequence_pairs, min_length):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(alignment_sequence_pairs), min_length):
        chunks.append(alignment_sequence_pairs[i:i + min_length])
    return chunks


def find_candidates(contig, small_chunk_keys, p_threshold):
    # for chunk_key in small_chunk_keys:
    candidate_variants = defaultdict(list)
    reference_dict = defaultdict()

    for file_name, chunk_name in small_chunk_keys:

        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}

        smaller_chunks = sorted(smaller_chunks)

        for chunk in smaller_chunks:
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['bases'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]
                ref_seq = hdf5_file['predictions'][contig][chunk_name][chunk]['ref_seq'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions = np.array(bases, dtype=np.int)

            for pos, indx, base_pred, ref_code in zip(positions, indices, base_predictions, ref_seq):
                if indx < 0 or pos < 0 or indx > 0:
                    continue
                # update the reference base
                reference_base = decode_ref_base(ref_code)

                # get the predictions and see if this site has a candidate
                base_pred_norm = base_pred / sum(base_pred)
                max_pred_value = max(base_pred)
                # first see if this position has any candidates
                has_candidate = False
                for pred_code, (pred_value_norm, pred_value) in enumerate(zip(base_pred_norm, base_pred)):
                    if pred_value_norm >= p_threshold or pred_value >= max_pred_value:
                        pred_bases = decode_bases(pred_code)
                        for pred_base in pred_bases:
                            if pred_base != reference_base:
                                has_candidate = True
                                break

                # if it has candidates, then update the dictionaries
                if has_candidate:
                    for pred_code, (pred_value_norm, pred_value) in enumerate(zip(base_pred_norm, base_pred)):
                        if pred_value_norm >= p_threshold or pred_value >= max_pred_value:
                            pred_bases = decode_bases(pred_code)
                            predicted_bases = pred_bases

                            # candidate bases
                            reference_dict[pos] = reference_base
                            if predicted_bases not in candidate_variants[pos]:
                                candidate_variants[pos].append(predicted_bases)

    return contig, reference_dict, candidate_variants


def find_SNP_candidates(contig, sequence_chunk_keys, threads, p_threshold):
    all_candidates = defaultdict(list)
    global_reference_dict = defaultdict()
    positions_with_candidates = set()

    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_keys, max(2,
                                                      int(len(sequence_chunk_keys) / threads) + 1))

        futures = [executor.submit(find_candidates, contig, file_chunk, p_threshold)
                   for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                contig, reference_dict, candidate_variants = fut.result()
                for pos, alt_alleles in candidate_variants.items():
                    if alt_alleles:
                        global_reference_dict[pos] = reference_dict[pos]
                        positions_with_candidates.add(pos)
                        all_candidates[pos].extend(alt_alleles)
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144
    positions_with_candidates = list(positions_with_candidates)

    return all_candidates, global_reference_dict, positions_with_candidates
