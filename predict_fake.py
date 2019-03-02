import argparse
import sys
import torch
import numpy as np
from torch.utils.data import DataLoader
from modules.python.models.dataloader_fake import SequenceDataset
from torchvision import transforms
import multiprocessing
from modules.python.TextColor import TextColor
from collections import defaultdict
from modules.python.vcf_writer import VCFWriter
from modules.python.FileManager import FileManager
import operator
import pickle
from tqdm import tqdm
import os
import time
"""
This script uses a trained model to call variants on a given set of images generated from the genome.
The process is:
- Create a prediction table/dictionary using a trained neural network
- Convert those predictions to a VCF file

INPUT:
- A trained model
- Set of images for prediction

Output:
- A VCF file containing all the variants.
"""
SNP = 1
IN = 2
DEL = 3
HOM = 0
HET = 1
HOM_ALT = 2

prediction_dict = defaultdict(list)
chromosome_list = set()
position_dict = defaultdict(set)


def predict(test_file, batch_size, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param test_file: File to predict on
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param gpu_mode: If true, predictions will be done over GPU
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # the prediction table/dictionary
    chromosome_name = ''
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    # data loader
    test_data = SequenceDataset(test_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)
    sys.stderr.write(TextColor.CYAN + 'Test data loaded\n')

    # TO HERE
    with torch.no_grad():
        for labels, chromosome_name, position, index in tqdm(test_loader, ncols=50):
            start_index = 0
            end_index = labels.size(1)

            for seq_index in range(start_index, end_index):
                batches = labels.size(0)

                for batch in range(batches):

                    true_label = labels[batch, seq_index]
                    fake_probs = [0.0] * 6
                    fake_probs[true_label] = 1.0
                    top_n, top_i = torch.FloatTensor(fake_probs).topk(1)
                    predicted_label = top_i[0].item()

                    chromosome = chromosome_name[batch]
                    genomic_pos = position[batch, seq_index].item()
                    genomic_index = index[batch, seq_index].item()

                    if chromosome not in chromosome_list:
                        chromosome_list.add(chromosome)
                    position_dict[chromosome].add((genomic_pos, genomic_index))
                    prediction_dict[(chromosome, genomic_pos, genomic_index)].append(predicted_label)


def get_record_from_prediction(pos, positional_record):
    predictions = prediction_dict[pos]
    genotype, qual, gq = VCFWriter.process_prediction(pos, predictions)
    return positional_record, genotype, qual, gq


def produce_vcf_records(chromosome_name, output_dir, thread_no, pos_list):
    """
    Convert prediction dictionary to a VCF file
    :param: arg_tuple: Tuple of arguments containing these values:
    - chromosome_name: Chromosome name
    - pos_list: List of positions where we will search for variants
    - prediction_dict: prediction dictionary containing predictions of each image records
    - reference_dict: Dictionary containing reference information
    - bam_file_path: Path to the BAM file
    - sample_name: Name of the sample in the BAM file
    - output_dir: Output directory
    - thread_id: Unique id assigned to each thread
    :return:
    """
    # object that can write and handle VCF
    # vcf_writer = VCFWriter(bam_file_path, sample_name, output_dir, thread_id)
    # collate multi-allelic records to a single record
    current_allele_dict = ''
    allele_dict = {}
    record_file = open(output_dir + chromosome_name + "_" + str(thread_no) + ".tsv", 'w')
    for pos in pos_list:
        allele_dict_path = reference_dict[pos]
        if allele_dict_path != current_allele_dict:
            allele_dict = pickle.load(open(allele_dict_path, 'rb'))
            current_allele_dict = allele_dict_path

        if pos not in allele_dict:
            continue
        positional_record = allele_dict[pos] if pos in allele_dict else None
        if positional_record is None:
            continue

        positional_record, genotype, qual, gq = get_record_from_prediction(pos, positional_record)

        if genotype == '0/0':
            continue
        # print('BEFORE', record)
        ref, alts, genotype = VCFWriter.get_proper_alleles(positional_record, genotype)
        # print('AFTER', record)
        if len(alts) == 1:
            alts.append('.')
        rec_end = int(pos + len(ref) - 1)
        record_string = chromosome_name + "\t" + str(pos) + "\t" + str(rec_end) + "\t" + ref + "\t" + '\t'.join(alts) \
                        + "\t" + genotype + "\t" + str(qual) + "\t" + str(gq) + "\t" + "\n"
        record_file.write(record_string)


def merge_call_files(vcf_file_directory):
    filemanager_object = FileManager()
    # get all bed file paths from the directory
    file_paths = filemanager_object.get_file_paths_from_directory(vcf_file_directory)
    all_records = []
    for file_path in file_paths:
        with open(file_path, 'r') as tsv:
            for line in tsv:
                chr_name, pos_st, pos_end, ref, alt1, alt2, genotype, qual, gq = line.strip().split('\t')
                alts = []
                pos_st, pos_end, qual, gq = int(pos_st), int(pos_end), float(qual), float(gq)
                if alt1 != '.':
                    alts.append(alt1)
                if alt2 != '.':
                    alts.append(alt2)
                all_records.append((chr_name, pos_st, pos_end, ref, alts, genotype, qual, gq))

    filemanager_object.delete_files(file_paths)
    os.rmdir(vcf_file_directory)

    return all_records


def label_to_base(label):
    if label == 2:
        return 'A'
    if label == 3:
        return 'C'
    if label == 4:
        return 'T'
    if label == 5:
        return 'G'
    return ''

def polish_genome(csv_file, batch_size, num_workers, bam_file_path, output_dir, max_threads):
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + output_dir + "\n")
    predict(csv_file, batch_size, num_workers)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "COMPILING PREDICTIONS TO CALL VARIANTS.\n")

    fasta_file = open(output_dir + "helen_polished.fa", 'w')

    for chromosome in chromosome_list:
        pos_list = list(position_dict[chromosome])
        pos_list = sorted(list(pos_list), key=lambda element: (element[0], element[1]))
        sequence = ''
        for pos, index in pos_list:
            prediction = max(prediction_dict[(chromosome, pos, index)],
                             key=prediction_dict[(chromosome, pos, index)].count)
            sequence += label_to_base(prediction)

        if len(sequence) > 0:
            fasta_file.write('>'+chromosome+"\n")
            fasta_file.write(sequence+"\n")


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv_file",
        type=str,
        required=True,
        help="CSV file containing all image segments for prediction."
    )
    parser.add_argument(
        "--bam_file",
        type=str,
        required=True,
        help="Path to the BAM file."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default='vcf_output',
        help="Output directory."
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=8,
        help="Number of maximum threads for this region."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)
    polish_genome(FLAGS.csv_file,
                  FLAGS.batch_size,
                  FLAGS.num_workers,
                  FLAGS.bam_file,
                  FLAGS.output_dir,
                  FLAGS.max_threads)

