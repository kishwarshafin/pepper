import sys
import torch
import torch.onnx
import pickle
from datetime import datetime
from torch.utils.data import DataLoader
from pepper_variant.build import PEPPER_VARIANT
from pepper_variant.modules.python.models.dataloader import SequenceDatasetFake
from pepper_variant.modules.python.DataStorePredict import DataStore


def predict_pytorch_fake(input_filepath, output_filepath, batch_size, num_workers):
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: STARTING DRY PREDICTION." + ".\n")
    sys.stderr.flush()

    # create output file
    output_filename = output_filepath + "pepper_prediction_fake" + ".hdf"
    output_hdf_file = DataStore(output_filename, 'w')

    # data loader
    input_data = SequenceDatasetFake(input_filepath)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             collate_fn=SequenceDatasetFake.my_collate)

    batch_completed = 0
    total_batches = len(data_loader)

    with torch.no_grad():

        # for contig, depth, candidates, candidate_frequency, images, position, output_base, output_type in data_loader:
        for contigs, positions, depths, candidates, candidate_frequencies, images, base_predictions, type_predictions in data_loader:
            sys.stderr.flush()
            # all_candidate_predictions = []
            # for i in range(images.size(0)):
            #     predicted_candidate = PEPPER_VARIANT.CandidateImagePrediction(contigs[i],
            #                                                                   positions[i],
            #                                                                   depths[i],
            #                                                                   candidates[i],
            #                                                                   candidate_frequencies[i],
            #                                                                   base_predictions[i],
            #                                                                   type_predictions[i])
            #     all_candidate_predictions.append(predicted_candidate)
            output_hdf_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, base_predictions)
            # pickle.dump(all_candidate_predictions, prediction_data_file)

            batch_completed += 1
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
            sys.stderr.flush()
        # prediction_data_file.close()


def predict_distributed_cpu_fake(filepath, output_filepath, batch_size, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers to start
    :param threads_per_caller: Number of threads per caller.
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    predict_pytorch_fake(filepath,  output_filepath, batch_size, num_workers)


