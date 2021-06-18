import sys
import os
import torch
import time
import torch.onnx
import torch.nn as nn
import onnxruntime
from datetime import datetime
import concurrent.futures
from torch.utils.data import DataLoader
import numpy as np
from numpy import argmax

from pepper_variant.modules.python.models.dataloader_predict import SequenceDataset
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.Options import ImageSizeOptions, TrainOptions
from pepper_variant.modules.python.DataStorePredict import DataStore


def predict(input_filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, threads, thread_id):
    # session options
    sess_options = onnxruntime.SessionOptions()
    sess_options.intra_op_num_threads = threads
    sess_options.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
    sess_options.graph_optimization_level = onnxruntime.GraphOptimizationLevel.ORT_ENABLE_ALL

    ort_session = onnxruntime.InferenceSession(model_path + ".onnx", sess_options=sess_options)
    torch.set_num_threads(threads)

    # create output file
    output_filename = output_filepath + "pepper_prediction_" + str(thread_id) + ".hdf"
    # prediction_data_file = DataStore(output_filename, mode='w')

    print(thread_id, file_chunks)
    # data loader
    input_data = SequenceDataset(input_filepath, file_chunks)

    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             collate_fn=SequenceDataset.my_collate)

    batch_completed = 0
    total_batches = len(data_loader)
    with torch.no_grad():
        for contigs, positions, depths, candidates, candidate_frequencies, images in data_loader:
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            # run inference on onnx mode, which takes numpy inputs
            ort_inputs = {ort_session.get_inputs()[0].name: images.cpu().numpy(),
                          ort_session.get_inputs()[1].name: hidden.cpu().numpy(),
                          ort_session.get_inputs()[2].name: cell_state.cpu().numpy()}
            output_base = ort_session.run(None, ort_inputs)

            # prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_base)
            batch_completed += 1

            if thread_id == 0 and batch_completed % 5 == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()


def predict_pytorch(input_filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, threads):
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: SETTING THREADS TO: " + str(threads) + ".\n")
    torch.set_num_threads(threads)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: INTER OP-THREAD SET TO: " + str(torch.get_num_threads()) + ".\n")
    sys.stderr.flush()

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                    num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)
    transducer_model.eval()
    # create output file
    output_filename = output_filepath + "pepper_prediction" + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDataset(input_filepath)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             collate_fn=SequenceDataset.my_collate)

    transducer_model.eval()

    batch_completed = 0
    total_batches = len(data_loader)

    with torch.no_grad():
        for contigs, positions, depths, candidates, candidate_frequencies, images in data_loader:
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            # run inference
            output_base = transducer_model(images, hidden, cell_state, False)
            # output_base, output_type = transducer_model(images, hidden, cell_state, False)

            output_base = output_base.detach().cpu().numpy()
            # output_type = output_type.detach().cpu().numpy()

            prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_base)
            batch_completed += 1

            if batch_completed % 10 == 0 and batch_completed != 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()


def predict_distributed_cpu(filepath, file_chunks, output_filepath, model_path, batch_size, total_callers, threads_per_caller, num_workers):
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
    # predict_pytorch(filepath, file_chunks[0],  output_filepath, model_path, batch_size, num_workers, threads_per_caller)
    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                    num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)
    transducer_model.eval()

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADING TO ONNX\n")
    x = torch.zeros(1, ImageSizeOptions.CANDIDATE_WINDOW_SIZE + 1, ImageSizeOptions.IMAGE_HEIGHT)
    h = torch.zeros(1, 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
    c = torch.zeros(1, 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

    if not os.path.isfile(model_path + ".onnx"):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SAVING MODEL TO ONNX\n")
        torch.onnx.export(transducer_model, (x, h, c),
                          model_path + ".onnx",
                          training=False,
                          opset_version=10,
                          do_constant_folding=True,
                          input_names=['image', 'hidden', 'cell_state'],
                          output_names=['output_base'],
                          dynamic_axes={'image': {0: 'batch_size'},
                                        'hidden': {0: 'batch_size'},
                                        'cell_state': {0: 'batch_size'},
                                        'output_base': {0: 'batch_size'}})

    start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=total_callers) as executor:
        futures = [executor.submit(predict, filepath, file_chunks[thread_id], output_filepath, model_path, batch_size, num_workers, threads_per_caller, thread_id)
                   for thread_id in range(0, total_callers)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                thread_id = fut.result()
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREAD "
                                 + str(thread_id) + " FINISHED SUCCESSFULLY.\n")
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PREDICTION\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec\n")


