import sys
import os
import torch
import time
import torch.nn as nn
from torch.utils.data import DataLoader
import concurrent.futures
import multiprocessing
from torch.nn.parallel import DistributedDataParallel
from pepper_snp.modules.python.models.dataloader_predict import SequenceDataset
from datetime import datetime
from pepper_snp.modules.python.models.ModelHander import ModelHandler
from pepper_snp.modules.python.Options import ImageSizeOptions, TrainOptions
from pepper_snp.modules.python.DataStorePredict import DataStore
os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'


def predict(input_filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, theads_per_caller, device_id, rank):
    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)
    transducer_model.eval()
    transducer_model = transducer_model.eval()
    # create output file
    output_filename = output_filepath + "pepper_prediction_" + str(rank) + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDataset(input_filepath, file_chunks)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)
    torch.set_num_threads(theads_per_caller)

    torch.cuda.set_device(device_id)
    transducer_model.to(device_id)
    transducer_model.eval()
    transducer_model = DistributedDataParallel(transducer_model, device_ids=[device_id])

    batch_completed = 0
    total_batches = len(data_loader)
    with torch.no_grad():
        for contig, contig_start, contig_end, chunk_id, images, position, index, ref_seq in data_loader:
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            prediction_base_tensor = torch.zeros((images.size(0), images.size(1), ImageSizeOptions.TOTAL_LABELS))

            images = images.to(device_id)
            hidden = hidden.to(device_id)
            prediction_base_tensor = prediction_base_tensor.to(device_id)

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break
                chunk_start = i
                chunk_end = i + TrainOptions.TRAIN_WINDOW
                # chunk all the data
                image_chunk = images[:, chunk_start:chunk_end]

                # run inference
                output_base, hidden = transducer_model(image_chunk, hidden)

                # now calculate how much padding is on the top and bottom of this chunk so we can do a simple
                # add operation
                top_zeros = chunk_start
                bottom_zeros = ImageSizeOptions.SEQ_LENGTH - chunk_end

                # do softmax and get prediction
                # we run a softmax a padding to make the output tensor compatible for adding
                inference_layers = nn.Sequential(
                    nn.Softmax(dim=2),
                    nn.ZeroPad2d((0, 0, top_zeros, bottom_zeros))
                )
                inference_layers = inference_layers.to(device_id)

                # run the softmax and padding layers
                base_prediction = (inference_layers(output_base) * 10).type(torch.IntTensor).to(device_id)

                # now simply add the tensor to the global counter
                prediction_base_tensor = torch.add(prediction_base_tensor, base_prediction)

                del inference_layers
                torch.cuda.empty_cache()

            prediction_base_tensor = prediction_base_tensor.cpu().numpy().astype(int)

            for i in range(images.size(0)):
                prediction_data_file.write_prediction(contig[i], contig_start[i], contig_end[i], chunk_id[i],
                                                      position[i], index[i], prediction_base_tensor[i], ref_seq[i])
            batch_completed += 1

            if rank == 0 and batch_completed % 5 == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")


def predict_distributed_gpu(filepath, file_chunks, output_filepath, model_path, batch_size, total_callers, threads_per_caller, device_ids, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers to spawn
    :param threads_per_caller: Number of threads to use per caller
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # print("TOTAL CALLERS: ", total_callers)
    # print("DEVICE IDs: ", device_ids)
    # exit()
    # setup args
    start_time = time.time()
    predict_args = []
    for i in range(total_callers):
        predict_args.append([(filepath, file_chunks[i], output_filepath, model_path, batch_size, num_workers, threads_per_caller, device_ids[i], i)])

    multiprocessing.set_start_method('spawn')
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: STARTING PROCESS POOL\n")
    # for arg in predict_args:
    #     print(arg)
    # exit()
    with multiprocessing.Pool(processes=total_callers) as pool:
        pool.starmap(predict, predict_args)

    # with concurrent.futures.ProcessPoolExecutor(max_workers=total_callers, mp_context=multiprocessing.get_context('spawn')) as executor:
    #     futures = [executor.submit(predict, filepath, file_chunks[thread_id], output_filepath, model_path, batch_size, num_workers, threads_per_caller, device_ids[thread_id], thread_id)
    #                for thread_id in range(0, total_callers)]
    #
    #     for fut in concurrent.futures.as_completed(futures):
    #         if fut.exception() is None:
    #             # get the results
    #             thread_id = fut.result()
    #             sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREAD "
    #                              + str(thread_id) + " FINISHED SUCCESSFULLY.\n")
    #         else:
    #             sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
    #         fut._result = None  # python issue 27144

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PREDICTION\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec\n")
