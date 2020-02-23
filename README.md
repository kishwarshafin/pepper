# P.E.P.P.E.R.
###### Program for Evaluating Patterns in the Pileups of Erroneous Reads

[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`P.E.P.P.E.R.` is a deep neural network based polisher designed to work with Oxford Nanopore Sequencing technology. `P.E.P.P.E.R.` uses a Recurrent Neural Network (RNN) based encoder-decoder model to call a consensus sequence from the summary statistics of each genomic position. The local realignment process using [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) library enables the module to be standalone and not require any prior polishing with other tools.

© 2020 Kishwar Shafin, Trevor Pesout, Miten Jain, Benedict Paten. <br/>
Computational Genomics Lab (CGL), University of California, Santa Cruz.

## Workflow
 * Sequence a genome and get a basecalled reads file (`reads.fastq`).
 * Use an assembler to get an assembly from the basecalled data (`assembly.fa`).
 * Use [minimap2](https://github.com/lh3/minimap2) to map `reads.fastq` to `assembly.fa` and get a bam file (`reads_2_assembly.bam`).
 * Use `pepper polish` to polish a genome.
 <p align="center">
 <img src="img/PEPPER_pipeline.svg" alt="pipeline.svg" height="640p">
 </p>


## Installation
We recommend using `Linux` environment to run `PEPPER`.
```bash
sudo apt-get -y install cmake make git gcc g++ autoconf bzip2 lzma-dev zlib1g-dev \
libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev \
liblzma-dev libhdf5-dev
```
We also require `CMake>=3.11`. To install `CMake 3.14` on Ubuntu follow this:
```bash
# https://cmake.org/install/
wget https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3.tar.gz
tar -xvf cmake-3.14.3.tar.gz
cd cmake-3.14.3
./bootstrap
make
make install # You may require sudo for this. Or you can use the cmake executable from this directory
```

Now download and install `PEPPER`:
```bash
git clone https://github.com/kishwarshafin/pepper.git
cd pepper
./build.sh
```
`PEPPER` uses `pytorch` and few other python libraries. Please install those using:
```bash
python3 -m pip install h5py tqdm torchnet numpy pyyaml
```

#### Install pytorch
We have not included `torch` in the requirement as it needs to be carefully installed to enable `CUDA`. Please refer to the official [Pytorch Website](https://pytorch.org/) to install the correct version with or without `CUDA`. We highly recommend using version `1.3` or later for `pytorch`.

## Usage

#### Step 1: Generate Images
```bash
python3 1_pepper_make_images.py \
-b </path/to/reads_2_draft_assembly.bam> \
-d <path/to/draft_assembly.fasta> \
-o <path/to/output_image_dir/> \
-t <number_of_threads>

ARGUMENT DETAILS:
1_pepper_make_images.py generates summary statistics from the aligment of reads to the draft assembly.
  -h, --help            
                        show this help message and exit
  -b BAM, --bam BAM     
                        BAM file containing mapping between reads and the draft assembly.
  -d DRAFT, --draft DRAFT
                        FASTA file containing the draft assembly.
  -r REGION, --region REGION
                        Region in [chr_name:start-end] format
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory, if it does not exist it will  be created.
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 5.
```

#### Step 2: Inference
```bash
python3 2_pepper_call_consensus.py \
-i <path/to/output_image_dir/> \
-m <path/to/pepper/models/XXX.pkl> \
-b <batch_size> \
-w <number_of_workers> \
-o <path/to/output_polished_sequence/> \
-g

ARGUMENT DETAILS:
2_pepper_call_consensus.py performs inference on images using a trained model.

  -h, --help            
                        show this help message and exit
  -i IMAGE_DIR, --image_dir IMAGE_DIR
                        Path to directory containing all HDF5 images.
  -m MODEL_PATH, --model_path MODEL_PATH
                        Path to a trained model.
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Batch size for testing, default is 100. Suggested values: 256/512/1024.
  -w NUM_WORKERS, --num_workers NUM_WORKERS
                        Number of workers for loading images. Default is 4.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to the output directory.
  -g, --gpu_mode        
                        If set then PyTorch will use GPUs for inference [CUDA REQUIRED]
```

#### Step 3: Stitch
```bash
python3 3_pepper_stitch.py \
-i <path/to/output_polished_sequence/pepper_predictions.hdf> \
-o <path/to/output_polished_sequence/> \
-t <number_of_threads>

ARGUMENT DETAILS:
3_pepper_stitch.py performs the final stitching to generate the polished sequences.

  -i INPUT_HDF, --input_hdf INPUT_HDF
                        Input hdf prediction file.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory.
  -t THREADS, --threads THREADS
                        Number of threads.

```

## Results

#### PEPPER achieves lower error rate than ONT suggested pipeline.
We compared `PEPPER` against `Racon-Medaka` pipeline and we demonstrate significantly better results for microbial genomes. We used Staphylococcus Aureus samples to evaluate these two pipelines. The PEPPER microbial model was trained on samples excluding Staphylococcus Aureus. We used `r941_prom_high` model to run `Medaka`.
<p align="center">
<img src="img/PEPPER_error_rate.png" alt="PEPPER_error_rate.png" height="420">
</p>

#### New R10 chemistry shows further improvement in polishing results
The new `R10` data is now available for `MinION` and we polished the assembly generated with `R9` data using the `R10` reads. The R10 data provides significant improvement in overall quality of the genome.
<p align="center">
<img src="img/PEPPER_chemistry.png" alt="PEPPER_chemistry.png" height="420p">
</p>

## Acknowledgement
We are thankful to the developers of these packages: </br>
* [Medaka](https://github.com/nanoporetech/medaka)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [htslib & samtools](http://www.htslib.org/)
* [ssw library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
* [hdf5 python (h5py)](https://www.h5py.org/)
* [pytorch](https://pytorch.org/)

## Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="guppy235" width="240p"> <br/>

The name "P.E.P.P.E.R." is also inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616). PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.


© 2020 Kishwar Shafin, Trevor Pesout, Miten Jain, Benedict Paten.
