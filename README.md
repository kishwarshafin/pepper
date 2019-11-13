# P.E.P.P.E.R.
###### Program for Evaluating Patterns in the Pileups of Erroneous Reads

[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a deep neural network based polisher designed to work with Oxford Nanopore Sequencing technology. `PEPPER` uses a Recurrent Neural Network (RNN) based encoder-decoder model to call a consensus sequence from the summary statistics of each genomic position. The local realignment process enabled by [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) library enables the module to be standalone and not require any prior polishing with other tools.

© 2019 Kishwar Shafin, Trevor Pesout, Benedict Paten. <br/>
Computational Genomics Lab (CGL), University of California, Santa Cruz.

## Workflow
 * Sequence a genome and get a basecalled reads file (`reads.fastq`).
 * Use Shasta assembler to get an assembly from the basecalled data (`shasta_assembly.fa`).
 * Use [Minimap2](https://github.com/lh3/minimap2) to map `reads.fastq` to `shasta_assembly.fa` and get a bam file (`reads_2_shasta.bam`).
 * Use `1_pepper_make_images.py` to generate pileup summaries.
 * Use `2_pepper_call_consensus.py` to generate a consensus sequence.
 * Use `3_pepper_stitch.py` to stitch chunks and get the polished sequence.
<p align="center">
<img src="img/PEPPER_pipeline.svg" alt="pipeline.svg" height="640p">
</p>


## Installation
We recommend using `Linux` environment to run `PEPPER`.
```bash
sudo apt-get -y install cmake make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev \
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
```

#### Step 2: Inference
```bash
python3 2_pepper_call_consensus.py \
-i <path/to/output_image_dir/> \
-m <path/to/pepper/models/XXX.pkl> \
-b <256/512/1024> \
-w <number_of_workers> \
-o <path/to/output_polished_sequence/> \
-g
```

#### Step 3: Inference
```bash
python3 3_pepper_stitch.py \
-i <path/to/output_polished_sequence/pepper_predictions.hdf> \
-o <path/to/output_polished_sequence/> \
-t <number_of_threads>
```



## Acknowledgement
We are thankful to the developers of these packages: </br>
* [Medaka](https://github.com/nanoporetech/medaka)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [htslib & samtools](http://www.htslib.org/)
* [ssw library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
* [hdf5 python (h5py)](https://www.h5py.org/)
* [pytorch](https://pytorch.org/)


© 2019 Kishwar Shafin, Trevor Pesout, Benedict Paten.
