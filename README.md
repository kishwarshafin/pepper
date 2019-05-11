# P.E.P.P.E.R.
###### Program for Evaluating Patterns in the Pileups of Erroneous Reads

[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a deep neural network based polisher designed to work with Oxford Nanopore Sequencing technology. `PEPPER` uses a Recurrent Neural Network (RNN) based encoder-decoder model to call a consensus sequence from the summary statistics of each genomic position. We are developing a pileup based approach where we encode read attributes in multiple channels and use a hybrid CNN-RNN based network to call consensus. The `pytorch` backend support provides easy GPU distribution making easier scaling on GPU while providing competitive results to other deep neural network based approaches.

© 2019 Kishwar Shafin, Trevor Pesout, Benedict Paten. <br/>
Computational Genomics Lab (CGL), University of California, Santa Cruz.

Workflow
--------
 * Sequence a genome and get a basecalled reads file (`reads.fa` or `reads.fastq`).
 * Use Shasta assembler to get an assembly from the basecalled data (`shasta_assembly.fa`).
 * Use `minimap2` to map `reads.fastq` to `shasta_assembly.fa` and get a bam file (`reads_2_shasta.bam`).
 * Use `1_pepper_make_images.py` to generate pileup summaries.
 * Use `2_pepper_call_consensus.py` to generate a consensus sequence.
 * Use `3_pepper_stitch.py` to stitch chunks and get the polished sequence.


Installation
------------
The pre-release can only be installed from the source. We plan to release binaries with our first release.

As we use `htslib`, we require few libraries to be installed before installing `PEPPER`. On a Linux machine, run:
```bash
sudo apt-get install gcc \
zlib1g-dev libbz2-dev liblzma-dev \
libncurses5-dev libffi-dev libssl-dev \
make wget python3-all-dev
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
python3 -m pip -r install requirements.txt
```
We have not included `torch` in the requirement as it needs to be carefully installed to enable `CUDA`. Please refer to the official [Pytorch Website](https://pytorch.org/) to install the correct version with or without `CUDA`. We tested `PEPPER` on `pytorch 1.0` and `1.1`.

USAGE
-----

Write down usage here. Each script and the expected output.

HELP
-----
Email address of someone helpful.


## Acknowledgement
We are thankful to the developers of these packages: </br>
* [Medaka](https://github.com/nanoporetech/medaka)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [htslib & samtools](http://www.htslib.org/)
* [ssw library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
* [hdf5 python (h5py)](https://www.h5py.org/)
* [pytorch](https://pytorch.org/)


© 2019 Kishwar Shafin, Trevor Pesout, Benedict Paten.
