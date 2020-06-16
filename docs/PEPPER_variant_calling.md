## Oxford Nanopore variant calling workflow (Limited support release)
In collaboration with the [DeepVariant](https://github.com/google/deepvariant) group, we are developing a haplotype-aware variant calling pipeline for Oxford Nanopore sequencing technology. We are releasing the pipeline with limited support until we finalize the pipeline.

<img src="../img/PEPPER_DeepVariant.png" alt="PEPPER DeepVariant workflow">


### Limited support statement
We are releasing this pipeline with limited support as we are actively developing the pipeline. This pipeline is developed outside the main repository of [DeepVariant](https://github.com/google/deepvariant). If you have issues with running this pipeline, please open a github issue in this repository (PEPPER). We will update this page as we evaluate the pipeline across all platforms.

### Expected input
We expect the input files meet the following requirements to produce results that we have reported.
```bash
Alignment file (BAM):
      Coverage:   ~50-80x
      Pore:       R9.4.1
      Basecaller: Guppy 3.6.0
      Platform:   PromethION or MinION
```

We have not yet tested the pipeline on different pores and devices. We will update this document as we validate the pipeline.

##### Disk-space requirements
As we use WhatsHap to phase the BAM file, please expect disk-space usage of ~4x the size of the BAM file.

### Benchmarking
You can run this pipeline on run-time efficient mode or high-accuracy mode (high-accuracy mode enabled by setting parameter `-x 1`). The difference between these two modes is the haplotype-aware model used for DeepVariant. The high-accuracy mode uses a better feature representation known as `rows` in DeepVariant that results in higher accuracy but slower run-time.

#### Expected run-time
We have tested the high-accuracy and run-time efficient mode of this pipeline on `~50x HG002 chr20` data.

Resources available:
```bash
CPUs:
  CPU(s):                          32
  CPU MHz:                         1979.549
  CPU max MHz:                     3400.0000
  CPU min MHz:                     2200.0000
GPU:
  GPU(s):                          1
  GPU Model:                       GeForce GTX 1080ti
  GPU Memory:                      11177MiB
```

|      Mode     | Platform | Sample | Region | Coverage | Run-time (wall-clock) |
|:-------------:|:--------:|:------:|:------:|:--------:|:---------------------:|
|    Run-time   |    CPU   |  HG002 |  chr20 |   ~50x   |                       |
|    Run-time   |    GPU   |  HG002 |  chr20 |   ~50x   |       71.03 mins      |
| High-accuracy |    CPU   |  HG002 |  chr20 |   ~50x   |                       |
| High-accuracy |    GPU   |  HG002 |  chr20 |   ~50x   |      164.71 mins      |

### How to run
We have combined all of the steps to run sequentially using one script. You can run the pipeline on your CPU-only or GPU machines.

#### Running on a CPU-only machine (Whole genome run)
Running the variant calling pipeline is as simple as:
On a **CPU**-based machine:
```bash
# For Run-time efficient mode:
docker run --ipc=host \
-v /input_dir:/input_dir \
-v /output_dir:/output_dir \
kishwars/pepper_deepvariant_cpu:latest \
/opt/run_pepper_deepvariant.sh \
-b </input/bam_file.bam> \
-f </input/reference_file.fasta> \
-o </output/output_dir/> \
-t <number_of_threads>

# For High-Accuracy mode:
docker run --ipc=host \
-v /input_dir:/input_dir \
-v /output_dir:/output_dir \
kishwars/pepper_deepvariant_cpu:latest \
/opt/run_pepper_deepvariant.sh \
-b </input/bam_file.bam> \
-f </input/reference_file.fasta> \
-o </output/output_dir/> \
-t <number_of_threads> \
-x 1
```

#### Running on a GPU machine (Whole genome run)
On a **GPU** based machine (make sure to install [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)):
```bash
# For Run-time efficient mode:
docker run --ipc=host \
-v /input_dir:/input_dir \
-v /output_dir:/output_dir \
--gpus all \
kishwars/pepper_deepvariant_gpu:latest \
/opt/run_pepper_deepvariant.sh \
-b </input/bam_file.bam> \
-f </input/reference_file.fasta> \
-o </output/output_dir/> \
-t <number_of_threads>

# For High-Accuracy mode:
docker run --ipc=host \
-v /input_dir:/input_dir \
-v /output_dir:/output_dir \
--gpus all \
kishwars/pepper_deepvariant_gpu:latest \
/opt/run_pepper_deepvariant.sh \
-b </input/bam_file.bam> \
-f </input/reference_file.fasta> \
-o </output/output_dir/> \
-t <number_of_threads> \
-x 1
```

## Quickstart
Here is an example on how to run the pipeline on a small example.

#### Download data
```bash
INPUT_DIR=${PWD}/pepper_deepvariant_input_data
mkdir -p ${INPUT_DIR}

DATA_DIR=https://storage.googleapis.com/kishwar-helen/ont-quickstart
wget -P ${INPUT_DIR} "${DATA_DIR}"/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam
wget -P ${INPUT_DIR} "${DATA_DIR}"/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam.bai
wget -P ${INPUT_DIR} "${DATA_DIR}"/GRCh38_no_alt_chr20.fa
wget -P ${INPUT_DIR} "${DATA_DIR}"/GRCh38_no_alt_chr20.fa.fai
wget -P ${INPUT_DIR} "${DATA_DIR}"/HG002_GRCh38_GIABv4.1.TRUTH.TEST.bed
wget -P ${INPUT_DIR} "${DATA_DIR}"/HG002_GRCh38_GIABv4.1.TRUTH.TEST.vcf.gz
wget -P ${INPUT_DIR} "${DATA_DIR}"/HG002_GRCh38_GIABv4.1.TRUTH.TEST.vcf.gz.tbi

OUTPUT_DIR=${PWD}/pepper_deepvariant_output_data
mkdir -p ${OUTPUT_DIR}
```

#### Run the pipeline
You can run the pipeline on a CPU-only or GPU machine.

##### CPU-only Machine:
```bash
# Runtime efficient mode
sudo docker run --ipc=host \
-v ${INPUT_DIR}:/input \
-v ${OUTPUT_DIR}:/output \
kishwars/pepper_deepvariant_cpu:latest \
/opt/run_pepper_deepvariant.sh \
-b /input/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam \
-f /input/GRCh38_no_alt_chr20.fa \
-o /output/ \
-t 1 \
-r chr20:1000000-1020000

# High-accuracy mode (Used for pFDA). USE:
sudo docker run --ipc=host \
-v ${INPUT_DIR}:/input \
-v ${OUTPUT_DIR}:/output \
kishwars/pepper_deepvariant_cpu:latest \
/opt/run_pepper_deepvariant.sh \
-b /input/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam \
-f /input/GRCh38_no_alt_chr20.fa \
-o /output/ \
-t 1 \
-r chr20:1000000-1020000 \
-x 1 # this enables high-accuracy version
```

##### GPU-based Machine:
On a GPU machine with CUDA installed, you can run the following commands to run the pipeline. Please make sure you have [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) installed.
```bash
# Runtime efficient mode
sudo docker run --ipc=host \
--gpus all \
-v ${INPUT_DIR}:/input \
-v ${OUTPUT_DIR}:/output \
kishwars/pepper_deepvariant_gpu:latest \
/opt/run_pepper_deepvariant.sh \
-b /input/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam \
-f /input/GRCh38_no_alt_chr20.fa \
-o /output/ \
-t 1 \
-r chr20:1000000-1020000

# High-accuracy mode (Used for pFDA). USE:
sudo docker run --ipc=host \
--gpus all \
-v ${INPUT_DIR}:/input \
-v ${OUTPUT_DIR}:/output \
kishwars/pepper_deepvariant_gpu:latest \
/opt/run_pepper_deepvariant.sh \
-b /input/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam \
-f /input/GRCh38_no_alt_chr20.fa \
-o /output/ \
-t 1 \
-r chr20:1000000-1020000 \
-x 1 # this enables high-accuracy version
```



#### Evaluate the variants:
You can use hap.py to evaluate the variants.
```bash
sudo docker run -it \
-v "${INPUT_DIR}":"/input" \
-v "${OUTPUT_DIR}:/output" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
/input/HG002_GRCh38_GIABv4.1.TRUTH.TEST.vcf.gz \
/output/PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
-f /input/HG002_GRCh38_GIABv4.1.TRUTH.TEST.bed \
-r /input/GRCh38_no_alt_chr20.fa \
-o /output/happy.output \
--pass-only \
--engine=vcfeval \
--threads=1 \
-l chr20:1000000-1020000
```

Expected output:

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
|-------|--------|-------------|----------|---------------|------------------|-----------------|
| INDEL | ALL    | 3           | 2        | 0.666667      | 1.0              | 0.8             |
| SNP   | ALL    | 45          | 45       | 1.0           | 1.0              | 1.0             |

### Authors:
This pipeline is a collaboration between UCSC genomics institute and the genomics team at Google health.

#### UCSC Genomics Institute:
* Kishwar Shafin
* Trevor Pesout
* Miten Jain
* Benedict Paten

#### Genomics team at Google Health:
* Pi-Chuan Chang
* Alexey Kolesnikov
* Maria Nattestad
* Gunjan Baid
* Sidharth Goel
* Howard Yang
* Andrew Carroll
