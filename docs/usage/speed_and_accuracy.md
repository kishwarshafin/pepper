# How to tune parameters for improved runtime.

We have tuned our pipeline to provide a balance between speed vs accuracy. However, accuracy gets priority over speed in most cases. We understand the need for speed in many rapid applications. Here we discuss and provide some guidance on how to tune parameters for your use-case that will provide you the best runtime.

Before we discuss parameters, we want to remind a few basic ways to improve runtime:
* **Use of GPU:** We can use an instance with NVIDIA-GPU where CUDA is installed to get accelerated inference time in our pipeline.
* **Using multiple instances**: Our pipeline takes a `-r` parameter where we can provide a contig/chrosome name. Instead of using one node, we can use multiple instances where we process contigs in parallel. We can confirm that contig-level parallelization does not affect performance and is the cleanest way to improve runtime. We have designed all of the components in our pipeline for maximum resource usage, so we will not see any periodic drop of utilization due to any bottleneck.

### Experimental setup
To assess how our parameters are changing runtime, we will use the following setup and data:
```
Instance type:                   n2-standard-80
Architecture:                    x86_64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
Address sizes:                   46 bits physical, 48 bits virtual
CPU(s):                          80
Model name:                      Intel(R) Xeon(R) CPU @ 2.80GHz
Stepping:                        7
CPU MHz:                         2800.214
BogoMIPS:                        5600.42
```
We will use the following dataset for the runtime:
```
Sample:                   HG003 (chr1)
Coverage:                 ~60x
Chemistry:                R9.4.1
Basecaller:               Guppy 5.0.7 "Sup"
```

### Download data
```bash
# You can change base to a directory of your choice
BASE="${HOME}/ont-speed-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_guppy_507_2_GRCh38_pass.60x.chr1.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_guppy_507_2_GRCh38_pass.60x.chr1.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr1.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr1.fa.fai

# Download truth VCFs
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

### Download docker
```bash
docker pull kishwars/pepper_deepvariant:r0.7
```

### 1. Baseline performance
First we run `PEPPER-Margin-DeepVariant` with basic command to derive the baseline performance in accuracy and runtime.
```bash
# Set up input variables
OUTPUT_PREFIX="HG003_ONT_30x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1"
OUTPUT_VCF="HG003_ONT_30x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1.vcf.gz"
# Set the number of CPUs to use
THREADS="80"

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:test-r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup
```

```bash
# Set up input data
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
```
