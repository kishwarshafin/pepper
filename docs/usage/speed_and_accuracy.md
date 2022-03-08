# How to tune parameters for improved runtime.

##### NOTE: This analysis has been demonstrated on r0.7. The overall concept of the analysis remains the same for r0.8.

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
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG003_GRCh38_1_22_v4.2.1_benchmark.bed
```

### Download docker
```bash
sudo docker pull kishwars/pepper_deepvariant:r0.7
sudo docker pull jmcdani20/hap.py:v0.3.12
```

### Baseline performance
First we run `PEPPER-Margin-DeepVariant` with basic command to derive the baseline performance in accuracy and runtime.
```bash
# Set up input variables
BAM="HG003_guppy_507_2_GRCh38_pass.60x.chr1.bam"
REF="GRCh38_no_alt.chr1.fa"
OUTPUT_PREFIX="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1"
OUTPUT_VCF="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1.vcf.gz"
PMDV_OUTPUT_DIR="${OUTPUT_DIR}"/BASELINE_PMDV
# Set the number of CPUs to use
THREADS="80"

time sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:test-r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${PMDV_OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup
```
To understand how we split the variant calling we can look at the log of `PEPPER HP` using the following commands:
```bash
cat "${OUTPUT_DIR}"/BASELINE_PMDV/logs/3_pepper_hp.log
# You will see the following lines:
# INFO: FINISHED PROCESSING, TOTAL CANDIDATES FOUND: 454884
# INFO: FINISHED PROCESSING, TOTAL VARIANTS IN PEPPER: 264691
# INFO: FINISHED PROCESSING, TOTAL VARIANTS SELECTED FOR RE-GENOTYPING: 190193
```
This means:
* We have found total 454884 candidates (`PEPPER_HP_VARIANT_FULL.vcf.gz`).
* Out of 454884 variants 264691 has QV above the set threshold and we will use those as-is (`PEPPER_HP_VARIANT_OUTPUT_PEPPER.vcf.gz`).
* Finally, we have 190193 candidates that need to be re-genotyped with DeepVariant (`PEPPER_HP_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz`).

After re-genotyping the variants with DeepVariant we merge the variants together.

```bash
# Set up input data
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark.bed"
QUERY_VCF=${PMDV_OUTPUT_DIR}/${OUTPUT_VCF}
HAPPY_OUTPUT_DIR=${OUTPUT_DIR}/happy_output
THREADS="80"

mkdir -p ${HAPPY_OUTPUT_DIR}

# Run hap.py
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${QUERY_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${HAPPY_OUTPUT_DIR}/happy.output.baseline" \
--pass-only \
-l chr1 \
--engine=vcfeval \
--threads="${THREADS}"
```

### Approach 1: Lowering QV cutoffs (affects precision)
There are two parameters that we can use to change the behavior of our pipeline so less variants get re-genotyped by DeepVariant:
```bash
--pepper_snp_q_cutoff
--pepper_indel_q_cutoff
```
The default values of these two parameters are: `--pepper_snp_q_cutoff 15, --pepper_indel_q_cutoff 10`. If a SNP variant has QV less than `--pepper_snp_q_cutoff` then the variant gets re-genotyped. If we set the values lower then smaller number of variants will be re-genotyped by DeepVariant. To improve runtime, we will now set both of these values to 10. QV10 is equivalent to 90% confidence.

**Note:** We do not recommend setting these two parameters lower than 10 for ONT.
```bash
# Set up input variables
BAM="HG003_guppy_507_2_GRCh38_pass.60x.chr1.bam"
REF="GRCh38_no_alt.chr1.fa"
OUTPUT_PREFIX="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1"
OUTPUT_VCF="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1.vcf.gz"
PMDV_OUTPUT_DIR="${OUTPUT_DIR}"/QV_CUTOFF_SPEED_UP_PMDV
# Set the number of CPUs to use
THREADS="80"

time sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:test-r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${PMDV_OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup \
--pepper_snp_q_cutoff 10 \
--pepper_indel_q_cutoff 10
```

```bash
cat "${OUTPUT_DIR}"/BASELINE_PMDV/logs/3_pepper_hp.log
# You will see the following lines:
# INFO: FINISHED PROCESSING, TOTAL CANDIDATES FOUND: 454884
# INFO: FINISHED PROCESSING, TOTAL VARIANTS IN PEPPER: 264691
# INFO: FINISHED PROCESSING, TOTAL VARIANTS SELECTED FOR RE-GENOTYPING: 190193

cat "${OUTPUT_DIR}"/QV_CUTOFF_SPEED_UP_PMDV/logs/3_pepper_hp.log
# INFO: FINISHED PROCESSING, TOTAL CANDIDATES FOUND: 454220
# INFO: FINISHED PROCESSING, TOTAL VARIANTS IN PEPPER: 307379
# INFO: FINISHED PROCESSING, TOTAL VARIANTS SELECTED FOR RE-GENOTYPING: 146841
```
With approach 1, we reduced the total number of variants to be re-genotyped by DeepVariant from 190193 to 146841.

```bash
# Set up input data
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark.bed"
QUERY_VCF=${PMDV_OUTPUT_DIR}/${OUTPUT_VCF}
THREADS="80"

# Run hap.py
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${QUERY_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${HAPPY_OUTPUT_DIR}/happy.output.approach1" \
--pass-only \
-l chr1 \
--engine=vcfeval \
--threads="${THREADS}"
```

### Approach 2: Increasing candidate predicted value threshold (affects sensitivity)
There are three parameters that we can use to limit the number of candidates the PEPPER picks for genotyping:
```bash
--pepper_snp_p_value
--pepper_insert_p_value
--pepper_delete_p_value
```
The neural network model of PEPPER predicts three classes `P = [hom (0/0), het (0/1), hom-alt (1/1)]`. For each variant type (SNP, Insert or Delete), if `max(P[hom(0/1)], P[hom(1/1)])` is higher than the set associated threshold then we pick that candidate to be genotyped by DeepVariant as a potential true variant. The default values of these parameters are: `--pepper_snp_p_value 0.1, --pepper_insert_p_value 0.25 --pepper_delete_p_value 0.25`. To improve runtime, we can set the SNP threshold to 0.2 and insert and delete to 0.3. This will affect the sensitivity as more variants will be excluded from the candidate set.
```bash
# Set up input variables
BAM="HG003_guppy_507_2_GRCh38_pass.60x.chr1.bam"
REF="GRCh38_no_alt.chr1.fa"
OUTPUT_PREFIX="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1"
OUTPUT_VCF="HG003_ONT_60x_2_GRCh38_PEPPER_Margin_DeepVariant.chr1.vcf.gz"
PMDV_OUTPUT_DIR="${OUTPUT_DIR}"/QV_CUTOFF_SPEED_UP_PMDV
# Set the number of CPUs to use
THREADS="80"

time sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:test-r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${PMDV_OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup \
--pepper_snp_p_value 0.2 \
--pepper_insert_p_value 0.3 \
--pepper_delete_p_value 0.3
```

```bash
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark.bed"
QUERY_VCF=${PMDV_OUTPUT_DIR}/${OUTPUT_VCF}
HAPPY_OUTPUT_DIR=${OUTPUT_DIR}/happy_output
THREADS="80"

mkdir -p ${HAPPY_OUTPUT_DIR}

# Run hap.py
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${QUERY_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${HAPPY_OUTPUT_DIR}/happy.output.approach2" \
--pass-only \
-l chr1 \
--engine=vcfeval \
--threads="${THREADS}"
```
