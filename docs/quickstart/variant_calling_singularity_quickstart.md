## Variant calling quickstart (Using singularity)
Here we show a quickstart for the Oxford Nanopore and PacBio HiFi variant calling.

### Install docker
Please see our [ONT](../pipeline_singularity/ONT_variant_calling_singularity.md) case-study documentation to see instructions on how to install docker.

## Quickstart: Nanopore variant calling
```bash
BASE="${HOME}/ont-quickstart"

# Set up input data
INPUT_DIR="${BASE}/input/data"
REF="GRCh38_no_alt.chr20.fa"
BAM="HG002_ONT_2_GRCh38.chr20.quickstart.bam"

# Set the number of CPUs to use
THREADS="1"

# Set up output directory
OUTPUT_DIR="${BASE}/output"
OUTPUT_PREFIX="HG002_ONT_2_GRCh38_PEPPER_Margin_DeepVariant.chr20"
OUTPUT_VCF="HG002_ONT_2_GRCh38_PEPPER_Margin_DeepVariant.chr20.vcf.gz"
TRUTH_VCF="HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz"
TRUTH_BED="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/HG002_ONT_2_GRCh38.chr20.quickstart.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/HG002_ONT_2_GRCh38.chr20.quickstart.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/GRCh38_no_alt.chr20.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/GRCh38_no_alt.chr20.fa.fai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/quickstart_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed

# Pull the docker images
singularity pull docker://jmcdani20/hap.py:v0.3.12
singularity pull docker://kishwars/pepper_deepvariant:r0.5

# Run PEPPER-Margin-DeepVariant
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.5.sif \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t ${THREADS} \
-r chr20:1000000-1020000 \
--ont

# Run hap.py
singularity exec --bind /usr/lib/locale/ \
hap.py_v0.3.12.sif \
/opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${OUTPUT_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}/happy.output" \
--pass-only \
-l chr20:1000000-1020000 \
--engine=vcfeval \
--threads="${THREADS}"
```

**Expected output:**

|  Type | Truth<br>total | True<br>positives | False<br>negatives | False<br>positives | Recall | Precision | F1-Score |
|:-----:|:--------------:|:-----------------:|:------------------:|:------------------:|:------:|:---------:|:--------:|
| INDEL |        2       |         2         |          0         |          1         |   1.0  |    0.6    |    0.8   |
|  SNP  |       39       |         39        |         39         |         39         |   1.0  |    1.0    |    1.0   |
### Authors:
This pipeline is developed in a collaboration between UCSC genomics institute and the genomics team at Google health.
