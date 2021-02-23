## PacBio-HiFi variant calling workflow
PEPPER-Margin-DeepVariant is a haplotype-aware variant calling pipeline for long reads.

<img src="../img/PMDV_variant_calling_HiFi.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="820p">

----
### PacBio-HiFi HG002 chr20 case-study
We evaluated this pipeline on `~35x` HG002 data. The data is publicly available, please feel free to download, run and evaluate the pipeline.
```bash
Sample:     HG002
Coverage:   ~35x
Region:     chr20
Reference:  GRCh38_no_alt
```

#### Command-line instructions
##### Step 1: Install docker
Please install docker and wget if you don't have it installed already. You can install docker for other distros from here:
* [CentOS](https://docs.docker.com/engine/install/centos/) docker installation guide
* [Debian/Raspbian](https://docs.docker.com/engine/install/debian/) docker installation guide
* [Fedora](https://docs.docker.com/engine/install/fedora/) installation guide
* [Ubuntu](https://docs.docker.com/engine/install/ubuntu/) installation guide

We show the installation instructions for Ubuntu here:
```bash
# Install wget to download data files.
sudo apt-get -qq -y update
sudo apt-get -qq -y install wget

# Install docker using instructions on:
# https://docs.docker.com/install/linux/docker-ce/ubuntu/
sudo apt-get -qq -y install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

sudo add-apt-repository \
"deb [arch=amd64] https://download.docker.com/linux/ubuntu \
$(lsb_release -cs) \
stable"

sudo apt-get -qq -y update
sudo apt-get -qq -y install docker-ce
docker --version
```

##### Step 2: Download and prepare input data
```bash
BASE="${HOME}/hifi-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
REF="GRCh38_no_alt.chr20.fa"
BAM="HG002_PacBio_HiFi_35x_2_GRCh38_no_alt.chr20.bam"

# Set the number of CPUs to use
THREADS="64"

# Set up output directory
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_PacBio_HiFi_35x_2_GRCh38_no_alt.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_PacBio_HiFi_35x_2_GRCh38_no_alt.chr20.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa.fai
```

##### Step 3: Run PEPPER-Margin to generate a phased bam
```bash
## Pull the docker image.
sudo docker pull kishwars/pepper_deepvariant:r0.4

# Run PEPPER-Margin-DeepVariant
sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-t ${THREADS} \
--ccs

# This generates a Phased bam in the output directory: MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam
```
##### Step 4: Run DeepVariant
```bash
PHASED_BAM=MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam
OUTPUT_PREFIX="HG002_HiFi_35x_2_GRCh38_PEPPER_Margin_DeepVariant.chr20"
OUTPUT_VCF="HG002_HiFi_35x_2_GRCh38_PEPPER_Margin_DeepVariant.chr20.vcf.gz"

sudo docker pull google/deepvariant:1.1.0

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
google/deepvariant:1.1.0 \
/opt/deepvariant/bin/run_deepvariant \
--model_type=PACBIO \
--ref="${INPUT_DIR}/${REF}" \
--reads="${OUTPUT_DIR}/${PHASED_BAM}" \
--output_vcf="${OUTPUT_DIR}/${OUTPUT_VCF}" \
--num_shards=${THREADS} \
--use_hp_information
```

##### Step 5: Phase the output VCF with Margin (Optional)
```bash
sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
margin phase \
"${INPUT_DIR}/${BAM}" \
"${INPUT_DIR}/${REF}" \
"${OUTPUT_DIR}/${OUTPUT_VCF}" \
/opt/margin_dir/params/misc/allParams.phase_vcf.json \
-t ${THREADS} \
-M \
-o "${OUTPUT_DIR}/${OUTPUT_PREFIX}"

OUTPUT_PHASED_VCF=${OUTPUT_PREFIX}.phased.vcf
```

###### Evaluation using hap.py (Optional)
You can evaluate the variants using `hap.py`.
Download benchmarking data:
```bash
# Set up input data
TRUTH_VCF="HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# Download truth VCFs
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

Run hap.py:
```bash
# Pull the docker image
sudo docker pull jmcdani20/hap.py:v0.3.12

# Run hap.py
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${OUTPUT_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}/happy.output" \
--pass-only \
-l chr20 \
--engine=vcfeval \
--threads="${THREADS}"
```

**Expected output:**

|  Type | Truth<br>total | True<br>positives | False<br>negatives | False<br>positives |  Recall  | Precision | F1-Score |
|:-----:|:--------------:|:-----------------:|:------------------:|:------------------:|:--------:|:---------:|:--------:|
| INDEL |      11256     |       11175       |         81         |         94         | 0.992804 |  0.991953 | 0.992378 |
|  SNP  |      71333     |       71277       |         56         |          8         | 0.999215 |  0.999888 | 0.999551 |

### Authors:
This pipeline is developed in a collaboration between UCSC genomics institute and the genomics team at Google health.
