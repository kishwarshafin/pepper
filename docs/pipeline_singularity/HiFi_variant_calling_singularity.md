## PacBio-HiFi variant calling workflow [Using singularity]
PEPPER-Margin-DeepVariant is a haplotype-aware variant calling pipeline for long reads.

<img src="../../img/PMDV_variant_calling_HiFi.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="820p">

----
### PacBio-HiFi HG002 chr20 case-study
We evaluated this pipeline on `~35x` HG002 data.
```bash
Sample:     HG002
Coverage:   ~35x
Region:     chr20
Reference:  GRCh38_no_alt
```

#### Command-line instructions
##### Step 1: Install singularity [must be installed by root]
Please install [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps). This must be installed by the system admin.

Follow these installation [instructions](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps) to install Singularity 3.7, if you want to install a newer version please follow instructions from the [singulaity website](https://sylabs.io/).
```bash
# Install dependencies
sudo apt-get update && sudo apt-get install -y \
build-essential \
libssl-dev \
uuid-dev \
libgpgme11-dev \
squashfs-tools \
libseccomp-dev \
wget \
pkg-config \
git \
cryptsetup


# Install Go on linux: https://golang.org/doc/install
export VERSION=1.14.12 OS=linux ARCH=amd64 && \
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
rm go$VERSION.$OS-$ARCH.tar.gz

# Set environment variable
echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# Download and install singularity
export VERSION=3.7.0 && # adjust this as necessary \
wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
tar -xzf singularity-${VERSION}.tar.gz && \
cd singularity

# install sigularity
./mconfig && \
make -C builddir && \
sudo make -C builddir install  

# After installation is complete log out and log in
sigularity help
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
## Pull the docker image to sigularity, this is a 6.6GB download
singularity pull docker://kishwars/pepper_deepvariant:r0.4

# The pull command creates pepper_deepvariant_r0.4.sif file locally

# Run PEPPER-Margin-DeepVariant
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
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

singularity pull docker://google/deepvariant:1.1.0

# The pull command creates deepvariant_1.1.0.sif file locally

singularity exec --bind /usr/lib/locale/ \
deepvariant_1.1.0.sif \
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
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
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
singularity pull docker://jmcdani20/hap.py:v0.3.12

# The pull command creates hap.py_v0.3.12.sif file locally.

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
