## Oxford Nanopore assembly polishing with PacBio HiFi reads [Using singularity]
PEPPER-Margin-DeepVariant can be used to polish nanopore-based assemblies in a diploid manner.

<img src="../../img/PMDV_polishing.png" alt="PEPPER-Margin-DeepVariant Polishing Workflow">

----

### HG002 chr20 Shasta assembly polishing case-study
Here we evaluate our pipeline on Shasta assembly and polish it with `~35x` PacBio HiFi data.
```bash
Sample:     HG002
Assembler:  Shasta
Data:       PacBio HiFi
Coverage:   ~35x
Region:     chr20
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
singularity help
```

##### Step 2: Download and prepare input data
```bash
BASE="${HOME}/hifi-polishing-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
ASM="HG002_Shasta_run1.chr20.fa"
BAM="HG002_HiFi_35x_2_Shasta_assembly.chr20.bam"
SAMPLE_NAME="HG002"
# Set the number of CPUs to use
THREADS="64"

# Set up output directory
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_HiFi_35x_2_Shasta_assembly.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_HiFi_35x_2_Shasta_assembly.chr20.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_Shasta_run1.chr20.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_Shasta_run1.chr20.fa.fai
```

##### Step 3: Run PEPPER-Margin-DeepVariant
```bash
## Pull the docker image to sigularity, this is a 6.6GB download
singularity pull docker://kishwars/pepper_deepvariant:r0.4

# The pull command creates pepper_deepvariant_r0.4.sif file locally

# Run PEPPER-Margin-DeepVariant
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
run_pepper_margin_deepvariant polish_assembly \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${ASM}" \
-o "${OUTPUT_DIR}" \
-t ${THREADS} \
-s ${SAMPLE_NAME} \
--ccs

# this generates 2 VCFs, one per haplotype
HAP1_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP1.vcf.gz
HAP2_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP2.vcf.gz

POLISHED_ASM_HAP1=HG002_Shasta_run1.PMDV.HiFi_polished.HAP1.fasta
POLISHED_ASM_HAP2=HG002_Shasta_run1.PMDV.HiFi_polished.HAP2.fasta

# Apply the VCF to the assembly
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
bcftools consensus \
-f "${INPUT_DIR}/${ASM}" \
-H 2 \
-s "${SAMPLE_NAME}" \
-o "${OUTPUT_DIR}/${POLISHED_ASM_HAP1}" \
"${OUTPUT_DIR}/${HAP1_VCF}"

singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
bcftools consensus \
-f "${INPUT_DIR}/${ASM}" \
-H 2 \
-s "${SAMPLE_NAME}" \
-o "${OUTPUT_DIR}/${POLISHED_ASM_HAP2}" \
"${OUTPUT_DIR}/${HAP2_VCF}"

# HG002_Shasta_run1.PMDV.HiFi_polished.HAP1.fasta and HG002_Shasta_run1.PMDV.HiFi_polished.HAP2.fasta are the polished assemblies.
```

### Evaluation of the polished assemblies using hap.py (Optional)
You can evaluate the variants using `hap.py`.

Download benchmarking data:
```bash
# Set up input data
REF="GRCh38_no_alt.chr20.fa"
TRUTH_VCF="HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
OUTPUT_DIR_ASM_TO_REF=${OUTPUT_DIR}/asm_to_ref_vcf_output/

# Create output directory
mkdir -p ${OUTPUT_DIR_ASM_TO_REF}

# Download truth VCFs
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa.fai
```

Apply the assembly back to the GRCh38 reference
```bash
# pull the docker
singularity pull docker://kishwars/asm_to_ref_vcf:latest

# apply the assembly to the reference using dipcall
singularity exec --bind /usr/lib/locale/ \
asm_to_ref_vcf_latest.sif \
bash /opt/asm_to_ref_vcf.sh \
-a "${INPUT_DIR}/${ASM}" \
-m "${OUTPUT_DIR}/${HAP1_VCF}" \
-p "${OUTPUT_DIR}/${HAP2_VCF}" \
-r "${INPUT_DIR}/${REF}" \
-s ${SAMPLE_NAME} \
-o "${OUTPUT_DIR_ASM_TO_REF}"

OUTPUT_REF_TO_VCF=ref_dipcall_output.vcf.gz

# Pull the docker image
singularity pull docker://jmcdani20/hap.py:v0.3.12

# Run hap.py
singularity exec --bind /usr/lib/locale/ \
hap.py_v0.3.12.sif \
/opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR_ASM_TO_REF}/${OUTPUT_REF_TO_VCF} \
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
| INDEL |      11256     |       10223       |        1033        |        1080        | 0.908227 |  0.906631 | 0.907428 |
|  SNP  |      71333     |       70058       |        1275        |         249        | 0.982126 |  0.996463 | 0.989243 |

### Authors:
This pipeline is developed in a collaboration between UCSC genomics institute and the genomics team at Google health.
