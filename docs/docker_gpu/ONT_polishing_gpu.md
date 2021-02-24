## Oxford Nanopore assembly polishing with Nanopore reads
PEPPER-Margin-DeepVariant can be used to polish nanopore-based assemblies in a diploid manner.

<img src="../../img/PMDV_polishing.png" alt="PEPPER-Margin-DeepVariant Polishing Workflow" width="820p">

----

### HG002 chr20 Shasta assembly polishing case-study
Here we evaluate our pipeline on Shasta assembly of `~50x` HG002 nanopore data.
```bash
Sample:     HG002
Assembler:  Shasta
Coverage:   ~50x
Region:     chr20
```

#### Command-line instructions
##### Preprocessing: Install CUDA [must be installed by root]
Install CUDA toolkit 11.0 from the [CUDA archive](https://developer.nvidia.com/cuda-toolkit-archive).

Here are the instructions to install CUDA 11.0 on Ubuntu 20.04LTS:
```bash
# Verify you have CUDA capable GPUs:
lspci | grep -i nvidia

# Verify Linux version
uname -m && cat /etc/*release
# Expected output: x86_64

sudo apt-get -qq -y update
sudo apt-get -qq -y install gcc wget make

# Install proper kernel headers: This is for ubuntu
# Details: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
sudo apt-get install linux-headers-$(uname -r)

wget https://developer.download.nvidia.com/compute/cuda/11.2.0/local_installers/cuda_11.2.0_460.27.04_linux.run
sudo sh cuda_11.2.0_460.27.04_linux.run
```
##### Step 1.1: Install docker
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

# To add the user to avoid running docker with sudo:
# Details: https://docs.docker.com/engine/install/linux-postinstall/

sudo groupadd docker
sudo usermod -aG docker $USER

# Log out and log back in so that your group membership is re-evaluated.

# After logging back in.
docker run hello-world

# If you can run docker without sudo then change the following commands accordingly.
```

##### Step 1.1: Install nvidia-docker
Install nvidia docker following these [instructions](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#getting-started).

```bash
distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
&& curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add - \
&& curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update

sudo apt-get install -y nvidia-docker2

sudo systemctl restart docker

sudo docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
# The output show show all your GPUs, if you enabled Docker for users then you should be able to run nvidia-docker without sudo
```

##### Step 2: Download and prepare input data
```bash
BASE="${HOME}/ont-polishing-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
ASM="HG002_Shasta_run1.chr20.fa"
BAM="HG002_ONT_50x_2_Shasta_assembly.chr20.bam"
SAMPLE_NAME="HG002"
# Set the number of CPUs to use
THREADS="64"

# Set up output directory
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_ONT_50x_2_Shasta_assembly.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_ONT_50x_2_Shasta_assembly.chr20.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_Shasta_run1.chr20.fa
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_Shasta_run1.chr20.fa.fai
```

##### Step 3: Run PEPPER-Margin-DeepVariant
```bash
## Pull the docker image.
sudo docker pull kishwars/pepper_deepvariant:r0.4

# Run PEPPER-Margin-DeepVariant
sudo docker run --ipc=host \
--gpus all \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant polish_assembly \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${ASM}" \
-o "${OUTPUT_DIR}" \
-t ${THREADS} \
-s ${SAMPLE_NAME} \
-g \
--ont

# this generates 2 VCFs, one per haplotype
HAP1_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP1.vcf.gz
HAP2_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP2.vcf.gz

POLISHED_ASM_HAP1=HG002_Shasta_run1.PMDV.HAP1.fasta
POLISHED_ASM_HAP2=HG002_Shasta_run1.PMDV.HAP2.fasta

# Apply the VCF to the assembly
sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools consensus \
-f "${INPUT_DIR}/${ASM}" \
-H 2 \
-s "${SAMPLE_NAME}" \
-o "${OUTPUT_DIR}/${POLISHED_ASM_HAP1}" \
"${OUTPUT_DIR}/${HAP1_VCF}"

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools consensus \
-f "${INPUT_DIR}/${ASM}" \
-H 2 \
-s "${SAMPLE_NAME}" \
-o "${OUTPUT_DIR}/${POLISHED_ASM_HAP2}" \
"${OUTPUT_DIR}/${HAP2_VCF}"

# HG002_Shasta_run1.PMDV.HAP1.fasta and HG002_Shasta_run1.PMDV.HAP2.fasta are the polished assemblies.
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
sudo docker pull kishwars/asm_to_ref_vcf:latest

# apply the assembly to the reference using dipcall
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/asm_to_ref_vcf:latest /opt/asm_to_ref_vcf.sh \
-a "${INPUT_DIR}/${ASM}" \
-m "${OUTPUT_DIR}/${HAP1_VCF}" \
-p "${OUTPUT_DIR}/${HAP2_VCF}" \
-r "${INPUT_DIR}/${REF}" \
-s ${SAMPLE_NAME} \
-o "${OUTPUT_DIR_ASM_TO_REF}"

OUTPUT_REF_TO_VCF=ref_dipcall_output.vcf.gz

# Pull the docker image
sudo docker pull jmcdani20/hap.py:v0.3.12

# Run hap.py
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
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
| INDEL |      11256     |        7034       |        4222        |        21706       | 0.624911 |  0.248537 | 0.355633 |
|  SNP  |      71334     |       69078       |        2256        |        2028        | 0.968374 |  0.971502 | 0.969935 |

### Authors:
This pipeline is developed in a collaboration between UCSC genomics institute and the genomics team at Google health.
