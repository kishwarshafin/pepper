# How to train PEPPER (Advanced- requires GPU)
##### Written by: Kishwar Shafin, Verified by: Jimin Park
##### We thank [Guillaume Holley](https://github.com/GuillaumeHolley) for the kind suggestions and contributions to the training pipeline.

In this walkthrough, we will see how to train `PEPPER` and replace default model with custom models. In this excercise we will train a model on `Guppy 4.2.2` data which is currently not supported by `PEPPER-Margin-DeepVariant` as the error-rate of the basecaller is too high.

## Training PEPPER
The first step is to train `PEPPER` that we use to find initial variants for phasing.

### Step 1: Install PEPPER locally
We recommend installing PEPPER locally for training PEPPER as it requires GPU for training.
```bash
git clone https://github.com/kishwarshafin/pepper.git
cd pepper
make install
. ./venv/bin/activate
```
This will activate pepper:
```bash
(pepper) bash-4.4$ pepper_variant --version
```

We will use the path where we installed `PEPPER` as `/path/to/PEPPER` for the rest of the workflow.

### Step 2: Download training data

For training `PEPPER` we need:
* A bam file where reads of a sample are aligned to a reference,
* Associated reference file to which the reads are aligned (FASTA).
* Truth VCF file ideally from Genome-In-a-Bottle for the sample.

In this case we are going to use 50x HG002 data basecalled with Guppy 4.2.2 basecaller.
```bash
# Cutomize this directory to where you want to store your training data
TRAIN_DIR=/data/PEPPER_TRAINING/
mkdir -p "${TRAIN_DIR}"

INPUT_DIR="${TRAIN_DIR}"/TRAIN_INPUTS/
OUTPUT_DIR="${TRAIN_DIR}"/TRAIN_OUTPUTS/

mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Download alignment files
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_guppy422_2_GRCh38_no_alt.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_guppy422_2_GRCh38_no_alt.bam.bai

# Download reference files
wget -P ${INPUT_DIR} https://storage.googleapis.com/kishwar-helen/variant_calling_data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wget -P ${INPUT_DIR} https://storage.googleapis.com/kishwar-helen/variant_calling_data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

# Download GIAB truth files
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
```

### Step 3: Downsample bam at different coverages
We want to train PEPPER for different coverages. We start by calculating the average coverage of the alignment file:
```bash
samtools depth -r chr20 ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.bam | awk '{sum+=$3} END { print "Average Coverage= ",sum/NR}'
# Expected OUTPUT:
# Average Coverage=  50.7078
```
Given that we have 50x coverage, we can train on `30x, 40x, 50x` coverage for the alignment file. We now downsample the bam to the expected coverages we want to train on:
```bash
THREADS=$(nproc)

for coverage in 30 40 50
do
  total_coverage=51
  downsample_fraction=0.$((coverage * 100 / total_coverage))
  echo "Coverage= ${coverage}, Downsample fraction = ${downsample_fraction}"

  samtools view -s $downsample_fraction -b -@${THREADS} ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.bam > ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.bam
  samtools index -@${THREADS} ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.bam
done
```

This will generate three alignment files. We will use these to generate training examples.
```bash
ls -lha ${INPUT_DIR}
# HG002_guppy422_2_GRCh38_no_alt.30x.bam
# HG002_guppy422_2_GRCh38_no_alt.40x.bam
# HG002_guppy422_2_GRCh38_no_alt.50x.bam
```

### Step 4: Train PEPPER

The first step of training PEPPER is to make training images. First we would activate PEPPER:
```bash
cd /path/to/PEPPER
. ./venv/bin/activate
pepper_variant_train make_train_images --help
```

#### Make images
Then we can start generating the training images:
```bash
TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT=${INPUT_DIR}/PEPPER_TRAIN_IMAGES
TEST_OUTPUT=${INPUT_DIR}/PEPPER_TEST_IMAGES

for coverage in 30 40 50
do
  BAM=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.bam

  echo ${BAM}

  time pepper_variant_train make_train_images \
  -b ${BAM} \
  -f ${REF} \
  -tv ${TRUTH_VCF} \
  -r chr1-19 \
  -rb ${TRUTH_BED} \
  -t ${THREADS} \
  -o ${TRAIN_OUTPUT} \
  -d 1.0 \
  -p 0.3 \
  --ont_r9_guppy5_sup

  time pepper_variant_train make_train_images \
  -b ${BAM} \
  -f ${REF} \
  -tv ${TRUTH_VCF} \
  -r chr20 \
  -rb ${TRUTH_BED} \
  -t ${THREADS} \
  -o ${TEST_OUTPUT} \
  -d 1.0 \
  -p 1.0 \
  --ont_r9_guppy5_sup   
done
```

Parameters:
* **-b** : Input BAM file.
* **-f** : Input reference fasta file.
* **-tv** : Input truth VCF file.
* **-r** : Regions (chr1-19 for training and chr20 for test).
* **-rb** : GIAB high-confidence regions bed file.
* **-o** : Output directory.
* **-d** : Downsample fraction.
* **-p** : Probability of a homozygous site being selected for training. This is to maintain class balance.
* **--ont_r9_guppy5_sup** : Use ONT presets for candidate finding.



#### Train model

Then we train a model for `PEPPER`:
```bash
MODEL_OUTPUT_DIR=${INPUT_DIR}/PEPPER_MODEL_OUTPUT

pepper_variant_train train_model \
-train ${TRAIN_OUTPUT} \
-test ${TEST_OUTPUT} \
-o ${MODEL_OUTPUT_DIR} \
-bs 128 \
--test_batch_size 512 \
-lr 0.001 \
-wd 0.00001 \
-s 10000 \
-e 1000 \
-w 0 \
-g
```

Parameters:
* **-train** : Path to train images directory.
* **-test** : Path to test images directory.
* **-o** : Output directory where models will be saved.
* **-bs** : Batch size during training.
* **--test_batch_size** : Batch size during test (does not affect training).
* **-lr** : Learning rate value.
* **-wd** : Weight decay value.
* **-s** : Step size.
* **-e** : Total epochs.
* **-w** : Number of workers.
* **-g** : Set to use GPUs.

The models will be saved in :
```bash
ls -lha ${MODEL_OUTPUT_DIR}
# ${MODEL_OUTPUT_DIR}/trained_models_********_******
```
In the `trained_models*` directory we will see the saved model files (`.pkl`) and a `stats` directory.
```bash
ls -lha ${MODEL_OUTPUT_DIR}/trained_models_********_******
# pkl files
# stats_********_******/
```
In the stats directory, we will find three files:
```bash
base_confusion_matrix.txt
test_loss.csv
train_loss.csv
```
We can use `test_loss.csv` to find which models has the lowest loss that we can evaluate. We would pick top 5 models with lowest loss values and then run the following evaluation pipeline and pick the model that performs best.


### Evaluating a trained model
First we run `PEPPER` with a model we select:
```bash
BAM=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.30x.bam
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# Select the model you want to test
MODEL=${MODEL_OUTPUT_DIR}/trained_models_****_****/PEPPER_VARIANT_STEP_****_checkpoint.pkl
EVAL_OUTPUT_DIR=$OUTPUT_DIR/pepper_output/pepper_model_****

pepper_variant call_variant \
-b $BAM \
-f $REF \
-m $MODEL \
-o $EVAL_OUTPUT_DIR \
-s HG003 \
-t $THREADS \
-w 0 \
-bs 512 \
-r chr20 \
--ont_r9_guppy5_sup
```

Then we benchmark the output of PEPPER using `hap.py`:
```bash
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

OUTPUT_VCF=${EVAL_OUTPUT_DIR}/PEPPER_VARIANT_FULL.vcf
HAPPY_OUTPUT_DIR=${OUTPUT_DIR}/happy_outputs
# Put the step number instead of **** so we can keep track of the performance of each model
HAPPY_OUTPUT_FILE=${HAPPY_OUTPUT_DIR}/HG002_30x_pepper_model_****

mkdir -p ${HAPPY_OUTPUT_DIR}

time docker run -it -v /data:/data -v /data2:/data2 \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${TRUTH_VCF} \
${OUTPUT_VCF} \
-f ${TRUTH_BED} \
-r ${REF} \
-o ${HAPPY_OUTPUT_FILE} \
-l chr20 \
--pass-only \
--no-roc \
--no-json \
--engine=vcfeval \
--threads=71
```


## Replace the default PEPPER model
You can use the following command to run `PEPPER-Margin-DeepVariant` with your trained model instead of the default model using `--pepper_model`:
```bash
BAM=HG002_guppy422_2_GRCh38_no_alt.30x.bam
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
PEPPER_OUTPUT_DIR="${OUTPUT_DIR}"/pepper_deepvariant_output

docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${PEPPER_OUTPUT_DIR}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup \
--pepper_model /PATH/TO/PEPPER_MODEL.pkl
```
