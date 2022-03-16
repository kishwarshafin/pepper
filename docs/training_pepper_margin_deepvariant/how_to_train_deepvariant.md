# How to train DeepVariant (with PEPPER) (Advanced- Requires GPU)
##### Prepared by: Kishwar Shafin, Verified by: Jimin Park
##### We thank [Guillaume Holley](https://github.com/GuillaumeHolley) for the kind suggestions and contributions to the training pipeline.

In this walkthrough, we will see how to train `DeepVariant` and replace default model with custom models. In this exercise, we will train a model on `Guppy 4.2.2` data which is currently not supported by `PEPPER-Margin-DeepVariant` as the error-rate of the basecaller is too high.

**NOTE:** This tutorial goes through the basic steps required for training DeepVariant. However, there are many more options available to ease the training process of DeepVariant. [This documentation](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-training-case-study.md) has a much more detailed explanation of the training process. If you have followed the official tutorial then you can use the `make_examples` stage from this tutorial and the rest of the steps should be exactly the same.

## Setup
In this setup we train one model with `--alt_aligned_pileup "diff_channels"` for calling SNPs and INDELs. However, you can train a `rows` model for INDEL calling and `none` model for SNP calling. 

To achieve two model setup, please repeat the DeepVariant training step with `--alt_aligned_pileup "none"` and `--alt_aligned_pileup "rows"` respectively.
Once you have the models, you can then use `--dv_model_snp <model_path>` and `--dv_model_indel <model_path>` to specify the model path for SNP and INDEL
independently. You can also select the alt alignment by `--dv_alt_aligned_pileup_snp` or `--dv_alt_aligned_pileup_indel` where `alt_aligned_pileup` can be `diff_channels`, `none` or `rows`.

## Training DeepVariant
We will now train `DeepVariant` that we use to genotype the candidates proposed by `PEPPER`.

### Step 1: Download training data

For training `DeepVariant` we need:
* A bam file where reads of a sample are aligned to a reference,
* Associated reference file to which the reads are aligned (FASTA).
* Truth VCF file ideally from Genome-In-a-Bottle for the sample.
* Annotated candidates from PEPPER.

In this case we are going to use 50x HG002 data basecalled with Guppy 4.2.2 basecaller.
```bash
# Cutomize this directory to where you want to store your training data
TRAIN_DIR=/data/PEPPER_TRAINING
mkdir -p "${TRAIN_DIR}"

INPUT_DIR="${TRAIN_DIR}"/TRAIN_INPUTS
OUTPUT_DIR="${TRAIN_DIR}"/TRAIN_OUTPUTS

mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Download alignment files
# If you have gsutil installed, you can use the following command for faster download of the bam file:
# gsutil -o GSUtil:parallel_thread_count=1 -o GSUtil:sliced_object_download_max_components=8 cp gs://pepper-deepvariant-public/pepper_deepvariant_training/HG002_guppy422_2_GRCh38_no_alt.bam ${INPUT_DIR}
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

### Step 2: Downsample bam at different coverages
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
### Step 3: Generate haplotagged BAM files and PEPPER candidates
Training DeepVariant requires haplotagged alignment files and candidates from PEPPER. The following command can be used to generate the haplotagged bam file and candidates using `PEPPER-Margin`:

```bash
for coverage in 30 40 50
do
  BAM=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.bam
  REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  PEPPER_OUTPUT_DIR=${OUTPUT_DIR}/PEPPER_SNP_HAPLOTAG_${coverage}x_output

  docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -u `id -u`:`id -g` \
  kishwars/pepper_deepvariant:r0.8 \
  run_pepper_margin_deepvariant call_variant \
  -b $BAM \
  -f $REF \
  -o $PEPPER_OUTPUT_DIR \
  -t $THREADS \
  -s HG002 \
  --ont_r9_guppy5_sup \
  --only_haplotag \
  -k
  # If you trained a custom PEPPER model, then use the following parameter
  # to use the custom model to generate the haplotagged bams:
  # --pepper_model </path/to/PEPPER_MODEL.pkl>


  mv ${OUTPUT_DIR}/PEPPER_SNP_HAPLOTAG_${coverage}x_output/intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.haplotagged.bam
  mv ${OUTPUT_DIR}/PEPPER_SNP_HAPLOTAG_${coverage}x_output/intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam.bai ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.haplotagged.bam.bai
  mv ${OUTPUT_DIR}/PEPPER_SNP_HAPLOTAG_${coverage}x_output/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.vcf.gz
  mv ${OUTPUT_DIR}/PEPPER_SNP_HAPLOTAG_${coverage}x_output/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz.tbi ${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.vcf.gz.tbi
done
```

This will generate three alignment files. We will use these to generate training examples.
```bash
ls -lha ${INPUT_DIR}
# Haplotagged bam files:
# HG002_guppy422_2_GRCh38_no_alt.30x.haplotagged.bam
# HG002_guppy422_2_GRCh38_no_alt.40x.haplotagged.bam
# HG002_guppy422_2_GRCh38_no_alt.50x.haplotagged.bam

# PEPPER candidate files
# HG002_guppy422_2_GRCh38_no_alt.30x.candidates.vcf.gz
# HG002_guppy422_2_GRCh38_no_alt.40x.candidates.vcf.gz
# HG002_guppy422_2_GRCh38_no_alt.50x.candidates.vcf.gz
```
### Step 4: Annotate candidate variants with truth
```bash
for coverage in 30 40 50
do
  CANDIDATE_VCF=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.vcf.gz
  TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
  TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  OUTPUT_TRUTH_VCF=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.vcf.gz
  OUTPUT_TRUTH_INTERSECTED=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.INTERSECTED.vcf
  OUTPUT_TRUTH_INTERSECTED_GZ=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.INTERSECTED.vcf.gz

  docker run -it \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -u `id -u`:`id -g` \
  kishwars/pepper_annotation:v0.1 \
  python3 /opt/label_candidates_positional.py \
  --truth_vcf $TRUTH_VCF \
  --candidate_vcf $CANDIDATE_VCF \
  --output_vcf $OUTPUT_TRUTH_VCF

  docker run -it \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -u `id -u`:`id -g` \
  kishwars/pepper_annotation:v0.1 \
  bedtools intersect \
  -header \
  -a $OUTPUT_TRUTH_VCF \
  -b $TRUTH_BED > "$OUTPUT_TRUTH_INTERSECTED"; \
  bgzip ${OUTPUT_TRUTH_INTERSECTED}; \
  tabix -p vcf ${OUTPUT_TRUTH_INTERSECTED_GZ}
done
```
#### Optional: Look at annotation accuracy and class distribution
As we annotate our candidates with the truth, the precision of the truth-annotated candidates are expected to be high. Also, this gives us a chance to look at the recall or sensitivity of our candidate finding approach. We can also look at the class distribution of each candidate set.
```bash
HAPPY_OUTPUT_DIR=${OUTPUT_DIR}/happy_outputs

mkdir -p ${HAPPY_OUTPUT_DIR}

for coverage in 30 40 50
do
  TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
  TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  OUTPUT_TRUTH_INTERSECTED_GZ=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.INTERSECTED.vcf.gz

  HAPPY_OUTPUT_FILE=${HAPPY_OUTPUT_DIR}/HG002_${i}x_truth_candidates

  time docker run -it -v /data:/data -v /data2:/data2 \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  ${TRUTH_VCF} \
  ${OUTPUT_TRUTH_INTERSECTED_GZ} \
  -f ${TRUTH_BED} \
  -r ${REF} \
  -o ${HAPPY_OUTPUT_FILE} \
  --pass-only \
  --no-roc \
  --no-json \
  --engine=vcfeval \
  --threads=71
done
```

To look at the class distribution:
```bash
for coverage in 30 40 50
do
  OUTPUT_TRUTH_INTERSECTED_GZ=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.INTERSECTED.vcf.gz
  echo "CLASS DISTRIBUTION: " ${OUTPUT_TRUTH_INTERSECTED_GZ}
  docker run -it \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -u `id -u`:`id -g` \
  kishwars/pepper_annotation:v0.1 \
  python3 /opt/vcf_print_GT_stat.py \
  -v ${OUTPUT_TRUTH_INTERSECTED_GZ}
done
```

### Step 5: DeepVariant training

#### Make examples
Then we can start generating the training images:
```bash
TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT_DIR=${INPUT_DIR}/deepvariant_train_examples
TEST_OUTPUT_DIR=${INPUT_DIR}/deepvariant_test_examples

mkdir -p ${TRAIN_OUTPUT_DIR}
mkdir -p ${TEST_OUTPUT_DIR}

for coverage in 30 40 50
do
  BAM=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.bam
  PEPPER_TRUTH_CANDIDATES=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.${coverage}x.candidates.TRUTH.INTERSECTED.vcf.gz

  # Train samples from chr1-19
  (time seq 0 $((THREADS-1)) | \
  parallel --halt 2 --line-buffer \
      docker run \
      -v "${INPUT_DIR}":"${INPUT_DIR}" \
      -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
      -u `id -u`:`id -g` \
      kishwars/pepper_deepvariant:r0.8 \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref $REF \
      --reads $BAM \
      --examples "${TRAIN_OUTPUT_DIR}/training_set.with_label.${coverage}x.tfrecord@${THREADS}.gz" \
      --alt_aligned_pileup "diff_channels" \
      --norealign_reads \
      --min_mapping_quality 5 \
      --min_base_quality 1 \
      --add_hp_channel \
      --parse_sam_aux_fields \
      --sort_by_haplotypes \
      --variant_caller "vcf_candidate_importer" \
      --truth_variants ${PEPPER_TRUTH_CANDIDATES} \
      --task {} \
      --exclude_regions "'chr20 chr21 chr22'" )

  # Test samples from chr20
  (time seq 0 $((THREADS-1)) | \
  parallel --halt 2 --line-buffer \
      docker run \
      -v "${INPUT_DIR}":"${INPUT_DIR}" \
      -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
      -u `id -u`:`id -g` \
      kishwars/pepper_deepvariant:r0.8 \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref $REF \
      --reads $BAM \
      --examples "${TEST_OUTPUT_DIR}/test_set.with_label.${coverage}x.tfrecord@${THREADS}.gz" \
      --alt_aligned_pileup "diff_channels" \
      --norealign_reads \
      --min_mapping_quality 5 \
      --min_base_quality 1 \
      --add_hp_channel \
      --parse_sam_aux_fields \
      --sort_by_haplotypes \
      --variant_caller "vcf_candidate_importer" \
      --truth_variants ${PEPPER_TRUTH_CANDIDATES} \
      --task {} \
      --regions "'chr20'" )
done
```

#### Shuffle the examples

**PLEASE READ:** If you have a memory limited machine then a better way to shuffle the examples without running out of memory has been implemented by [Guillaume Holley](https://github.com/GuillaumeHolley) in this repository: https://github.com/GuillaumeHolley/TFrecordShuffler.



Next step is to shuffle the models. This is an important step in training deepvariant. There are several ways to do it as explained in the [documentation](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-training-case-study.md). However, we will use the local version for simplicity. Please choose one that works best for your use-case:

```bash
SHUFFLE_SCRIPT_DIR="${OUTPUT_DIR}/shuffle_script"
mkdir -p ${SHUFFLE_SCRIPT_DIR}
wget https://raw.githubusercontent.com/google/deepvariant/r1.3/tools/shuffle_tfrecords_beam.py -O ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py

# Depending on your setup, you may need to install the following modules, we suggest using virtual environment for this
virtualenv venv --python=python3 --prompt "(shuffle-tf) "
. ./venv/bin/activate
pip3 install setuptools --upgrade
pip3 install apache_beam[gcp]
pip3 install tensorflow
```
Shuffle training and test set:
##### Note on input_pattern_list parameter (Please read)
You do not need to replace ?? in `input_pattern_list` parameter. The `"?"` values in this parameter is used as a wildcard.

```bash
# Shuffle training dataset
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
--input_pattern_list="${TRAIN_OUTPUT_DIR}"/training_set.with_label.??x.tfrecord-?????-of-000${THREADS}.gz \
--output_pattern_prefix="${TRAIN_OUTPUT_DIR}/training_set.with_label.shuffled" \
--output_dataset_name="HG002" \
--output_dataset_config_pbtxt="${OUTPUT_DIR}/training_set.dataset_config.pbtxt" \
--job_name=shuffle-tfrecords \
--runner=DirectRunner \
--direct_num_workers=0

# Shuffle test dataset
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
--input_pattern_list="${TEST_OUTPUT_DIR}"/test_set.with_label.??x.tfrecord-?????-of-000${THREADS}.gz \
--output_pattern_prefix="${TEST_OUTPUT_DIR}/test_set.with_label.shuffled" \
--output_dataset_name="HG002" \
--output_dataset_config_pbtxt="${OUTPUT_DIR}/test_set.dataset_config.pbtxt" \
--job_name=shuffle-tfrecords \
--runner=DirectRunner \
--direct_num_workers=0
```

Now we can look at the statistics of the training data:
```bash
cat $OUTPUT_DIR/test_set.dataset_config.pbtxt
cat $OUTPUT_DIR/training_set.dataset_config.pbtxt
```

#### Train model
Now we are ready to start training a model. It is highly recommended to use `--start_from_checkpoint` option as that reduces the total consumed time for convergence. We can first copy the DeepVariant ONT model provided within the docker:
```bash
ONT_DV_MODEL_DIR=${INPUT_DIR}/ont_dv_model

docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-u `id -u`:`id -g` \
kishwars/pepper_deepvariant:r0.8 \
cp -r /opt/dv_models/ont_deepvariant_vc/ ${ONT_DV_MODEL_DIR}
```
Now we can see the models in the directory:
```bash
ls ${ONT_DV_MODEL_DIR}
# dv_ont_r9_guppy5_sup_vc_model.data-00000-of-00001
# dv_ont_r9_guppy5_sup_vc_model.index
# dv_ont_r9_guppy5_sup_vc_model.meta
```
Next we fetch the gpu docker of PEPPER-DeepVariant:
```bash
docker pull kishwars/pepper_deepvariant:r0.8-gpu
```
And then we can start training and evaluating concurrently:
```bash
DEEPVARIANT_MODEL_OUTPUT_DIR=${OUTPUT_DIR}/deepvariant_models
mkdir -p ${DEEPVARIANT_MODEL_OUTPUT_DIR}

( time docker run --ipc=host \
  --gpus all \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -u `id -u`:`id -g` \
  kishwars/pepper_deepvariant:r0.8-gpu \
  /opt/deepvariant/bin/model_train \
  --dataset_config_pbtxt="${OUTPUT_DIR}/training_set.dataset_config.pbtxt" \
  --train_dir="${DEEPVARIANT_MODEL_OUTPUT_DIR}" \
  --model_name="inception_v3" \
  --number_of_steps=50000 \
  --save_interval_secs=300 \
  --batch_size=32 \
  --learning_rate=0.0005 \
  --start_from_checkpoint="${ONT_DV_MODEL_DIR}/dv_ont_r9_guppy5_sup_vc_model" \
) > "${OUTPUT_DIR}/deepvariant_train.log" 2>&1 &
```
As we launch the training in background we won't see any errors if there are any issues with the run.

Please keep an eye on the output log for a bit to make sure that training has started and then move to the next step.

We highly recommend using `top` or `htop` to see if the processes launched successfully and are running. You can also inspect the log by running:
```bash
# This will print the training log every 1.0 sec
watch -n 1 tail -n 20 "${OUTPUT_DIR}/deepvariant_train.log"

# Once you are sure that training is going smoothly, proceed to the next step.
# You can cancel this watch by running ctrl+c
```

Now we can also launch the evaluation process.
```bash
# We can use the pepper-deepvariant docker for evaluation
docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-u `id -u`:`id -g` \
kishwars/pepper_deepvariant:r0.8 \
/opt/deepvariant/bin/model_eval \
--dataset_config_pbtxt="${OUTPUT_DIR}/test_set.dataset_config.pbtxt" \
--checkpoint_dir="${DEEPVARIANT_MODEL_OUTPUT_DIR}" \
--batch_size=512 > "${OUTPUT_DIR}/deepvariant_eval.log" 2>&1 &
```

After the model starts generating checkpoints, you will see them appear in the output directory:
```bash
ls ${DEEPVARIANT_MODEL_OUTPUT_DIR}
# model.ckpt-****
```
The evaluation process starts generating the best checkpoint it has observed so far:
```bash
cat ${DEEPVARIANT_MODEL_OUTPUT_DIR}/best_checkpoint.txt
# /path/to/model.ckpt-***
```
Once qw want think the training has converged and want to test the best checkpoint, we can move to the next step.
### Evaluating a trained deepvariant model
We can use the `PEPPER-Margin-DeepVariant` pipeline and replace the default DeepVariant model with the model we trained:
```bash
BAM=${INPUT_DIR}/HG002_guppy422_2_GRCh38_no_alt.30x.haplotagged.bam
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# Select the model you want to test
# Put the step number instead of **** so we can keep track of the performance of each model
DEEPVARIANT_MODEL=${DEEPVARIANT_MODEL_OUTPUT_DIR}/model.ckpt-****
EVAL_OUTPUT_DIR=$OUTPUT_DIR/pepper_margin_deepvariant_eval

docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b "${BAM}" \
-f "${REF}" \
-o "${EVAL_OUTPUT_DIR}" \
-t "${THREADS}" \
-r "chr20" \
--ont_r9_guppy5_sup \
--dv_model "${DEEPVARIANT_MODEL}"

# Use the following parameters if you also trained different PEPPER models:
# --pepper_model /PATH/TO/PEPPER_MODEL.pkl
```

Then we benchmark the output of PEPPER using `hap.py`:
```bash
REF=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
TRUTH_BED=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed
TRUTH_VCF=${INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

OUTPUT_VCF=${EVAL_OUTPUT_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz
HAPPY_OUTPUT_DIR=${OUTPUT_DIR}/happy_outputs
HAPPY_OUTPUT_FILE=${HAPPY_OUTPUT_DIR}/HG002_30x_pepper_margin_deepvariant

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

This concludes our training tutorial. Happy training!
