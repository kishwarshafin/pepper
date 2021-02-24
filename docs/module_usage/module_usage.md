# Use PEPPER-Margin locally

## Using docker
The PEPPER or Margin module can be used independently in the following way:
```bash
docker pull kishwars/pepper_deepvariant:r0.4

# PEPPER SNP
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.4 \
pepper_snp --help

# PEPPER SNP models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl
# /opt/pepper_models/PEPPER_SNP_HiFi_V4.pkl

# PEPPER HP
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.4 \
pepper_hp --help

# PEPPER HP models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_HP_R941_ONT_V4.pkl

# Margin phase
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.4 \
margin phase --help

# Margin models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.phase_vcf.json

# Margin haplotag
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.4 \
margin halotag --help

# Margin halotag models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.ont_haplotag.json
# /opt/margin_dir/params/misc/allParams.ccs_haplotag.json
```

## Using Singularity
The PEPPER or Margin module can be used independently in the following way:
```bash
# This is a 6.6GB download
singularity pull docker://kishwars/pepper_deepvariant:r0.4
# This will download pepper_deepvariant_r0.4.sif

# PEPPER SNP
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
pepper_snp --help

# PEPPER SNP models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl
# /opt/pepper_models/PEPPER_SNP_HiFi_V4.pkl

# PEPPER HP
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
pepper_hp --help

# PEPPER HP models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_HP_R941_ONT_V4.pkl

# Margin phase
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
margin phase --help

# Margin models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.phase_vcf.json

# Margin haplotag
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.4.sif \
margin halotag --help

# Margin halotag models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.ont_haplotag.json
# /opt/margin_dir/params/misc/allParams.ccs_haplotag.json
```
