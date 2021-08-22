# Use PEPPER-Margin locally

## Using docker
The PEPPER or Margin module can be used independently in the following way:
```bash
docker pull kishwars/pepper_deepvariant:r0.5

# PEPPER
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.5 \
pepper_variant --help

# PEPPER models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_VARIANT_R941_ONT_V5.pkl

# Margin phase
docker run --ipc=host \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
kishwars/pepper_deepvariant:r0.5 \
margin phase --help

# Margin models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.phase_vcf.json

# Margin halotag models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.ont_haplotag.guppy507.json
# /opt/margin_dir/params/misc/allParams.ccs_haplotag.json
```

## Using Singularity
The PEPPER or Margin module can be used independently in the following way:
```bash
# This is a 6.6GB download
singularity pull docker://kishwars/pepper_deepvariant:r0.5
# This will download pepper_deepvariant_r0.5.sif

# PEPPER
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.5.sif \
pepper_variant --help

# PEPPER models can be invoked from inside the docker:
# /opt/pepper_models/PEPPER_VARIANT_R941_ONT_V5.pkl

# Margin phase
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.5.sif \
margin phase --help

# Margin models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.phase_vcf.json

# Margin halotag models can be invoked from inside the docker
# /opt/margin_dir/params/misc/allParams.ont_haplotag.guppy507.json
```
