## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the [PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline. This pipeline enables nanopore-based variant calling with [DeepVariant](https://github.com/google/deepvariant).

<p align="center">
<img src="./img/PMDV_variant_calling_ONT_v5.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="720p"></img>
</p>

---
### Version r0.7 update

The new r0.7 variant calling results report:
* Detailed performance evaluation on ONT and PacBio-HiFi data.
* Included training documentation for PEPPER-Margin-DeepVariant.
* Detailed explanation of method update.
* Examples on how to tune parameters to balance speed and accuracy.

### Variant calling performance

Detailed variant calling performance evaluation:
* [ONT R9.4.1 Guppy 5 "Sup" HG003 whole genome (PEPPER vs Clair3).](./docs/performance_evaluation/Oxford_nanopore_r9_whole_genome.md)
* [ONT R9.4.1 Guppy 5 "Sup" HG002-HG006 chr20 (PEPPER vs Clair3).](./docs/performance_evaluation/Oxford_nanopore_r9_whole_genome.md)
* [PacBio HiFi (PEPPER-Margin-DeepVariant vs DeepVariant-WhatsHap-DeepVariant).](./docs/performance_evaluation)

---

### Useful links to documentations
* [Description of PEPPER-Margin-DeepVariant method.](./docs/misc/pepper_v0.7_method_update.md)
* How to train PEPPER-DeepVariant:
    * [How to train PEPPER-SNP.](./docs/training_pepper_margin_deepvariant/how_to_train_pepper_snp.md)
    * [How to train PEPPER-HP.](./docs/training_pepper_margin_deepvariant/how_to_train_pepper_hp.md)
    * [How to train DeepVariant (using candidates from PEPPER).](./docs/training_pepper_margin_deepvariant/how_to_train_deepvariant.md)
* [How to install PEPPER locally.](./docs/local_install/install_pepper_locally.md)
* [List of parameters and description.](./docs/usage/usage_and_parameters.md)
* [How to tune parameters to balance speed and accuracy.](./docs/usage/speed_and_accuracy.md)


### How to cite
Please cite the following manuscript if you are using `PEPPER-Margin-DeepVariant`:


<details>
<summary><a href="https://www.nature.com/articles/s41592-021-01299-w"><b>Nature Methods:</b> Haplotype-aware variant calling enables high accuracy in nanopore long-reads using deep neural networks.</a></summary>
Authors: Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov, Sidharth Goel, <br/> Gunjan Baid, Mikhail Kolmogorov, Jordan M. Eizenga, Karen H. Miga, Paolo Carnevali, Miten Jain, Andrew Carroll & Benedict Paten.
</details>

---
### How to run
PEPPER-Margin-DeepVariant can be run using **Docker** or **Singularity**. A simple docker command looks like:
```bash
sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup

# --ont_r9_guppy5_sup is preset for ONT R9.4.1 Guppy 5 "Sup" basecaller
# for ONT R10.4 Q20 reads: --ont_r10_q20
# for PacBio-HiFi reads: --hifi
```

### Case studies

The variant calling pipeline can be run on [Docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps). The case studies are designed on `chr20` of `HG002` sample for ONT and `HG003` for PacBio-HiFi.

#### Oxford Nanopore Variant calling
The case-studies include input data and benchmarking of the run:
* Nanopore variant calling using **Docker**: [Link](./docs/pipeline_docker/ONT_variant_calling.md)
* Nanopore variant calling using **Singularity**: [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity.md)
* **Nanopore R10.4 Q20** variant calling: [Link](./docs/pipeline_docker/ONT_variant_calling_r10_q20.md)

#### PacBio-HiFi variant calling
* PacBio-HiFi variant calling using **Docker**: [Link](./docs/pipeline_docker/HiFi_variant_calling.md)
* PacBio-HiFi variant calling using **Singularity**: [Link](./docs/pipeline_singularity/HiFi_variant_calling_singularity.md)

### License
[PEPPER license](./LICENSE), [Margin License](https://github.com/UCSC-nanopore-cgl/margin/blob/master/LICENSE.txt) and [DeepVariant License](https://github.com/google/deepvariant/blob/r1.1/LICENSE) extend to the trained models (PEPPER, Margin and DeepVariant) and container environment (Docker and Singularity).

### Acknowledgement
We are thankful to the developers of these packages:
* [htslib & samtools](http://www.htslib.org/)
* [pytorch](https://pytorch.org/)
* [ONNX](https://onnx.ai/)
* [hdf5 python (h5py)](https://www.h5py.org/)

### Authors
[PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline is developed in a collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/) and the [Genomics team in Google Health](https://health.google/health-research/genomics/).


### Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="Iron-Man" width="240p"> <br/>

The name "P.E.P.P.E.R." is inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616).

PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.
