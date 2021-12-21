## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the [PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline. This pipeline enables nanopore-based variant calling with [DeepVariant](https://github.com/google/deepvariant).

<p align="center">
<img src="./img/PMDV_variant_calling_ONT_v5.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="720p"></img>
</p>

### How to cite
Please cite the following manuscript if you are using `PEPPER-Margin-DeepVariant`:

<details>
<summary><a href="https://www.nature.com/articles/s41592-021-01299-w"><b>Nature Methods:</b> Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads. </a>
</summary>

Authors: Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov, Sidharth Goel, <br/> Gunjan Baid, Mikhail Kolmogorov, Jordan M. Eizenga, Karen H. Miga, Paolo Carnevali, Miten Jain, Andrew Carroll & Benedict Paten.

**Free access link to the manuscript:** [https://rdcu.be/cABfv]()
</details>



---

### Long read variant calling performance evaluation

Detailed variant calling performance evaluation:
* **ONT R9.4.1 Guppy 5.0.6 Sup:**
  * [HG003 whole genome (PEPPER vs Clair3).](./docs/performance_evaluation/Oxford_nanopore_r9_whole_genome.md)
  * [HG002-HG006 chr20 (PEPPER vs Clair3).](./docs/performance_evaluation/Oxford_nanopore_r9_chr20.md)
* **ONT R10.4 Q20 Guppy 5.0.15 Sup:**
  * [HG002 chr20 (PEPPER vs Clair3).](./docs/performance_evaluation/Oxford_nanopore_r10.4_q20_chr20.md)
* **PacBio HiFi:**
  * [HG003 chr20 (PEPPER-Margin-DeepVariant vs DeepVariant-WhatsHap-DeepVariant).](./docs/performance_evaluation/PacBio_hifi_HG003_chr20.md)

---

### Useful links to documentations
* [Quickstarts to check system configuration.](#quickstarts-small-runs-to-test-system-configuration)
* [Case-studies to reproduce performance.](#case-studies-chromosome-20-runs-for-performance-reproducibility)
* [Description of PEPPER-Margin-DeepVariant method.](./docs/misc/pepper_v0.7_method_update.md)
* How to train PEPPER-DeepVariant:
    * [How to train PEPPER-SNP.](./docs/training_pepper_margin_deepvariant/how_to_train_pepper_snp.md)
    * [How to train PEPPER-HP.](./docs/training_pepper_margin_deepvariant/how_to_train_pepper_hp.md)
    * [How to train DeepVariant (using candidates from PEPPER).](./docs/training_pepper_margin_deepvariant/how_to_train_deepvariant.md)
* [How to install PEPPER locally.](./docs/local_install/install_pepper_locally.md)
* [List of parameters and description.](./docs/usage/usage_and_parameters.md)
* [How to tune parameters to balance speed and accuracy.](./docs/usage/speed_and_accuracy.md)

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

### Quickstarts (small runs to test system configuration)
|                  Test type                    |                          Links                                             |
|:---------------------------------------------:|:--------------------------------------------------------------------------:|
|              Docker quickstart                |  [Link](./docs/quickstart/variant_calling_docker_quickstart.md)            |
|              Singularity quickstart           |  [Link](./docs/quickstart/variant_calling_singularity_quickstart.md)       |
|              Docker-gpu quickstart            |  [Link](./docs/quickstart/variant_calling_docker_gpu_quickstart.md)        |

### Case studies (chromosome 20 runs for performance reproducibility)

|                       Pipeline                                 |                         Docker                                 |                               Singularity                                       |                     NVIDIA-docker<br>(GPU)                             |
|:--------------------------------------------------------------:|:--------------------------------------------------------------:|:-------------------------------------------------------------------------------:|:----------------------------------------------------------------------:|
|              Nanopore R9.4.1<br>variant calling                |  [Link](./docs/pipeline_docker/ONT_variant_calling.md)         |  [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity.md)         |  [Link](./docs/pipeline_docker_gpu/ONT_variant_calling_gpu.md)         |
|              Nanopore R10.4 Q20<br>variant calling             |  [Link](./docs/pipeline_docker/ONT_variant_calling_r10_q20.md) |  [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity_r10_q20.md) |  [Link](./docs/pipeline_docker_gpu/ONT_variant_calling_r10_q20_gpu.md) |
|                  PacBio HiFi<br>variant calling                |  [Link](./docs/pipeline_docker/HiFi_variant_calling.md)        |  [Link](./docs/pipeline_singularity/HiFi_variant_calling_singularity.md)        | [Link](./docs/pipeline_docker_gpu/HiFi_variant_calling_gpu.md)         |


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
