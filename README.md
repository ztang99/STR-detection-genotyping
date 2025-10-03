# str-detection-genotyping

## Introduction

This repository holds scripts of a short tandem repeat (STR) expansion detection pipeline.

This pipeline successfully detects *RFC1* STR expansions using short-read whole-genome sequencing (WGS) data in a idiopathic Peripheral Neuropathy (iPN) patient cohort. 
The pipeline includes modules for quality control, WGS alignment, STR calling, gene annotation, and STR genotyping prediction through unsupervised clustering.

## Table of Contents
- [Introduction](#introduction)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
    - [Prepare Run](#prepare-run)
        - [Prep Softwares](#prep-software-environment-docker-used)
        - [Clone Repo](#clone-repository)
        - [Prep Inputs](#prep-inputs)
    - [Run Pipeline](#run-pipeline)
- [Citation](#citation)
- [License](#license)

## Repository Structure
```markdown
├── Alignment/
    └── parabricksFsq2Bams.sh
├── QC/
    ├── ethnicity_pred_gnomad_cont.py
    ├── ExtractInfo4Table.sh
    ├── IPNSamtoolsCheck.sh
    ├── snvstory_ethnicity_check.sh
    └── snvstory_selectVariant.sh
├── STR_detection_pipeline/
    ├── python_scripts/
    ├── 00_RunAll.sh
    ├── 1_EHdn_GenerateStrProfile.sh
    ├── 2_EHdn_GenerateManifestFile.sh
    ├── 3_EHdn_RunAnnotEHdn.sh
    ├── 4_EH_RunEH.sh
    ├── 5_CombineEHResult.sh
    ├── 6_RunBLAT.sh
    ├── 7_BuildDatabase.sh
    ├── 8_QueryDatabase.sh
    └── 9_UnsupervisedIPNFinalStep.R
└── Gene_annotation/
    ├── processRawVariant.py
    └── runSIFT.sh
```

## Getting Started

### Prepare Run

#### Prep Software Environment (Docker used)
- NVIDIA Parabricks 4.0.0-1 (`nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1`)
- SIFT4G Annotator (`staphb/snpeff`)
- samtools (`elle72/basic:vszt`)
- SNVstory v1.1 (`mgibio/snvstory:v1.1-buster`)
- ExpansionHunterDenovo (`ztang301/exph:v1.2`)
- ExpansionHunter v2.1 (`ztang301/exph:v2.1`)
- BLAT (`ztang301/basic:vsamtools_samblaster_blat_2`)
- R with libraries: tidyverse, cluster, factoextra, mclust, NbClust, dendextend, readxl ([RStudio](https://posit.co/downloads/))

#### Clone Repository

First, `cd` into your working directory, and clone this repository:

```bash
git clone https://github.com/ztang99/STR-detection-genotyping.git
```

#### Prep Inputs

To prepare running the pipeline, place the following three files in `input/` folder:
- `cases.txt`: text file with each line denoting full path to a case's BAM file
- `controls.txt`: text file with each line denoting full path to a control's BAM file
- `roi.bed`: regions-of-interest BED file with each line denoting a STR expansion region of interest

Your folder structure should look like this when done:

```markdown
├── STR_detection_pipeline/
├── input/
    ├── cases.txt
    ├── controls.txt
    └── roi.bed
└── output/
```

### Run Pipeline

To run the entire STR detection pipeline, simply do:

```bash
cd STR_detection_pipeline
bash 00_RunAll.sh $JOB_NAME
```
where `$JOB_NAME` denotes the name to the current detection job. We recommend you use a different job name every time you run this pipeline. It will create a different subfolder to store the results.

You should expect the following folder structure after the pipeline finishes running (only subdirectories and key files are included below):

```markdown
├── STR_detection_pipeline/
├── input/
└── output/
    ├── BLAT/$JOB_NAME/
        ├── FASTAs/
        ├── PSLs/
        └── SAMs/
    ├── EH/$JOB_NAME/
        ├── cases_results/
        ├── controls_results/
        └── EH_combined.vcf
    ├── EHdn/$JOB_NAME/
        ├── EHdn_cases_str-profiles/
        ├── EHdn_controls_str-profiles/
        └── EHdn_combined_results.csv
    ├── QueryResults/$JOB_NAME/
        └── [One .txt file per gene+motif queried]
    └── STR_DBs/$JOB_NAME/
        └── [One .db file per gene+motif established]
```

## Citation
If you use this pipeline in your research, please cite our paper:

> Tang Z, Ovunc SS, Mehinovic E, et al. Heterozygous and Homozygous RFC1 AAGGG Repeat Expansions are Common in Idiopathic Peripheral Neuropathy. Preprint. medRxiv. 2025;2025.04.18.25325809. Published 2025 Apr 23. doi:10.1101/2025.04.18.25325809

## License
See the LICENSE file for details.
