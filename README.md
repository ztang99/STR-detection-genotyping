# str-detection-genotyping

This repository contains the bioinformatics pipeline used for the analysis of whole genome sequencing data in the study of Inherited Peripheral Neuropathy (IPN). The pipeline includes modules for alignment, variant calling, quality control, gene annotation, and Short Tandem Repeat (STR) detection.

## Repository Structure
```markdown
├── Alignment\
│   └── parabricksFsq2Bams.sh            # FASTQ to BAM conversion using Parabricks\
├── Gene_annotation\
│   ├── processRawVariant.py             # Extract deleterious variants from annotations\
│   └── runSIFT.sh                       # Run SIFT for variant annotation\
├── QC\
│   ├── ethnicity_pred_gnomad_cont.py    # Ethnicity prediction using gnomAD\
│   ├── ExtractInfo4Table.sh             # Extract QC metrics from BAM files\
│   ├── IPNSamtoolsCheck.sh              # BAM file quality check using samtools\
│   ├── snvstory_ethnicity_check.sh      # Ethnicity check using SNVstory\
│   └── snvstory_selectVariant.sh        # Select variants from gVCF files\
└── STR_detection_pipeline\
    ├── python_scripts/                  # Helper Python scripts\
    ├── 00_RunAll.sh                     # Master script to run the entire STR pipeline\
    ├── 1_EHdn_GenerateStrProfile.sh     # Generate STR profiles using ExpansionHunterDenovo\
    ├── 2_EHdn_GenerateManifestFile.sh   # Create manifest file for EHdn\
    ├── 3_EHdn_RunAnnotEHdn.sh           # Run and annotate EHdn results\
    ├── 4_EH_RunEH.sh                    # Run ExpansionHunter\
    ├── 5_CombineEHResult.sh             # Combine EH results\
    ├── 6_RunBLAT.sh                     # Run BLAT alignment of STR regions\
    ├── 7_BuildDatabase.sh               # Build STR sequence database\
    ├── 8_QueryDatabase.sh               # Query STR database\
    └── 9_UnsupervisedIPNFinalStep.R     # Unsupervised clustering for RFC1 repeat analysis
```

## Pipeline Components

### Alignment
- **parabricksFsq2Bams.sh**: GPU-accelerated FASTQ to BAM alignment using NVIDIA Parabricks

### Variant Calling and Annotation
- **runSIFT.sh**: Run SIFT for deleterious gene annotation
- **processRawVariant.py**: Extract deleterious variants from SIFT annotation files

### Quality Control
- **IPNSamtoolsCheck.sh**: Quality check of BAM files using samtools
- **ExtractInfo4Table.sh**: Extract QC metrics from samtools output
- **ethnicity_pred_gnomad_cont.py**: Predict sample ancestry using gnomAD continental probabilities
- **snvstory_ethnicity_check.sh**: Run SNVstory for ancestry determination
- **snvstory_selectVariant.sh**: Select variants from gVCF files for SNVstory

### STR Detection Pipeline
1. **1_EHdn_GenerateStrProfile.sh**: Generate STR profiles using ExpansionHunterDenovo
2. **2_EHdn_GenerateManifestFile.sh**: Generate manifest file for ExpansionHunterDenovo
3. **3_EHdn_RunAnnotEHdn.sh**: Run and annotate ExpansionHunterDenovo results
4. **4_EH_RunEH.sh**: Run ExpansionHunter on detected STR regions
5. **5_CombineEHResult.sh**: Combine ExpansionHunter results
6. **6_RunBLAT.sh**: Run BLAT alignment of STR regions
7. **7_BuildDatabase.sh**: Build STR sequence database
8. **8_QueryDatabase.sh**: Query STR database
9. **9_UnsupervisedIPNFinalStep.R**: R script for unsupervised clustering of repeat expansions
10. All python scripts being called by bash scripts are stored inside `python_scripts/`.

## Usage

To run the entire STR detection pipeline:

```bash
cd STR_detection_pipeline
bash 00_RunAll.sh
```
Individual components can be run separately as needed.


## Requirements
- NVIDIA Parabricks 4.0.0-1
- SIFT4G Annotator
- samtools
- SNVstory v1.1
- ExpansionHunterDenovo
- ExpansionHunter v2.1
- BLAT
- R with libraries: tidyverse, cluster, factoextra, mclust, NbClust, dendextend, readxl

## Citation
If you use this pipeline in your research, please cite our paper:

[Paper citation information to be added upon publication]

## License
See the LICENSE file for details.
