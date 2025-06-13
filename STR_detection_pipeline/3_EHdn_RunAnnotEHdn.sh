#!/bin/bash

# Check if required arguments are provided
if [ -z "$2" ]; then
    echo "Error: Please provide 1) manifest file path"
    exit 1
fi

PROJECT_NAME="$1"
SUBNAME="$2"

OUTPUT_DIR="${PROJECT_NAME}/output"
REF="Homo_sapiens_assembly38.fasta"

ANNOVAR_VARIATION="annovar_20191024/annotate_variation.pl"
ANNOVAR_HUMANDB="annovar_20191024/humandb"

EHDN_OUTPUT_PREFIX="${OUTPUT_DIR}/EHdn/${SUBNAME}"
MANIFEST_FILE="${EHDN_OUTPUT_PREFIX}/EHdn_manifest.tsv"
EHDN_MULTI_PROFILE="${EHDN_OUTPUT_PREFIX}/${SUBNAME}.multisample_profile.json"
EHDN_OTL_LOCUS="${OUTPUT_DIR}/EHdn/${SUBNAME}/outliers_locus.tsv"
EHDN_CACO_LOCUS="${OUTPUT_DIR}/EHdn/${SUBNAME}/casecontrol_locus.tsv"
EHDN_OTL_LOCUS_ANNOT="${OUTPUT_DIR}/EHdn/${SUBNAME}/outliers_locus_hg38_annotated.tsv"
EHDN_CACO_LOCUS_ANNOT="${OUTPUT_DIR}/EHdn/${SUBNAME}/casecontrol_locus_hg38_annotated.tsv"


# Run EHdn merge step
echo "Running ExpansionHunterDenovo merge..."
bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=4GB]' -a 'docker(ztang301/exph:v1.2)' \
    /ExpansionHunterDenovo/build/ExpansionHunterDenovo merge \
    --reference "${REF}" \
    --manifest "${MANIFEST_FILE}" \
    --output-prefix "${EHDN_OUTPUT_PREFIX}/${SUBNAME}"

# Check if merge was successful
if [ ! -f "${EHDN_MULTI_PROFILE}" ]; then
    echo "Error: Merge output file not found: ${EHDN_MULTI_PROFILE}"
    exit 1
fi

# Run outlier analysis
echo "Running outlier analysis..."
bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=4GB]' -a 'docker(ztang301/exph:v1.2)' \
    /usr/bin/python3 /ExpansionHunterDenovo/scripts/outlier.py locus \
    --manifest "${MANIFEST_FILE}" \
    --multisample-profile "${EHDN_MULTI_PROFILE}" \
    --output "${EHDN_OTL_LOCUS}"

if [ $? -ne 0 ]; then
    echo "Error: Outlier analysis failed"
    exit 1
fi

# Run case-control analysis
echo "Running case-control analysis..."
bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=4GB]' -a 'docker(ztang301/exph:v1.2)' \
    /usr/bin/python3 /ExpansionHunterDenovo/scripts/casecontrol.py locus \
    --manifest "${MANIFEST_FILE}" \
    --multisample-profile "${EHDN_MULTI_PROFILE}" \
    --output "${EHDN_CACO_LOCUS}"

if [ $? -ne 0 ]; then
    echo "Error: Case control analysis failed"
    exit 1
fi


# Run gene annotation
echo "Running OTL gene annotation..."
bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=8GB]' -a 'docker(ztang301/exph:v1.2)' \
    bash /ExpansionHunterDenovo/scripts/annotate_ehdn.sh \
        --ehdn-results "${EHDN_OTL_LOCUS}" \
        --ehdn-annotated-results "${EHDN_OTL_LOCUS_ANNOT}" \
        --annovar-annotate-variation "${ANNOVAR_VARIATION}" \
        --annovar-humandb "${ANNOVAR_HUMANDB}" \
        --annovar-buildver hg38

echo "Running CACO gene annotation..."
bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=8GB]' -a 'docker(ztang301/exph:v1.2)' \
    bash /ExpansionHunterDenovo/scripts/annotate_ehdn.sh \
        --ehdn-results "${EHDN_CACO_LOCUS}" \
        --ehdn-annotated-results "${EHDN_CACO_LOCUS_ANNOT}" \
        --annovar-annotate-variation "${ANNOVAR_VARIATION}" \
        --annovar-humandb "${ANNOVAR_HUMANDB}" \
        --annovar-buildver hg38

if [ $? -ne 0 ]; then
    echo "Error: Gene annotation failed"
    exit 1
fi

echo "Running EHdn gene-based annotation..."
REPEAT_MASKER="Samplemaps/hg38RM_simple_repeats.bed"
REF_PN_GENES="Samplemaps/PNPAN_PNrelated.csv"
REF_FUNC_GENES="Samplemaps/NIHGene_CellFunc.csv"
REF_STR_GENES="Samplemaps/Malik_STR_genes.csv"
REF_HIGH_DRG_EXP="Samplemaps/high_DRGexp_genes.csv"

FILTERED_RESULT_DIR="${OUTPUT_DIR}/EHdn/${SUBNAME}"
mkdir -p "${FILTERED_RESULT_DIR}"

bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=20GB]' -a 'docker(elle72/basic:vszt)' \
    /opt/conda/bin/python AllScripts/python_scripts/wdl_filter_ehdn_results.py \
        --outlier-locus "${EHDN_OTL_LOCUS_ANNOT}" \
        --casecontrol-locus "${EHDN_CACO_LOCUS_ANNOT}" \
        --output-dir "${FILTERED_RESULT_DIR}" \
        --output-file "${FILTERED_RESULT_DIR}/EHdn_combined_results.csv" \
        --repeatmasker-file "${REPEAT_MASKER}" \
        --case-count 788 \
        --control-count 879 \
        --gene-list-files "${REF_PN_GENES},${REF_FUNC_GENES},${REF_STR_GENES},${REF_HIGH_DRG_EXP}"

if [ $? -ne 0 ]; then
    echo "Error: EHdn results annotation failed"
    exit 1
fi



# Check if output files were created
if [ -f "${EHDN_MULTI_PROFILE}" ] && \
   [ -f "${EHDN_OTL_LOCUS}" ] && \
   [ -f "${EHDN_CACO_LOCUS}" ] && \
   [ -f "${EHDN_OTL_LOCUS_ANNOT}" ] && \
   [ -f "${EHDN_CACO_LOCUS_ANNOT}" ]; then
    echo "Processing completed successfully"
    echo "Multisample profile: ${EHDN_MULTI_PROFILE}"
    echo "Outliers locus file: ${EHDN_OTL_LOCUS}"
    echo "Case control locus file: ${EHDN_CACO_LOCUS}"
    echo "Annotated outliers results: ${EHDN_OTL_LOCUS_ANNOT}"
    echo "Annotated case-control results: ${EHDN_CACO_LOCUS_ANNOT}"
else
    echo "Error: Some output files are missing"
    exit 1
fi