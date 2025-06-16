#!/bin/bash

##############################################################################

# Step 6: Combine ExpansionHunterDenovo and ExpansionHunter results.

## Docker: ztang301/exph:v2.1

## author: Zitian Tang
## contact: tang.zitian@wustl.edu

##############################################################################

# Usage: bash 5_RunBLAT.sh <ehdn_results> <eh_results> <case_bams_list> <control_bams_list> [roi_bed]

# Check minimum required arguments (5 arguments, with roi_bed being optional)
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <ehdn_results> <eh_results> <case_bams_list> <control_bams_list> [roi_bed]"
    exit 1
fi

PROJECT_NAME=$1
SUBNAME=$2
CASE_BAM_PATHS=$3
CONTROL_BAM_PATHS=$4
ROI_BED=$5  # Optional parameter, can be empty

OUTPUT_DIR="/storage1/fs1/jin810/Active/testing/ztang/code/TRE_IPN/${PROJECT_NAME}/output"
REF="/storage1/fs1/jin810/Active/References/GRCh38/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta"

EHDN_RESULTS="${OUTPUT_DIR}/EHdn/${SUBNAME}/EHdn_combined_results.csv"
EH_RESULTS="${OUTPUT_DIR}/EH/${SUBNAME}/EH_combined_all4DRG20genes.vcf"
WORKDIR="${OUTPUT_DIR}/BLAT/${SUBNAME}"

mkdir -p ${WORKDIR}

ALL_BAMS_LIST="${WORKDIR}/all_bams.txt"
cat ${CASE_BAM_PATHS} ${CONTROL_BAM_PATHS} > ${ALL_BAMS_LIST}

COMBINED_JSON="${WORKDIR}/ConsensusSTRMotifs.json"

# CombineEHdnEH
echo "Combining EHdn and EH results..."

## 20250512 - passed in skipRM to bypass repeatmasker check and keep all shared motifs
CMD="/opt/conda/bin/python python_scripts/wdl_combine_ehdn_eh.py \
    --ehdn-results ${EHDN_RESULTS} \
    --eh-results ${EH_RESULTS} \
    --bams ${ALL_BAMS_LIST} \
    --min-overlap-percent 10 \
    --output-file ${COMBINED_JSON} \
    --skipRM"

# Add ROI_BED parameter only if it's provided
if [ ! -z "${ROI_BED}" ]; then
    CMD="${CMD} --roi-bed ${ROI_BED}"
fi

# # Execute the command using bsub
${CMD}

if [ ! -f "${COMBINED_JSON}" ]; then
    echo "Error: Combined results file not found: ${COMBINED_JSON}"
    exit 1
fi