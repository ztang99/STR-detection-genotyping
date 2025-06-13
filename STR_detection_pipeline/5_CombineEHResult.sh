#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Error: Please provide required arguments: 1) EHdn combined result; 2) BAM txt for cases; 3) BAM txt for controls."
    exit 1
fi

PROJECT_NAME="$1"
SUBNAME="$2"

OUTPUT_DIR="${PROJECT_NAME}/output"
REF="Homo_sapiens_assembly38.fasta"

MANIFEST_FILE="${OUTPUT_DIR}/EHdn/${SUBNAME}/EHdn_manifest.tsv"

WORKDIR="${OUTPUT_DIR}/EH/${SUBNAME}"
PATH_TO_CASE_RESULTS="${WORKDIR}/cases_results"
PATH_TO_CONTROL_RESULTS="${WORKDIR}/controls_results"


EH_COMBINED_IR="${WORKDIR}/EH_combined_ir4all.vcf"
EH_COMBINED_RFC1="${WORKDIR}/EH_combined_all4rfc1.vcf"


echo "Combining ExpansionHunter results in INREPEAT mode..."

vcf_files=()
for vcf in ${PATH_TO_CASE_RESULTS}/*.vcf ${PATH_TO_CONTROL_RESULTS}/*.vcf; do
    if [ -f "$vcf" ]; then
        sample_name=$(basename "$vcf" .vcf)
        if grep -q "$sample_name" "$MANIFEST_FILE"; then
            vcf_files+=("$vcf")
        fi
    fi
done

if [ ! -f "${EH_COMBINED_IR}" ]; then
    /opt/conda/bin/python AllScripts/python_scripts/wdl_filter_eh_vcfs.py \
        --output-file "${EH_COMBINED_IR}" \
        --vcf-files "${vcf_files[@]}" \
        --mode inrepeat4all
fi

if [ ! -f "${EH_COMBINED_RFC1}" ]; then
    /opt/conda/bin/python AllScripts/python_scripts/wdl_filter_eh_vcfs.py \
        --output-file "${EH_COMBINED_RFC1}" \
        --vcf-files "${vcf_files[@]}" \
        --mode all4gene \
        --gene-name RFC1
fi


if [ $? -ne 0 ]; then
    echo "Error: Failed to combine ExpansionHunter results"
    exit 1
fi

echo "All processes completed successfully!"
