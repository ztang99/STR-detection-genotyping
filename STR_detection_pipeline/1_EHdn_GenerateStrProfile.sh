#!/bin/bash

if [ -z "$3" ]; then
    echo "Error: Please provide 1) txt file with BAM paths; 2) Mode - cases or controls; 3) outdir."
    exit 1
fi

BAM_PATHS="$1"
MODE="$2"
OUTPUT_DIR="$3"

REF="Homo_sapiens_assembly38.fasta"

PROFILE_DIR="${OUTPUT_DIR}/EHdn/EHdn_${MODE}_str-profiles"
mkdir -p "${PROFILE_DIR}"

while IFS= read -r bam_path; do
    sample_name=$(basename "${bam_path}" .bam | sed 's/\.cram$//')

    if [ -f "${PROFILE_DIR}/${sample_name}.str_profile.json" ]; then
        echo "Skipping ${sample_name} ------------------------>>"
        continue
    fi
    
    echo "Processing sample: ${sample_name}"
    echo "BAM path: ${bam_path}"
    
    # bsub -G compute-jin810 -q general-interactive -n 2 -R 'rusage[mem=6GB]' -a 'docker(ztang301/exph:v1.2)' \
    /ExpansionHunterDenovo/build/ExpansionHunterDenovo profile \
        --reads "${bam_path}" \
        --reference "${REF}" \
        --output-prefix "${PROFILE_DIR}/${sample_name}" \
        --min-anchor-mapq 50 \
        --max-irr-mapq 40
    
    if [ $? -eq 0 ]; then
        echo "Successfully processed ${sample_name}"
    fi
    
    echo "----------------------------------------"
done < "${BAM_PATHS}"

echo "All samples processed. Output directory: ${PROFILE_DIR}"