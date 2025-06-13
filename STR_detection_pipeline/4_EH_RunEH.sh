#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Error: Please provide required arguments: 1) EHdn combined result; 2) BAM txt for cases; 3) BAM txt for controls."
    exit 1
fi

export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /scratch1/fs1/jin810:/scratch1/fs1/jin810 /storage2/fs1/epigenome/Active:/storage2/fs1/epigenome/Active /home/tang.zitian:/home/tang.zitian'

PROJECT_NAME="$1"
SUBNAME="$2"
CASE_BAM_PATHS="$3"
CONTROL_BAM_PATHS="$4"

OUTPUT_DIR="${PROJECT_NAME}/output"
REF="Homo_sapiens_assembly38.fasta"

EHDN_RESULTS="${OUTPUT_DIR}/EHdn/${SUBNAME}/EHdn_combined_results.csv"
WORKDIR="${OUTPUT_DIR}/EH/${SUBNAME}"
mkdir -p "${WORKDIR}"

EH_CATALOG_JSON="${WORKDIR}/EH_variant_catalog.json"
PATH_TO_CASE_RESULTS="${WORKDIR}/cases_results"
PATH_TO_CONTROL_RESULTS="${WORKDIR}/controls_results"

EH_COMBINED="${WORKDIR}/EH_combined_allRFC1.vcf"


run_expansion_hunter() {
    local bam_list="$1"
    local output_dir="$2"

    echo "Running ExpansionHunter..."
    while IFS= read -r bam_path; do
        sample_name_base=$(basename "${bam_path}" .bam | sed 's/\.cram$//')
        sample_name=$(echo "$sample_name_base" | sed 's/\^/_/g')
        out_prefix="${output_dir}/${sample_name}"
        out_vcf="${out_prefix}.vcf"

        if [ -f "${out_vcf}" ]; then
            echo "Skipping ${sample_name} - output file already exists"
            echo "----------------------------------------"
            continue
        fi

        echo "Processing sample: ${sample_name}"
        echo "BAM path: ${bam_path}"

        bsub -G compute-jin810 -q general-interactive -R 'rusage[mem=4GB]' -a 'docker(ztang301/exph:v2.1)' \
            /ExpansionHunter/build/install/bin/ExpansionHunter \
                --reads ${bam_path} \
                --reference ${REF} \
                --variant-catalog ${EH_CATALOG_JSON} \
                --output-prefix ${out_prefix}

        if [ $? -eq 0 ]; then
            echo "Successfully processed ${sample_name}"
        else
            echo "Error processing ${sample_name}"
        fi
        
        echo "----------------------------------------"
    done < "${bam_list}"
}

# Generate EH catalog
if [ ! -f "${EH_CATALOG_JSON}" ]; then
    echo "Generating ExpansionHunter catalog..."
    bsub -K -G compute-jin810 -q general-interactive -n 1 -R 'rusage[mem=4GB]' -a 'docker(elle72/basic:vszt)' \
        /opt/conda/bin/python AllScripts/python_scripts/wdl_IPN_generate_EHcatalog.py \
        "${EH_CATALOG_JSON}" "${EHDN_RESULTS}"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate EH catalog"
        exit 1
    fi
fi

# Run EH for cases and controls
run_expansion_hunter "${CASE_BAM_PATHS}" "${PATH_TO_CASE_RESULTS}"
run_expansion_hunter "${CONTROL_BAM_PATHS}" "${PATH_TO_CONTROL_RESULTS}"

# Wait for all EH jobs to complete
echo "Waiting for all ExpansionHunter jobs to complete..."
job_id_string=$(IFS=":"; echo "${job_ids[*]}")
if [ -n "$job_id_string" ]; then
    bwait -w "ended(${job_id_string})"
fi

# Check outputs
echo -e "\nChecking EH individual outputs..."

# Check EH catalog
if [ -f "${EH_CATALOG_JSON}" ]; then
    echo "EH catalog file generated successfully"
fi

# Count case VCF files
case_vcf_count=$(ls "${PATH_TO_CASE_RESULTS}"/*.vcf 2>/dev/null | wc -l)
total_cases=$(wc -l < "${CASE_BAM_PATHS}")
echo "Cases processed: ${case_vcf_count}/${total_cases}"

# Count control VCF files
control_vcf_count=$(ls "${PATH_TO_CONTROL_RESULTS}"/*.vcf 2>/dev/null | wc -l)
total_controls=$(wc -l < "${CONTROL_BAM_PATHS}")
echo "Controls processed: ${control_vcf_count}/${total_controls}"

# Summary status
echo -e "\nFinal Status:"
if [ -f "${EH_CATALOG_JSON}" ] && \
   [ ${case_vcf_count} -eq ${total_cases} ] && \
   [ ${control_vcf_count} -eq ${total_controls} ]; then
    echo "Workflow completed successfully"
    echo "  - EH catalog generated"
    echo "  - All ${total_cases} case samples processed"
    echo "  - All ${total_controls} control samples processed"
else
    echo "Workflow completed with missing outputs:"
    [ ! -f "${EH_CATALOG_JSON}" ] && echo "  - Missing EH catalog"
    [ ${case_vcf_count} -ne ${total_cases} ] && echo "  - Missing $(($total_cases - $case_vcf_count)) case VCF files"
    [ ${control_vcf_count} -ne ${total_controls} ] && echo "  - Missing $(($total_controls - $control_vcf_count)) control VCF files"
fi




