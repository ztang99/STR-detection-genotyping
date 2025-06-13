#!/bin/bash

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

OUTPUT_DIR="${PROJECT_NAME}/output"

EHDN_RESULTS="${OUTPUT_DIR}/EHdn/${SUBNAME}/EHdn_combined_results.csv"
EH_RESULTS="${OUTPUT_DIR}/EH/${SUBNAME}/EH_combined_all4rfc1.vcf"

WORKDIR="${OUTPUT_DIR}/BLAT/${SUBNAME}"
mkdir -p ${WORKDIR}

ALL_BAMS_LIST="${WORKDIR}/all_bams.txt"
cat ${CASE_BAM_PATHS} ${CONTROL_BAM_PATHS} > ${ALL_BAMS_LIST}

COMBINED_JSON="${WORKDIR}/str_motifs_RFC1.json"

# CombineEHdnEH
echo "Combining EHdn and EH results..."

CMD="/opt/conda/bin/python AllScripts/python_scripts/wdl_combine_ehdn_eh.py \
    --ehdn-results ${EHDN_RESULTS} \
    --eh-results ${EH_RESULTS} \
    --bams ${ALL_BAMS_LIST} \
    --min-overlap-percent 10 \
    --output-file ${COMBINED_JSON}"

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

# Step 2: Generate SAM files

cat > ${WORKDIR}/process_json.py << 'EOL'
import json
import sys

json_file = sys.argv[1]
with open(json_file) as f:
    data = json.load(f)
    
for entry in data:
    region = f"chr{entry['chrom']}:{entry['start']}-{entry['end']}"
    for carrier in entry['carriers']:
        print(f"{entry['gene']}\t{entry['motif']}\t{region}\t{carrier}")
EOL

# Process the JSON file using bsub to get motif-carrier pairs
/opt/conda/bin/python ${WORKDIR}/process_json.py ${COMBINED_JSON} > ${WORKDIR}/motif_carriers.txt

while IFS=$'\t' read -r gene motif region bam_path; do
    output_dir=${WORKDIR}/SAMs/${gene}_${motif}
    mkdir -p ${output_dir}
    sample_name=$(basename ${bam_path} .bam | sed 's/\.cram$//')
    output_sam=${output_dir}/${sample_name}.sam

    echo "Generating SAM for ${sample_name}, ${gene}, ${motif}..."
    samtools view ${bam_path} ${region} > ${output_sam}
    
    if [ ! -f "${output_sam}" ]; then
        echo "Error: SAM file not generated for ${sample_name}, ${gene}, ${motif}"
        exit 1
    fi
done < ${WORKDIR}/motif_carriers.txt

rm ${WORKDIR}/process_json.py ${WORKDIR}/motif_carriers.txt


# Step 3: Generate FASTA files
GENE_OF_INTEREST="RFC1"

for dir in ${WORKDIR}/SAMs/${GENE_OF_INTEREST}_*; do
    if [ ! -d "$dir" ]; then
        echo "No directories found matching ${GENE_OF_INTEREST}_*"
        continue
    fi

    gene_motif=$(basename ${dir})
    echo "Processing ${gene_motif}..."

    if [ "${gene_motif}" != "RFC1_AAGGG" ]; then # or AAGAC
        echo "Skipping ${gene_motif}"
        continue
    fi
    
    for sam_file in ${dir}/*.sam; do
        if [ ! -f "$sam_file" ]; then
            echo "No SAM files found in ${dir}"
            continue
        fi

        sample_name=$(basename ${sam_file} .sam)

        FASTA_DIR="${WORKDIR}/FASTAs/${gene_motif}"
        PSL_DIR="${WORKDIR}/PSLs/${gene_motif}"
        mkdir -p "${FASTA_DIR}"
        mkdir -p "${PSL_DIR}"

        fasta_file="${FASTA_DIR}/${sample_name}.fa"
        psl_file="${PSL_DIR}/${sample_name}.psl"

        echo "Processing SAM file: ${sample_name}"
        # bsub -G compute-jin810 -q general-interactive -R 'rusage[mem=4GB]' -a 'docker(elle72/basic:vszt)' \
        if [ ! -f "$fasta_file" ]; then
            /opt/conda/bin/python AllScripts/python_scripts/wdl_query_STR_db.py \
                filter_reads_to_fasta \
                --sam-file ${sam_file} \
                --output-file ${fasta_file} \
                --mapq-threshold 1 \
                --append false
        fi

        if [ ! -f "${fasta_file}" ]; then
            echo "Error: FASTA file not generated for ${sample_name}"
            continue
        fi

        echo "Submitted BLAT job for ${sample_name}..."
        bsub -g /tang.zitian/misc -G compute-jin810-t3 -q subscription -sla jin810_t3 -R 'rusage[mem=6GB]' -a 'docker(ztang301/basic:vsamtools_samblaster_blat_2)' \
            blat ${REF} ${fasta_file} ${psl_file} -t=dna -q=dna -repMatch=1000000
    done
done
