#!/bin/bash

##############################################################################

# Step 7: Run BLAT in defined STR regions.

## Docker: ztang301/exph:v2.1
## Docker for BLAT: ztang301/basic:vsamtools_samblaster_blat_2

## author: Zitian Tang
## contact: tang.zitian@wustl.edu

##############################################################################

# Check minimum required arguments (5 arguments, with roi_bed being optional)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <project_name> <subname>"
    exit 1
fi

PROJECT_NAME=$1
SUBNAME=$2

REF="Homo_sapiens_assembly38.fasta"

OUTPUT_DIR="${PROJECT_NAME}/output"

WORKDIR="${OUTPUT_DIR}/BLAT/${SUBNAME}"
mkdir -p ${WORKDIR}

COMBINED_JSON="${WORKDIR}/str_motifs_RFC1.json"


# Generate SAM files

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


# Generate FASTA files

for dir in ${WORKDIR}/SAMs/*; do
    if [ ! -d "$dir" ]; then
        echo "No directories found"
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
        blat ${REF} ${fasta_file} ${psl_file} -t=dna -q=dna -repMatch=1000000
    done
done
