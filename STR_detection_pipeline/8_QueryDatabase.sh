#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <project_name> <subname> [roi_bed]"
    exit 1
fi

PROJECT_NAME=$1
SUBNAME=$2
ROI_BED=$3 # optional

OUTPUT_DIR="${PROJECT_NAME}/output"
DB_DIR="${OUTPUT_DIR}/STR_DBs"

BLAT_DIR="${OUTPUT_DIR}/BLAT/${SUBNAME}/PSLs"
JSON_FILE="${OUTPUT_DIR}/BLAT/${SUBNAME}/str_motifs_RFC1.json"

echo "Querying STR sequences..."

QUERY_OUTPUT_DIR="${OUTPUT_DIR}/QueryResults"
mkdir -p ${QUERY_OUTPUT_DIR}


process_query() {
    local gene_motif=$1
    local gene=$(echo ${gene_motif} | cut -d'_' -f1)
    local motif=$(echo ${gene_motif} | cut -d'_' -f2)
    local db_path="${DB_DIR}/${SUBNAME}_${gene_motif}.db"
    local output_file="${QUERY_OUTPUT_DIR}/${SUBNAME}_${gene_motif}_results_try2.txt"

    if [ -f "${db_path}" ]; then
        echo "Querying sequences for ${gene_motif}..."
        
        local roi_args=""
        if [ -n "${ROI_BED}" ]; then
            roi_args="--roi-bed ${ROI_BED}"
        else
            roi_args="--json-file ${JSON_FILE}"
        fi

        /opt/conda/bin/python AllScripts/python_scripts/wdl_query_STR_db.py \
            --db-path ${db_path} \
            --gene ${gene} \
            --motif ${motif} \
            ${roi_args} \
            --output-file ${output_file} \
            --allowed-patterns AAAAG ACGGG
    else
        echo "Database ${db_path} not found for ${gene_motif}"
    fi
}

if [ -n "${ROI_BED}" ]; then
    while IFS=$'\t' read -r chr start end gene motif; do
        process_query "${gene}_${motif}"
    done < ${ROI_BED}
else
    for psl_dir in ${BLAT_DIR}/*AAGGG*/; do
        if [ -d "${psl_dir}" ]; then
            gene_motif=$(basename ${psl_dir})
            process_query "${gene_motif}"
        fi
    done
fi

echo "Querying complete!"
