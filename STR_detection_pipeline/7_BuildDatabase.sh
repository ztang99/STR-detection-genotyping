#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <project_name> <subname> [roi_bed]"
    exit 1
fi

PROJECT_NAME=$1
SUBNAME=$2
ROI_BED=$3 #optional

OUTPUT_DIR="${PROJECT_NAME}/output"
DB_DIR="${OUTPUT_DIR}/STR_DBs"
mkdir -p ${DB_DIR}

BLAT_DIR="${OUTPUT_DIR}/BLAT/${SUBNAME}/PSLs"
SAM_DIR="${OUTPUT_DIR}/BLAT/${SUBNAME}/SAMs"

process_gene_motif() {
    local gene_motif=$1
    local gene=$(echo ${gene_motif} | cut -d'_' -f1)
    local motif=$(echo ${gene_motif} | cut -d'_' -f2)

    echo "Processing ${gene_motif}..."
    PSL_DIR="${BLAT_DIR}/${gene_motif}"
    SAM_SUBDIR="${SAM_DIR}/${gene_motif}"
    
    if [ -d "${PSL_DIR}" ]; then
        psl_count=$(ls ${PSL_DIR}/*.psl 2>/dev/null | wc -l)
        
        if [ "${psl_count}" -gt 0 ]; then
            echo "Found ${psl_count} PSL files under ${PSL_DIR}"

            DB_NAME="${SUBNAME}_${gene_motif}.db"

            DB_PATH="${DB_DIR}/${DB_NAME}"
            
            if [ ! -f "${DB_PATH}" ]; then
                echo "Processing SAM files and initializing database: ${DB_NAME}"
                /opt/conda/bin/python AllScripts/python_scripts/wdl_addBlatResult2db.py \
                    --mode init \
                    --db-path ${DB_PATH} \
                    --gene ${gene} \
                    --sam-dir ${SAM_SUBDIR}
            else
                echo "Database ${DB_NAME} already exists, skipping initialization"
            fi
            
            echo "Processing BLAT results for ${gene}_${motif}..."
            
            /opt/conda/bin/python /storage1/fs1/jin810/Active/testing/ztang/code/TRE_IPN/pipeline_scripts/AllScripts/python_scripts/wdl_addBlatResult2db.py \
                --mode blat \
                --db-path ${DB_PATH} \
                --psl-files ${PSL_DIR}/*.psl
                
            echo "Processing complete for ${gene}_${motif}. Database at: ${DB_PATH}"
        else
            echo "No PSL files found in directory ${PSL_DIR}"
        fi
    else
        echo "${gene}_${motif} was not detected in the current ${PROJECT_NAME}_${SUBNAME} by EHdn and EH"
    fi
    
    echo "----------------------------------------"
}

if [ -n "${ROI_BED}" ]; then
    # Process from ROI bed file
    while IFS=$'\t' read -r chr start end gene motif; do
        process_gene_motif "${gene}_${motif}"
    done < ${ROI_BED}
else
    # Process all subdirectories in BLAT_DIR
    for psl_dir in ${BLAT_DIR}/*/; do
        if [ -d "${psl_dir}" ]; then
            gene_motif=$(basename ${psl_dir})
            process_gene_motif "${gene_motif}"
        fi
    done
fi

echo "Finished building database!"


