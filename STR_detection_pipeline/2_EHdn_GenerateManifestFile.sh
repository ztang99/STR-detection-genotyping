#!/bin/bash

if [ -z "$4" ]; then
    echo "Error: Please provide at least 4 arguments."
    exit 1
fi

PROJECT_NAME="$1"
SUBNAME="$2"
CASE_BAMs="$3"
CASES_PROFILE_DIR="$4"
CONTROL_BAMs="$5"
CONTROLS_PROFILE_DIR="$6"

OUTPUT_DIR="${PROJECT_NAME}/output"
EHDN_OUTPUT_PREFIX="${OUTPUT_DIR}/EHdn/${SUBNAME}"
mkdir -p "${OUTPUT_DIR}" "${EHDN_OUTPUT_PREFIX}"

MANIFEST_FILE="${EHDN_OUTPUT_PREFIX}/EHdn_manifest.tsv"
> "${MANIFEST_FILE}"

# Process cases
echo "Processing cases..."
for filepath in "${CASES_PROFILE_DIR}"/*.str_profile.json; do
    if [ -f "$filepath" ]; then
        sample_name=$(basename "$filepath" .str_profile.json)
        if grep -q "${sample_name}" "$CASE_BAMs"; then  # only process when samplename also exist in BAMpath txt
            full_path=$(realpath "$filepath")
            echo -e "${sample_name}\tcase\t${full_path}" >> "${MANIFEST_FILE}"
        fi
    fi
done

# Process controls
echo "Processing controls..."
for filepath in "${CONTROLS_PROFILE_DIR}"/*.str_profile.json; do
    if [ -f "$filepath" ]; then
        sample_name=$(basename "$filepath" .str_profile.json)
        if grep -q "/${sample_name}[/_]" "$CONTROL_BAMs"; then
           full_path=$(realpath "$filepath")
           echo -e "${sample_name}\tcontrol\t${full_path}" >> "${MANIFEST_FILE}"
       fi
    fi
done



# Check if manifest file was created successfully
if [ -f "${MANIFEST_FILE}" ]; then
    echo "Manifest file created successfully at: ${MANIFEST_FILE}"
    echo "Number of samples processed:"
    echo "Cases: $(grep -c "case" "${MANIFEST_FILE}")"
    echo "Controls: $(grep -c "control" "${MANIFEST_FILE}")"
else
    echo "Error: Failed to create manifest file"
    exit 1
fi


