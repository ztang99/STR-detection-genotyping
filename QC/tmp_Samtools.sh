#!/bin/bash

##############################################################################

# Helper to compute average coverage based on `samtools depth`.

## author: Zitian Tang
## contact: tang.zitian@wustl.edu

##############################################################################

set -e  # Exit immediately if a command fails

# Check if a BAM file is provided as input
if [ $# -ne 2 ]; then
    echo "Usage: $0 <bam_file> <output>"
    exit 1
fi

# Get the BAM file path
bam_file="$1"
output_file="$2"

if [ ! -f "$bam_file" ]; then
    echo "Error: BAM file '$bam_file' not found"
    exit 1
fi

# Calculate average coverage and write to output file
samtools depth "$bam_file" | awk '{sum+=$3} END {print "Average coverage = " sum/NR}' > "$output_file"
