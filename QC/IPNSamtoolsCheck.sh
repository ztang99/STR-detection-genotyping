#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file_list> <output_directory>"
    exit 1
fi

input_list="$1"
out_dir="$2"

if [ ! -f "$input_list" ]; then
    echo "Error: Input file list '$input_list' does not exist"
    exit 1
fi

mkdir -p "$out_dir" "${out_dir}/stats" "${out_dir}/depths" 

# Process each BAM file
while IFS= read -r input_bam; do
    [ -z "$input_bam" ] && continue
    
    if [ ! -f "$input_bam" ]; then
        echo "Warning: BAM file '$input_bam' does not exist, skipping..."
        continue
    fi
    
    sample_name=$(basename "${input_bam}" .bam | sed 's/\.cram$//')
    
    echo "Processing sample: $sample_name"
    
    stats_output="${out_dir}/stats/${sample_name}_samtools.stats.txt"
    depth_output="${out_dir}/depths/${sample_name}_samtools.depth.txt"
    
    if [ -s "$stats_output" ] && [ -s "$depth_output" ]; then
        echo ">>>Both output files already exist and are non-empty for $sample_name."
        continue
    fi
    
    if [ ! -s "$stats_output" ]; then
        echo "Running samtools stats..."
        bsub -R "rusage[mem=4GB]" -G compute-jin810 -q general -o /dev/null -e logs/20250204_IPNQC/IPN_%J.log -a 'docker(elle72/basic:vszt)' /bin/bash -c '(samtools stats -@ 4 '${input_bam}' > '${stats_output}')'
    else
        echo ">>>stats output already exists and is non-empty."
    fi
    
    if [ ! -s "$depth_output" ]; then
        echo "Running samtools depth..."
        bsub -g /tang.zitian/TRE -R "rusage[mem=4GB]" -G compute-jin810 -q general -oo logs/20250205_IPNQC_depth/%J.log -a "docker(elle72/basic:vszt)" bash code/QC/run_qc/tmp_Samtools.sh "${input_bam}" "${depth_output}"
    else
        echo ">>>depth output already exists and is non-empty."
    fi
    
    echo "Completed processing for $sample_name"
    echo "----------------------------------------"

done < "$input_list"

echo "All processing complete!"