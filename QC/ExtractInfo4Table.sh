#!/bin/bash


# Function to get value after a specific pattern from samtools stats
get_stats_value() {
    local file=$1
    local pattern=$2
    if [ "$pattern" = "bases mapped (cigar)" ]; then
        awk -F'\t' '$2 ~ /bases mapped \(cigar\)/ {print $3}' "$file"
    else
        awk -F'\t' -v pat="$pattern" '$2 ~ pat {print $3}' "$file"
    fi
}

# Function to get depth value
get_depth_value() {
    local file=$1
    awk -F'=' '/Average coverage/{print $2}' "$file" | tr -d ' '
}

# Main processing function
process_metrics() {
    local work_dir=$1
    local stats_dir="${work_dir}/stats"
    local depths_dir="${work_dir}/depths"
    local output_csv="${work_dir}/QCmetricsCombined.csv"

    echo "Sample,Read Length,Number of reads (G),Mean_Coverage (X),Reads mapped and paired (%),Reads duplicated (%),Bases mapped (%),Bases duplicated (%),Insert size,,Mean error rate (%)" > "$output_csv"

    # Process each stats file
    for stats_file in "${stats_dir}"/*_samtools.stats.txt; do
        # Extract sample name
        local sample_name=$(basename "$stats_file" _samtools.stats.txt)
        local depth_file="${depths_dir}/${sample_name}_samtools.depth.txt"

        echo "Processing sample: $sample_name"

        # Get values from samtools stats
        local read_length=$(get_stats_value "$stats_file" "average length")
        local total_reads=$(get_stats_value "$stats_file" "raw total sequences")
        local mapped_paired=$(get_stats_value "$stats_file" "reads mapped and paired")
        local duplicated_reads=$(get_stats_value "$stats_file" "reads duplicated")
        local total_length=$(get_stats_value "$stats_file" "total length")
        local bases_mapped=$(get_stats_value "$stats_file" "bases mapped (cigar)")
        local bases_duplicated=$(get_stats_value "$stats_file" "bases duplicated")
        local error_rate=$(get_stats_value "$stats_file" "error rate")
        local insert_size=$(get_stats_value "$stats_file" "insert size average")
        local insert_std=$(get_stats_value "$stats_file" "insert size standard deviation")

        # Get mean coverage from depth file
        local mean_coverage="NA"
        if [ -f "$depth_file" ]; then
            mean_coverage=$(get_depth_value "$depth_file")
        fi

        # Calculate metrics
        local num_reads_g=$(awk "BEGIN {printf \"%.2f\", $total_reads/1000000000}")
        local pct_mapped_paired=$(awk "BEGIN {printf \"%.2f\", $mapped_paired/$total_reads * 100}")
        local pct_duplicated=$(awk "BEGIN {printf \"%.2f\", $duplicated_reads/$total_reads * 100}")
        local pct_bases_mapped=$(awk "BEGIN {printf \"%.2f\", $bases_mapped/$total_length * 100}")
        local pct_bases_duplicated=$(awk "BEGIN {printf \"%.2f\", $bases_duplicated/$total_length * 100}")
        local error_rate_pct=$(awk "BEGIN {printf \"%.2f\", $error_rate * 100}")
        local insert_size_combined="${insert_size} ± ${insert_std}"

        # Write to CSV
        echo "${sample_name},${read_length},${num_reads_g},${mean_coverage},${pct_mapped_paired},${pct_duplicated},${pct_bases_mapped},${pct_bases_duplicated},\"${insert_size_combined}\",${error_rate_pct}" >> "$output_csv"
    done

    echo "Processing complete! Results written to: $output_csv"
}

# Check if path argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <working_directory>"
    echo "Note: working_directory should contain 'stats' and 'depths' subdirectories"
    exit 1
fi

# Check if the required directories exist
if [ ! -d "$1/stats" ] || [ ! -d "$1/depths" ]; then
    echo "Error: 'stats' and 'depths' subdirectories must exist in the working directory"
    exit 1
fi

process_metrics "$1"

