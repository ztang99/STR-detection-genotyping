#!/bin/bash

##############################################################################

# This file is for aligning combined FASTQs to BAMs. Each sample need to have only 1 FASTQ for R1 and R2, respectively.
# Script uses parabricks implementation of GATK Germline pipeline with 2 GPUs.
# Script takes in samplenames either through manual entering or a .txt file with each line denote a sample.
# Script can skip samples which has BAMs already generated in out_dir (i.e., resubmit only a proportion of samples).

##############################################################################

### EXPORTS ###

export PATH="/opt/miniconda/bin:$PATH"
export LSF_DOCKER_NETWORK=host
#export LSF_DOCKER_RUN_LOGLEVEL=DEBUG
# Use entry point because the parabricks container has other entrypoints but our cluster, by default, requires /bin/sh
export LSF_DOCKER_ENTRYPOINT=/bin/sh

### REF DO NOT MODIFY ###
PARA_REF="Homo_sapiens_assembly38.fasta"
pKS4="Homo_sapiens_assembly38.known_indels.vcf.gz"
T2T_REF="GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

txt_file="20240906_batch3_allsamples.txt"
sample_names=()
while IFS= read -r line; do
    line=$(echo "$line" | sed ''s/[^[:alnum:]-]//g'')
    sample_names+=("$line")
done < "$txt_file"

## File counts & output paths ##
count=0
current_date=$(date +"%Y%m%d")
jobname="rerun_batch3_3samples"
work_dir="INPUT_FILES"


fsq_dir="${work_dir}/temp_combined_fastqs/IPN_batch3_combined_fsqs"
out_bams="${work_dir}/BAM/batch3_bams"
out_gvcfs="${work_dir}/BAM/batch3_gvcfs"
out_metrics="${work_dir}/BAM/batch3_metrics"

### Rerun_batch3 3 samples ###
out_bams="${work_dir}/BAM/${jobname}"
out_gvcfs="${work_dir}/BAM/${jobname}"
out_metrics="${work_dir}/BAM/${jobname}"
#####

[ ! -d $out_bams ] && mkdir -p $out_bams
[ ! -d $out_gvcfs ] && mkdir -p $out_gvcfs
[ ! -d $out_metrics ] && mkdir -p $out_metrics
log_dir="logs/${current_date}_${jobname}"
[ ! -d $log_dir ] && mkdir -p $log_dir
####################################################################################################

for sample_name in "${sample_names[@]}"; do
    fastq_r1="${fsq_dir}/${sample_name}_R1.fastq.gz"
    fastq_r2="${fsq_dir}/${sample_name}_R2.fastq.gz"
    out_bam="${out_bams}/${sample_name}.bam"
    out_gvcf="${out_gvcfs}/${sample_name}.g.vcf"
    out_recal="${out_metrics}/${sample_name}.report.txt"
    out_duplicate_metrics="${out_metrics}/${sample_name}.dup_metrics.txt"
    
    if [ ! -f "${out_bams}/${sample_name}.bam" ]; then
        RGTAG="@RG\tID:${sample_name}\tLB:lib1\tPL:Illumina\tSM:${sample_name}\tPU:${sample_name}"
        mkdir -p "${work_dir}/BAM/TMP/${sample_name}"
        
        bsub -o "${log_dir}/%J_${sample_name}.out" -e "${log_dir}/%J_${sample_name}.err" -M 500GB -R 'gpuhost rusage[mem=350GB, tmp=40GB]' -G compute-jin810 -q general -gpu "num=2:j_exclusive=no" -a 'docker(nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1)' TMPDIR="${work_dir}/BAM/TMP/${sample_name}" pbrun germline --ref "${PARA_REF}" \
            --in-fq "${fastq_r1}" "${fastq_r2}" "${RGTAG}" \
            --knownSites "${pKS4}" \
            --out-bam "${out_bam}" \
            --out-variants "${out_gvcf}" \
            --out-recal-file "${out_recal}" \
            --out-duplicate-metrics "${out_duplicate_metrics}" \
            --bwa-options='-Y -K 100000000' \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --gvcf \
            --num-gpus 2 \
            --tmp-dir "${work_dir}/BAM/TMP/${sample_name}"
    else
        echo "Skipping ${sample_name} as .bam and .g.vcf files already exist."
        count=$((count+1))
    fi
done
echo "Number of samples already processed: ${count}"


