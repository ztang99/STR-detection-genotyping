#!/bin/bash
# SNVStory Runner - Consolidated Script

export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /scratch1/fs1/jin810:/scratch1/fs1/jin810 /storage2/fs1/epigenome/Active:/storage2/fs1/epigenome/Active /home/tang.zitian:/home/tang.zitian'

### Select variant (i.e., non 0/0 ones) from gvcf files; bgzip and tabix ###
# cd /storage1/fs1/jin810/Active/Projects/TRE_IPN/INPUT_FILES/UMD_BAM/parents_abv40_gvcfs || exit 1

# for gvcf in *.g.vcf; do
#     bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -R 'rusage[mem=5GB]' -a 'docker(elle72/basic:vszt)' \
#         bash /storage1/fs1/jin810/Active/testing/ztang/code/QC/run_qc/snvstory_selectVariant.sh "$gvcf"
# done

cd /storage1/fs1/jin810/Active/Projects/TRE_IPN/INPUT_FILES/UMD_BAM/parents_abv40_gvcfs/test_var_only || exit 1

### Multi-allelic split ###
# for curr_vcf in *.vcf.gz; do
#     base=$(basename "$curr_vcf" .vcf.gz)
#     bsub -g /tang.zitian/misc -G compute-jin810-t3 -q subscription -sla jin810_t3 -R 'rusage[mem=10GB]'  -a 'docker(elle72/basic:vszt)' \
#         -oo "/storage1/fs1/jin810/Active/testing/ztang/logs/20250221_snvstory/MS_%J.log" \
#         bcftools norm --check-ref s -f /storage1/fs1/jin810/Active/References/GRCh38/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
#         -m -both "${curr_vcf}" -o "${base}_MS.vcf.gz" -O z
# done

### Remove alt sites with NON-REF
# for curr_vcf in *_MS.vcf.gz; do
#     base=$(basename "$curr_vcf" _MS.vcf.gz)
#     bcftools view -i 'ALT!~"<NON_REF>"' -Oz -o "${base}_MS_filtered.vcf.gz" "${base}_MS.vcf.gz"
# done

### SNVstory run ###
# -G compute-jin810-t3 -q subscription -sla jin810_t3

out_dir="/storage1/fs1/jin810/Active/Projects/TRE_IPN/INPUT_FILES/BAM/QC_post_alignment/QC_STRpaper_20250204/UMD_controls/snvstory"

for curr_vcf in *_MS_filtered.vcf.gz; do
    base=$(basename "$curr_vcf" _MS_filtered.vcf.gz)

    out_sample_dir="${out_dir}/${base}_try3"
    mkdir -p $out_sample_dir

    bsub -g /tang.zitian/misc -G compute-jin810 -q general -R 'select[mem>32MB && tmp>32MB] rusage[mem=40GB,tmp=40GB]' \
        -a 'docker(mgibio/snvstory:v1.1-buster)' -oo "/storage1/fs1/jin810/Active/testing/ztang/logs/20250224_snvstory/${base}_%J.log" \
        /opt/conda/bin/python3 -m igm_churchill_ancestry \
        --path "/storage1/fs1/jin810/Active/Projects/TRE_IPN/INPUT_FILES/UMD_BAM/parents_abv40_gvcfs/test_var_only/${curr_vcf}" \
        --resource /storage1/fs1/jin810/Active/testing/elle/Lab_Members/Zitian/2025-02-19/resource/ \
        --output-dir "${out_sample_dir}" \
        --genome-ver 38 \
        --mode WGS
done

