#!/bin/bash

### Running SIFT for deleterious gene annotation.

# Setting up snpEff environment
export SNPEFF_HOME=/tmp

# Downloading reference genome databases
java -jar /snpEff/snpEff.jar download -dataDir /tmp hg38
java -jar /snpEff/snpEff.jar download -dataDir /tmp GRCh38.86

# Annotating VCF files with snpEff
java -jar /snpEff/snpEff.jar -v -dataDir /tmp hg38 chr4_IPN_Batch1-3.vcf > Annotated.vcf
java -jar /snpEff/SnpSift.jar varType Annotated.vcf > a2.vcf
java -jar /snpEff/snpEff.jar -v -dataDir /tmp GRCh38.86 chr4_IPN_Batch1-3.vcf > Annotated2.vcf
java -jar /snpEff/SnpSift.jar varType Annotated2.vcf > b2.vcf

# Using a different directory for data
export SNPEFF_HOME=2025-01-21/tmp
java -jar /snpEff/snpEff.jar download -dataDir 2025-01-21/tmp GRCh38.86
java -jar /snpEff/snpEff.jar -v -dataDir 2025-01-21/tmp GRCh38.86 chr4_IPN_Batch1-3.vcf > Annotated2.vcf
java -Xms2g -Xmx8g -jar /snpEff/snpEff.jar -v -dataDir 2025-01-21/tmp GRCh38.86 chr4_IPN_Batch1-3.vcf > Annotated2.vcf

# Filtering variants
java -jar /snpEff/SnpSift.jar filter "(SIFT < 0.05)" Annotated.vcf > a3.vcf
java -jar /snpEff/SnpSift.jar filter "(VARTYPE=DEL)" Annotated2.vcf > a3.vcf

# Checking BCF file
bcftools view -H -i 'N_ALT>1 && TYPE="snp"' chr1_IPN_Batch1-3.bcf | head -n1

# Launching interactive sessions with LSF
bsub -Is -G compute-jin810 -q general-interactive -a 'docker(staphb/snpeff)' /bin/bash
bsub -Is -G compute-jin810 -q general-interactive -a 'docker(elle72/basic:vs5)' /bin/bash

# File processing with awk
awk '{split($2, a, /[pq]/); print a[1], $3, $4}' genes.txt | sort -k1 -V > output_gene_pos.txt
awk '{print $0 > $1"_output.txt"}' output_gene_pos.txt

# Creating BED files
awk '{print $1 "\t" $2 "\t" $3}' gene_list/1_output.txt > BED_FILES/1_regions.bed
cd gene_list
aa=$(find . -name "*output.txt" |  cut -f2 -d'/' | cut -f1 -d'_') 
for i in ${aa[@]}; do awk '{print $1 "\t" $2 "\t" $3}' ${i}_output.txt > ../BED_FILES/${i}_regions.bed ;done

# Creating BED files for all chromosomes
for chr in {1..22} X Y; do
    awk '{print "chr"$1 "\t" $2 "\t" $3}' 2025-01-21/gene_list/${chr}_output.txt > 2025-01-21/BED_FILES/${chr}_regions.bed
done

# Processing VCF files for all chromosomes
for chr in {1..22} X Y; do
    mkdir -p 2025-01-21/VCFs/chr${chr}
    cp GLNEXUS/IPN_Batch1-3/BCFs/chr${chr}/chr${chr}_IPN_Batch1-3.bcf 2025-01-21/TEMP_VCFs/TEMP_chr${chr}_IPN_Batch1-3.bcf
    bcftools index -c 2025-01-21/TEMP_VCFs/TEMP_chr${chr}_IPN_Batch1-3.bcf
    bcftools view -R 2025-01-21/BED_FILES/${chr}_regions.bed -O z -o 2025-01-21/VCFs/chr${chr}/chr${chr}_filtered.vcf.gz 2025-01-21/TEMP_VCFs/TEMP_chr${chr}_IPN_Batch1-3.bcf
done

# Indexing VCF files
for chr in {1..22} X Y; do
    tabix -p vcf 2025-01-21/VCFs/chr${chr}/chr${chr}_filtered.vcf.gz
done

# Creating result directories
for chr in {1..22} X Y; do
    mkdir -p RESULTS/chr${chr}
done

# Converting compressed VCF to uncompressed
for chr in {1..22} X Y; do
    bcftools view -O v -o 2025-01-21/VCFs/chr${chr}/chr${chr}_filtered.vcf 2025-01-21/VCFs/chr${chr}/chr${chr}_filtered.vcf.gz
done

# Running SIFT4G annotation on all chromosomes
for chr in {1..22} X Y; do
    java -jar SIFT4G_Annotator.jar -c -i 2025-01-21/VCFs/chr${chr}/chr${chr}_filtered.vcf -d 2025-01-21/tmp/GRCh38.83.chr -r 2025-01-21/RESULTS/chr${chr} -t
done