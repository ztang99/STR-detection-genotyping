#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.g.vcf>"
    exit 1
fi

gvcf="$1"

if [ ! -f "$gvcf" ]; then
    echo "Error: Input file $gvcf does not exist"
    exit 1
fi

if [[ ! "$gvcf" =~ \.g\.vcf$ ]]; then
    echo "Error: Input file must be a .g.vcf file"
    exit 1
fi

mkdir -p test_var_only

base=$(basename "$gvcf" .g.vcf)

bcftools view -i 'GT!="./." && GT!=".|." && GT!="0/0" && GT!="0|0"' "$gvcf" \
    | bgzip -c > "test_var_only/${base}.vcf.gz"

tabix -p vcf "test_var_only/${base}.vcf.gz"

echo "Done with ${base}!"