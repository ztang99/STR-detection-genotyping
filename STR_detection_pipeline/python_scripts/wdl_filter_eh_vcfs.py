
##############################################################################

# Helper script for filtering INREPEAT calls from ExpansionHunter results.
# Called by 5_CombineEHResult.sh

## author: Zitian Tang
## contact: tang.zitian@wustl.edu

##############################################################################

"""
Usage:
# Process multiple VCF files for a single gene
python vcf_processor.py \
  --gene ATXN1 \
  --vcf-files sample1.vcf sample2.vcf sample3.vcf \
  --output-prefix /path/to/results/EH_combined_result

# Process multiple VCF files for genes listed in a CSV file
python vcf_processor.py \
  --gene-list genes_of_interest.csv \
  --vcf-files *.vcf \
  --output-prefix /home/user/analysis/EH_combined_result
"""

import sys
import os
import csv
import pandas as pd
import argparse

def read_gene_list(gene_csv):
    """Read gene names from a CSV file, one gene per line."""
    genes = []
    with open(gene_csv, 'r') as csv_file:
        for line in csv_file:
            gene = line.strip()
            if gene:
                genes.append(gene)
    return genes

def process_vcf_file_for_gene(vcf_path, gene_name):
    """Process a single VCF file for a specific gene and categorize results."""
    results = {'IRRonly': [], 'SP': [], 'FL': [], 'ALL': []}
    sample_id = os.path.basename(vcf_path).replace('.vcf', '')
    
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith('#') and 'INREPEAT' in line:
                fields = line.strip().split('\t')
                chrom = fields[0]
                info = fields[7]
                
                if f"VARID={gene_name}" in info:
                    pos = fields[1]
                    id_ = fields[2]
                    ref = fields[3]
                    alt = fields[4]
                    format_ = fields[8]
                    sample_data = fields[9]
                    
                    record = [sample_id, chrom, pos, id_, ref, alt, info, format_, sample_data]
                    results['ALL'].append(record)

                    if 'SPANNING' in line:
                        results['SP'].append(record)
                    elif 'FLANKING' in line:
                        results['FL'].append(record)
                    else:
                        results['IRRonly'].append(record)
    
    return results

def combine_vcfs(vcf_files, output_prefix, gene_list):
    """Combine multiple VCF files and generate three categorized output files."""
    header = ['SampleID', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'FORMAT', 'VARIANTS']
    
    all_results = {'IRRonly': [], 'SP': [], 'FL': [], 'ALL': []}
    
    for gene_name in gene_list:
        print(f"Processing gene: {gene_name}")
        for vcf_file in vcf_files:
            if os.path.exists(vcf_file):
                results = process_vcf_file_for_gene(vcf_file, gene_name)
                for category in all_results:
                    all_results[category].extend(results[category])
            else:
                print(f"Warning: VCF file not found: {vcf_file}")
    
    output_files = {
        'IRRonly': f"{output_prefix}_IRRonly.csv",
        'SP': f"{output_prefix}_SP.csv",
        'FL': f"{output_prefix}_FL.csv",
        'ALL': f"{output_prefix}.csv"
    }
    
    for category, filename in output_files.items():
        if os.path.exists(filename):
            os.remove(filename)
        
        with open(filename, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(header)
            writer.writerows(all_results[category])
    
    return output_files

def main():
    parser = argparse.ArgumentParser(description='Process VCF files and categorize INREPEAT variants.')
    parser.add_argument('--output-prefix', required=True,
                      help='Output prefix for the three result files (path and prefix)')
    parser.add_argument('--vcf-files', required=True, nargs='+',
                      help='One or more VCF files to process')
    
    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument('--gene',
                          help='Single gene name to process')
    gene_group.add_argument('--gene-list',
                          help='CSV file with gene names, one per line')

    args = parser.parse_args()

    if not args.vcf_files:
        print("Error: No VCF files provided")
        sys.exit(1)
    
    if args.gene:
        gene_list = [args.gene]
    else:
        if not os.path.exists(args.gene_list):
            parser.error(f"Gene list file not found: {args.gene_list}")
        gene_list = read_gene_list(args.gene_list)
        if not gene_list:
            parser.error("No genes found in the provided CSV file")
    
    output_files = combine_vcfs(args.vcf_files, args.output_prefix, gene_list)
    
    print("VCF processing completed. Output files created:")
    for category, filename in output_files.items():
        print(f"  {category}: {filename}")

if __name__ == '__main__':
    main()