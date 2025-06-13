import sys
import os
import csv
import pandas as pd
import argparse

def process_vcf_file(vcf_path):
    """Process a single VCF file and return the relevant data."""
    results = []
    sample_id = os.path.basename(vcf_path).replace('.vcf', '')
    
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith('#') and 'INREPEAT' in line:
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = fields[1]
                id_ = fields[2]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]
                format_ = fields[8]
                sample_data = fields[9]
                
                results.append([
                    sample_id, chrom, pos, id_, ref, alt,
                    info, format_, sample_data
                ])
    return results

def process_vcf_file_all4gene(vcf_path, gene_name):
    """Process a single VCF file and return the relevant data."""
    results = []
    sample_id = os.path.basename(vcf_path).replace('.vcf', '')
    
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith('#'):
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
                    
                    results.append([
                        sample_id, chrom, pos, id_, ref, alt,
                        info, format_, sample_data
                    ])
    return results

def combine_vcfs(vcf_files, output_file, mode, gene_name):
    """Combine multiple VCF files into a single output file."""
    
    if os.path.exists(output_file):
        os.remove(output_file)
        
    header = ['SampleID', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 
              'INFO', 'FORMAT', 'VARIANTS']
    
    all_results = []
    for vcf_file in vcf_files:
        if os.path.exists(vcf_file):
            if mode == "inrepeat4all": # INREPEAT calls for all loci
                results = process_vcf_file(vcf_file)
            elif mode == "all4gene": # ALL calls for a specific gene loci
                results = process_vcf_file_all4gene(vcf_file, gene_name)
            all_results.extend(results)
        else:
            print(f"Warning: VCF file not found: {vcf_file}")
    
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(header)
        writer.writerows(all_results)

def main():
    parser = argparse.ArgumentParser(description='Process VCF files based on specified mode.')
    parser.add_argument('--output-file', required=True,
                      help='Output file path for the combined results')
    parser.add_argument('--vcf-files', required=True, nargs='+',
                      help='One or more VCF files to process')
    parser.add_argument('--mode', required=True, choices=['inrepeat4all', 'all4gene'],
                      help='Processing mode: "inrepeat4all" or "all4gene"')
    parser.add_argument('--gene-name', default="RFC1",
                      help='Gene name to search for in all4gene mode (default: RFC1)')

    args = parser.parse_args()

    if not args.vcf_files:
        print("Error: No VCF files provided")
        sys.exit(1)
        
    if args.mode == 'all4gene' and not args.gene_name:
        parser.error("--gene-name is required when mode is all4gene")

    combine_vcfs(args.vcf_files, args.output_file, args.mode, args.gene_name)
    print(f"Combined VCF data written to: {args.output_file}")


if __name__ == '__main__':
    main()








# archived 01212025
# def main():
#     if len(sys.argv) < 3:
#         print("Usage: python filter_eh_vcfs.py <output_file> <vcf_files...>")
#         sys.exit(1)

#     output_file = sys.argv[1]
#     vcf_files = sys.argv[2:]

#     if not vcf_files:
#         print("Error: No VCF files provided")
#         sys.exit(1)

#     combine_vcfs(vcf_files, output_file)
#     print(f"Combined VCF data written to: {output_file}")