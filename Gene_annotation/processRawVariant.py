#!/usr/bin/env python3

import os
import csv
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Extract deleterious variants from SIFT annotation files.')
    parser.add_argument('--root-dir', required=True, help='Root directory containing chromosome subfolders')
    parser.add_argument('--output', default='deleterious_variants.csv', help='Output CSV file (default: deleterious_variants.csv)')
    parser.add_argument('--sift-threshold', type=float, default=0.05, help='SIFT score threshold (default: 0.05)')
    parser.add_argument('--include-borderline', action='store_true', help='Include borderline variants (within 0.01 of threshold)')
    parser.add_argument('--gene-list', help='CSV file containing list of genes to filter by (one gene per line)')
    return parser.parse_args()

def is_deleterious(sift_score, sift_prediction, threshold, include_borderline):
    """Determine if a variant is deleterious based on SIFT score and prediction."""
    if sift_score == 'NA' or sift_prediction == 'NA':
        return False
    
    try:
        score = float(sift_score)
        # Strictly deleterious
        if score <= threshold:
            return True
        # Borderline deleterious (optional)
        if include_borderline and score <= (threshold + 0.01):
            return True
    except ValueError:
        return False
    
    return False

def parse_xls_file(xls_file, threshold, include_borderline):
    """Parse XLS file and extract deleterious variants."""
    deleterious_variants = []
    
    with open(xls_file, 'r') as f:
        # Skip header
        header = next(f).strip().split('\t')
        
        # Find important column indices
        chrom_idx = header.index('CHROM')
        pos_idx = header.index('POS')
        ref_idx = header.index('REF_ALLELE')
        alt_idx = header.index('ALT_ALLELE')
        gene_idx = header.index('GENE_NAME')
        region_idx = header.index('REGION')
        variant_type_idx = header.index('VARIANT_TYPE')
        sift_score_idx = header.index('SIFT_SCORE')
        sift_pred_idx = header.index('SIFT_PREDICTION')
        
        for line in f:
            fields = line.strip().split('\t')
            
            # Skip lines with insufficient fields
            if len(fields) <= max(chrom_idx, pos_idx, ref_idx, alt_idx, gene_idx, region_idx, 
                               variant_type_idx, sift_score_idx, sift_pred_idx):
                continue
            
            sift_score = fields[sift_score_idx]
            sift_prediction = fields[sift_pred_idx]
            
            if is_deleterious(sift_score, sift_prediction, threshold, include_borderline):
                variant = {
                    'chrom': fields[chrom_idx],
                    'pos': fields[pos_idx],
                    'ref': fields[ref_idx],
                    'alt': fields[alt_idx],
                    'gene': fields[gene_idx],
                    'region': fields[region_idx],
                    'variant_type': fields[variant_type_idx],
                    'sift_score': sift_score,
                    'sift_prediction': sift_prediction
                }
                
                # Add more fields if you need them from the XLS file
                for i, field in enumerate(fields):
                    if i not in [chrom_idx, pos_idx, ref_idx, alt_idx, gene_idx, region_idx, 
                               variant_type_idx, sift_score_idx, sift_pred_idx]:
                        variant[header[i]] = field
                
                deleterious_variants.append(variant)
    
    return deleterious_variants

def get_sample_genotypes(vcf_file, variant):
    """Extract genotype information for a specific variant from VCF file."""
    chrom = variant['chrom']
    pos = variant['pos']
    ref = variant['ref']
    alt = variant['alt']
    
    genotypes = {}
    
    with open(vcf_file, 'r') as f:
        # Skip header lines
        header_line = ""
        sample_names = []
        
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header_line = line.strip()
                sample_names = header_line.split('\t')[9:]  # Samples start at column 10
                break
        
        if not header_line:
            return genotypes
        
        # Search for the variant
        for line in f:
            fields = line.strip().split('\t')
            
            # Skip lines with insufficient fields
            if len(fields) < 10:  # Need at least FORMAT field + 1 sample
                continue
            
            vcf_chrom = fields[0]
            vcf_pos = fields[1]
            vcf_ref = fields[3]
            vcf_alt = fields[4]
            
            # Match the variant
            if (vcf_chrom == chrom and vcf_pos == pos and 
                vcf_ref == ref and vcf_alt == alt):
                
                format_field = fields[8].split(':')
                gt_idx = format_field.index('GT') if 'GT' in format_field else -1
                
                if gt_idx >= 0:
                    for i, sample in enumerate(sample_names):
                        if i + 9 < len(fields):  # +9 to account for fixed VCF fields
                            sample_data = fields[i + 9].split(':')
                            if gt_idx < len(sample_data):
                                genotype = sample_data[gt_idx]
                                if genotype not in ['./.', '0/0']:  # Only include non-reference and non-missing
                                    genotypes[sample] = genotype
                
                break  # Found the variant, no need to continue
    
    return genotypes

def load_gene_list(gene_list_file):
    """Load gene list from a CSV file."""
    genes = set()
    if gene_list_file and os.path.exists(gene_list_file):
        with open(gene_list_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:  # Skip empty lines
                    genes.add(gene)
        print(f"Loaded {len(genes)} genes from {gene_list_file}")
    return genes

def main():
    args = parse_args()
    
    # Load gene list if provided
    gene_set = set()
    if args.gene_list:
        gene_set = load_gene_list(args.gene_list)
    
    # Open output CSV file
    with open(args.output, 'w', newline='') as outfile:
        csv_writer = csv.writer(outfile)
        
        # Write header
        header = ['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'REGION', 'VARIANT_TYPE', 
                 'SIFT_SCORE', 'SIFT_PREDICTION', 'SAMPLES_GENOTYPES']
        csv_writer.writerow(header)
        
        # Process each chromosome directory
        for chrom_dir in sorted(os.listdir(args.root_dir)):
            chrom_path = os.path.join(args.root_dir, chrom_dir)
            
            # Skip if not a directory or not a chromosome directory
            if not os.path.isdir(chrom_path) or not chrom_dir.startswith('chr'):
                continue
            
            print(f"Processing {chrom_dir}...")
            
            # Find XLS annotation file
            xls_files = [f for f in os.listdir(chrom_path) 
                        if f.endswith('_SIFTannotations.xls')]
            
            if not xls_files:
                print(f"  No SIFT annotation file found in {chrom_dir}")
                continue
            
            xls_file = os.path.join(chrom_path, xls_files[0])
            
            # Find VCF prediction file
            vcf_files = [f for f in os.listdir(chrom_path) 
                        if f.endswith('_SIFTpredictions.vcf')]
            
            if not vcf_files:
                print(f"  No SIFT prediction VCF file found in {chrom_dir}")
                continue
            
            vcf_file = os.path.join(chrom_path, vcf_files[0])
            
            # Parse the XLS file to find deleterious variants
            deleterious_variants = parse_xls_file(
                xls_file, args.sift_threshold, args.include_borderline)
            
            print(f"  Found {len(deleterious_variants)} deleterious variants")
            
            # Process each deleterious variant
            for variant in deleterious_variants:
                # Skip variants if gene list is provided and gene is not in the list
                if gene_set and variant['gene'] not in gene_set:
                    continue
                    
                # Get genotype information from VCF
                genotypes = get_sample_genotypes(vcf_file, variant)
                
                # Format all samples with genotypes as "sample:genotype;sample:genotype;..."
                if genotypes:
                    genotype_str = ";".join([f"{sample}:{genotype}" for sample, genotype in genotypes.items()])
                else:
                    genotype_str = "NA"
                
                # Write a single row for this variant with all samples
                row = [
                    variant['chrom'],
                    variant['pos'],
                    variant['ref'],
                    variant['alt'],
                    variant['gene'],
                    variant['region'],
                    variant['variant_type'],
                    variant['sift_score'],
                    variant['sift_prediction'],
                    genotype_str
                ]
                csv_writer.writerow(row)
    
    print(f"Results written to {args.output}")
    
    if args.gene_list:
        print(f"Results filtered to only include variants in genes from {args.gene_list}")

if __name__ == "__main__":
    main()