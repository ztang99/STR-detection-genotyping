import os
import sys
import numpy as np
import pandas as pd
import argparse
from scipy.stats import fisher_exact

def compute_fisher_pvalue(caco_raw_data, case_total, control_total):
    """
    Compute Fisher's exact test p-value for case-control comparison.
    """
    case_count = sum("PNRR" in item for item in caco_raw_data.split(','))
    control_count = len(caco_raw_data.split(',')) - case_count
    table = [[case_count, control_count], 
             [case_total - case_count, control_total - control_count]]
    _, p_value = fisher_exact(table)
    return p_value

def read_exdn_otl(otl_data_path):
    unfiltered = pd.read_csv(otl_data_path, sep='\t', header=0, 
                            names=['chr', 'start', 'end', 'motif', 'gene', 'region', 
                                  'top_zscore', 'raw_data', 'all_counts'])
    print(f"EXDN initial output count: {len(unfiltered)}.")
    return unfiltered

def merge_exdn_caco_output(otl_data, exdn_caco_path):
    exdn_caco_out = pd.read_csv(exdn_caco_path, sep='\t', header=0, 
                               names=['chr', 'start', 'end', 'motif', 'gene', 'region', 
                                     'p_val', 'bonf_pval', 'caco_raw_data'])
    exdn_caco_out['chr'] = exdn_caco_out['chr'].str[3:]
    exdn_caco_out = exdn_caco_out[
        (exdn_caco_out['chr'].str.isnumeric()) | 
        (exdn_caco_out['chr'].isin(['X', 'Y']))
    ].copy()
    
    merged = otl_data.merge(exdn_caco_out, on=['chr', 'start', 'end', 'motif', 'gene', 'region'])
    merged = merged.dropna(subset=['caco_raw_data'])
    print(f"Counts after incorporating controls: {len(merged)}.")
    return merged

def filter_chromosomes(in_file):
    in_file['chr'] = in_file['chr'].str[3:]
    out_file = in_file[
        (in_file['chr'].str.isnumeric()) | 
        (in_file['chr'].isin(['X', 'Y']))
    ].copy()
    print(f"Counts after chromosome filtering: {len(out_file)}.")
    return out_file

def filter_motif_lengths(in_data, min_motif_len=2, max_motif_len=6):
    out_data = in_data[
        (in_data['motif'].str.len() <= max_motif_len) & 
        (in_data['motif'].str.len() >= min_motif_len)
    ]
    print(f"Counts after motif length filtering: {len(out_data)}.")
    return out_data

def normalize_motif(motif):
    def reverse_complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement[base] for base in reversed(seq))
    
    def get_cyclic_permutations(seq):
        return {seq[i:] + seq[:i] for i in range(len(seq))}
    
    return get_cyclic_permutations(motif) | get_cyclic_permutations(reverse_complement(motif))

def extract_genes(gene_str):
    """Extract clean gene names from gene string."""
    extracted = []
    genes = gene_str.split(',')
    for g in genes:
        ext_g = g.split('(')[0]
        extracted.append(ext_g)
    return extracted

def gene_in_list(gene_str, gene_list):
    """Check if any gene from the gene string is in gene list."""
    extracted = extract_genes(gene_str)
    return any(gene in gene_list for gene in extracted)

def filter_genes(in_data, gene_list):
    """Filter data to include only genes in the provided list."""
    return in_data[in_data['gene'].apply(lambda x: gene_in_list(x, gene_list))]

def annotate_with_repeatmasker(output_data, repeat_masker_path):
    repeatmasker_data = pd.read_csv(repeat_masker_path, sep='\t', header=None, 
        names=['chr', 'begin', 'end', 'score', 'div', 'del', 'ins', 'sequence', 
               'repeat_begin', 'repeat_end', 'left', 'strand', 'repeat', 
               'class/family', 'repeat_start', 'repeat_finish', 'left2', 'ID'])
    
    # Convert numeric columns
    output_data['start'] = pd.to_numeric(output_data['start'], errors='coerce')
    output_data['end'] = pd.to_numeric(output_data['end'], errors='coerce')
    repeatmasker_data['begin'] = pd.to_numeric(repeatmasker_data['begin'], errors='coerce')
    repeatmasker_data['end'] = pd.to_numeric(repeatmasker_data['end'], errors='coerce')
    
    repeatmasker_data['repeat_cleaned'] = repeatmasker_data['repeat'].str.extract(r'\((.*?)\)n')[0]
    output_data['RepeatMasker_ID'] = None

    for index, row in output_data.iterrows():
        chr_mask = repeatmasker_data['sequence'] == f'chr{row["chr"]}'
        chr_repeatmasker = repeatmasker_data[chr_mask]
        
        normalized_motifs = normalize_motif(row['motif'])
        overlap_mask = (
            ((row['start'] <= chr_repeatmasker['begin']) & 
             (chr_repeatmasker['begin'] <= row['end'])) | 
            ((row['start'] <= chr_repeatmasker['end']) & 
             (chr_repeatmasker['end'] <= row['end']))
        )
        
        matches = chr_repeatmasker['repeat_cleaned'].apply(
            lambda x: any(nm in str(x) for nm in normalized_motifs)
        )
        
        matching_rows = chr_repeatmasker[overlap_mask & matches]
        if not matching_rows.empty:
            annotations = []
            for _, matching_row in matching_rows.iterrows():
                annotation = f"{matching_row['repeat']}; {matching_row['chr']}:{matching_row['begin']}-{matching_row['end']}"
                annotations.append(annotation)
            output_data.at[index, 'RepeatMasker_ID'] = " | ".join(annotations)

    print(f"Counts after RepeatMasker annotation: {len(output_data)}")
    return output_data

def filter_by_gene_list(out_data, gene_file, case_total, control_total, save_path):
    """Filter and save results for a specific gene list."""
    gene_list = pd.read_csv(gene_file, header=None)[0].tolist()
    gene_list_name = os.path.splitext(os.path.basename(gene_file))[0]
    
    filtered = filter_genes(out_data, gene_list)
    # Add Fisher's exact test p-values
    filtered['fisher_p_value'] = filtered['caco_raw_data'].apply(
        lambda x: compute_fisher_pvalue(x, case_total, control_total)
    )
    
    output_file = f'{save_path}/{gene_list_name}.csv'
    filtered.to_csv(output_file, index=False)
    print(f"Saved results for {gene_list_name}: {len(filtered)} entries")
    return output_file, filtered

def main():
    ## 1212 new ver. ##
    parser = argparse.ArgumentParser(description='Combine and filter EHDN results')
    parser.add_argument('--outlier-locus', required=True,
                      help='Path to EHdn outlier locus file')
    parser.add_argument('--casecontrol-locus', required=True,
                      help='Path to EHdn case-control locus file')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for individual gene list files')
    parser.add_argument('--output-file', required=True,
                        help='Path to combined output file')
    parser.add_argument('--repeatmasker-file', required=True,
                      help='Path to RepeatMasker annotation file')
    parser.add_argument('--case-count', type=int, required=True,
                      help='Total number of cases')
    parser.add_argument('--control-count', type=int, required=True,
                      help='Total number of controls')
    parser.add_argument('--gene-list-files', required=True,
                      help='Comma-separated list of gene list files')

    args = parser.parse_args()
    gene_list_files = args.gene_list_files.split(',')
    
    
    motif_len_min = 2
    motif_len_max = 6


    # Process data
    print("Reading EXDN outlier data...")
    otl_data = read_exdn_otl(args.outlier_locus)
    
    print("Filtering chromosomes...")
    filtered_data = filter_chromosomes(otl_data)
    
    print("Filtering motif lengths...")
    filtered_data = filter_motif_lengths(filtered_data, motif_len_min, motif_len_max)
    
    print("Merging with case-control data...")
    merged_data = merge_exdn_caco_output(filtered_data, args.casecontrol_locus)
    
    print("Annotating with RepeatMasker...")
    annotated_data = annotate_with_repeatmasker(merged_data, args.repeatmasker_file)

    # Process each gene list
    output_files = []
    all_dfs = []
    for gene_file in gene_list_files:
        print(f"\nProcessing gene list: {gene_file}")
        output_file, df = filter_by_gene_list(
            annotated_data, gene_file, args.case_count, args.control_count, args.output_dir
        )
        output_files.append(output_file)
        all_dfs.append(df)

    combined_results = pd.concat(all_dfs).drop_duplicates().reset_index(drop=True)

    # Save combined results
    # combined_output = f'{args.output_dir}/EHdn_combined_results.csv'
    # combined_results.to_csv(combined_output, index=False)
    # print(f"\nSaved combined results to: {combined_output}")
    combined_results.to_csv(args.output_file, index=False)
    print(f"\nSaved combined results to: {args.output_file}")

    print("\nAll processing complete!")

if __name__ == '__main__':
    main()
    
    
    













## old ver archived ##
# if len(sys.argv) != 5:
#         print("Usage: python filter_exdn_results.py <outlier_locus> <casecontrol_locus> " +
#               "<output_dir> <repeatmasker_file> <gene_list_files>")
#         sys.exit(1)

#     exdn_otl_path = sys.argv[1]
#     exdn_caco_path = sys.argv[2]
#     output_dir = sys.argv[3]
#     repeat_masker_path = sys.argv[4]
#     gene_list_files = sys.argv[5].split(',')

#     # Constants for Fisher's exact test
#     case_total = 798        #need to get rid of this
#     control_total = 887     #need to get rid of this (or input from outside)

    # os.makedirs(output_dir, exist_ok=True)