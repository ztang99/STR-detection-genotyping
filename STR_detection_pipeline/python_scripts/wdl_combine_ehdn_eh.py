import argparse
import os
import re
import pandas as pd
import json
from wdl_str_motif import STRMotif

def clean_sample_name(sample_name):
    """Standardize and clean sample names by removing file extensions and invalid chars."""
    name = os.path.basename(sample_name)
    name = name.replace('.bam', '').replace('.cram', '')
    return name.replace('^', '_')

def load_bam_paths(bams_file):
    """Load BAM paths and create a mapping of cleaned sample names to BAM paths."""
    bam_mapping = {}
    with open(bams_file, 'r') as f:
        for line in f:
            bam_path = line.strip()
            if bam_path:
                cleaned_name = clean_sample_name(bam_path)
                bam_mapping[cleaned_name] = bam_path
    return bam_mapping

def get_sample_sets(ehdn_data, eh_data, gene, motif):
    """Extract and clean sample sets from EHdn and EH data for a given gene/motif."""
    # Get EHdn samples
    ehdn_matches = ehdn_data[
        (ehdn_data['gene'] == gene) & 
        (ehdn_data['motif'] == motif)
    ]
    ehdn_samples = set()
    for raw_data in ehdn_matches['raw_data']:
        samples = [clean_sample_name(s.split(':')[0]) for s in str(raw_data).split(',')]
        ehdn_samples.update(samples)
    
    # Get EH samples
    eh_samples = set(eh_data[
        (eh_data['REPID'] == gene) & 
        (eh_data['RU'] == motif)
    ]['SampleID'])

    return ehdn_samples, eh_samples

def check_sample_overlap(ehdn_samples, eh_samples, min_overlap_percent):
    """Check if sample overlap meets the minimum percentage requirement."""
    overlap_samples = ehdn_samples & eh_samples
    total_samples = ehdn_samples | eh_samples
    
    if not total_samples:
        return False, set()
    
    overlap_percent = len(overlap_samples) * 100 / len(total_samples)
    return overlap_percent >= min_overlap_percent, total_samples

def create_str_motif(gene, motif, row, total_samples, bam_mapping):
    """Create and populate a STR motif object."""
    str_motif = STRMotif(
        gene=gene,
        motif=motif,
        start=int(row['start']),
        end=int(row['end']),
        chrom=str(row['chr'])
    )

    missing_samples = []
    for sample in total_samples:
        if sample in bam_mapping:
            str_motif.add_carrier(bam_mapping[sample])
        else:
            missing_samples.append(sample)
    
    if missing_samples:
        print(f"Warning: could not find BAM for samples in {gene}: {', '.join(missing_samples)}")
    
    return str_motif

def load_ehdn_results(ehdn_file):
    """Load EHdn results from CSV file."""
    return pd.read_csv(ehdn_file)

def process_info_field(df, info_str):
    """Process INFO field from VCF-like format into separate columns."""
    info_fields = set()
    for info in info_str.str.split(';'):
        for field in info:
            field_name = field.split('=')[0]
            info_fields.add(field_name)
    
    info_data = {field: [] for field in info_fields}
    
    for info in info_str:
        current_values = {field: '' for field in info_fields}
        for item in info.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                current_values[key] = value
        for field in info_fields:
            info_data[field].append(current_values[field])
    
    return info_data

def load_eh_results(eh_results_file):
    """Load and process EH results file."""
    df = pd.read_csv(eh_results_file)
    
    # Process INFO field
    info_data = process_info_field(df, df['INFO'])
    
    # Process FORMAT field
    format_cols = df.iloc[0]['FORMAT'].split(':')
    variant_data = df['VARIANTS'].str.split(':', expand=True)
    
    # Construct output dataframe
    new_df = pd.DataFrame()
    
    # Add basic columns
    basic_cols = ['SampleID', 'CHROM', 'POS', 'ID', 'REF', 'ALT']
    for col in basic_cols:
        new_df[col] = df[col].astype(str)
    
    # Add INFO columns
    for field, values in info_data.items():
        new_df[field] = values
    
    # Add FORMAT columns
    for i, col in enumerate(format_cols):
        new_df[col] = variant_data[i].astype(str)
    
    # Clean up data
    new_df = new_df.replace('nan', '')
    
    return new_df

def filter_motifs_by_evidence(motifs, ehdn_data, eh_data, bam_mapping, min_overlap_percent):
    """Filter motifs by requiring evidence from both EHdn and EH."""
    filtered_motifs = []
    
    for motif in motifs:
        ehdn_samples, eh_samples = get_sample_sets(ehdn_data, eh_data, motif.gene, motif.motif)
        passes_overlap, total_samples = check_sample_overlap(ehdn_samples, eh_samples, min_overlap_percent)
        
        if passes_overlap:
            # Use the region info from original motif but update carriers based on evidence
            str_motif = STRMotif(
                gene=motif.gene,
                motif=motif.motif,
                start=motif.start,
                end=motif.end,
                chrom=motif.chrom
            )
            
            for sample in total_samples:
                if sample in bam_mapping:
                    str_motif.add_carrier(bam_mapping[sample])
            
            filtered_motifs.append(str_motif)
    
    return filtered_motifs

def identify_motifs_from_results(ehdn_data, eh_data, bam_mapping, min_overlap_percent):
    """Identify and create STR motifs from EHdn and EH results."""
    # Filter by RepeatMasker
    passing_mask = ehdn_data.apply(check_repeatmasker_motif, axis=1)
    filtered_data = ehdn_data[passing_mask].copy()

    motifs = []
    for _, row in filtered_data.iterrows():
        gene = row['gene']
        motif_seq = row['motif']
        
        ehdn_samples, eh_samples = get_sample_sets(ehdn_data, eh_data, gene, motif_seq)
        passes_overlap, total_samples = check_sample_overlap(ehdn_samples, eh_samples, min_overlap_percent)
        
        if passes_overlap:
            str_motif = create_str_motif(gene, motif_seq, row, total_samples, bam_mapping)
            motifs.append(str_motif)
    
    return motifs

def load_motifs_from_bed(bed_file):
    """Load STR motifs from input ROI bed file."""
    motifs = []
    if bed_file and os.path.exists(bed_file):
        with open(bed_file, 'r') as f:
            for line in f:
                if line.strip():  # Skip empty lines
                    try:
                        motif = STRMotif.from_bed_line(line)
                        motifs.append(motif)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Skipping malformed bed line: {line.strip()}, Error: {e}")
    return motifs

def check_repeatmasker_motif(row):
    """Check if the motif matches RepeatMasker annotation."""
    if pd.isna(row['RepeatMasker_ID']):
        return True
        
    motif = row['motif']
    normalized_motifs = set(motif[i:] + motif[:i] for i in range(len(motif)))
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp = ''.join(complement[base] for base in reversed(motif))
    normalized_motifs.update(rev_comp[i:] + rev_comp[:i] for i in range(len(rev_comp)))
    
    rm_annotation = str(row['RepeatMasker_ID'])
    return not any(nm in rm_annotation for nm in normalized_motifs)

def main():
    parser = argparse.ArgumentParser(description='Combine EHdn and EH results')
    parser.add_argument('--ehdn-results', required=True, help='EHdn result file')
    parser.add_argument('--eh-results', required=True, help='EH results file')
    parser.add_argument('--roi-bed', help='Optional BED file with regions of interest')
    parser.add_argument('--bams', required=True, help='File storing BAM files paths')
    parser.add_argument('--min-overlap-percent', type=float, default=10,
                       help='Minimum percentage of samples detected by both EHdn and EH')
    parser.add_argument('--output-file', required=True, help='JSON output file for the motifs')
    args = parser.parse_args()

    # Load all required data
    bam_mapping = load_bam_paths(args.bams)
    ehdn_data = load_ehdn_results(args.ehdn_results)
    eh_data = load_eh_results(args.eh_results)
    
    # Get motifs either from ROI bed or from results
    if args.roi_bed:
        initial_motifs = load_motifs_from_bed(args.roi_bed)
        motifs = filter_motifs_by_evidence(
            initial_motifs, ehdn_data, eh_data, bam_mapping, args.min_overlap_percent
        )
    else:
        motifs = identify_motifs_from_results(
            ehdn_data, eh_data, bam_mapping, args.min_overlap_percent
        )

    # Print summary
    print(f"\nFound {len(motifs)} valid STR motifs")
    for motif in motifs:
        print(f"\nMotif: {motif}")
        print(f"Carriers: {len(motif.carriers)}")

    # Save results
    if args.output_file:
        with open(args.output_file, 'w') as f:
            json.dump([m.to_dict() for m in motifs], f, indent=4)

if __name__ == '__main__':
    main()
    
    
    

### archived 12/17 ###
# def find_bam_path(sample_id, bams):
#     """Search for BAM file matching the sampleID in all 3 batches sub-directories."""
#     batch_dirs = [d for d in os.listdir(bam_dir) if os.path.isdir(os.path.join(bam_dir, d)) and 'bams' in d]
    
#     for bd in batch_dirs:
#         bd_path = os.path.join(bam_dir, bd)

#         patterns = [
#             f"TWHJ-PNRR-{sample_id}.bam", 
#             f"TWHJ-PNRR-{sample_id}-{sample_id}.bam",
#             f"PNRR-{sample_id}.bam"
#         ]

#         for pp in patterns:
#             potential_path = os.path.join(bd_path, pp)
#             if os.path.exists(potential_path):
#                 return potential_path
    
#     return None

# def load_ehdn_results(ehdn_files):
#     """Load and combine EHdn results from multiple files."""
#     dfs = []
#     for file in ehdn_files:
#         if os.path.exists(file):
#             df = pd.read_csv(file)
#             dfs.append(df)
#     return pd.concat(dfs).drop_duplicates().reset_index(drop=True) if dfs else pd.DataFrame()

# def identify_motifs_from_results(ehdn_data, eh_data, bam_dir, min_overlap_percent=50):
#     """Identify and create STR motifs from EHdn and EH results."""
#     passing_mask = ehdn_data.apply(check_repeatmasker_motif, axis=1)
#     filtered_data = ehdn_data[passing_mask].copy()

#     motifs = []
#     for _, row in filtered_data.iterrows():
#         gene = row['gene']
#         motif_seq = row['motif']
#         print(gene, motif_seq)
        
#         # ehdn_samples = set([s.split(':')[0] for s in str(row['raw_data']).split(',')])
#         ehdn_samples = set([re.sub(r'[^A-Za-z0-9_]', '_', s.split(':')[0]) for s in str(row['raw_data']).split(',')])
#         # print(ehdn_samples)
#         eh_samples = set(eh_data[(eh_data['REPID'] == gene) & (eh_data['RU'] == motif_seq)]['SampleID'])
#         # print(eh_samples)

#         overlap_samples = ehdn_samples & eh_samples
#         print(overlap_samples)
#         total_samples = ehdn_samples | eh_samples
#         print(total_samples)

#         if total_samples:
#             overlap_percent = len(overlap_samples) * 100 / len(total_samples)
#             if overlap_percent >= min_overlap_percent:
#                 str_motif = STRMotif(
#                     gene=gene,
#                     motif=motif_seq,
#                     start=int(row['start']),
#                     end=int(row['end']),
#                     chrom=str(row['chr'])
#                 )

#                 missing_samples = []
#                 for sample in total_samples:
#                     if sample in bam_mapping:
#                         str_motif.add_carrier(bam_mapping[sample])
#                     else:
#                         missing_samples.append(sample)
                
#                 if missing_samples:
#                     print(f"Warning: could not find BAM for samples in {gene}: {', '.join(missing_samples)}")
                
                
#                 # for ss in total_samples:
#                 #     bam_path = find_bam_path(ss, bam_dir)
#                 #     if bam_path:
#                 #         str_motif.add_carrier(bam_path)
#                 #     else:
#                 #         missing_samples.append(ss)
                
#                 # if missing_samples:
#                 #     print(f"Warning: could not find BAM for the following samples in {gene}: ")
#                 #     print(", ".join(missing_samples))

#                 motifs.append(str_motif)
    
#     return motifs

# def load_eh_results(eh_results_file):
#     """Process EH results file to split FORMAT column into separate columns. """
#     df = pd.read_csv(eh_results_file)
    
#     info_fields = set()
#     for info in df['INFO'].str.split(';'):
#         for field in info:
#             field_name = field.split('=')[0]
#             info_fields.add(field_name)
#     info_fields = sorted(list(info_fields)) # ['END', 'REF', 'REPID', 'RL', 'RU', 'VARID']
    
#     info_data = {}
#     for field in info_fields:
#         info_data[field] = []
        
#     for info in df['INFO']:
#         current_values = {field: '' for field in info_fields}
        
#         for item in info.split(';'):
#             if '=' in item:
#                 key, value = item.split('=', 1)
#                 current_values[key] = value
        
#         for field in info_fields:
#             info_data[field].append(current_values[field])
    
#     format_cols = df.iloc[0]['FORMAT'].split(':')
#     variant_data = df['VARIANTS'].str.split(':', expand=True)
    
#     new_df = pd.DataFrame()
    
#     basic_cols = ['SampleID', 'CHROM', 'POS', 'ID', 'REF', 'ALT']
#     for col in basic_cols:
#         new_df[col] = df[col].astype(str)
    
#     for field in info_fields:
#         new_df[field] = info_data[field]
    
#     for i, col in enumerate(format_cols):
#         new_df[col] = variant_data[i].astype(str)
   
#     for col in new_df.columns:
#         new_df[col] = new_df[col].astype(str)
#         new_df[col] = new_df[col].replace('nan', '')
    
#     # df = df.drop(['FORMAT', 'VARIANTS'], axis=1)

#     # for i, col in enumerate(format_cols):
#     #     df[col] = variant_data[i]
    
#     return new_df
