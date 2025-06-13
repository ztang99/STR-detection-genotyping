

import os
import sys
import sqlite3
import argparse
import re
import json
import math

def filter_reads_to_fasta(sam_file, output_file, mapq_threshold=1, append=False):
    """Filter reads from SAM file with MAPQ >= threshold and save to FASTA"""
    mode = 'a' if append else 'w'
    with open(sam_file, 'r') as sam, open(output_file, mode) as fasta:
        for line in sam:
            if line.startswith('@'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
                
            qname, flag, _, _, mapq, _, _, _, _, seq = fields[:10]
            
            if int(mapq) >= mapq_threshold:
                fasta.write(f">{qname}\n{seq}\n")

def get_roi_coordinates(roi_bed=None, json_file=None, gene=None, motif=None):
    """Get ROI coordinates and chromosome from either bed file or json file"""
    if roi_bed:
        with open(roi_bed, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chr_, start, end, bed_gene, bed_motif = fields[:5]
                    if bed_gene == gene and bed_motif == motif:
                        return chr_, int(start), int(end)
    elif json_file:
        with open(json_file, 'r') as f:
            data = json.load(f)
            for entry in data:
                if entry['gene'] == gene and entry['motif'] == motif:
                    return entry['chrom'], entry['start'], entry['end']
    
    raise ValueError(f"Could not find ROI coordinates for {gene}_{motif}")

def get_read_length_from_cigar(cigar):
    """Calculate read length from CIGAR string"""
    if not cigar:
        return 0
    
    numbers = re.split('[MIDNSHP=X]', cigar)
    numbers = [int(x) for x in numbers if x]
    
    return sum(numbers)


def normalize_pattern(pattern):
    """Get canonical form of a pattern to avoid rotational duplicates"""
    # Check if pattern is single character repeating
    if len(set(pattern)) == 1:
        return None
    rotations = {pattern[i:] + pattern[:i] for i in range(len(pattern))}
    return min(rotations)

def find_other_repeats(sequence, motif_len, target_motif, allowed_patterns):
    """Find other repetitive patterns of length motif_len ± 1"""
    patterns = set()
    normalized_target = normalize_pattern(target_motif)
    allowed_patterns = set(allowed_patterns or [])
    
    for pattern_len in [motif_len-1, motif_len, motif_len+1]:
        if pattern_len < 3:
            continue
            
        for i in range(len(sequence) - 3*pattern_len + 1):
            potential_pattern = sequence[i:i+pattern_len]
            if (sequence[i:i+pattern_len] == potential_pattern and
                sequence[i+pattern_len:i+2*pattern_len] == potential_pattern and
                sequence[i+2*pattern_len:i+3*pattern_len] == potential_pattern):
                
                normalized = normalize_pattern(potential_pattern)
                if normalized and normalized != normalized_target and normalized not in allowed_patterns:
                    patterns.add(normalized)
    
    return list(patterns)

def find_max_str_length(sequence, motif, allowed_patterns=None):
    """Find the maximum STR length with 1bp tolerance and detect other repeats"""
    motif_len = len(motif)
    motif_variants = {motif[i:] + motif[:i] for i in range(motif_len)}
    
    if allowed_patterns is None:
        allowed_patterns = []
    other_patterns = find_other_repeats(sequence, motif_len, motif, allowed_patterns)
   
    def hamming_distance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
   
    def find_first_occurrence():
        """Find the first occurrence of any motif variant"""
        for i in range(len(sequence) - motif_len + 1):
            current = sequence[i:i+motif_len]
            if any(hamming_distance(current, variant) == 0 for variant in motif_variants):
                return i
        return -1
    
    start_pos = find_first_occurrence()
    if start_pos == -1:
        return 0
   
    current_pos = start_pos
    max_length = 0
    
    while current_pos <= len(sequence) - motif_len:
        window = sequence[current_pos:current_pos + motif_len]
        if any(hamming_distance(window, variant) <= 1 for variant in motif_variants):
            current_length = motif_len
            next_pos = current_pos + motif_len
            
            while next_pos <= len(sequence) - motif_len:
                next_window = sequence[next_pos:next_pos + motif_len]
                if any(hamming_distance(next_window, variant) <= 1 for variant in motif_variants):
                    current_length += motif_len
                    next_pos += motif_len
                else:
                    break
            
            max_length = max(max_length, current_length)
            current_pos = next_pos
        else:
            current_pos += 1
           
    return max_length//2 if other_patterns else max_length, other_patterns

def query_sequences(db_path, sample_name, motif, threshold_len, roi_chr, roi_start, roi_end, allowed_patterns=None):
    """Query STR sequences to further validate carriers"""
    try:
        table_name = os.path.splitext(os.path.basename(db_path))[0]
        conn = sqlite3.connect(db_path)
        curr = conn.cursor()

        # Debug: Print ROI info and check database content
        print(f"Debug: ROI coordinates - {roi_chr}:{roi_start}-{roi_end}")
        minn = math.floor(roi_start / 1000) * 1000
        maxx = math.ceil(roi_end / 1000) * 1000
        print(f"Debug: Query range - {minn}-{maxx}")
        
        # Verify database has relevant data
        curr.execute(f"SELECT COUNT(*) FROM {table_name} WHERE sample_name = ?", (sample_name,))
        total_samples = curr.fetchone()[0]
        print(f"Debug: Records for sample {sample_name}: {total_samples}")
        
        curr.execute(f"SELECT COUNT(*) FROM {table_name} WHERE sample_name = ? AND top_N_blat_results IS NOT NULL", (sample_name,))
        blat_samples = curr.fetchone()[0]
        print(f"Debug: Records with BLAT results for sample {sample_name}: {blat_samples}")
        
        # If no BLAT results, early exit
        if blat_samples == 0:
            print(f"Warning: No BLAT results found for sample {sample_name}. Check if BLAT results were properly loaded.")
            return sample_name, 0, 0, 0, 0, 0, 'None'

        # Get a few examples to verify data format
        curr.execute(f"SELECT qname, top_N_blat_results FROM {table_name} WHERE sample_name = ? AND top_N_blat_results IS NOT NULL LIMIT 3", (sample_name,))
        examples = curr.fetchall()
        for ex in examples:
            print(f"Debug: Example BLAT result - {ex[0]}: {ex[1]}")

        # Original query
        query = f'''
        SELECT seq, top_N_blat_results, pos, cigar FROM {table_name}
        WHERE sample_name = ?
        AND mapq >= 1
        AND pos BETWEEN ? AND ?;
        '''
        
        curr.execute(query, (sample_name, minn, maxx))
        results = curr.fetchall()
        print(f"Debug: Query returned {len(results)} results")
        
        total_reads = 0
        read_count = 0
        max_str_length = 0
        total_str_length = 0 
        all_other_patterns = set()
        
        # Debug counters
        no_blat_count = 0
        format_error_count = 0
        region_mismatch_count = 0
        
        for row in results:
            sequence, top_N_blat, pos, cigar = row
            
            # Skip if no BLAT results
            if not top_N_blat:
                no_blat_count += 1
                continue
                
            try:
                blat_results = top_N_blat.split(';')
                if not blat_results:
                    continue
                    
                first_match = blat_results[0].split(':')
                if len(first_match) != 6:
                    format_error_count += 1
                    continue
                    
                index, score, chr_, start, end, strand = first_match
                
                # Safe conversion
                try:
                    start = int(float(start))
                    end = int(float(end))
                except (ValueError, TypeError):
                    format_error_count += 1
                    continue
                
                # Region filter
                if (index == '1' and chr_ == roi_chr and 
                    minn <= start <= maxx and minn <= end <= maxx):
                    total_reads += 1
                    
                    read_len = get_read_length_from_cigar(cigar)
                    if threshold_len is None:
                        threshold_len = max(read_len * 0.1, 1)
                    
                    # Find STR length
                    if allowed_patterns is None:
                        allowed_patterns = []
                    result = find_max_str_length(sequence, motif, allowed_patterns)
                    if isinstance(result, tuple):
                        str_length, other_patterns = result
                        all_other_patterns.update(other_patterns)
                    else:
                        str_length = result
                        other_patterns = []
                    
                    max_str_length = max(max_str_length, str_length)
                    total_str_length += str_length
                    
                    if str_length >= threshold_len:
                        read_count += 1
                else:
                    region_mismatch_count += 1
            except Exception as e:
                print(f"Error processing row: {e}")
                continue
        
        # Print debugging stats
        print(f"Debug: Reads with no BLAT results: {no_blat_count}")
        print(f"Debug: Reads with format errors: {format_error_count}")
        print(f"Debug: Reads outside region: {region_mismatch_count}")
        print(f"Debug: Total reads passing filters: {total_reads}")
        print(f"Debug: Reads above threshold: {read_count}")
        
        # Prevent division by zero
        percentage_reads = (read_count / total_reads * 100) if total_reads > 0 else 0
        mean_str_length = (total_str_length / total_reads) if total_reads > 0 else 0
        other_patterns_str = ';'.join(sorted(all_other_patterns)) if all_other_patterns else 'None'
        
        return sample_name, total_reads, read_count, percentage_reads, max_str_length, mean_str_length, other_patterns_str

    except sqlite3.Error as e:
        print(f"Database error: {e}")
        return sample_name, 0, 0, 0, 0, 0, 'None'  
    
    finally:
        if conn:
            conn.close()


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "filter_reads_to_fasta":
        # Create a parser for filter_reads_to_fasta
        filter_parser = argparse.ArgumentParser(description='Filter reads from SAM file to FASTA')
        filter_parser.add_argument('--sam-file', required=True, help='Path to SAM file')
        filter_parser.add_argument('--output-file', required=True, help='Output FASTA file path')
        filter_parser.add_argument('--mapq-threshold', type=int, default=1, help='MAPQ threshold (default: 1)')
        filter_parser.add_argument('--append', type=lambda x: (x.lower() == 'true'), default=False, 
                                  help='Append to existing file if true, overwrite if false (default: false)')
        
        # Parse only the remaining arguments (excluding "filter_reads_to_fasta")
        filter_args = filter_parser.parse_args(sys.argv[2:])
        
        # Call the filter_reads_to_fasta function
        filter_reads_to_fasta(
            sam_file=filter_args.sam_file,
            output_file=filter_args.output_file,
            mapq_threshold=filter_args.mapq_threshold,
            append=filter_args.append
        )
        return

    parser = argparse.ArgumentParser(description='Query and analyze STR sequences')
    parser.add_argument('--db-path', help='Path to SQLite database')
    parser.add_argument('--gene', help='Gene of interest')
    parser.add_argument('--motif', help='Query motif')
    parser.add_argument('--threshold-len', type=int, help='Minimum length threshold (default: 10% of read length)')
    parser.add_argument('--allowed-patterns', nargs='+', help='List of patterns to ignore when detecting other repeats') ## added
    parser.add_argument('--roi-bed', help='Path to ROI bed file (optional)')
    parser.add_argument('--json-file', help='Path to json file with ROI coordinates (optional)')
    parser.add_argument('--output-file', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    if not args.roi_bed and not args.json_file:
        parser.error("Either --roi-bed or --json-file must be provided")
    
    roi_chr, roi_start, roi_end = get_roi_coordinates(
        roi_bed=args.roi_bed,
        json_file=args.json_file,
        gene=args.gene,
        motif=args.motif
    )
    
    table_name = os.path.splitext(os.path.basename(args.db_path))[0]
    conn = sqlite3.connect(args.db_path)
    curr = conn.cursor()
    curr.execute(f"SELECT DISTINCT sample_name FROM {table_name}")
    sample_names = [row[0] for row in curr.fetchall()]
    conn.close()
    # print(args.allowed_patterns)

    all_results = []
    for sample_name in sample_names:
        result = query_sequences(args.db_path, sample_name, args.motif, args.threshold_len, roi_chr, roi_start, roi_end, args.allowed_patterns)
        all_results.append(result)
    
    # with open(args.output_file, 'w') as file:
    #     file.write("sample_name\ttotal_read_num\t"
    #             "num_reads_above_threshold\t"
    #             "percent_reads\t"
    #             "max_str_length\t"
    #             "mean_str_length\n")
    #     for result in all_results:
    #         file.write(f"{result[0]}\t{result[1]}\t{result[2]}\t"
    #                 f"{result[3]:.2f}\t{result[4]}\t"
    #                 f"{result[5]:.2f}\n")
    print(all_results[:10])
    with open(args.output_file, 'w') as file:
        file.write("sample_name\ttotal_read_num\t"
                "num_reads_above_threshold\t"
                "percent_reads\t"
                "max_str_length\t"
                "mean_str_length\t"
                "other_detected_patterns\n")
        for result in all_results:
            if result[2] == 0 and result[3] == 0:
                continue
            file.write(f"{result[0]}\t{result[1]}\t{result[2]}\t"
                    f"{result[3]:.2f}\t{result[4]}\t"
                    f"{result[5]:.2f}\t{result[6]}\n")


if __name__ == '__main__':
    main()
    
    
## 03072025 archived ##
# def query_sequences(db_path, sample_name, motif, threshold_len, roi_chr, roi_start, roi_end, allowed_patterns=None):
#     """Query STR sequences to further validate carriers"""
#     try:
#         table_name = os.path.splitext(os.path.basename(db_path))[0]
#         conn = sqlite3.connect(db_path)
#         curr = conn.cursor()

#         minn = math.floor(roi_start / 1000) * 1000
#         maxx = math.ceil(roi_end / 1000) * 1000

#         query = f'''
#         SELECT seq, top_N_blat_results, pos, cigar FROM {table_name}
#         WHERE sample_name = ?
#         AND mapq >= 1
#         AND pos BETWEEN ? AND ?;
#         '''
        
#         curr.execute(query, (sample_name, minn, maxx))
#         results = curr.fetchall()
        
#         total_reads = 0
#         read_count = 0
#         max_str_length = 0
#         total_str_length = 0 
        
#         all_other_patterns = set() ## added
        
#         for row in results:
#             sequence, top_N_blat, pos, cigar = row
            
#             if top_N_blat:
#                 blat_results = top_N_blat.split(';')
#                 if blat_results:
#                     first_match = blat_results[0].split(':')
#                     if len(first_match) == 6:
#                         index, score, chr_, start, end, strand = first_match
#                         # start, end = int(start), int(end)
#                         start = int(float(start))
#                         end = int(float(end))
                        
#                         if (index == '1' and chr_ == roi_chr and 
#                             minn <= start <= maxx and minn <= end <= maxx):
#                             total_reads += 1
                            
#                             read_len = get_read_length_from_cigar(cigar)
#                             if threshold_len is None:
#                                 threshold_len = max(read_len * 0.1, 1)
                            
#                             # str_length = find_max_str_length(sequence, motif)
#                             # if allowed_patterns is None:
#                             #     allowed_patterns = []
#                             # str_length, other_patterns = find_max_str_length(sequence, motif, allowed_patterns)
#                             # all_other_patterns.update(other_patterns)
                            
#                             if allowed_patterns is None:
#                                 allowed_patterns = []
#                             result = find_max_str_length(sequence, motif, allowed_patterns)
#                             if isinstance(result, tuple):
#                                 str_length, other_patterns = result
#                                 all_other_patterns.update(other_patterns)
#                             else:
#                                 str_length = result
#                                 other_patterns = []
#                             # print(">>>")
#                             # print(sequence)
#                             # print(str_length)
#                             # print(other_patterns)
                            
#                             max_str_length = max(max_str_length, str_length)
#                             total_str_length += str_length
                            
#                             if str_length >= threshold_len:
#                                 read_count += 1

#         percentage_reads = (read_count / total_reads * 100) if total_reads > 0 else 0
#         mean_str_length = (total_str_length / total_reads) if total_reads > 0 else 0
#         other_patterns_str = ';'.join(sorted(all_other_patterns)) if all_other_patterns else 'None'
#         print(sample_name, total_reads, read_count, percentage_reads, max_str_length, mean_str_length, other_patterns_str)

#         return sample_name, total_reads, read_count, percentage_reads, max_str_length, mean_str_length, other_patterns_str

#     except sqlite3.Error as e:
#         print(f"Database error: {e}")
#         return sample_name, 0, 0, 0, 0, 0, 'None'  # Add missing values for new columns
    
#     finally:
#         if conn:
#             conn.close()

    

## 02072025 archived ##

# def find_max_str_length(sequence, motif):
#     """Find the maximum STR length with 1bp tolerance per extension"""
#     motif_len = len(motif)
#     motif_variants = {motif[i:] + motif[:i] for i in range(motif_len)}
    
#     def hamming_distance(s1, s2):
#         """Calculate Hamming distance between two strings of equal length"""
#         return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
#     def find_first_occurrence():
#         """Find the first occurrence of any motif variant"""
#         for i in range(len(sequence) - motif_len + 1):
#             current = sequence[i:i+motif_len]
#             if any(hamming_distance(current, variant) == 0 for variant in motif_variants):
#                 return i
#         return -1
    
#     start_pos = find_first_occurrence()
#     if start_pos == -1:
#         return 0
        
#     current_pos = start_pos
#     total_length = 0
    
#     while current_pos <= len(sequence) - motif_len:
#         window = sequence[current_pos:(current_pos + motif_len)]
#         min_distance = min(hamming_distance(window, variant) for variant in motif_variants)
        
#         if min_distance <= 1:
#             total_length = current_pos + motif_len - start_pos
#             current_pos += 1
#         else:
#             break
            
#     return total_length


## 01082025 archived ##
# def get_roi_coordinates(roi_bed=None, json_file=None, gene=None, motif=None):
#     """Get ROI coordinates from either bed file or json file"""
#     if roi_bed:
#         with open(roi_bed, 'r') as f:
#             for line in f:
#                 fields = line.strip().split('\t')
#                 if len(fields) >= 5:
#                     _, start, end, bed_gene, bed_motif = fields[:5]
#                     if bed_gene == gene and bed_motif == motif:
#                         return int(start), int(end)
#     elif json_file:
#         with open(json_file, 'r') as f:
#             data = json.load(f)
#             for entry in data:
#                 if entry['gene'] == gene and entry['motif'] == motif:
#                     return entry['start'], entry['end']
    
#     raise ValueError(f"Could not find ROI coordinates for {gene}_{motif}")

# def query_sequences(db_path, sample_name, subsequence, N, roi_start, roi_end):
#     """
#     Queries the database for sequences with MAPQ cutoff and do NOT have secondary alignments
#     that belong to the given sample. Only considers reads within the ROI.
    
#     Parameters:
#     - db_path (str): Path to SQLite database
#     - sample_name (str): The name of the sample to filter by
#     - subsequence (str): The subsequence to look for within the sequence
#     - N (int): The minimum number of occurrences of the subsequence
#     - roi_start (int): Start position of the ROI
#     - roi_end (int): End position of the ROI
    
#     Returns:
#     - tuple: (sample_name, total_reads, reads_with_N_occurrences, percentage, total_occurrences)
#     """
#     try:
#         table_name = os.path.splitext(os.path.basename(db_path))[0]
#         conn = sqlite3.connect(db_path)
#         curr = conn.cursor()

#         pos_margin = 151 - len(subsequence)
#         pos_lower_bound = roi_start - pos_margin
#         pos_upper_bound = roi_end + pos_margin

#         query = f'''
#         SELECT seq, top_N_blat_results, pos FROM {table_name}
#         WHERE sample_name = ?
#         AND mapq >= 1
#         AND pos BETWEEN ? AND ?;
#         '''
        
#         curr.execute(query, (sample_name, pos_lower_bound, pos_upper_bound))
#         results = curr.fetchall()
        
#         total_reads = len(results)
#         read_count = 0
#         subsequence_count = 0
        
#         for row in results:
#             sequence = row[0]
#             pos = row[2]

#             top_N_blat = row[1]
#             if top_N_blat:
#                 blat_results = top_N_blat.split(';')
#                 if blat_results:
#                     first_match = blat_results[0].split(':')
#                     if len(first_match) == 6:
#                         index, score, chr_, start, end, strand = first_match
                        
#                         if index == '1' and chr_ == 'chr4' and 39348000 <= int(start) <= 39349000 and 39348000 <= int(end) <= 39349000:
#                             total_reads += 1
#                             occurrences = len(re.findall(re.escape(subsequence), sequence))
#                             subsequence_count += occurrences
#                             if occurrences >= N:
#                                 read_count += 1

#         percentage_reads = (read_count / total_reads) * 100 if total_reads > 0 else 0
#         return sample_name, total_reads, read_count, percentage_reads, subsequence_count

#     except sqlite3.Error as e:
#         print(f"Database error: {e}")
#         return sample_name, 0, 0, 0, 0
    
#     finally:
#         if conn:
#             conn.close()