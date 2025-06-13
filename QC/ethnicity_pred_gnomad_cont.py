#!/usr/bin/env python3

import os
import glob
import pandas as pd
import numpy as np
import ast

def predict_ancestry(gnomad_continental_probs):
    populations = ['afr', 'amr', 'asj', 'eas', 'eur', 'sas']
    
    # Convert string representation of list to actual list
    if isinstance(gnomad_continental_probs, str):
        try:
            # Using ast.literal_eval to safely evaluate the string representation of list
            probs = np.array(ast.literal_eval(gnomad_continental_probs))
        except:
            print(f"Error parsing probabilities: {gnomad_continental_probs}")
            return None, None
    else:
        probs = np.array(gnomad_continental_probs)
    
    max_prob_idx = np.argmax(probs)
    predicted_pop = populations[max_prob_idx]
    confidence = probs[max_prob_idx]
    
    return predicted_pop, confidence

def process_samples(parent_dir):
    results = []
    
    # Walk through all directories in parent_dir
    for base in os.listdir(parent_dir):
        sample_dir = os.path.join(parent_dir, base)
        if not os.path.isdir(sample_dir):
            continue
            
        output_dir = os.path.join(sample_dir, 'output')
        if not os.path.isdir(output_dir):
            print(f"Warning: No output directory found for {base}")
            continue
            
        # Find CSV file
        csv_files = glob.glob(os.path.join(output_dir, '*.csv'))
        if not csv_files:
            print(f"Warning: No CSV file found for {base}")
            continue
            
        csv_file = csv_files[0]
        # print(f"Processing {base} from file: {csv_file}")
        
        try:
            # Read the CSV file using pandas
            df = pd.read_csv(csv_file)
            
            # Get the gnomAD_continental values
            if 'gnomAD_continental' in df.columns:
                gnomad_values = df['gnomAD_continental'].iloc[0]  # Get first row
                ancestry, confidence = predict_ancestry(gnomad_values)
                
                if ancestry is not None and confidence is not None:
                    results.append({
                        'SampleID': base,
                        'Ancestry': ancestry,
                        'Confidence': confidence
                    })
                    print(f"Successfully processed {base}: {ancestry} ({confidence:.4f})")
                
        except Exception as e:
            print(f"Error processing {base}: {str(e)}")
            continue
    
    # Create DataFrame and save to CSV in the parent directory
    if results:
        df = pd.DataFrame(results)
        output_file = os.path.join(parent_dir, 'ancestry_combined_gnomAD_continental.csv')
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        print(f"Processed {len(results)} samples successfully")
    else:
        print("No results found to save")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python predict_ancestry.py /path/to/parent_dir")
        sys.exit(1)
        
    parent_dir = sys.argv[1]
    if not os.path.isdir(parent_dir):
        print(f"Error: Directory {parent_dir} does not exist")
        sys.exit(1)
        
    process_samples(parent_dir)