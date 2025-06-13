import os
import sys
import json
import pandas as pd

def create_eh_catalog_entry(row):
    return {
        "LocusId": row["gene"],
        "LocusStructure": f'({row["motif"]})*',
        "ReferenceRegion": f'{row["chr"]}:{row["start"]}-{row["end"]}',
        "VariantType": "Repeat"
    }

def generate_catalog(input_file, output_json):
    df = pd.read_csv(input_file)
    json_data = [create_eh_catalog_entry(row) for _, row in df.iterrows()]
    
    with open(output_json, 'w') as outfile:
        json.dump(json_data, outfile, indent=4)
    print(f"Generated catalog for {len(json_data)} variants in {output_json}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python IPN_generate_EHcatalog.py <output_json> <input_csv_files...>")
        sys.exit(1)

    output_json = sys.argv[1]
    input_file = sys.argv[2]
    
    if os.path.exists(input_file):
        df = pd.read_csv(input_file)
    else:
        print(f"Warning: Input file {input_file} not found")
        sys.exit(1)
    
    generate_catalog(input_file, output_json)
    

if __name__ == "__main__":
    main()