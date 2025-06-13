import os
import sys
import numpy as np
import sqlite3
import csv
import re
import glob
import pandas as pd
from wdl_addBlatResult2db import parse_psl_file, compute_blat_score

def initialize_database(db_name='IPNstrs.db'):
    """
    Initializes the SQLite database and creates the IPNstrs table if it does not exist.
    """
    conn = sqlite3.connect(db_name)
    curr = conn.cursor()

    curr.execute('''
    CREATE TABLE IF NOT EXISTS IPNstrs (
        qname TEXT,
        flag INTEGER,
        rname TEXT,
        pos INTEGER,
        mapq INTEGER,
        cigar TEXT,
        rnext TEXT,
        pnext INTEGER,
        tlen INTEGER,
        seq TEXT,
        qual TEXT,
        sample_name TEXT,
        gene TEXT,
        case_control TEXT,
        secondAlign INTEGER,
        top_N_blat_results TEXT
    )
    ''')
    
    return conn, curr


def read_exists(curr, qname, flag):
    """
    Checks if a read with the given qname and flag already exists in the database.

    Parameters:
    - curr (sqlite3.Cursor): Cursor object for the SQLite connection.
    - qname (str): The qname of the read to check.
    - flag (int): The flag of the read to check.

    Returns:
    - bool: True if the read with the same qname and flag exists, False otherwise.
    """
    curr.execute("SELECT 1 FROM IPNstrs WHERE qname = ? AND flag = ?", (qname, flag))
    return curr.fetchone() is not None

def sample_exists(curr, sample_name):
    """
    Checks if a sample with the given sample_name already exists in the database.

    Parameters:
    - curr (sqlite3.Cursor): Cursor object for the SQLite connection.
    - sample_name (str): The sample name to check.

    Returns:
    - bool: True if the sample already exists, False otherwise.
    """
    curr.execute("SELECT 1 FROM IPNstrs WHERE sample_name = ?", (sample_name,))
    return curr.fetchone() is not None

def read_in_sample_exists(curr, qname, flag, sample_name):
    """
    Checks if a read with the given qname, flag, and sample_name already exists in the database.

    Parameters:
    - curr (sqlite3.Cursor): Cursor object for the SQLite connection.
    - qname (str): The qname of the read to check.
    - flag (int): The flag of the read to check.
    - sample_name (str): The sample name to check.

    Returns:
    - bool: True if the read with the same qname, flag, and sample_name exists, False otherwise.
    """
    curr.execute("SELECT 1 FROM IPNstrs WHERE qname = ? AND flag = ? AND sample_name = ?", 
                 (qname, flag, sample_name))
    return curr.fetchone() is not None


def process_optional_fields(optional_fields):
    """
    Extracts the sample name from optional fields of a SAM file.
    
    Parameters:
    - optional_fields (list): List of optional field strings from a SAM file.
    Returns:
    - sample_name (str): The extracted sample name or None if not found.
    """
    sample_name = None
    secondAlign = False
    for field in optional_fields:
        if field.startswith("RG:Z:"):
            sample_name = field.split(":")[-1]
            sample_name = sample_name.split("-")[-1]
        if field.startswith("SA:Z:"):
            secondAlign = True
    return sample_name, secondAlign

def parse_sam_file(file_path, curr, gene_of_interest):
    """
    Parses a SAM file and inserts data into the SQLite database.

    Parameters:
    - file_path (str): The path to the SAM file to be parsed.
    - curr (sqlite3.Cursor): Cursor object for the SQLite connection.
    """
    with open(file_path, 'r') as sam_file:
        for line in sam_file:
            if line.startswith('@'):
                continue
                
            fields = line.strip().split('\t')
            optional_fields = fields[11:]
            sample_name, secondAlign = process_optional_fields(optional_fields)
            
            # RG tag is not unique for MAP controls; so reassigning samplename acc. to filename
            cram_file = os.path.basename(file_path)
            if "LINDSEY" in cram_file:
                sample_name = f"{cram_file.split('_')[-1][:-9]}_{cram_file.split('_')[2]}"

            if read_in_sample_exists(curr, fields[0], int(fields[1]), sample_name):
                continue

            # Determine case/control status based on sample name length
            case_control = 'Case' if len(sample_name) == 5 else 'Control'

            curr.execute('''
                INSERT INTO IPNstrs (
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, 
                    seq, qual, sample_name, gene, case_control, secondAlign
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                fields[0], int(fields[1]), fields[2], int(fields[3]), 
                int(fields[4]), fields[5], fields[6], int(fields[7]), 
                int(fields[8]), fields[9], fields[10], sample_name, 
                gene_of_interest, case_control, secondAlign
            ))

    sam_file.close()
    
def parse_directory(directory, curr, gene_of_interest):
    """
    Parses all SAM files within the specified directory and stores the data in the database.

    Parameters:
    - directory (str): The path to the directory containing SAM files.
    - gene_of_interest: gene that we're importing data from.
    """
    sam_files = glob.glob(os.path.join(directory, "*.sam"))
    
    for sfile in sam_files:
        print(f"Importing information from {sfile}...")
        parse_sam_file(sfile, curr, gene_of_interest)

def print_current_db(conn):
    df = pd.read_sql_query("SELECT * FROM IPNstrs", conn)
    print(df)


def main():
    if len(sys.argv) != 5:
        print("Usage: python build_STR_db.py <gene> <sam_files_dir> <blat_results> <output_db>")
        sys.exit(1)
    
    gene = sys.argv[1]
    sam_dir = sys.argv[2]
    blat_results = sys.argv[3]
    output_db = sys.argv[4]
    
    # Initialize database
    conn, curr = initialize_database(output_db)
    
    # Process SAM files
    parse_directory(sam_dir, curr, gene)
    
    # Update case/control status
    update_case_control_based_on_sample_name(conn)
    
    # Process BLAT results if available
    if os.path.exists(blat_results):
        print(f"Processing BLAT results from {blat_results}...")
        parse_psl_file(blat_results, output_db, N=3)
    
    # Print final database state
    print_current_db(conn)
    
    conn.commit()
    conn.close()


    
if __name__ == '__main__':
    main()
