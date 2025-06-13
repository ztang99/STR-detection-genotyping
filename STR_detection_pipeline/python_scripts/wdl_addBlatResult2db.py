
import os
import sys
import sqlite3
import argparse
import pandas as pd
import numpy as np
import csv
import re
import glob
from contextlib import closing

def get_table_name(db_path):
    """Extract table name from database filename"""
    return os.path.splitext(os.path.basename(db_path))[0]

def connect_to_db(db_path):
    """Connect to SQLite database with proper settings for reliability"""
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")  # Use Write-Ahead Logging for better reliability
    conn.execute("PRAGMA synchronous=NORMAL")  # Better performance while maintaining reliability
    return conn

def initialize_database(db_path):
    """Initialize database with table name matching db filename"""
    table_name = get_table_name(db_path)
    
    # Use 'with' statement to ensure proper closing of connection
    with closing(connect_to_db(db_path)) as conn:
        with closing(conn.cursor()) as curr:
            curr.execute(f'''
            CREATE TABLE IF NOT EXISTS {table_name} (
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
                top_N_blat_results TEXT
            )
            ''')
            conn.commit()
    
    return table_name

def read_exists(curr, table_name, qname, flag, sample_name):
    """Check if read exists in sample"""
    curr.execute(f"SELECT 1 FROM {table_name} WHERE qname = ? AND flag = ? AND sample_name = ? LIMIT 1", 
                 (qname, flag, sample_name))
    return curr.fetchone() is not None

def parse_sam_file(file_path, conn, table_name, gene_of_interest):
    """Parses a SAM file and inserts data into the SQLite database"""
    sample_name = os.path.basename(file_path).split('.')[0]
    
    with closing(conn.cursor()) as curr:
        with open(file_path, 'r') as sam_file:
            # Use batch inserts for better performance
            batch_size = 1000
            batch = []
            
            for line in sam_file:
                if line.startswith('@'):
                    continue
                    
                fields = line.strip().split('\t')
                
                qname = fields[0]
                flag = int(fields[1])

                if read_exists(curr, table_name, qname, flag, sample_name):
                    continue

                rname = fields[2]
                pos = int(fields[3])
                mapq = int(fields[4])
                cigar = fields[5]
                rnext = fields[6]
                pnext = int(fields[7])
                tlen = int(fields[8])
                seq = fields[9]
                qual = fields[10]

                batch.append((qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, sample_name, gene_of_interest))
                
                if len(batch) >= batch_size:
                    curr.executemany(f'''
                        INSERT INTO {table_name} (
                            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, sample_name, gene
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', batch)
                    conn.commit()
                    batch = []
            
            # Insert any remaining records
            if batch:
                curr.executemany(f'''
                    INSERT INTO {table_name} (
                        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, sample_name, gene
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', batch)
                conn.commit()
            
def parse_directory(directory, conn, table_name, gene_of_interest):
    """Parses all SAM files within the specified directory and stores the data in the database"""
    sam_files = glob.glob(os.path.join(directory, "*.sam"))
    
    for sfile in sam_files:
        print(f"Importing information from {sfile}...")
        parse_sam_file(sfile, conn, table_name, gene_of_interest)

def compute_blat_score(rows, N):
    """Compute the BLAT score for each row and return the top N rows as a formatted string."""
    score_cols = ['match', 'mis-match', 'rep.match', 'Q gap bases', 'T gap bases']
    
    # Convert columns to numeric, with better error handling
    for col in score_cols:
        rows[col] = pd.to_numeric(rows[col], errors='coerce').fillna(0)

    # Calculate BLAT score
    rows['blat_score'] = rows['match'] - rows['mis-match'] - rows['rep.match'] - rows['Q gap bases'] - rows['T gap bases']
    
    # Get top N rows
    top_rows = rows.nlargest(N, 'blat_score')
    
    # Format results
    results = []
    for i, row in enumerate(top_rows.itertuples(), 1):
        result = f"{i}:{row.blat_score}:{row.T_name}:{row.T_start}:{row.T_end}:{row.strand}"
        results.append(result)
    
    return ';'.join(results)

def add_column_if_not_exists(conn, table_name, column_name, column_type):
    """Add a column to the table if it doesn't already exist"""
    with closing(conn.cursor()) as curr:
        curr.execute(f"PRAGMA table_info({table_name})")
        columns = [info[1] for info in curr.fetchall()]
        
        if column_name not in columns:
            curr.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {column_type}")
            conn.commit()
            print(f"Added column '{column_name}' to table '{table_name}'")
        else:
            print(f"Column '{column_name}' already exists in table '{table_name}'")

def parse_single_psl_file(psl_file, db_path, table_name, N=3):
    """Parse a single PSL file and update the SQLite database with the top N BLAT results for this sample."""
    columns = ['match', 'mis-match', 'rep.match', 'N\'s', 'Q gap count', 'Q gap bases', 
              'T gap count', 'T gap bases', 'strand', 'Q_name', 'Q size', 'Q start', 
              'Q end', 'T_name', 'T size', 'T_start', 'T_end', 'block count', 
              'blockSizes', 'qStarts', 'tStarts']
    
    # Extract sample name from PSL file name
    sample_name = os.path.basename(psl_file).split('.')[0]
    print(f"Processing {sample_name} from {psl_file}...")
    
    try:
        # Read PSL file
        df = pd.read_csv(psl_file, sep='\t', header=None, names=columns, 
                        comment='#', skiprows=5)
        
        if df.empty:
            print(f"Warning: No data found in {psl_file}")
            return
            
        # Group by read name to compute BLAT scores
        unique_qnames = df['Q_name'].unique()
        print(f"Found {len(unique_qnames)} unique read names in {sample_name}")
        
        # Connect to database
        with closing(connect_to_db(db_path)) as conn:
            # Add the top_N_blat_results column if it doesn't exist
            add_column_if_not_exists(conn, table_name, "top_N_blat_results", "TEXT")
            
            # First, get list of reads for this sample to process only those
            with closing(conn.cursor()) as curr:
                curr.execute(f"SELECT qname FROM {table_name} WHERE sample_name = ?", (sample_name,))
                db_qnames = set(row[0] for row in curr.fetchall())
                
                # Find intersection of reads in PSL and database
                qnames_to_process = set(unique_qnames).intersection(db_qnames)
                print(f"Processing {len(qnames_to_process)} reads that exist in both PSL and database")
                
                # Process in batches to improve performance
                batch_size = 100
                batch = []
                processed_count = 0
                
                for qname in qnames_to_process:
                    rows = df[df['Q_name'] == qname].copy() 
                    top_blat_results = compute_blat_score(rows, N)
                    
                    batch.append((top_blat_results, qname, sample_name))
                    
                    if len(batch) >= batch_size:
                        curr.executemany(f"""
                            UPDATE {table_name} 
                            SET top_N_blat_results = ? 
                            WHERE qname = ? AND sample_name = ?
                        """, batch)
                        conn.commit()
                        processed_count += len(batch)
                        print(f"  Progress: {processed_count}/{len(qnames_to_process)} reads processed")
                        batch = []
                
                # Process any remaining reads
                if batch:
                    curr.executemany(f"""
                        UPDATE {table_name} 
                        SET top_N_blat_results = ? 
                        WHERE qname = ? AND sample_name = ?
                    """, batch)
                    conn.commit()
                    processed_count += len(batch)
                
                # Check how many reads were actually updated
                curr.execute(f"SELECT COUNT(*) FROM {table_name} WHERE sample_name = ? AND top_N_blat_results IS NOT NULL", 
                            (sample_name,))
                updated_count = curr.fetchone()[0]
                
                print(f"Sample {sample_name}: {processed_count} reads processed, {updated_count} reads updated with BLAT results")
                
                # Print a few examples for verification
                curr.execute(f"""
                    SELECT qname, top_N_blat_results 
                    FROM {table_name} 
                    WHERE sample_name = ? AND top_N_blat_results IS NOT NULL 
                    LIMIT 3
                """, (sample_name,))
                sample_rows = curr.fetchall()
                if sample_rows:
                    print("Sample updated entries:")
                    for row in sample_rows:
                        print(f"  {row[0]}: {row[1]}")
                
        print(f"{sample_name} processing completed")
        
    except Exception as e:
        print(f"Error processing {psl_file}: {e}")
        import traceback
        traceback.print_exc()

def process_psl_directory(psl_dir, db_path, table_name, N=3):
    """Process each PSL file in directory separately"""
    psl_files = glob.glob(os.path.join(psl_dir, "*.psl"))
    
    if not psl_files:
        print(f"No PSL files found in directory: {psl_dir}")
        return
        
    print(f"Found {len(psl_files)} PSL files to process")
    
    # Try to repair the database if needed before processing any files
    try:
        with closing(connect_to_db(db_path)) as conn:
            pass
    except sqlite3.DatabaseError as e:
        print(f"Attempting to recover corrupted database: {e}")
        try:
            # Create a temporary database to recover data
            temp_db = f"{db_path}.temp"
            os.system(f"echo .dump | sqlite3 {db_path} | sqlite3 {temp_db}")
            os.rename(temp_db, db_path)
            print(f"Database recovery attempted. Original database backed up.")
        except Exception as e:
            print(f"Database recovery failed: {e}")
            print("Please create a new database or restore from backup.")
            return
    
    # Process each PSL file separately
    for i, psl_file in enumerate(psl_files, 1):
        print(f"Processing file {i}/{len(psl_files)}: {psl_file}")
        parse_single_psl_file(psl_file, db_path, table_name, N)
        print(f"Completed file {i}/{len(psl_files)}")
        print("-" * 50)

def main():
    parser = argparse.ArgumentParser(description='Process SAM and PSL files for STR analysis')
    parser.add_argument('--mode', choices=['init', 'blat'], required=True,
                       help='Mode: initialize from SAM files or process BLAT results')
    parser.add_argument('--db-path', required=True,
                       help='Path to SQLite database')
    parser.add_argument('--gene', help='Gene of interest (for init mode)')
    parser.add_argument('--sam-dir', help='Directory containing SAM files (for init mode)')
    parser.add_argument('--psl-dir', help='Directory containing PSL files (for blat mode)')
    parser.add_argument('--top-n', type=int, default=3, help='Number of top BLAT results to save (default: 3)')
    
    args = parser.parse_args()
    
    if args.mode == 'init':
        if not args.gene or not args.sam_dir:
            parser.error("init mode requires --gene and --sam-dir")
            
        table_name = initialize_database(args.db_path)
        
        with closing(connect_to_db(args.db_path)) as conn:
            parse_directory(args.sam_dir, conn, table_name, args.gene)
            
            # Print some sample data for verification
            print("\nSample data after initialization:")
            with closing(conn.cursor()) as curr:
                curr.execute(f"SELECT qname, sample_name, gene FROM {table_name} LIMIT 5")
                sample_rows = curr.fetchall()
                for row in sample_rows:
                    print(row)
                    
    elif args.mode == 'blat':
        if not args.psl_dir:
            parser.error("blat mode requires --psl-dir")
            
        table_name = get_table_name(args.db_path)
        process_psl_directory(args.psl_dir, args.db_path, table_name, args.top_n)

if __name__ == '__main__':
    main()


# import os
# import sys
# import sqlite3
# import argparse
# import pandas as pd
# import numpy as np
# import csv
# import re
# import glob

# def initialize_database(db_path):
#     """Initialize database with table name matching db filename"""
#     table_name = os.path.splitext(os.path.basename(db_path))[0]
#     conn = sqlite3.connect(db_path)
#     curr = conn.cursor()

#     curr.execute(f'''
#     CREATE TABLE IF NOT EXISTS {table_name} (
#         qname TEXT,
#         flag INTEGER,
#         rname TEXT,
#         pos INTEGER,
#         mapq INTEGER,
#         cigar TEXT,
#         rnext TEXT,
#         pnext INTEGER,
#         tlen INTEGER,
#         seq TEXT,
#         qual TEXT,
#         sample_name TEXT,
#         gene TEXT,
#         case_control TEXT,
#         top_N_blat_results TEXT
#     )
#     ''')
    
#     return conn, curr, table_name

# def read_in_sample_exists(curr, table_name, qname, flag, sample_name):
#     """Check if read exists in sample"""
#     curr.execute(f"SELECT 1 FROM {table_name} WHERE qname = ? AND flag = ? AND sample_name = ?", 
#                  (qname, flag, sample_name))
#     return curr.fetchone() is not None

# def parse_sam_file(file_path, curr, table_name, gene_of_interest):
#     """Parses a SAM file and inserts data into the SQLite database"""
#     sample_name = os.path.basename(file_path).split('.')[0]
    
#     with open(file_path, 'r') as sam_file:
#         for line in sam_file:
#             if line.startswith('@'):
#                 continue
#             fields = line.strip().split('\t')
            
#             qname = fields[0]
#             flag = int(fields[1])

#             if read_in_sample_exists(curr, table_name, qname, flag, sample_name):
#                 continue

#             rname = fields[2]
#             pos = int(fields[3])
#             mapq = int(fields[4])
#             cigar = fields[5]
#             rnext = fields[6]
#             pnext = int(fields[7])
#             tlen = int(fields[8])
#             seq = fields[9]
#             qual = fields[10]

#             curr.execute(f'''
#                 INSERT INTO {table_name} (
#                     qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, sample_name, gene
#                 ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
#             ''', (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, sample_name, gene_of_interest))
            
#     sam_file.close()
    
# def parse_directory(directory, curr, table_name, gene_of_interest):
#     """Parses all SAM files within the specified directory and stores the data in the database"""
#     sam_files = glob.glob(os.path.join(directory, "*.sam"))
    
#     for sfile in sam_files:
#         print(f"Importing information from {sfile}...")
#         parse_sam_file(sfile, curr, table_name, gene_of_interest)

# def compute_blat_score(rows, N):
#     """Compute the BLAT score for each row and return the top N rows as a formatted string."""
#     rows['match'] = pd.to_numeric(rows['match'], errors='coerce')
#     rows['mis-match'] = pd.to_numeric(rows['mis-match'], errors='coerce')
#     rows['rep.match'] = pd.to_numeric(rows['rep.match'], errors='coerce')
#     rows['Q gap bases'] = pd.to_numeric(rows['Q gap bases'], errors='coerce')
#     rows['T gap bases'] = pd.to_numeric(rows['T gap bases'], errors='coerce')

#     rows['blat_score'] = rows['match'] - rows['mis-match'] - rows['rep.match'] - rows['Q gap bases'] - rows['T gap bases']
#     top_rows = rows.sort_values(by='blat_score', ascending=False).head(N)
    
#     results = []
#     for i, row in enumerate(top_rows.itertuples(), 1):
#         result = f"{i}:{row.blat_score}:{row.T_name}:{row.T_start}:{row.T_end}:{row.strand}"
#         results.append(result)
    
#     return ';'.join(results)

# def parse_psl_files(psl_files, sqlite_db, table_name, N=3):
#     """Parse one or more PSL files and update the SQLite database with the top N BLAT results."""
#     columns = ['match', 'mis-match', 'rep.match', 'N\'s', 'Q gap count', 'Q gap bases', 
#               'T gap count', 'T gap bases', 'strand', 'Q_name', 'Q size', 'Q start', 
#               'Q end', 'T_name', 'T size', 'T_start', 'T_end', 'block count', 
#               'blockSizes', 'qStarts', 'tStarts']
    
#     dfs = []
#     for psl_file in psl_files:
#         df = pd.read_csv(psl_file, sep='\t', header=None, names=columns, 
#                         comment='#', skiprows=5)
#         dfs.append(df)
    
#     # Combine all dataframes
#     combined_df = pd.concat(dfs, ignore_index=True)
#     unique_qnames = combined_df['Q_name'].unique()

#     conn = sqlite3.connect(sqlite_db)
#     curr = conn.cursor()

#     try:
#         curr.execute(f"ALTER TABLE {table_name} ADD COLUMN top_N_blat_results TEXT")
#     except sqlite3.OperationalError as e:
#         if "duplicate column" not in str(e).lower():
#             raise e
#         print("Column 'top_N_blat_results' already exists.")

#     curr.execute(f"PRAGMA table_info({table_name});")
#     columns = curr.fetchall()
#     print(f"Columns in {table_name}: {columns}")
    
#     curr.execute(f"SELECT COUNT(*) FROM {table_name}")
#     count = curr.fetchone()
#     print(f"Number of rows in {table_name}: {count[0]}")

#     for qname in unique_qnames:
#         rows = combined_df[combined_df['Q_name'] == qname].copy() 
#         top_blat_results = compute_blat_score(rows, N)
        
#         # print(table_name)
#         curr.execute(f"SELECT * FROM {table_name} WHERE qname=?", (qname,))
#         db_entry = curr.fetchone()
#         # print(db_entry)

#         if db_entry:
#             curr.execute(f"""
#                 UPDATE {table_name} 
#                 SET top_N_blat_results = ? 
#                 WHERE qname = ?
#             """, (top_blat_results, qname))
#         else:
#             print(f"No database entry found for {qname}")

#     print_current_db(conn, table_name)
#     conn.commit()
#     conn.close()

# def print_current_db(conn, table_name):
#     df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
#     print(df)

# def main():
#     parser = argparse.ArgumentParser(description='Process SAM and PSL files for STR analysis')
#     parser.add_argument('--mode', choices=['init', 'blat'], required=True,
#                        help='Mode: initialize from SAM files or process BLAT results')
#     parser.add_argument('--db-path', required=True,
#                        help='Path to SQLite database')
#     parser.add_argument('--gene', help='Gene of interest (for init mode)')
#     parser.add_argument('--sam-dir', help='Directory containing SAM files (for init mode)')
#     parser.add_argument('--psl-files', nargs='+', help='PSL files to process (for blat mode)')
    
#     args = parser.parse_args()
    
#     if args.mode == 'init':
#         if not args.gene or not args.sam_dir:
#             parser.error("init mode requires --gene and --sam-dir")
#         conn, curr, table_name = initialize_database(args.db_path)
#         parse_directory(args.sam_dir, curr, table_name, args.gene)
#         conn.commit()
#         print_current_db(conn, table_name)
#         conn.close()
        
#     elif args.mode == 'blat':
#         if not args.psl_files:
#             parser.error("blat mode requires --psl-files")
#         table_name = os.path.splitext(os.path.basename(args.db_path))[0]
#         parse_psl_files(args.psl_files, args.db_path, table_name)

# if __name__ == '__main__':
#     main()