#!/usr/bin/env python3
"""
Extract Taxids from Division Matches

This script reads the division_matches.csv file, extracts the taxid column,
and outputs just the valid taxids to taxid.txt.

Usage:
    python extract_phyla_taxids.py

Output:
    taxid.txt: A text file with one taxid per line
"""

import pandas as pd
import os

def extract_phyla_taxids():
    """
    Extract taxids from division_matches.csv and save to a text file.
    """
    # Input file
    input_file = "division_matches.csv"
    
    # Output file
    output_file = "taxid.txt"
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        return
    
    try:
        # Read the CSV file
        df = pd.read_csv(input_file)
        
        # Check if taxid column exists
        if 'taxid' not in df.columns:
            print(f"Error: 'taxid' column not found in {input_file}")
            print(f"Available columns: {df.columns.tolist()}")
            return
        
        # Extract taxids, filtering out "NA" values
        taxids = df['taxid'].tolist()
        valid_taxids = [str(taxid) for taxid in taxids if str(taxid) != "NA" and str(taxid).strip()]
        
        # Write to output file - just the taxids, one per line
        with open(output_file, 'w') as f:
            for taxid in valid_taxids:
                f.write(f"{taxid}\n")
        
        print(f"Successfully wrote {len(valid_taxids)} taxids to {output_file}")
        print(f"Filtered out {len(taxids) - len(valid_taxids)} invalid taxids")
        
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    extract_phyla_taxids()
