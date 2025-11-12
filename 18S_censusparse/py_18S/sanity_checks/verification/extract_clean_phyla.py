#!/usr/bin/env python3
"""
Extract Clean Phyla Names for Taxonkit

This script extracts phyla (division) names from the division_matches.csv file
and outputs them to a text file in a clean format suitable for taxonkit:
- One name per line
- No numbers or list formatting
- No headers or comments
- Alphabetically sorted

Usage:
    python extract_clean_phyla.py

Output:
    phyla_names_for_taxonkit.txt: A clean text file containing just phyla names
"""

import pandas as pd
import os

def extract_clean_phyla():
    """
    Extract phyla names from division_matches.csv and save to a clean text file for taxonkit.
    """
    # Input file
    input_file = "division_matches.csv"
    
    # Output file
    output_file = "phyla_names_for_taxonkit.txt"
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        return
    
    try:
        # Read the CSV file
        df = pd.read_csv(input_file)
        
        # Check if taxon_name column exists
        if 'taxon_name' not in df.columns:
            print(f"Error: 'taxon_name' column not found in {input_file}")
            print(f"Available columns: {df.columns.tolist()}")
            return
        
        # Extract phyla names
        phyla_names = df['taxon_name'].tolist()
        
        # Sort alphabetically
        phyla_names.sort()
        
        # Write to output file - just the names, one per line
        with open(output_file, 'w') as f:
            for name in phyla_names:
                f.write(f"{name}\n")
        
        print(f"Successfully wrote {len(phyla_names)} clean phyla names to {output_file}")
        
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    extract_clean_phyla()
