#!/usr/bin/env python3
"""
Extract Taxids from Phylum Taxid File

This script reads the phylum_taxid.txt file, which contains phylum names and taxids,
and extracts just the taxid column to a new file.

Usage:
    python extract_phylum_taxids.py

Input:
    phylum_taxid.txt: A tab-separated file with phylum names and taxids

Output:
    phylum_taxids_only.txt: A text file with one taxid per line
"""

def extract_taxids():
    """
    Extract taxids from phylum_taxid.txt and save to a new file.
    """
    # Input and output files
    input_file = "phylum_taxid.txt"
    output_file = "phylum_taxids_only.txt"
    
    try:
        # Read the input file
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Extract taxids
        taxids = []
        for line in lines:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            # Split by tab or multiple spaces
            parts = line.split('\t')
            if len(parts) < 2:  # If not tab-separated, try splitting by spaces
                parts = line.split()
                
            if len(parts) >= 2:
                taxid = parts[-1].strip()  # Get the last part as taxid
                if taxid.isdigit():  # Ensure it's a valid taxid
                    taxids.append(taxid)
        
        # Write taxids to output file
        with open(output_file, 'w') as f:
            for taxid in taxids:
                f.write(f"{taxid}\n")
        
        print(f"Successfully extracted {len(taxids)} taxids to {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    extract_taxids()
