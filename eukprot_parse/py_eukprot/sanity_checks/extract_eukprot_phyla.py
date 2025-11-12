#!/usr/bin/env python3
"""
Extract EukProt Phylum Names

This script parses the phylum_grouped.csv file and extracts only the phylum names.

Usage:
    python extract_eukprot_phyla.py
"""

import pandas as pd

def main():
    # Input file
    input_file = "phylum_grouped.csv"
    
    try:
        # Read the CSV file using pandas
        df = pd.read_csv(input_file)
        
        # Check which column contains the phylum names
        print(f"Columns in the file: {', '.join(df.columns)}")
        
        # Assuming the phylum column is named 'phylum'
        # If it's named differently, change this line
        if 'phylum' in df.columns:
            phylum_col = 'phylum'
        elif 'Phylum' in df.columns:
            phylum_col = 'Phylum'
        else:
            # Try to find a column that might contain phylum names
            potential_cols = [col for col in df.columns if 'phyl' in col.lower()]
            if potential_cols:
                phylum_col = potential_cols[0]
                print(f"Using column '{phylum_col}' for phylum names")
            else:
                # If we can't find a suitable column, print the first few rows
                # so the user can identify the correct column
                print("Could not identify phylum column. Here are the first few rows:")
                print(df.head())
                return
        
        # Extract phylum names
        phylum_names = df[phylum_col].unique().tolist()
        
        # Remove any NaN values
        phylum_names = [p for p in phylum_names if isinstance(p, str)]
        
        # Sort alphabetically
        phylum_names.sort()
        
        # Print phylum names
        print("\nEukProt Phylum Names:")
        for i, name in enumerate(phylum_names):
            print(f"{i+1}. {name}")
        
        # Write phylum names to a text file
        with open("eukprot_phylum_names.txt", "w") as f:
            for name in phylum_names:
                f.write(f"{name}\n")
        
        print(f"\nTotal phyla: {len(phylum_names)}")
        print("Phylum names written to eukprot_phylum_names.txt")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
