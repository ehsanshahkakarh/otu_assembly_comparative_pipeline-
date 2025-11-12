#!/usr/bin/env python3
"""
Extract Taxonomic Ranks from EukCensus 16S Clusters TSV File

This script parses the eukcensus_16S.clusters.97.tsv file (metadata) and extracts
the unique taxonomic ranks (phylum, family, genus) into separate text files.

No taxonkit or filtering is performed - just a simple extraction of the unique values.

Input file:
- eukcensus_16S.clusters.97.tsv

Output files:
- phylum_names.txt
- family_names.txt
- genus_names.txt
"""

import pandas as pd
from tqdm import tqdm

def main():
    # Input file path
    input_file = "eukcensus_16S.clusters.97.tsv"
    
    # Output file paths
    phylum_output = "phylum_names.txt"
    family_output = "family_names.txt"
    genus_output = "genus_names.txt"
    
    print(f"Reading input file: {input_file}")
    
    # Read the TSV file
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        return
    
    print(f"Successfully read {len(df)} rows from the input file")
    
    # Check if required columns exist
    required_columns = ['phylum', 'familiy', 'genus']  # Note the typo in 'familiy'
    for col in required_columns:
        if col not in df.columns:
            print(f"Error: Required column '{col}' not found in the input file")
            return
    
    # Extract unique values for each taxonomic rank
    print("Extracting unique taxonomic ranks...")
    
    # Get unique phyla
    phyla = df['phylum'].dropna().unique().tolist()
    phyla.sort()
    print(f"Found {len(phyla)} unique phyla")
    
    # Get unique families
    families = df['familiy'].dropna().unique().tolist()  # Note the typo in 'familiy'
    families.sort()
    print(f"Found {len(families)} unique families")
    
    # Get unique genera
    genera = df['genus'].dropna().unique().tolist()
    genera.sort()
    print(f"Found {len(genera)} unique genera")
    
    # Write phyla to file
    print(f"Writing phyla to {phylum_output}")
    with open(phylum_output, 'w') as f:
        for phylum in phyla:
            f.write(f"{phylum}\n")
    
    # Write families to file
    print(f"Writing families to {family_output}")
    with open(family_output, 'w') as f:
        for family in families:
            f.write(f"{family}\n")
    
    # Write genera to file
    print(f"Writing genera to {genus_output}")
    with open(genus_output, 'w') as f:
        for genus in genera:
            f.write(f"{genus}\n")
    
    print("Done! Generated the following files:")
    print(f"- {phylum_output}")
    print(f"- {family_output}")
    print(f"- {genus_output}")

if __name__ == "__main__":
    main()
