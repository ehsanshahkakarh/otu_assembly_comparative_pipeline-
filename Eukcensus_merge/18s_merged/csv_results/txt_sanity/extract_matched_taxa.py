#!/usr/bin/env python3
"""
Extract all matched taxa from 18S NCBI merged files.
Outputs separate .txt files for phyla, families, and genera with their counts.

Created: 2025-01-12
"""

import pandas as pd
import os
from datetime import datetime

def load_18s_ncbi_data():
    """Load the 18S NCBI merged data files."""
    base_dir = "."
    
    # Load the three taxonomic levels
    phylum_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_phylum.csv"))
    family_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_family.csv"))
    genus_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_genus.csv"))
    
    # Filter for matched entries only
    phylum_matched = phylum_df[phylum_df['match_status'] == 'matched'].copy()
    family_matched = family_df[family_df['match_status'] == 'matched'].copy()
    genus_matched = genus_df[genus_df['match_status'] == 'matched'].copy()
    
    return phylum_matched, family_matched, genus_matched

def extract_taxa_to_files(phylum_df, family_df, genus_df):
    """Extract matched taxa and save to separate .txt files."""
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Extract phyla
    phylum_data = []
    for _, row in phylum_df.iterrows():
        phylum_data.append({
            'taxon': row['phylum'],
            'census_otu_count': row['census_otu_count'],
            'ncbi_species_count': row['ncbi_species_count']
        })
    
    # Extract families
    family_data = []
    for _, row in family_df.iterrows():
        family_data.append({
            'taxon': row['family'],
            'census_otu_count': row['census_otu_count'],
            'ncbi_species_count': row['ncbi_species_count']
        })
    
    # Extract genera
    genus_data = []
    for _, row in genus_df.iterrows():
        genus_data.append({
            'taxon': row['genus'],
            'census_otu_count': row['census_otu_count'],
            'ncbi_species_count': row['ncbi_species_count']
        })
    
    # Sort by NCBI species count (descending), then by census OTU count (descending)
    phylum_data.sort(key=lambda x: (-x['ncbi_species_count'], -x['census_otu_count']))
    family_data.sort(key=lambda x: (-x['ncbi_species_count'], -x['census_otu_count']))
    genus_data.sort(key=lambda x: (-x['ncbi_species_count'], -x['census_otu_count']))
    
    # Write phyla file
    with open('18s_ncbi_matched_phyla.txt', 'w') as f:
        f.write(f"18S NCBI Matched Phyla\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"Total matched phyla: {len(phylum_data)}\n")
        f.write("=" * 60 + "\n\n")
        f.write("Phylum\tCensus_OTU_Count\tNCBI_Species_Count\n")
        f.write("-" * 60 + "\n")
        
        for item in phylum_data:
            f.write(f"{item['taxon']}\t{item['census_otu_count']}\t{item['ncbi_species_count']}\n")
    
    # Write families file
    with open('18s_ncbi_matched_families.txt', 'w') as f:
        f.write(f"18S NCBI Matched Families\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"Total matched families: {len(family_data)}\n")
        f.write("=" * 60 + "\n\n")
        f.write("Family\tCensus_OTU_Count\tNCBI_Species_Count\n")
        f.write("-" * 60 + "\n")
        
        for item in family_data:
            f.write(f"{item['taxon']}\t{item['census_otu_count']}\t{item['ncbi_species_count']}\n")
    
    # Write genera file
    with open('18s_ncbi_matched_genera.txt', 'w') as f:
        f.write(f"18S NCBI Matched Genera\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"Total matched genera: {len(genus_data)}\n")
        f.write("=" * 60 + "\n\n")
        f.write("Genus\tCensus_OTU_Count\tNCBI_Species_Count\n")
        f.write("-" * 60 + "\n")
        
        for item in genus_data:
            f.write(f"{item['taxon']}\t{item['census_otu_count']}\t{item['ncbi_species_count']}\n")
    
    return len(phylum_data), len(family_data), len(genus_data)

def main():
    """Main function to extract matched taxa."""
    print("Loading 18S NCBI merged data...")
    phylum_df, family_df, genus_df = load_18s_ncbi_data()
    
    print(f"Loaded data:")
    print(f"  - Matched Phyla: {len(phylum_df)}")
    print(f"  - Matched Families: {len(family_df)}")
    print(f"  - Matched Genera: {len(genus_df)}")
    print()
    
    print("Extracting matched taxa to separate files...")
    phyla_count, families_count, genera_count = extract_taxa_to_files(phylum_df, family_df, genus_df)
    
    print(f"Successfully extracted:")
    print(f"  - {phyla_count} phyla → 18s_ncbi_matched_phyla.txt")
    print(f"  - {families_count} families → 18s_ncbi_matched_families.txt")
    print(f"  - {genera_count} genera → 18s_ncbi_matched_genera.txt")
    print()
    print("All files saved with census_otu_count and ncbi_species_count columns.")

if __name__ == "__main__":
    main()
