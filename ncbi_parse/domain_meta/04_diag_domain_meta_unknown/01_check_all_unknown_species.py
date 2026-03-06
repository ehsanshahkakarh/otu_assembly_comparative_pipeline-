#!/usr/bin/env python3
"""
Sanity Check Script: Check if taxid == species_taxid for ALL unknown species.

This script reads the unknown_species.csv file and checks if all assemblies
for each unknown species have taxid == species_taxid in the raw assembly file.

Usage:
    python check_all_unknown_species.py
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

def find_assembly_file():
    """Find the NCBI assembly summary file."""
    script_dir = Path(__file__).parent
    
    possible_paths = [
        script_dir.parent.parent.parent / "00assembly_summary_genbank.txt",
        script_dir.parent.parent.parent / "old_pipeline" / "ncbi_parse_old" / "metadata" / "00assembly_summary_genbank.txt",
    ]
    
    for path in possible_paths:
        if path.exists():
            return path
    
    print("❌ Could not find 00assembly_summary_genbank.txt")
    return None

def check_all_unknown_species():
    """Check taxid vs species_taxid for all unknown species."""
    
    script_dir = Path(__file__).parent
    unknown_species_file = script_dir.parent / "output" / "unknown_species.csv"
    
    # Find the assembly file
    assembly_file = find_assembly_file()
    if assembly_file is None:
        return
    
    print("="*80)
    print("SANITY CHECK: Checking taxid vs species_taxid for ALL unknown species")
    print("="*80)
    print()
    
    # Read unknown species file
    print(f"Reading unknown species file: {unknown_species_file}")
    unknown_df = pd.read_csv(unknown_species_file)
    print(f"Total unknown species: {len(unknown_df):,}")
    print()
    
    # Get list of species_taxids to check
    species_taxids = unknown_df['species_taxid'].tolist()
    print(f"Checking {len(species_taxids):,} species_taxids...")
    print()
    
    # Read assembly file
    print(f"Reading assembly file: {assembly_file}")
    print(f"File size: {assembly_file.stat().st_size / (1024**3):.2f} GB")
    print("This may take a minute...")
    
    # Find header line
    with open(assembly_file, 'r') as f:
        lines = f.readlines()
    
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#assembly_accession'):
            header_line = line.strip().lstrip('#')
            data_start = i + 1
            break
    
    # Read assembly data
    assembly_df = pd.read_csv(assembly_file, sep='\t', skiprows=data_start, 
                              names=header_line.split('\t'), low_memory=False)
    
    print(f"Total assemblies in file: {len(assembly_df):,}")
    print()
    
    # Filter for unknown species
    print("Filtering for unknown species...")
    filtered_df = assembly_df[assembly_df['species_taxid'].isin(species_taxids)]
    print(f"Total assemblies for unknown species: {len(filtered_df):,}")
    print()
    
    # Check taxid vs species_taxid
    print("="*80)
    print("RESULTS: TaxID vs Species_TaxID Check")
    print("="*80)
    print()
    
    total_assemblies = len(filtered_df)
    same_count = (filtered_df['taxid'] == filtered_df['species_taxid']).sum()
    different_count = (filtered_df['taxid'] != filtered_df['species_taxid']).sum()
    
    print(f"Total assemblies checked: {total_assemblies:,}")
    print(f"  Same (taxid == species_taxid): {same_count:,} ({same_count/total_assemblies*100:.2f}%)")
    print(f"  Different (taxid != species_taxid): {different_count:,} ({different_count/total_assemblies*100:.2f}%)")
    print()
    
    if different_count > 0:
        print(f"⚠️  WARNING: Found {different_count:,} assemblies where taxid != species_taxid")
        print()
        
        # Group by species_taxid to see which species have mismatches
        print("Species with taxid != species_taxid:")
        print("-" * 80)
        
        diff_df = filtered_df[filtered_df['taxid'] != filtered_df['species_taxid']]
        species_with_diff = diff_df.groupby('species_taxid').size().sort_values(ascending=False)
        
        for species_taxid, count in species_with_diff.items():
            species_name = unknown_df[unknown_df['species_taxid'] == species_taxid]['species_name'].values[0]
            total_for_species = (filtered_df['species_taxid'] == species_taxid).sum()
            print(f"  species_taxid {species_taxid} ({species_name}):")
            print(f"    {count:,} / {total_for_species:,} assemblies have different taxid")
            
            # Show unique taxid values for this species
            unique_taxids = diff_df[diff_df['species_taxid'] == species_taxid]['taxid'].unique()
            print(f"    Unique taxid values: {sorted(unique_taxids)[:10]}")  # Show first 10
            print()
    else:
        print("✓ ALL assemblies have taxid == species_taxid")
        print()
        print("This means that for all unknown species, every assembly record")
        print("has the same taxid and species_taxid value.")
    
    print("="*80)
    
    # Save summary to .txt file only
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = script_dir / f"unknown_species_taxid_check_{timestamp}.txt"

    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SANITY CHECK: TaxID vs Species_TaxID for Unknown Species\n")
        f.write("="*80 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Total unknown species checked: {len(species_taxids):,}\n")
        f.write(f"Total assemblies checked: {total_assemblies:,}\n\n")
        f.write(f"Same (taxid == species_taxid): {same_count:,} ({same_count/total_assemblies*100:.2f}%)\n")
        f.write(f"Different (taxid != species_taxid): {different_count:,} ({different_count/total_assemblies*100:.2f}%)\n\n")

        if different_count == 0:
            f.write("RESULT: ✓ CONFIRMED\n")
            f.write("All assemblies for unknown species have taxid == species_taxid\n\n")
        else:
            f.write("RESULT: ✗ DENIED\n")
            f.write("Some assemblies have taxid != species_taxid\n\n")
            f.write("Species with mismatches:\n")
            f.write("-" * 80 + "\n")
            for species_taxid, count in species_with_diff.items():
                species_name = unknown_df[unknown_df['species_taxid'] == species_taxid]['species_name'].values[0]
                total_for_species = (filtered_df['species_taxid'] == species_taxid).sum()
                f.write(f"species_taxid {species_taxid} ({species_name}): {count:,} / {total_for_species:,}\n")

    print(f"\nResult saved to: {output_file}")

if __name__ == "__main__":
    check_all_unknown_species()
    print("\n" + "="*80)
    print("Sanity check complete!")
    print("="*80)

