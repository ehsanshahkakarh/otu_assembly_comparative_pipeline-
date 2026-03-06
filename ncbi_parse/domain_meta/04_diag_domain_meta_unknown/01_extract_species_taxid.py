#!/usr/bin/env python3
"""
Sanity Check Script: Extract rows for any species_taxid from the original 
00assembly_summary_genbank.txt file.

This script extracts all rows with a specified species_taxid from the NCBI 
assembly summary file and checks if taxid == species_taxid for all entries.

Usage:
    python extract_species_taxid.py <species_taxid>
    python extract_species_taxid.py 749906
"""

import pandas as pd
from pathlib import Path
from datetime import datetime
import sys

def find_assembly_file():
    """Find the NCBI assembly summary file."""
    script_dir = Path(__file__).parent
    
    possible_paths = [
        script_dir.parent.parent.parent / "00assembly_summary_genbank.txt",
        script_dir.parent.parent.parent / "old_pipeline" / "ncbi_parse_old" / "metadata" / "00assembly_summary_genbank.txt",
        Path("00assembly_summary_genbank.txt"),
        Path("../00assembly_summary_genbank.txt"),
        Path("../../00assembly_summary_genbank.txt"),
        Path("../../../00assembly_summary_genbank.txt"),
    ]
    
    for path in possible_paths:
        if path.exists():
            return path
    
    print("❌ Could not find 00assembly_summary_genbank.txt")
    print("   Searched in:", [str(p) for p in possible_paths])
    return None

def extract_species_taxid(target_taxid):
    """Extract all rows with the specified species_taxid from 00assembly_summary_genbank.txt"""
    
    # Define paths
    script_dir = Path(__file__).parent
    output_dir = script_dir
    
    # Find the assembly file
    input_file = find_assembly_file()
    if input_file is None:
        return
    
    # Create output filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = output_dir / f"taxid_{target_taxid}_raw_assembly_{timestamp}.tsv"
    
    print(f"Reading input file: {input_file}")
    print(f"File size: {input_file.stat().st_size / (1024**3):.2f} GB")
    
    # Read the assembly file with proper header handling
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line (starts with #assembly_accession)
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#assembly_accession'):
            header_line = line.strip().lstrip('#')
            data_start = i + 1
            break
    
    if header_line is None:
        print("ERROR: Could not find header line in assembly file")
        return
    
    print(f"Found header at line {data_start}")
    print(f"Reading data starting from line {data_start + 1}...")
    
    # Read the data using pandas, skipping the comment lines
    df = pd.read_csv(input_file, sep='\t', skiprows=data_start, names=header_line.split('\t'), low_memory=False)
    
    print(f"Total rows in assembly file: {len(df):,}")
    print(f"Columns: {len(df.columns)}")
    
    # Filter for the target species_taxid
    filtered_df = df[df['species_taxid'] == target_taxid]
    
    print(f"\nRows with species_taxid {target_taxid}: {len(filtered_df):,}")
    
    if len(filtered_df) == 0:
        print(f"WARNING: No rows found with species_taxid {target_taxid}")
        return
    
    # Save to output file
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nData saved to: {output_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS:")
    print("="*80)
    print(f"Total assemblies for species_taxid {target_taxid}: {len(filtered_df):,}")
    
    # Get organism name
    if 'organism_name' in filtered_df.columns:
        organism_names = filtered_df['organism_name'].unique()
        print(f"Organism name(s): {', '.join(map(str, organism_names))}")
    
    print(f"\nFirst 3 rows:")
    print(filtered_df.head(3).to_string(index=False))
    print("\n" + "="*80)
    
    # Check if taxid and species_taxid are the same
    print("\n" + "="*80)
    print("TAXID vs SPECIES_TAXID CHECK:")
    print("="*80)
    if 'taxid' in filtered_df.columns and 'species_taxid' in filtered_df.columns:
        same_count = (filtered_df['taxid'] == filtered_df['species_taxid']).sum()
        different_count = (filtered_df['taxid'] != filtered_df['species_taxid']).sum()
        print(f"  Same (taxid == species_taxid): {same_count:,}")
        print(f"  Different (taxid != species_taxid): {different_count:,}")
        
        if different_count > 0:
            print(f"\n  ⚠️  WARNING: Found {different_count} entries where taxid != species_taxid")
            print("\n  Entries with different taxid and species_taxid:")
            diff_entries = filtered_df[filtered_df['taxid'] != filtered_df['species_taxid']][['assembly_accession', 'taxid', 'species_taxid', 'organism_name']]
            print(diff_entries.to_string(index=False))
            
            # Show unique taxid values
            unique_taxids = filtered_df[filtered_df['taxid'] != filtered_df['species_taxid']]['taxid'].unique()
            print(f"\n  Unique taxid values (when different from species_taxid {target_taxid}):")
            for tid in sorted(unique_taxids):
                count = (filtered_df['taxid'] == tid).sum()
                print(f"    taxid {tid}: {count:,} entries")
        else:
            print(f"  ✓ All {same_count:,} entries have taxid == species_taxid ({target_taxid})")
    print("="*80)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_species_taxid.py <species_taxid>")
        print("Example: python extract_species_taxid.py 749906")
        sys.exit(1)
    
    try:
        target_taxid = int(sys.argv[1])
    except ValueError:
        print(f"ERROR: Invalid species_taxid '{sys.argv[1]}'. Must be an integer.")
        sys.exit(1)
    
    print("="*80)
    print(f"SANITY CHECK: Extracting species_taxid {target_taxid}")
    print("FROM: 00assembly_summary_genbank.txt (raw NCBI assembly data)")
    print("="*80)
    print()
    
    extract_species_taxid(target_taxid)
    
    print("\n" + "="*80)
    print("Sanity check complete!")
    print("="*80)

