#!/usr/bin/env python3
"""
Sanity Check Script: Extract rows for species_taxid 749906 (gut metagenome)
from the original 00assembly_summary_genbank.txt file.

This script extracts all rows with species_taxid 749906 from the NCBI assembly
summary file to verify all the raw entries for this specific taxid.
"""

import pandas as pd
import os
from pathlib import Path
from datetime import datetime

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

def extract_taxid_749906():
    """Extract all rows with species_taxid 749906 from 00assembly_summary_genbank.txt"""

    # Define paths
    script_dir = Path(__file__).parent
    output_dir = script_dir

    # Find the assembly file
    input_file = find_assembly_file()
    if input_file is None:
        return

    # Create output filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = output_dir / f"taxid_749906_raw_assembly_{timestamp}.tsv"

    print(f"Reading input file: {input_file}")
    print(f"File size: {input_file.stat().st_size / (1024**3):.2f} GB")

    # Read the assembly file with proper header handling
    # The file has comment lines starting with ## and a header line starting with #
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
    df = pd.read_csv(input_file, sep='\t', skiprows=data_start, names=header_line.split('\t'))

    print(f"Total rows in assembly file: {len(df):,}")
    print(f"Columns: {len(df.columns)}")

    # Filter for species_taxid 749906
    # species_taxid is the 7th column (index 6)
    filtered_df = df[df['species_taxid'] == 749906]

    print(f"\nRows with species_taxid 749906: {len(filtered_df):,}")

    if len(filtered_df) == 0:
        print("WARNING: No rows found with species_taxid 749906")
        return

    # Save to output file
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nData saved to: {output_file}")

    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS:")
    print("="*80)
    print(f"Total assemblies for species_taxid 749906: {len(filtered_df):,}")
    print(f"\nFirst 5 rows:")
    print(filtered_df.head().to_string(index=False))
    print("\n" + "="*80)

    # Additional statistics
    print("\n" + "="*80)
    print("DETAILED STATISTICS:")
    print("="*80)

    # Check if taxid and species_taxid are the same
    print("\nTaxID vs Species_TaxID Check:")
    if 'taxid' in filtered_df.columns and 'species_taxid' in filtered_df.columns:
        same_count = (filtered_df['taxid'] == filtered_df['species_taxid']).sum()
        different_count = (filtered_df['taxid'] != filtered_df['species_taxid']).sum()
        print(f"  Same (taxid == species_taxid): {same_count:,}")
        print(f"  Different (taxid != species_taxid): {different_count:,}")

        if different_count > 0:
            print(f"\n  ⚠️  WARNING: Found {different_count} entries where taxid != species_taxid")
            print("\n  Sample of entries with different taxid and species_taxid:")
            diff_entries = filtered_df[filtered_df['taxid'] != filtered_df['species_taxid']][['assembly_accession', 'taxid', 'species_taxid', 'organism_name']]
            print(diff_entries.head(10).to_string(index=False))
        else:
            print(f"  ✓ All {same_count:,} entries have taxid == species_taxid (749906)")

    # Count by assembly level
    if 'assembly_level' in filtered_df.columns:
        print("\nAssembly Levels:")
        assembly_counts = filtered_df['assembly_level'].value_counts()
        for level, count in assembly_counts.items():
            print(f"  {level}: {count:,}")

    # Count by refseq_category
    if 'refseq_category' in filtered_df.columns:
        print("\nRefSeq Categories:")
        refseq_counts = filtered_df['refseq_category'].value_counts()
        for cat, count in refseq_counts.items():
            print(f"  {cat}: {count:,}")

    # Count by genome_rep
    if 'genome_rep' in filtered_df.columns:
        print("\nGenome Representation:")
        genome_rep_counts = filtered_df['genome_rep'].value_counts()
        for rep, count in genome_rep_counts.items():
            print(f"  {rep}: {count:,}")

    # Check excluded_from_refseq
    if 'excluded_from_refseq' in filtered_df.columns:
        excluded_count = filtered_df['excluded_from_refseq'].notna().sum()
        print(f"\nExcluded from RefSeq: {excluded_count:,} assemblies")
        if excluded_count > 0:
            print("  Exclusion reasons:")
            exclusion_reasons = filtered_df[filtered_df['excluded_from_refseq'].notna()]['excluded_from_refseq'].value_counts()
            for reason, count in exclusion_reasons.head(10).items():
                print(f"    {reason}: {count:,}")

    print("="*80)

if __name__ == "__main__":
    print("="*80)
    print("SANITY CHECK: Extracting species_taxid 749906 (gut metagenome)")
    print("FROM: 00assembly_summary_genbank.txt (raw NCBI assembly data)")
    print("="*80)
    print()

    extract_taxid_749906()

    print("\n" + "="*80)
    print("Sanity check complete!")
    print("="*80)

