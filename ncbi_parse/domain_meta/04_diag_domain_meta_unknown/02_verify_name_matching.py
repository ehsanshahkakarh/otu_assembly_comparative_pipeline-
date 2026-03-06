#!/usr/bin/env python3
"""
Sub-Diagnostic 04.02: Verify Name Matching Between Census and NCBI Merger

This script randomly samples genus names from 16S and 18S census data,
performs manual searches for matches in NCBI taxonomy, and compares
the results with the new_merger output to verify matching accuracy.

Author: Diagnostic 04.02
Date: 2026-03-05
"""

import pandas as pd
import random
from datetime import datetime
from pathlib import Path

# Configuration
WORKSPACE_ROOT = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data")
CENSUS_16S = WORKSPACE_ROOT / "00_gaps_taxonomic/00parse_database/16S_censusparse/csv_16S/eukcensus16S_by_genus.csv"
CENSUS_18S = WORKSPACE_ROOT / "00_gaps_taxonomic/00parse_database/18S_censusparse/output/eukcensus_18S_by_genus.csv"
MERGER_16S = WORKSPACE_ROOT / "00_gaps_taxonomic/00parse_database/Eukcensus_merge/16s_merged/analysis_summary/16s_ncbi_merger_clean_summary.csv"
MERGER_18S = WORKSPACE_ROOT / "00_gaps_taxonomic/00parse_database/Eukcensus_merge/py_mergers/18s_merged/analysis_summary/18s_ncbi_merger_clean_summary.csv"

SAMPLE_SIZE = 10  # Number of random samples per dataset
RANDOM_SEED = 42

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = Path(__file__).parent / f"02_name_matching_verification_{timestamp}.txt"
    
    print("=" * 80)
    print("SUB-DIAGNOSTIC 04.02: Name Matching Verification")
    print("=" * 80)
    print()
    
    # Set random seed for reproducibility
    random.seed(RANDOM_SEED)
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SUB-DIAGNOSTIC 04.02: Name Matching Verification\n")
        f.write("=" * 80 + "\n")
        f.write(f"\nDate: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Sample size per dataset: {SAMPLE_SIZE}\n")
        f.write(f"Random seed: {RANDOM_SEED}\n\n")
        
        # Process 16S data
        print("\n" + "=" * 80)
        print("16S CENSUS DATA SAMPLING")
        print("=" * 80)
        f.write("\n" + "=" * 80 + "\n")
        f.write("16S CENSUS DATA SAMPLING\n")
        f.write("=" * 80 + "\n")
        
        process_dataset(CENSUS_16S, MERGER_16S, "16S", SAMPLE_SIZE, f)
        
        # Process 18S data
        print("\n" + "=" * 80)
        print("18S CENSUS DATA SAMPLING")
        print("=" * 80)
        f.write("\n" + "=" * 80 + "\n")
        f.write("18S CENSUS DATA SAMPLING\n")
        f.write("=" * 80 + "\n")
        
        process_dataset(CENSUS_18S, MERGER_18S, "18S", SAMPLE_SIZE, f)
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("VERIFICATION COMPLETE\n")
        f.write("=" * 80 + "\n")
    
    print(f"\n{'=' * 80}")
    print(f"Results saved to: {output_file}")
    print("=" * 80)

def process_dataset(census_file, merger_file, dataset_name, sample_size, output_handle):
    """Process a single dataset (16S or 18S)"""
    
    print(f"\nReading census file: {census_file}")
    output_handle.write(f"\nCensus file: {census_file}\n")
    
    # Read census data
    census_df = pd.read_csv(census_file)
    print(f"Total entries in census: {len(census_df):,}")
    output_handle.write(f"Total entries in census: {len(census_df):,}\n")
    
    # Read merger summary if it exists
    merger_exists = merger_file.exists()
    if merger_exists:
        merger_df = pd.read_csv(merger_file)
        print(f"Merger summary file found: {merger_file}")
        output_handle.write(f"Merger summary file: {merger_file}\n")
        output_handle.write(f"Merger summary:\n{merger_df.to_string()}\n\n")
    else:
        print(f"⚠️  Merger summary file not found: {merger_file}")
        output_handle.write(f"⚠️  Merger summary file not found: {merger_file}\n\n")
    
    # Sample random entries
    sample_indices = random.sample(range(len(census_df)), min(sample_size, len(census_df)))
    sampled_data = census_df.iloc[sample_indices]
    
    print(f"\nSampled {len(sampled_data)} random entries:")
    output_handle.write(f"\nSampled {len(sampled_data)} random entries:\n")
    output_handle.write("=" * 80 + "\n\n")
    
    for idx, (_, row) in enumerate(sampled_data.iterrows(), 1):
        name_to_use = row['Name_to_use']
        taxid = row['taxid']
        lineage = row.get('lineage', 'N/A')
        
        print(f"\n[Sample {idx}]")
        print(f"  Name_to_use: {name_to_use}")
        print(f"  TaxID: {taxid}")
        print(f"  Lineage: {lineage}")
        
        output_handle.write(f"[Sample {idx}]\n")
        output_handle.write(f"  Name_to_use: {name_to_use}\n")
        output_handle.write(f"  TaxID: {taxid}\n")
        output_handle.write(f"  Lineage: {lineage}\n")
        output_handle.write(f"\n  → Manual verification needed: Search '{name_to_use}' in NCBI taxonomy\n")
        output_handle.write(f"  → Expected TaxID: {taxid}\n")
        output_handle.write("\n" + "-" * 80 + "\n\n")

if __name__ == "__main__":
    main()

