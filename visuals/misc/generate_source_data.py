#!/usr/bin/env python3
"""
Generate Source Data for Mega Visual Figures
Created: 2025-01-26
Purpose: Extract genome/isolate ratio data for all mega visual plots
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

def calculate_genome_isolate_ratio(data):
    """Calculate genome/isolate ratio with proper handling of zero isolates."""
    # Handle cases where isolate_count is 0
    max_ratio = data['ncbi_genome_count'].max() if len(data) > 0 else 1
    
    data['Genome_Isolate_Ratio'] = np.where(
        data['isolate_count'] > 0,
        data['ncbi_genome_count'] / data['isolate_count'],
        max_ratio
    )
    return data

def process_16s_data(level, domain_filter, data_dir, output_dir):
    """Process 16S data for a specific taxonomic level and domain."""
    file_path = data_dir / f"16s_ncbi_merged_clean_{level}.csv"
    
    if not file_path.exists():
        print(f"Warning: File not found: {file_path}")
        return None
    
    print(f"Processing 16S {level} data for {domain_filter}")
    
    # Load data
    data = pd.read_csv(file_path)
    
    # Filter by domain
    data = data[data['domain'] == domain_filter].copy()
    
    if len(data) == 0:
        print(f"No data found for domain {domain_filter}")
        return None
    
    # Calculate genome/isolate ratio
    data = calculate_genome_isolate_ratio(data)
    
    # Select relevant columns for source data
    source_data = data[[
        level, 'domain', 'census_otu_count', 'census_size_count',
        'ncbi_genome_count', 'ncbi_species_count', 'isolate_count',
        'isolate_percentage', 'Genome_Isolate_Ratio', 'coverage_percentage',
        'match_status'
    ]].copy()
    
    # Rename columns for clarity
    source_data.columns = [
        'Taxon', 'Domain', 'Census_OTU_Count', 'Census_Size_Count',
        'NCBI_Genome_Count', 'NCBI_Species_Count', 'Isolate_Count',
        'Isolate_Percentage', 'Genome_Isolate_Ratio', 'Coverage_Percentage',
        'Match_Status'
    ]
    
    # Sort by genome/isolate ratio (descending)
    source_data = source_data.sort_values('Genome_Isolate_Ratio', ascending=False)
    
    return source_data

def process_18s_data(level, data_dir, output_dir):
    """Process 18S data for a specific taxonomic level."""
    file_path = data_dir / f"18s_ncbi_merged_clean_{level}.csv"
    
    if not file_path.exists():
        print(f"Warning: File not found: {file_path}")
        return None
    
    print(f"Processing 18S {level} data")
    
    # Load data
    data = pd.read_csv(file_path)
    
    if len(data) == 0:
        print(f"No data found in {file_path}")
        return None
    
    # Calculate genome/isolate ratio
    data = calculate_genome_isolate_ratio(data)
    
    # Select relevant columns for source data
    source_data = data[[
        level, 'domain', 'census_otu_count', 'census_size_count',
        'ncbi_genome_count', 'ncbi_species_count', 'isolate_count',
        'isolate_percentage', 'Genome_Isolate_Ratio', 'coverage_percentage',
        'match_status'
    ]].copy()
    
    # Rename columns for clarity
    source_data.columns = [
        'Taxon', 'Domain', 'Census_OTU_Count', 'Census_Size_Count',
        'NCBI_Genome_Count', 'NCBI_Species_Count', 'Isolate_Count',
        'Isolate_Percentage', 'Genome_Isolate_Ratio', 'Coverage_Percentage',
        'Match_Status'
    ]
    
    # Sort by genome/isolate ratio (descending)
    source_data = source_data.sort_values('Genome_Isolate_Ratio', ascending=False)
    
    return source_data

def main():
    """Generate all source data files."""
    # Configuration - adjust paths for running from visuals directory
    base_dir = Path(".")
    data_16s = base_dir / "../Eukcensus_merge/16s_merged/csv_results"
    data_18s = base_dir / "../Eukcensus_merge/18s_merged/csv_results"
    output_dir = base_dir / "final_visuals"
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    print("=== Generating Source Data Files ===\n")
    
    # Define datasets to process
    datasets = [
        # 16S datasets
        ("16S Bacteria Phylum", lambda: process_16s_data("phylum", "Bacteria", data_16s, output_dir), "source_data_16s_bacteria_phylum.csv"),
        ("16S Bacteria Family", lambda: process_16s_data("family", "Bacteria", data_16s, output_dir), "source_data_16s_bacteria_family.csv"),
        ("16S Archaea Phylum", lambda: process_16s_data("phylum", "Archaea", data_16s, output_dir), "source_data_16s_archaea_phylum.csv"),
        ("16S Archaea Family", lambda: process_16s_data("family", "Archaea", data_16s, output_dir), "source_data_16s_archaea_family.csv"),
        # 18S datasets (use phylum for divisions since that's the actual column name)
        ("18S Divisions", lambda: process_18s_data("phylum", data_18s, output_dir), "source_data_18s_divisions.csv"),
        ("18S Family", lambda: process_18s_data("family", data_18s, output_dir), "source_data_18s_family.csv"),
    ]
    
    # Process each dataset
    for name, processor, filename in datasets:
        try:
            data = processor()
            if data is not None:
                output_path = output_dir / filename
                data.to_csv(output_path, index=False)
                print(f"✓ Saved: {output_path} ({len(data)} taxa)")
            else:
                print(f"✗ Failed to process: {name}")
        except Exception as e:
            print(f"✗ Error processing {name}: {e}")
    
    print(f"\n=== Source Data Generation Complete ===")
    print(f"All files saved to: {output_dir}")
    print("Files contain genome/isolate ratios and all relevant metrics for each taxon.")

if __name__ == "__main__":
    main()
