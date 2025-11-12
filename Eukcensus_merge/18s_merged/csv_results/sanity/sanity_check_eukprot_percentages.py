#!/usr/bin/env python3
"""
Sanity Check Script for EukProt Census Percentage Calculations
Created: 2025-01-13

This script validates that the census percentages in the EukProt merged files
match exactly with the original census data.
"""

import pandas as pd
from pathlib import Path

def load_original_census_data():
    """Load original census data for comparison."""
    # Use relative path to 18S_censusparse directory
    script_dir = Path(__file__).resolve().parent
    # Go up from sanity -> csv_results -> 18s_merged -> Eukcensus_merge -> OTU_97_eukcensus_merger -> parse_repaa_table
    census_dir = script_dir.parent.parent.parent.parent.parent.parent / "18S_censusparse" / "csv_outputs"
    
    census_data = {}
    
    for level in ['family', 'genus', 'division']:
        census_file = census_dir / f"eukcensus_18S_by_{level}.csv"
        if census_file.exists():
            df = pd.read_csv(census_file)
            # Create lookup dictionary by taxon name
            census_lookup = {}
            for _, row in df.iterrows():
                taxon_name = row['Name_to_use']
                census_lookup[taxon_name] = {
                    'otu_count': row['otu_count'],
                    'otu_percentage': row['otu_percentage'],
                    'size_count': row['size_count'],
                    'size_percentage': row['size_percentage']
                }
            census_data[level] = census_lookup
            print(f"Loaded original census {level}: {len(census_lookup)} taxa")
    
    return census_data

def check_eukprot_percentages(merged_file, original_census, level):
    """Check that EukProt merged file percentages match original census data."""
    if not merged_file.exists():
        return f"File not found: {merged_file}"
    
    df = pd.read_csv(merged_file)
    
    # Target taxa to check
    target_taxa = ['Insecta', 'Sordariomycetes', 'Mammalia', 'Teleostei']
    
    results = []
    results.append(f"\n=== {level.upper()} LEVEL CENSUS PERCENTAGE VALIDATION ===")
    results.append(f"File: {merged_file.name}")
    results.append("")
    
    for taxon in target_taxa:
        # Find the taxon in the merged file
        taxon_row = df[df[level] == taxon]
        
        if taxon_row.empty:
            results.append(f"{taxon}: NOT FOUND in merged file")
            continue
        
        row = taxon_row.iloc[0]
        
        # Extract values from merged file
        merged_otu_count = row.get('census_otu_count', 0)
        merged_otu_percentage = row.get('otu_percentage', 0)
        merged_size_count = row.get('census_size_count', 0)
        merged_size_percentage = row.get('size_percentage', 0)
        
        # Get original census data
        if taxon in original_census[level]:
            original = original_census[level][taxon]
            
            results.append(f"{taxon}:")
            results.append(f"  OTU Count:")
            results.append(f"    Merged file: {merged_otu_count:,}")
            results.append(f"    Original:    {original['otu_count']:,}")
            results.append(f"    Match: {'✓' if merged_otu_count == original['otu_count'] else '✗'}")
            results.append(f"  ")
            results.append(f"  OTU Percentage:")
            results.append(f"    Merged file: {merged_otu_percentage:.2f}%")
            results.append(f"    Original:    {original['otu_percentage']:.2f}%")
            results.append(f"    Match: {'✓' if abs(merged_otu_percentage - original['otu_percentage']) < 0.01 else '✗'}")
            results.append(f"  ")
            results.append(f"  Size Count:")
            results.append(f"    Merged file: {merged_size_count:,}")
            results.append(f"    Original:    {original['size_count']:,}")
            results.append(f"    Match: {'✓' if merged_size_count == original['size_count'] else '✗'}")
            results.append(f"  ")
            results.append(f"  Size Percentage:")
            results.append(f"    Merged file: {merged_size_percentage:.2f}%")
            results.append(f"    Original:    {original['size_percentage']:.2f}%")
            results.append(f"    Match: {'✓' if abs(merged_size_percentage - original['size_percentage']) < 0.01 else '✗'}")
            results.append("")
        else:
            results.append(f"{taxon}: NOT FOUND in original census data")
            results.append("")
    
    return "\n".join(results)

def main():
    """Main function to run EukProt sanity checks."""
    print("🔍 EukProt Census Percentage Sanity Check")
    print("=" * 50)
    
    # Load original census data
    print("Loading original census data...")
    original_census = load_original_census_data()
    
    # Check merged files
    results_dir = Path(__file__).parent.parent  # Go up from sanity to csv_results
    output_lines = []
    
    output_lines.append("EUKPROT CENSUS PERCENTAGE SANITY CHECK RESULTS")
    output_lines.append("=" * 50)
    output_lines.append(f"Generated: {pd.Timestamp.now()}")
    output_lines.append("")
    
    # Check each taxonomic level
    for level in ['family', 'genus', 'division']:
        merged_file = results_dir / f"18s_eukprot_merged_{level}.csv"
        
        if level in original_census:
            result = check_eukprot_percentages(merged_file, original_census, level)
            output_lines.append(result)
            print(result)
        else:
            print(f"No original census data available for {level} level")
    
    # Save results to file
    output_file = results_dir / "sanity_check_eukprot_percentages.txt"
    with open(output_file, 'w') as f:
        f.write("\n".join(output_lines))
    
    print(f"\n📄 Results saved to: {output_file}")
    print("EukProt sanity check complete!")

if __name__ == "__main__":
    main()
