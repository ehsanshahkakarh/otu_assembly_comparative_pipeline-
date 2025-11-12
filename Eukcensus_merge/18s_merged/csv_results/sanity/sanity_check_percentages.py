#!/usr/bin/env python3
"""
Sanity Check Script for NCBI Percentage Calculations
Created: 2025-01-13

This script validates the percentage calculations in the merged output files
by checking specific taxa against the original NCBI database totals.
"""

import pandas as pd
from pathlib import Path
import sys

def load_ncbi_totals():
    """Load NCBI database totals for percentage calculations."""
    # Use relative path to ncbi_parse directory
    script_dir = Path(__file__).resolve().parent
    # Go up from sanity -> csv_results -> 18s_merged -> Eukcensus_merge -> OTU_97_eukcensus_merger -> parse_repaa_table
    ncbi_dir = script_dir.parent.parent.parent.parent.parent.parent / "ncbi_parse" / "csv_ncbi"

    print(f"Looking for NCBI data in: {ncbi_dir}")
    if not ncbi_dir.exists():
        print(f"NCBI directory not found: {ncbi_dir}")
        return {}
    
    totals = {}
    
    for level in ['family', 'genus', 'phylum']:
        counts_file = ncbi_dir / f"ncbi_{level}_counts.csv"
        if counts_file.exists():
            df = pd.read_csv(counts_file)
            total_genomes = df[f'{level}_genome_count'].sum()
            total_species = df[f'{level}_species_count'].sum()
            
            totals[level] = {
                'total_genomes': total_genomes,
                'total_species': total_species
            }
            print(f"NCBI {level} database totals: {total_genomes:,} genomes, {total_species:,} species")
    
    return totals

def check_merged_percentages(merged_file, ncbi_totals, level):
    """Check percentage calculations in merged file."""
    if not merged_file.exists():
        return f"File not found: {merged_file}"
    
    df = pd.read_csv(merged_file)
    
    # Target taxa to check
    target_taxa = ['Insecta', 'Sordariomycetes', 'Mammalia', 'Teleostei']
    
    results = []
    results.append(f"\n=== {level.upper()} LEVEL ANALYSIS ===")
    results.append(f"File: {merged_file.name}")
    results.append(f"NCBI Database Totals: {ncbi_totals[level]['total_genomes']:,} genomes, {ncbi_totals[level]['total_species']:,} species")
    results.append("")
    
    for taxon in target_taxa:
        # Find the taxon in the merged file
        taxon_row = df[df[level] == taxon]
        
        if taxon_row.empty:
            results.append(f"{taxon}: NOT FOUND in merged file")
            continue
        
        row = taxon_row.iloc[0]
        
        # Extract values from merged file
        ncbi_genome_count = row.get('ncbi_genome_count', 0)
        ncbi_species_count = row.get('ncbi_species_count', 0)
        isolate_count = row.get('isolate_count', 0)
        
        genome_pct_db = row.get('genome_pct_db', 0)
        species_pct = row.get('species_pct', 0)
        isolate_percentage = row.get('isolate_percentage', 0)
        
        # Calculate expected percentages
        expected_genome_pct = (ncbi_genome_count / ncbi_totals[level]['total_genomes'] * 100) if ncbi_totals[level]['total_genomes'] > 0 else 0
        expected_species_pct = (ncbi_species_count / ncbi_totals[level]['total_species'] * 100) if ncbi_totals[level]['total_species'] > 0 else 0
        
        # Calculate isolate percentage (isolate_count / ncbi_genome_count * 100)
        expected_isolate_pct = (isolate_count / ncbi_genome_count * 100) if ncbi_genome_count > 0 else 0
        
        results.append(f"{taxon}:")
        results.append(f"  NCBI Genome Count: {ncbi_genome_count:,}")
        results.append(f"  NCBI Species Count: {ncbi_species_count:,}")
        results.append(f"  Isolate Count: {isolate_count:,}")
        results.append(f"  ")
        results.append(f"  Genome Percentage:")
        results.append(f"    Merged file: {genome_pct_db:.4f}%")
        results.append(f"    Expected:    {expected_genome_pct:.4f}%")
        results.append(f"    Match: {'✓' if abs(genome_pct_db - expected_genome_pct) < 0.01 else '✗'}")
        results.append(f"  ")
        results.append(f"  Species Percentage:")
        results.append(f"    Merged file: {species_pct:.4f}%")
        results.append(f"    Expected:    {expected_species_pct:.4f}%")
        results.append(f"    Match: {'✓' if abs(species_pct - expected_species_pct) < 0.01 else '✗'}")
        results.append(f"  ")
        results.append(f"  Isolate Percentage:")
        results.append(f"    Merged file: {isolate_percentage:.2f}%")
        results.append(f"    Expected:    {expected_isolate_pct:.2f}%")
        results.append(f"    Match: {'✓' if abs(isolate_percentage - expected_isolate_pct) < 0.01 else '✗'}")
        results.append("")
    
    return "\n".join(results)

def main():
    """Main function to run sanity checks."""
    print("🔍 NCBI Percentage Sanity Check")
    print("=" * 50)
    
    # Load NCBI database totals
    print("Loading NCBI database totals...")
    ncbi_totals = load_ncbi_totals()
    
    # Check merged files
    results_dir = Path(__file__).parent.parent  # Go up from sanity to csv_results
    output_lines = []
    
    output_lines.append("NCBI PERCENTAGE SANITY CHECK RESULTS")
    output_lines.append("=" * 50)
    output_lines.append(f"Generated: {pd.Timestamp.now()}")
    output_lines.append("")
    
    # Check each taxonomic level
    for level in ['family', 'genus', 'phylum']:
        merged_file = results_dir / f"18s_ncbi_merged_clean_{level}.csv"
        
        if level in ncbi_totals:
            result = check_merged_percentages(merged_file, ncbi_totals, level)
            output_lines.append(result)
            print(result)
        else:
            print(f"No NCBI totals available for {level} level")
    
    # Save results to file
    output_file = results_dir / "sanity_check_percentages.txt"
    with open(output_file, 'w') as f:
        f.write("\n".join(output_lines))
    
    print(f"\n📄 Results saved to: {output_file}")
    print("Sanity check complete!")

if __name__ == "__main__":
    main()
