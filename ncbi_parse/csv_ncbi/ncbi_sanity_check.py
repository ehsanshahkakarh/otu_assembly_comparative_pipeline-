#!/usr/bin/env python3
"""
NCBI CSV Sanity Check Script

Analyzes the NCBI accession files to provide statistics on:
1. Family groupings with counts and percentages
2. Isolate vs uncultured breakdown with percentages
3. Overall data quality metrics

Outputs results to a .txt file for easy review.
"""

import pandas as pd
from pathlib import Path
import sys

def find_ncbi_csv_files(search_dir=None):
    """Find NCBI CSV files in the specified directory or common locations."""
    if search_dir:
        search_path = Path(search_dir)
    else:
        # Try current directory first
        search_path = Path('.')

    # Look for CSV files in the specified directory
    csv_files = list(search_path.glob('ncbi_*_with_accessions.csv'))

    if not csv_files and not search_dir:
        # If not found in current directory, try common subdirectories
        common_paths = [
            Path('./csv_ncbi'),
            Path('../csv_ncbi'),
            Path('./ncbi_parse/py_ncbi/csv_ncbi'),
            Path('../ncbi_parse/py_ncbi/csv_ncbi'),
            Path('../../csv_ncbi')
        ]

        for path in common_paths:
            if path.exists():
                csv_files = list(path.glob('ncbi_*_with_accessions.csv'))
                if csv_files:
                    search_path = path
                    break

    return csv_files, search_path

def detect_taxonomic_level(df):
    """Detect which taxonomic level this file contains."""
    taxonomic_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']

    for level in taxonomic_levels:
        if level in df.columns:
            return level
    return None

def analyze_ncbi_files(search_dir=None):
    """Analyze NCBI CSV files and generate statistics."""

    print("NCBI CSV SANITY CHECK")
    print("=" * 50)

    # Find CSV files
    csv_files, search_path = find_ncbi_csv_files(search_dir)

    if not csv_files:
        print(f"No NCBI accession CSV files found in search directory: {search_path}")
        return

    print(f"Found {len(csv_files)} CSV files in: {search_path}")
    print()

    results = []
    results.append("NCBI CSV SANITY CHECK RESULTS")
    results.append("=" * 50)
    results.append("")

    for csv_file in sorted(csv_files):
        print(f"Analyzing: {csv_file.name}")
        results.append(f"FILE: {csv_file.name}")
        results.append("-" * 40)

        try:
            # Load the CSV file
            df = pd.read_csv(csv_file)
            total_entries = len(df)
            results.append(f"Total entries: {total_entries:,}")

            # Detect taxonomic level
            taxonomic_level = detect_taxonomic_level(df)
            if not taxonomic_level:
                results.append("ERROR: No recognized taxonomic level column found")
                results.append("")
                continue

            results.append(f"Taxonomic level: {taxonomic_level.upper()}")

            # Check required columns
            required_cols = [taxonomic_level, 'genome_source']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                results.append(f"ERROR: Missing columns: {missing_cols}")
                results.append("")
                continue
            
            # Taxonomic analysis
            results.append("")
            results.append(f"{taxonomic_level.upper()} BREAKDOWN:")
            taxon_counts = df[taxonomic_level].value_counts()

            # Top 20 taxa
            results.append(f"Top 20 {taxonomic_level}s:")
            for i, (taxon, count) in enumerate(taxon_counts.head(20).items(), 1):
                percentage = (count / total_entries) * 100
                results.append(f"  {i:2d}. {taxon}: {count:,} ({percentage:.2f}%)")

            # Taxonomic summary stats
            results.append("")
            results.append(f"{taxonomic_level.capitalize()} summary:")
            results.append(f"  Total unique {taxonomic_level}s: {len(taxon_counts):,}")
            results.append(f"  Largest {taxonomic_level}: {taxon_counts.iloc[0]:,} entries ({taxon_counts.iloc[0]/total_entries*100:.2f}%)")
            results.append(f"  Smallest {taxonomic_level}: {taxon_counts.iloc[-1]:,} entries ({taxon_counts.iloc[-1]/total_entries*100:.2f}%)")
            results.append(f"  Average entries per {taxonomic_level}: {taxon_counts.mean():.1f}")
            
            # Isolate vs Uncultured analysis
            results.append("")
            results.append("GENOME SOURCE BREAKDOWN:")
            source_counts = df['genome_source'].value_counts()
            
            for source, count in source_counts.items():
                percentage = (count / total_entries) * 100
                results.append(f"  {source.capitalize()}: {count:,} ({percentage:.2f}%)")
            
            # Domain-specific breakdowns
            if 'domain' in df.columns:
                results.append("")
                results.append(f"DOMAIN-SPECIFIC {taxonomic_level.upper()} BREAKDOWN:")

                domains = df['domain'].value_counts().index
                for domain in domains:
                    domain_df = df[df['domain'] == domain]
                    domain_total = len(domain_df)
                    domain_taxon_counts = domain_df[taxonomic_level].value_counts()

                    results.append("")
                    results.append(f"=== {domain.upper()} ({domain_total:,} entries) ===")
                    results.append(f"Top 10 {taxonomic_level}s in {domain}:")

                    for i, (taxon, count) in enumerate(domain_taxon_counts.head(10).items(), 1):
                        percentage = (count / domain_total) * 100
                        results.append(f"  {i:2d}. {taxon}: {count:,} ({percentage:.2f}%)")

                    # Isolate vs uncultured for this domain
                    domain_sources = domain_df['genome_source'].value_counts()
                    isolate_count = domain_sources.get('isolate', 0)
                    uncultured_count = domain_sources.get('uncultured', 0)
                    isolate_pct = (isolate_count / domain_total) * 100 if domain_total > 0 else 0
                    uncultured_pct = (uncultured_count / domain_total) * 100 if domain_total > 0 else 0

                    results.append(f"  Genome sources in {domain}:")
                    results.append(f"    Isolate: {isolate_count:,} ({isolate_pct:.1f}%)")
                    results.append(f"    Uncultured: {uncultured_count:,} ({uncultured_pct:.1f}%)")

            # Overall Cross-tabulation: Taxonomic level vs Genome Source (top 10 taxa)
            results.append("")
            results.append(f"OVERALL TOP 10 {taxonomic_level.upper()}S - ISOLATE vs UNCULTURED:")
            top_taxa = taxon_counts.head(10).index

            for taxon in top_taxa:
                taxon_df = df[df[taxonomic_level] == taxon]
                taxon_total = len(taxon_df)
                taxon_sources = taxon_df['genome_source'].value_counts()

                isolate_count = taxon_sources.get('isolate', 0)
                uncultured_count = taxon_sources.get('uncultured', 0)

                isolate_pct = (isolate_count / taxon_total) * 100 if taxon_total > 0 else 0
                uncultured_pct = (uncultured_count / taxon_total) * 100 if taxon_total > 0 else 0

                results.append(f"  {taxon}:")
                results.append(f"    Total: {taxon_total:,}")
                results.append(f"    Isolate: {isolate_count:,} ({isolate_pct:.1f}%)")
                results.append(f"    Uncultured: {uncultured_count:,} ({uncultured_pct:.1f}%)")
            
            # Data quality checks
            results.append("")
            results.append("DATA QUALITY:")

            # Check for missing values
            missing_taxon = df[taxonomic_level].isna().sum()
            missing_source = df['genome_source'].isna().sum()

            results.append(f"  Missing {taxonomic_level} names: {missing_taxon:,} ({missing_taxon/total_entries*100:.2f}%)")
            results.append(f"  Missing genome sources: {missing_source:,} ({missing_source/total_entries*100:.2f}%)")
            
            # Check for unusual values
            if 'organism_name' in df.columns:
                unique_organisms = df['organism_name'].nunique()
                results.append(f"  Unique organism names: {unique_organisms:,}")
            
            if 'domain' in df.columns:
                domain_counts = df['domain'].value_counts()
                results.append("  Domain distribution:")
                for domain, count in domain_counts.items():
                    percentage = (count / total_entries) * 100
                    results.append(f"    {domain}: {count:,} ({percentage:.2f}%)")
            
        except Exception as e:
            results.append(f"ERROR analyzing {csv_file.name}: {str(e)}")
        
        results.append("")
        results.append("")
    
    # Write results to file in the same directory as the CSV files
    output_file = search_path / 'ncbi_sanity_check_results.txt'
    with open(output_file, 'w') as f:
        f.write('\n'.join(results))

    print(f"Results saved to: {output_file}")
    print("Analysis complete!")

if __name__ == "__main__":
    # Allow command line argument for search directory
    search_directory = sys.argv[1] if len(sys.argv) > 1 else None
    analyze_ncbi_files(search_directory)
