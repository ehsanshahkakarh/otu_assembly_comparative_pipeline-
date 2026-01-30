#!/usr/bin/env python3
"""
Domain Parser
=============

Parse the species-grouped dataset by domain (Bacteria, Archaea, Eukaryota, Viruses).
Extracts domain information from lineage and generates statistics.
"""

import pandas as pd
import sys
from pathlib import Path
from datetime import datetime


def extract_domain(row):
    """
    Extract domain from lineage information.
    
    Returns one of: 'Bacteria', 'Archaea', 'Eukaryota', 'Viruses', or 'Unclassified'
    """
    lineage = str(row['lineage'])
    
    if pd.isna(lineage) or lineage == 'nan':
        return 'Unclassified'
    
    lineage_lower = lineage.lower()
    
    # Check for each domain
    if 'bacteria' in lineage_lower and 'archaea' not in lineage_lower:
        return 'Bacteria'
    elif 'archaea' in lineage_lower:
        return 'Archaea'
    elif 'eukaryota' in lineage_lower:
        return 'Eukaryota'
    elif 'virus' in lineage_lower or 'viridae' in lineage_lower:
        return 'Viruses'
    else:
        return 'Unclassified'


def parse_by_domain(input_file, output_dir):
    """
    Parse species file by domain and generate statistics.
    
    Args:
        input_file: Path to species_grouped CSV file
        output_dir: Directory to save output files
    """
    print("="*70)
    print("DOMAIN PARSER - Species Dataset Analysis")
    print("="*70)
    
    # Load data
    print(f"\n📁 Loading data from: {input_file}")
    df = pd.read_csv(input_file)
    print(f"   Total species: {len(df):,}")
    print(f"   Total genomes: {df['total_genome_count'].sum():,}")
    
    # Extract domain
    print("\n🔍 Extracting domain information...")
    df['domain'] = df.apply(extract_domain, axis=1)
    
    # Generate statistics
    print("\n📊 Generating domain statistics...")
    domain_stats = df.groupby('domain').agg({
        'species_taxid': 'count',
        'total_genome_count': 'sum',
        'isolate_genome_count': 'sum',
        'uncultured_genome_count': 'sum'
    }).rename(columns={
        'species_taxid': 'species_count',
        'total_genome_count': 'total_genomes',
        'isolate_genome_count': 'isolate_genomes',
        'uncultured_genome_count': 'uncultured_genomes'
    })
    
    # Calculate percentages
    total_genomes = domain_stats['total_genomes'].sum()
    domain_stats['genome_percentage'] = (domain_stats['total_genomes'] / total_genomes * 100).round(2)
    domain_stats['isolate_percentage'] = (domain_stats['isolate_genomes'] / domain_stats['total_genomes'] * 100).round(2)
    
    # Sort by genome count
    domain_stats = domain_stats.sort_values('total_genomes', ascending=False)
    
    # Display results
    print("\n" + "="*70)
    print("DOMAIN STATISTICS")
    print("="*70)
    print(domain_stats.to_string())
    
    # Save domain statistics
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stats_file = output_dir / f'domain_statistics_{timestamp}.csv'
    domain_stats.to_csv(stats_file)
    print(f"\n💾 Saved domain statistics: {stats_file}")
    
    # Save full dataset with domain column
    full_output = output_dir / f'species_with_domain_{timestamp}.csv'
    df.to_csv(full_output, index=False)
    print(f"💾 Saved full dataset with domain: {full_output}")
    
    # Generate domain-specific files
    print("\n📂 Generating domain-specific files...")
    for domain in ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']:
        domain_df = df[df['domain'] == domain]
        if len(domain_df) > 0:
            domain_file = output_dir / f'{domain.lower()}_species_{timestamp}.csv'
            domain_df.to_csv(domain_file, index=False)
            print(f"   💾 {domain}: {len(domain_df):,} species, {domain_df['total_genome_count'].sum():,} genomes → {domain_file.name}")
    
    return domain_stats, df


def main():
    # Find the most recent species_grouped file
    output_dir = Path(__file__).parent.parent / 'output'
    species_files = list(output_dir.glob('species_grouped_*.csv'))
    
    if not species_files:
        print("❌ No species_grouped files found in output directory!")
        sys.exit(1)
    
    # Use the most recent file
    input_file = max(species_files, key=lambda p: p.stat().st_mtime)
    
    # Create sanity_check output directory
    sanity_output = Path(__file__).parent / 'output'
    
    # Run analysis
    domain_stats, df_with_domain = parse_by_domain(input_file, sanity_output)
    
    print("\n✅ Domain parsing complete!")
    print(f"📁 Results saved to: {sanity_output}")


if __name__ == '__main__':
    main()

