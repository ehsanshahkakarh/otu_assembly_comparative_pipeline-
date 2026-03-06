#!/usr/bin/env python3
"""
Validate Domain Splitter Output
================================

Quick validation script to verify domain splitting worked correctly.
"""

import pandas as pd
from pathlib import Path


def validate_domain_split():
    """Validate that domain splitting worked correctly."""
    
    output_dir = Path(__file__).parent / 'output'
    
    print("=" * 80)
    print("DOMAIN SPLITTER OUTPUT VALIDATION")
    print("=" * 80)
    
    # Load all domain files
    domains = {
        'Bacteria': output_dir / 'bacteria_species.csv',
        'Archaea': output_dir / 'archaea_species.csv',
        'Eukaryota': output_dir / 'eukaryota_species.csv',
        'Viruses': output_dir / 'viruses_species.csv',
        'Unknown': output_dir / 'unknown_species.csv'
    }
    
    total_species = 0
    total_genomes = 0
    
    print("\n✅ Domain File Validation:\n")
    
    for domain_name, filepath in domains.items():
        if not filepath.exists():
            print(f"❌ {domain_name}: File not found!")
            continue
        
        df = pd.read_csv(filepath)
        species_count = len(df)
        genome_count = df['total_genome_count'].sum()
        
        total_species += species_count
        total_genomes += genome_count
        
        # Verify domain column
        if 'domain' in df.columns:
            unique_domains = df['domain'].unique()
            if len(unique_domains) == 1 and unique_domains[0] == domain_name:
                domain_check = "✅"
            else:
                domain_check = f"⚠️  (found: {unique_domains})"
        else:
            domain_check = "❌ No domain column"
        
        print(f"{domain_name:12} {domain_check}")
        print(f"  Species: {species_count:,}")
        print(f"  Genomes: {genome_count:,}")
        print(f"  File size: {filepath.stat().st_size / 1024 / 1024:.2f} MB")
        print()
    
    print("=" * 80)
    print(f"TOTAL: {total_species:,} species, {total_genomes:,} genomes")
    print("=" * 80)
    
    # Load summary file
    summary_files = list(output_dir.glob('domain_summary_*.csv'))
    if summary_files:
        latest_summary = max(summary_files, key=lambda p: p.stat().st_mtime)
        print(f"\n✅ Summary file: {latest_summary.name}")
        
        summary_df = pd.read_csv(latest_summary)
        print("\nDomain Summary:")
        print(summary_df.to_string(index=False))
    else:
        print("\n⚠️  No summary file found")
    
    print("\n" + "=" * 80)
    print("✅ Validation complete!")
    print("=" * 80)


if __name__ == '__main__':
    validate_domain_split()

