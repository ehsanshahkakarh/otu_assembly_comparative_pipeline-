#!/usr/bin/env python3
"""
Bacteria Tree Parser
====================

Parse bacterial species and extract taxonomic tree structure.
Generates phylum, class, order, family, genus level statistics.
"""

import pandas as pd
import sys
from pathlib import Path
from datetime import datetime


def extract_taxonomic_level(row, target_rank):
    """Extract taxon name and taxid for a specific rank."""
    lineage = str(row['lineage'])
    lineage_ranks = str(row['lineage_ranks'])
    lineage_taxids = str(row['lineage_taxids'])
    
    if pd.isna(lineage) or lineage == 'nan':
        return None, None
    
    taxa = lineage.split(';')
    ranks = lineage_ranks.split(';')
    taxids = lineage_taxids.split(';')
    
    for i, rank in enumerate(ranks):
        if rank.strip() == target_rank:
            return taxa[i].strip(), taxids[i].strip()
    
    return None, None


def parse_bacteria_tree(input_file, output_dir):
    """
    Parse bacterial species and generate taxonomic tree statistics.
    
    Args:
        input_file: Path to bacteria_species CSV file
        output_dir: Directory to save output files
    """
    print("="*70)
    print("BACTERIA TREE PARSER")
    print("="*70)
    
    # Load bacterial species
    print(f"\n📁 Loading bacterial species from: {input_file.name}")
    df = pd.read_csv(input_file)
    print(f"   Total bacterial species: {len(df):,}")
    print(f"   Total bacterial genomes: {df['total_genome_count'].sum():,}")
    
    # Extract taxonomic levels
    print("\n🌳 Extracting taxonomic tree levels...")
    
    taxonomic_levels = ['phylum', 'class', 'order', 'family', 'genus']
    
    for level in taxonomic_levels:
        print(f"   Extracting {level}...")
        df[f'{level}_name'] = None
        df[f'{level}_taxid'] = None
        
        for idx, row in df.iterrows():
            name, taxid = extract_taxonomic_level(row, level)
            df.at[idx, f'{level}_name'] = name
            df.at[idx, f'{level}_taxid'] = taxid
    
    # Generate statistics for each level
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    print("\n📊 Generating statistics for each taxonomic level...")
    
    for level in taxonomic_levels:
        level_df = df[df[f'{level}_name'].notna()].copy()
        
        stats = level_df.groupby(f'{level}_name').agg({
            'species_taxid': 'count',
            'total_genome_count': 'sum',
            'isolate_genome_count': 'sum',
            'uncultured_genome_count': 'sum',
            f'{level}_taxid': 'first'
        }).rename(columns={
            'species_taxid': 'species_count',
            'total_genome_count': 'total_genomes',
            'isolate_genome_count': 'isolate_genomes',
            'uncultured_genome_count': 'uncultured_genomes',
            f'{level}_taxid': 'taxid'
        })
        
        # Calculate percentages
        total_genomes = stats['total_genomes'].sum()
        stats['genome_percentage'] = (stats['total_genomes'] / total_genomes * 100).round(2)
        stats['isolate_percentage'] = (stats['isolate_genomes'] / stats['total_genomes'] * 100).round(2)
        
        # Sort by genome count
        stats = stats.sort_values('total_genomes', ascending=False)
        
        # Save
        output_file = output_dir / f'bacteria_{level}_{timestamp}.csv'
        stats.to_csv(output_file)
        
        print(f"   ✅ {level.capitalize():8s}: {len(stats):,} unique, {stats['total_genomes'].sum():,} genomes → {output_file.name}")
    
    # Save full dataset with all taxonomic levels
    full_output = output_dir / f'bacteria_full_taxonomy_{timestamp}.csv'
    df.to_csv(full_output, index=False)
    print(f"\n💾 Saved full bacterial dataset with taxonomy: {full_output.name}")
    
    # Print summary
    print("\n" + "="*70)
    print("BACTERIA TAXONOMIC SUMMARY")
    print("="*70)
    
    for level in taxonomic_levels:
        count = df[f'{level}_name'].nunique()
        with_assignment = df[f'{level}_name'].notna().sum()
        print(f"   {level.capitalize():8s}: {count:,} unique ({with_assignment:,} species assigned)")
    
    return df


def main():
    # Find the most recent bacteria_species file
    sanity_output = Path(__file__).parent / 'output'
    bacteria_files = list(sanity_output.glob('bacteria_species_*.csv'))
    
    if not bacteria_files:
        print("❌ No bacteria_species files found!")
        print("   Run domain_parser.py first to generate domain-specific files.")
        sys.exit(1)
    
    # Use the most recent file
    input_file = max(bacteria_files, key=lambda p: p.stat().st_mtime)
    
    # Run analysis
    df_bacteria = parse_bacteria_tree(input_file, sanity_output)
    
    print("\n✅ Bacteria tree parsing complete!")
    print(f"📁 Results saved to: {sanity_output}")


if __name__ == '__main__':
    main()

