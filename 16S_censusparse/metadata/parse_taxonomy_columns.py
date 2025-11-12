#!/usr/bin/env python3
"""
Parse taxonomy columns (phylum/divisions, family, genus) from EukCensus 16S metadata file
Created: 2025-08-18
"""

import pandas as pd
import os
from collections import Counter
from tqdm import tqdm
import argparse

def parse_taxonomy_columns(input_file, output_dir=None):
    """
    Parse phylum, family, and genus columns from the EukCensus 16S metadata file
    
    Args:
        input_file (str): Path to the input TSV file
        output_dir (str): Output directory for results (default: same as input file)
    """
    
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    
    print(f"📊 Reading metadata file: {input_file}")
    
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    print(f"✅ Loaded {len(df):,} records")
    print(f"📋 Columns: {list(df.columns)}")
    
    # Extract taxonomy columns (note: 'familiy' has a typo in the original file)
    taxonomy_cols = ['phylum', 'familiy', 'genus']

    print(f"\n🔍 Analyzing taxonomy columns...")

    # Create summary statistics for each taxonomic level
    results = {}

    for col in taxonomy_cols:
        if col in df.columns:
            print(f"\n--- {col.upper()} ---")

            # Count occurrences and sum sizes
            taxon_stats = df.groupby(col).agg({
                col: 'count',  # Number of clusters
                'size': 'sum'  # Total sequence count
            }).rename(columns={col: 'cluster_count', 'size': 'size_count'})

            # Sort by cluster count
            taxon_stats = taxon_stats.sort_values('cluster_count', ascending=False)

            total_clusters = len(df[col].dropna())
            total_sequences = df['size'].sum()
            unique_taxa = len(taxon_stats)

            print(f"Total clusters: {total_clusters:,}")
            print(f"Total sequences: {total_sequences:,}")
            print(f"Unique {col}: {unique_taxa:,}")
            print(f"Missing values: {df[col].isna().sum():,}")

            # Show top 10 most common
            print(f"\nTop 10 most common {col} (by cluster count):")
            for i, (taxon, row) in enumerate(taxon_stats.head(10).iterrows(), 1):
                cluster_pct = (row['cluster_count'] / total_clusters) * 100
                size_pct = (row['size_count'] / total_sequences) * 100
                print(f"  {i:2d}. {taxon:<30} {row['cluster_count']:>8,} clusters ({cluster_pct:5.1f}%) | {row['size_count']:>10,} seqs ({size_pct:5.1f}%)")

            # Store results
            results[col] = {
                'stats': taxon_stats,
                'total_clusters': total_clusters,
                'total_sequences': total_sequences,
                'unique_taxa': unique_taxa,
                'missing_values': df[col].isna().sum()
            }

            # Save detailed counts to file with both cluster and size counts
            output_file = os.path.join(output_dir, f"{col}_counts.csv")
            counts_df = pd.DataFrame({
                f'{col}_name': taxon_stats.index,
                f'{col}_cluster_count': taxon_stats['cluster_count'].values,
                f'{col}_cluster_percentage': (taxon_stats['cluster_count'].values / total_clusters) * 100,
                f'{col}_size_count': taxon_stats['size_count'].values,
                f'{col}_size_percentage': (taxon_stats['size_count'].values / total_sequences) * 100
            })
            counts_df.to_csv(output_file, index=False)
            print(f"💾 Saved detailed counts to: {output_file}")
    
    # Create a combined taxonomy summary
    print(f"\n📈 Creating combined taxonomy summary...")

    # Extract unique combinations with size information
    taxonomy_df = df[taxonomy_cols + ['size']].copy()
    taxonomy_df.columns = ['phylum', 'family', 'genus', 'size']  # Fix the typo

    # Count unique combinations and sum sizes
    combo_stats = taxonomy_df.groupby(['phylum', 'family', 'genus']).agg({
        'phylum': 'count',  # Number of clusters
        'size': 'sum'       # Total sequence count
    }).rename(columns={'phylum': 'cluster_count', 'size': 'size_count'})

    combo_stats = combo_stats.sort_values('cluster_count', ascending=False).reset_index()

    # Add percentages
    total_clusters = len(df)
    total_sequences = df['size'].sum()
    combo_stats['cluster_percentage'] = (combo_stats['cluster_count'] / total_clusters) * 100
    combo_stats['size_percentage'] = (combo_stats['size_count'] / total_sequences) * 100

    # Save combined taxonomy
    combined_output = os.path.join(output_dir, "taxonomy_combinations.csv")
    combo_stats.to_csv(combined_output, index=False)
    print(f"💾 Saved taxonomy combinations to: {combined_output}")

    print(f"\nTop 10 taxonomy combinations (by cluster count):")
    for i, row in combo_stats.head(10).iterrows():
        print(f"  {i+1:2d}. {row['phylum']:<20} | {row['family']:<25} | {row['genus']:<20} | {row['cluster_count']:>6,} clusters ({row['cluster_percentage']:5.1f}%) | {row['size_count']:>8,} seqs ({row['size_percentage']:5.1f}%)")
    
    # Create a summary report
    summary_output = os.path.join(output_dir, "taxonomy_summary.txt")
    with open(summary_output, 'w') as f:
        f.write("EukCensus 16S Taxonomy Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total records: {len(df):,}\n")
        f.write(f"Total sequences: {total_sequences:,}\n")
        f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for col in taxonomy_cols:
            if col in results:
                r = results[col]
                f.write(f"{col.upper()} Statistics:\n")
                f.write(f"  - Total clusters: {r['total_clusters']:,}\n")
                f.write(f"  - Total sequences: {r['total_sequences']:,}\n")
                f.write(f"  - Unique taxa: {r['unique_taxa']:,}\n")
                f.write(f"  - Missing values: {r['missing_values']:,}\n")
                f.write(f"  - Coverage: {((r['total_clusters'] - r['missing_values']) / len(df) * 100):.1f}%\n\n")

        f.write(f"Unique taxonomy combinations: {len(combo_stats):,}\n")

    print(f"💾 Saved summary report to: {summary_output}")

    return results, combo_stats

def main():
    parser = argparse.ArgumentParser(description='Parse taxonomy columns from EukCensus 16S metadata')
    parser.add_argument('--input', '-i', 
                       default='/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/16S_censusparse/metadata/eukcensus_16S.clusters.97.tsv',
                       help='Input TSV file path')
    parser.add_argument('--output', '-o', 
                       help='Output directory (default: same as input file)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file not found: {args.input}")
        return
    
    print(f"🚀 Starting taxonomy analysis...")
    results, combo_stats = parse_taxonomy_columns(args.input, args.output)
    print(f"✅ Analysis complete!")

if __name__ == "__main__":
    main()
