#!/usr/bin/env python3
"""
18S EukCensus Eukaryota.U.division Sanity Check Parser
Created: 2025-08-18

This script parses the 18S EukCensus metadata file to identify and group all entries
with "Eukaryota.U.division" in the division column, along with their family and genus.

Groups entries by exact matches of division, family, and genus combinations and 
provides counts for each unique combination.

Usage:
    python eukaryota_u_division_sanity_check.py [--input-file FILE] [--output-dir DIR]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys
from tqdm import tqdm

class EukaryotaUDivisionParser:
    """Parser for Eukaryota.U.division entries in 18S EukCensus data."""
    
    def __init__(self, input_file=None, output_dir=None):
        self.script_dir = Path(__file__).resolve().parent
        self.input_file = Path(input_file) if input_file else self._find_input_file()
        self.output_dir = Path(output_dir) if output_dir else self.script_dir
        
        if not self.input_file or not self.input_file.exists():
            raise FileNotFoundError(f"Could not find input file: {self.input_file}")
            
        print(f"📁 Input file: {self.input_file}")
        print(f"📁 Output directory: {self.output_dir}")
    
    def _find_input_file(self):
        """Find the 18S EukCensus TSV file."""
        possible_files = [
            self.script_dir / "eukcensus_18S.clusters.97.tsv",
            Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/18S_censusparse/metadata/eukcensus_18S.clusters.97.tsv")
        ]
        
        for file_path in possible_files:
            if file_path.exists():
                return file_path
        return None
    
    def load_data(self):
        """Load the 18S EukCensus data."""
        print("📂 Loading 18S EukCensus data...")
        
        # Load the TSV file
        df = pd.read_csv(self.input_file, sep='\t', low_memory=False)
        
        print(f"✅ Loaded {len(df):,} total records")
        print(f"📋 Columns: {list(df.columns)}")
        
        # Check for required columns
        required_columns = ['division', 'family', 'genus', 'size']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        return df
    
    def filter_eukaryota_u_division(self, df):
        """Filter for entries with 'Eukaryota.U.division' in the division column."""
        print("🔍 Filtering for 'Eukaryota.U.division' entries...")
        
        # Filter for exact match of "Eukaryota.U.division"
        euk_u_mask = df['division'] == 'Eukaryota.U.division'
        euk_u_df = df[euk_u_mask].copy()
        
        print(f"✅ Found {len(euk_u_df):,} entries with 'Eukaryota.U.division'")
        print(f"📊 This represents {len(euk_u_df)/len(df)*100:.2f}% of total records")
        
        if len(euk_u_df) == 0:
            print("❌ No entries found with 'Eukaryota.U.division'")
            return euk_u_df
        
        # Show some basic statistics
        total_sequences = euk_u_df['size'].sum()
        print(f"📈 Total sequences in Eukaryota.U.division entries: {total_sequences:,}")
        print(f"📈 Average sequences per entry: {total_sequences/len(euk_u_df):.1f}")
        
        return euk_u_df
    
    def group_by_taxonomy(self, euk_u_df):
        """Group entries by exact division, family, genus combinations."""
        print("🔄 Grouping by division, family, and genus combinations...")
        
        # Group by the three taxonomic columns and aggregate
        grouped = euk_u_df.groupby(['division', 'family', 'genus']).agg({
            'centroid': 'count',  # Number of clusters/entries
            'size': 'sum'         # Total sequence count
        }).reset_index()
        
        # Rename columns for clarity
        grouped.columns = ['division', 'family', 'genus', 'cluster_count', 'total_sequences']
        
        # Sort by cluster count (descending) to show most common combinations first
        grouped = grouped.sort_values('cluster_count', ascending=False)
        
        print(f"✅ Found {len(grouped):,} unique division-family-genus combinations")
        
        # Show summary statistics
        print(f"\n📊 Grouping statistics:")
        print(f"   Total unique combinations: {len(grouped):,}")
        print(f"   Total clusters: {grouped['cluster_count'].sum():,}")
        print(f"   Total sequences: {grouped['total_sequences'].sum():,}")
        print(f"   Average clusters per combination: {grouped['cluster_count'].mean():.1f}")
        print(f"   Average sequences per combination: {grouped['total_sequences'].mean():.1f}")
        
        return grouped
    
    def analyze_patterns(self, grouped_df):
        """Analyze patterns in the grouped data."""
        print("\n🔍 Analyzing taxonomic patterns...")
        
        # Check family patterns
        family_counts = grouped_df['family'].value_counts()
        print(f"\n📋 Family distribution:")
        print(f"   Unique families: {len(family_counts)}")
        print(f"   Top 10 families:")
        for i, (family, count) in enumerate(family_counts.head(10).items(), 1):
            total_clusters = grouped_df[grouped_df['family'] == family]['cluster_count'].sum()
            total_seqs = grouped_df[grouped_df['family'] == family]['total_sequences'].sum()
            print(f"     {i:2d}. {family}: {count} combinations, {total_clusters:,} clusters, {total_seqs:,} sequences")
        
        # Check genus patterns
        genus_counts = grouped_df['genus'].value_counts()
        print(f"\n📋 Genus distribution:")
        print(f"   Unique genera: {len(genus_counts)}")
        print(f"   Top 10 genera:")
        for i, (genus, count) in enumerate(genus_counts.head(10).items(), 1):
            total_clusters = grouped_df[grouped_df['genus'] == genus]['cluster_count'].sum()
            total_seqs = grouped_df[grouped_df['genus'] == genus]['total_sequences'].sum()
            print(f"     {i:2d}. {genus}: {count} combinations, {total_clusters:,} clusters, {total_seqs:,} sequences")
        
        # Show top combinations by cluster count
        print(f"\n🔝 Top 15 combinations by cluster count:")
        for i, row in grouped_df.head(15).iterrows():
            print(f"   {i+1:2d}. Division: {row['division']}")
            print(f"       Family: {row['family']}")
            print(f"       Genus: {row['genus']}")
            print(f"       Clusters: {row['cluster_count']:,}, Sequences: {row['total_sequences']:,}")
            print()
    
    def save_results(self, grouped_df, euk_u_df):
        """Save the analysis results to files."""
        print("💾 Saving results...")
        
        # Save grouped results
        grouped_output = self.output_dir / "eukaryota_u_division_grouped.csv"
        grouped_df.to_csv(grouped_output, index=False)
        print(f"✅ Saved grouped results: {grouped_output}")
        
        # Save all Eukaryota.U.division entries for reference
        all_entries_output = self.output_dir / "eukaryota_u_division_all_entries.csv"
        # Select relevant columns for the output
        output_columns = ['centroid', 'members', 'size', 'division', 'family', 'genus']
        euk_u_df[output_columns].to_csv(all_entries_output, index=False)
        print(f"✅ Saved all entries: {all_entries_output}")
        
        # Create summary report
        summary_output = self.output_dir / "eukaryota_u_division_summary.txt"
        with open(summary_output, 'w') as f:
            f.write("18S EukCensus Eukaryota.U.division Analysis Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("OVERVIEW:\n")
            f.write(f"  Total Eukaryota.U.division entries: {len(euk_u_df):,}\n")
            f.write(f"  Total sequences: {euk_u_df['size'].sum():,}\n")
            f.write(f"  Unique taxonomic combinations: {len(grouped_df):,}\n")
            f.write(f"  Average sequences per entry: {euk_u_df['size'].mean():.1f}\n\n")
            
            f.write("TAXONOMIC DIVERSITY:\n")
            f.write(f"  Unique families: {grouped_df['family'].nunique():,}\n")
            f.write(f"  Unique genera: {grouped_df['genus'].nunique():,}\n\n")
            
            f.write("TOP 10 COMBINATIONS BY CLUSTER COUNT:\n")
            for i, row in grouped_df.head(10).iterrows():
                f.write(f"  {i+1:2d}. {row['family']} | {row['genus']} | {row['cluster_count']:,} clusters | {row['total_sequences']:,} sequences\n")
        
        print(f"✅ Saved summary report: {summary_output}")
        
        return grouped_output, all_entries_output, summary_output
    
    def run_analysis(self):
        """Run the complete analysis."""
        try:
            print("🚀 Starting Eukaryota.U.division analysis...")
            
            # Load data
            df = self.load_data()
            
            # Filter for Eukaryota.U.division entries
            euk_u_df = self.filter_eukaryota_u_division(df)
            
            if len(euk_u_df) == 0:
                print("❌ No data to analyze")
                return
            
            # Group by taxonomy
            grouped_df = self.group_by_taxonomy(euk_u_df)
            
            # Analyze patterns
            self.analyze_patterns(grouped_df)
            
            # Save results
            output_files = self.save_results(grouped_df, euk_u_df)
            
            print(f"\n" + "="*60)
            print("ANALYSIS COMPLETE")
            print("="*60)
            print(f"📁 Output files:")
            for file_path in output_files:
                print(f"   - {file_path}")
            
            print(f"\n🎯 Key findings:")
            print(f"   - {len(euk_u_df):,} entries with 'Eukaryota.U.division'")
            print(f"   - {len(grouped_df):,} unique taxonomic combinations")
            print(f"   - {grouped_df['family'].nunique():,} unique families")
            print(f"   - {grouped_df['genus'].nunique():,} unique genera")
            
        except Exception as e:
            print(f"❌ Error during analysis: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description='18S EukCensus Eukaryota.U.division sanity check')
    parser.add_argument('--input-file', help='Input TSV file path')
    parser.add_argument('--output-dir', help='Output directory for results')
    
    args = parser.parse_args()
    
    analyzer = EukaryotaUDivisionParser(
        input_file=args.input_file,
        output_dir=args.output_dir
    )
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
