#!/usr/bin/env python3
"""
Taxonomic Combination Parser for 18S EukCensus Data

This script parses the exact taxonomic combinations (division|family|genus) from the 
18S EukCensus metadata file and groups entries by identical taxonomic strings.
For each unique combination, it provides:
- Number of rows (clusters) with that exact combination
- Total aggregated size (sequence count) for that combination
- Example centroid IDs

Created: 2025-01-15
"""

import pandas as pd
import os
from collections import defaultdict
from tqdm import tqdm
from datetime import datetime

def parse_taxonomic_combinations(input_file, output_dir):
    """
    Parse and group taxonomic combinations from TSV file
    
    Args:
        input_file (str): Path to input TSV file
        output_dir (str): Directory to save output files
    """
    
    print(f"📊 Reading TSV file: {input_file}")
    
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"✅ Loaded {len(df):,} records")
    print(f"📋 Columns: {list(df.columns)}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store combinations
    combinations = defaultdict(lambda: {
        'row_count': 0,
        'total_size': 0,
        'example_centroids': []
    })
    
    print("\n🔍 Processing taxonomic combinations...")
    
    # Process each row
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Grouping combinations"):
        division = row['division'] if pd.notna(row['division']) else 'NA'
        family = row['family'] if pd.notna(row['family']) else 'NA'
        genus = row['genus'] if pd.notna(row['genus']) else 'NA'
        size = row['size'] if pd.notna(row['size']) else 0
        centroid = row['centroid'] if pd.notna(row['centroid']) else 'NA'
        
        # Create combination key
        combination_key = f"{division}|{family}|{genus}"
        
        # Update statistics
        combinations[combination_key]['row_count'] += 1
        combinations[combination_key]['total_size'] += size
        
        # Store example centroids (limit to 3 examples)
        if len(combinations[combination_key]['example_centroids']) < 3:
            combinations[combination_key]['example_centroids'].append(centroid)
    
    print(f"✅ Found {len(combinations):,} unique taxonomic combinations")
    
    # Convert to DataFrame for easier handling
    results = []
    for combo_key, stats in combinations.items():
        division, family, genus = combo_key.split('|')
        results.append({
            'taxonomic_combination': combo_key,
            'division': division,
            'family': family,
            'genus': genus,
            'row_count': stats['row_count'],
            'total_size': stats['total_size'],
            'avg_size_per_cluster': round(stats['total_size'] / stats['row_count'], 2) if stats['row_count'] > 0 else 0,
            'example_centroids': '; '.join(stats['example_centroids'])
        })
    
    # Create DataFrame
    results_df = pd.DataFrame(results)

    # Custom hierarchical sorting: highest row_count first, then all from same division
    print("\n🔄 Applying hierarchical sorting by division groups...")

    # Create a list to store the final sorted order
    sorted_results = []
    processed_divisions = set()

    # Sort initially by row_count descending to get the order of divisions to process
    temp_sorted = results_df.sort_values('row_count', ascending=False)

    for _, row in temp_sorted.iterrows():
        division = row['division']

        # If we haven't processed this division yet
        if division not in processed_divisions:
            # Get all entries for this division, sorted by row_count descending
            division_entries = results_df[results_df['division'] == division].sort_values('row_count', ascending=False)

            # Add all entries from this division to our sorted results
            sorted_results.extend(division_entries.to_dict('records'))

            # Mark this division as processed
            processed_divisions.add(division)

    # Convert back to DataFrame with the new hierarchical order
    results_df = pd.DataFrame(sorted_results)
    results_df = results_df.reset_index(drop=True)
    
    # Add metadata row
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata_info = [
        f"# 18S EukCensus Taxonomic Combination Analysis",
        f"# Analysis date: {timestamp}",
        f"# Total original records: {len(df):,}",
        f"# Unique taxonomic combinations: {len(combinations):,}",
        f"# Total sequences: {df['size'].sum():,}",
        ""
    ]
    
    # Save detailed results
    output_file = os.path.join(output_dir, 'taxonomic_combinations_detailed.csv')
    with open(output_file, 'w') as f:
        # Write metadata
        for line in metadata_info:
            f.write(line + '\n')
        
        # Write CSV data
        results_df.to_csv(f, index=False)
    
    print(f"💾 Saved detailed results to: {output_file}")
    
    # Create summary statistics
    summary_stats = {
        'total_combinations': len(combinations),
        'total_records': len(df),
        'total_sequences': int(df['size'].sum()),
        'avg_clusters_per_combination': round(len(df) / len(combinations), 2),
        'avg_sequences_per_combination': round(df['size'].sum() / len(combinations), 2)
    }
    
    # Save summary
    summary_file = os.path.join(output_dir, 'taxonomic_combinations_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("18S EukCensus Taxonomic Combination Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis date: {timestamp}\n\n")
        
        f.write("OVERVIEW:\n")
        f.write(f"  Total unique combinations: {summary_stats['total_combinations']:,}\n")
        f.write(f"  Total records (clusters): {summary_stats['total_records']:,}\n")
        f.write(f"  Total sequences: {summary_stats['total_sequences']:,}\n")
        f.write(f"  Average clusters per combination: {summary_stats['avg_clusters_per_combination']}\n")
        f.write(f"  Average sequences per combination: {summary_stats['avg_sequences_per_combination']}\n\n")
        
        f.write("TOP 20 COMBINATIONS BY CLUSTER COUNT:\n")
        for i, row in results_df.head(20).iterrows():
            f.write(f"  {row.name + 1:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['row_count']:,} clusters, {row['total_size']:,} sequences)\n")
        
        f.write(f"\nTOP 20 COMBINATIONS BY SEQUENCE COUNT:\n")
        top_by_size = results_df.sort_values('total_size', ascending=False).head(20)
        for i, (idx, row) in enumerate(top_by_size.iterrows(), 1):
            f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['total_size']:,} sequences, {row['row_count']:,} clusters)\n")
    
    print(f"💾 Saved summary to: {summary_file}")
    
    # Show top 10 for immediate feedback
    print(f"\n📈 Top 10 combinations by cluster count:")
    for i, (idx, row) in enumerate(results_df.head(10).iterrows(), 1):
        print(f"  {i:2d}. {row['taxonomic_combination']:<50} "
              f"({row['row_count']:,} clusters, {row['total_size']:,} sequences)")
    
    return results_df, summary_stats

def main():
    """Main function"""
    
    # File paths - script is in metadata folder
    input_file = "eukcensus_18S.clusters.97.tsv"
    output_dir = "sanity_check"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file not found: {input_file}")
        return
    
    print("🚀 Starting 18S EukCensus taxonomic combination analysis...")
    print("=" * 70)
    
    # Parse taxonomic combinations
    results_df, summary_stats = parse_taxonomic_combinations(input_file, output_dir)
    
    print("\n" + "=" * 70)
    print("✅ Analysis completed successfully!")
    
    # Show output files
    print(f"\n📁 Output files created in '{output_dir}':")
    output_files = [
        'taxonomic_combinations_detailed.csv',
        'taxonomic_combinations_summary.txt'
    ]
    
    for filename in output_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            print(f"  • {filepath}")

if __name__ == "__main__":
    main()
