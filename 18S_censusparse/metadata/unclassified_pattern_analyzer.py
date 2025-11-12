#!/usr/bin/env python3
"""
Unclassified Pattern Analyzer for 18S EukCensus Data

This script analyzes patterns in unclassified (.U.) taxonomic entries from the 
taxonomic combinations analysis. It identifies:
- Combinations with .U. at different taxonomic levels
- Patterns of classification completeness
- Distribution of unclassified entries across major divisions

Created: 2025-01-15
"""

import pandas as pd
import os
from collections import defaultdict
from datetime import datetime

def analyze_unclassified_patterns(combinations_file, output_dir):
    """
    Analyze patterns in unclassified taxonomic entries
    
    Args:
        combinations_file (str): Path to taxonomic combinations CSV file
        output_dir (str): Directory to save output files
    """
    
    print(f"📊 Reading combinations file: {combinations_file}")
    
    # Read the combinations file, skipping metadata lines
    df = pd.read_csv(combinations_file, comment='#')
    
    print(f"✅ Loaded {len(df):,} taxonomic combinations")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze classification patterns
    patterns = {
        'fully_classified': [],
        'division_unclassified': [],
        'family_unclassified': [],
        'genus_unclassified': [],
        'multiple_unclassified': []
    }
    
    print("\n🔍 Analyzing classification patterns...")
    
    for idx, row in df.iterrows():
        division = row['division']
        family = row['family']
        genus = row['genus']
        
        # Check for .U. patterns
        div_unclass = '.U.' in division
        fam_unclass = '.U.' in family
        gen_unclass = '.U.' in genus
        
        # Categorize patterns
        if not (div_unclass or fam_unclass or gen_unclass):
            patterns['fully_classified'].append(row)
        elif div_unclass and not fam_unclass and not gen_unclass:
            patterns['division_unclassified'].append(row)
        elif not div_unclass and fam_unclass and not gen_unclass:
            patterns['family_unclassified'].append(row)
        elif not div_unclass and not fam_unclass and gen_unclass:
            patterns['genus_unclassified'].append(row)
        else:
            patterns['multiple_unclassified'].append(row)
    
    # Convert to DataFrames
    pattern_dfs = {}
    for pattern_name, pattern_data in patterns.items():
        if pattern_data:
            pattern_dfs[pattern_name] = pd.DataFrame(pattern_data)
        else:
            pattern_dfs[pattern_name] = pd.DataFrame()
    
    # Calculate statistics
    total_combinations = len(df)
    total_clusters = df['row_count'].sum()
    total_sequences = df['total_size'].sum()
    
    stats = {}
    for pattern_name, pattern_df in pattern_dfs.items():
        if not pattern_df.empty:
            stats[pattern_name] = {
                'combinations': len(pattern_df),
                'clusters': pattern_df['row_count'].sum(),
                'sequences': pattern_df['total_size'].sum(),
                'combinations_pct': round(len(pattern_df) / total_combinations * 100, 1),
                'clusters_pct': round(pattern_df['row_count'].sum() / total_clusters * 100, 1),
                'sequences_pct': round(pattern_df['total_size'].sum() / total_sequences * 100, 1)
            }
        else:
            stats[pattern_name] = {
                'combinations': 0, 'clusters': 0, 'sequences': 0,
                'combinations_pct': 0.0, 'clusters_pct': 0.0, 'sequences_pct': 0.0
            }
    
    # Create summary report
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    summary_file = os.path.join(output_dir, 'unclassified_patterns_summary.txt')
    
    with open(summary_file, 'w') as f:
        f.write("18S EukCensus Unclassified Pattern Analysis\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis date: {timestamp}\n\n")
        
        f.write("OVERVIEW:\n")
        f.write(f"  Total combinations analyzed: {total_combinations:,}\n")
        f.write(f"  Total clusters: {total_clusters:,}\n")
        f.write(f"  Total sequences: {total_sequences:,}\n\n")
        
        f.write("CLASSIFICATION PATTERNS:\n")
        for pattern_name, pattern_stats in stats.items():
            pattern_display = pattern_name.replace('_', ' ').title()
            f.write(f"  {pattern_display}:\n")
            f.write(f"    Combinations: {pattern_stats['combinations']:,} ({pattern_stats['combinations_pct']}%)\n")
            f.write(f"    Clusters: {pattern_stats['clusters']:,} ({pattern_stats['clusters_pct']}%)\n")
            f.write(f"    Sequences: {pattern_stats['sequences']:,} ({pattern_stats['sequences_pct']}%)\n\n")
        
        # Show top entries for each pattern
        for pattern_name, pattern_df in pattern_dfs.items():
            if not pattern_df.empty and len(pattern_df) > 0:
                pattern_display = pattern_name.replace('_', ' ').title()
                f.write(f"TOP 10 {pattern_display.upper()} BY CLUSTER COUNT:\n")
                top_entries = pattern_df.nlargest(10, 'row_count')
                for i, (idx, row) in enumerate(top_entries.iterrows(), 1):
                    f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                           f"({row['row_count']:,} clusters, {row['total_size']:,} sequences)\n")
                f.write("\n")
    
    print(f"💾 Saved summary to: {summary_file}")
    
    # Save detailed pattern files
    for pattern_name, pattern_df in pattern_dfs.items():
        if not pattern_df.empty:
            pattern_file = os.path.join(output_dir, f'{pattern_name}_combinations.csv')
            
            # Add metadata header
            with open(pattern_file, 'w') as f:
                f.write(f"# {pattern_name.replace('_', ' ').title()} Combinations\n")
                f.write(f"# Analysis date: {timestamp}\n")
                f.write(f"# Total combinations: {len(pattern_df):,}\n")
                f.write(f"# Total clusters: {pattern_df['row_count'].sum():,}\n")
                f.write(f"# Total sequences: {pattern_df['total_size'].sum():,}\n\n")
                
                # Write CSV data
                pattern_df.to_csv(f, index=False)
            
            print(f"💾 Saved {pattern_name} details to: {pattern_file}")
    
    # Analyze major division patterns
    print("\n🔍 Analyzing major division patterns...")
    
    major_divisions = {}
    for idx, row in df.iterrows():
        base_division = row['division'].split('.U.')[0]  # Get base division name
        if base_division not in major_divisions:
            major_divisions[base_division] = {
                'total_combinations': 0,
                'total_clusters': 0,
                'total_sequences': 0,
                'classified_combinations': 0,
                'classified_clusters': 0,
                'classified_sequences': 0
            }
        
        major_divisions[base_division]['total_combinations'] += 1
        major_divisions[base_division]['total_clusters'] += row['row_count']
        major_divisions[base_division]['total_sequences'] += row['total_size']
        
        # Check if fully classified
        if not any('.U.' in x for x in [row['division'], row['family'], row['genus']]):
            major_divisions[base_division]['classified_combinations'] += 1
            major_divisions[base_division]['classified_clusters'] += row['row_count']
            major_divisions[base_division]['classified_sequences'] += row['total_size']
    
    # Save division analysis
    division_file = os.path.join(output_dir, 'major_division_classification_rates.csv')
    division_data = []
    
    for division, stats in major_divisions.items():
        if stats['total_combinations'] > 0:
            classification_rate = round(stats['classified_combinations'] / stats['total_combinations'] * 100, 1)
            cluster_classification_rate = round(stats['classified_clusters'] / stats['total_clusters'] * 100, 1)
            sequence_classification_rate = round(stats['classified_sequences'] / stats['total_sequences'] * 100, 1)
            
            division_data.append({
                'division': division,
                'total_combinations': stats['total_combinations'],
                'classified_combinations': stats['classified_combinations'],
                'combination_classification_rate': classification_rate,
                'total_clusters': stats['total_clusters'],
                'classified_clusters': stats['classified_clusters'],
                'cluster_classification_rate': cluster_classification_rate,
                'total_sequences': stats['total_sequences'],
                'classified_sequences': stats['classified_sequences'],
                'sequence_classification_rate': sequence_classification_rate
            })
    
    division_df = pd.DataFrame(division_data)
    division_df = division_df.sort_values('total_clusters', ascending=False)
    
    with open(division_file, 'w') as f:
        f.write(f"# Major Division Classification Rates\n")
        f.write(f"# Analysis date: {timestamp}\n")
        f.write(f"# Shows classification completeness by major taxonomic division\n\n")
        division_df.to_csv(f, index=False)
    
    print(f"💾 Saved division analysis to: {division_file}")
    
    return stats, major_divisions

def main():
    """Main function"""
    
    # File paths
    combinations_file = "sanity_check/taxonomic_combinations_detailed.csv"
    output_dir = "sanity_check"
    
    # Check if input file exists
    if not os.path.exists(combinations_file):
        print(f"❌ Error: Input file not found: {combinations_file}")
        print("Please run taxonomic_combination_parser.py first.")
        return
    
    print("🚀 Starting unclassified pattern analysis...")
    print("=" * 60)
    
    # Analyze patterns
    stats, major_divisions = analyze_unclassified_patterns(combinations_file, output_dir)
    
    print("\n" + "=" * 60)
    print("✅ Analysis completed successfully!")
    
    # Show quick summary
    print(f"\n📈 Quick Summary:")
    for pattern_name, pattern_stats in stats.items():
        pattern_display = pattern_name.replace('_', ' ').title()
        print(f"  {pattern_display}: {pattern_stats['combinations']:,} combinations "
              f"({pattern_stats['combinations_pct']}%)")

if __name__ == "__main__":
    main()
