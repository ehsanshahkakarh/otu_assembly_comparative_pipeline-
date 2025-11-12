#!/usr/bin/env python3
"""
Size Distribution Analyzer for 18S EukCensus Data

This script analyzes the size distribution patterns in the taxonomic combinations
to identify potential outliers, unusual patterns, and provide statistics on
cluster sizes and sequence counts.

Created: 2025-01-15
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

def analyze_size_distributions(combinations_file, output_dir):
    """
    Analyze size distribution patterns in taxonomic combinations
    
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
    
    # Calculate basic statistics
    print("\n🔍 Calculating size distribution statistics...")
    
    # Cluster count statistics
    cluster_stats = {
        'mean': df['row_count'].mean(),
        'median': df['row_count'].median(),
        'std': df['row_count'].std(),
        'min': df['row_count'].min(),
        'max': df['row_count'].max(),
        'q25': df['row_count'].quantile(0.25),
        'q75': df['row_count'].quantile(0.75)
    }
    
    # Sequence count statistics
    sequence_stats = {
        'mean': df['total_size'].mean(),
        'median': df['total_size'].median(),
        'std': df['total_size'].std(),
        'min': df['total_size'].min(),
        'max': df['total_size'].max(),
        'q25': df['total_size'].quantile(0.25),
        'q75': df['total_size'].quantile(0.75)
    }
    
    # Average size per cluster statistics
    avg_size_stats = {
        'mean': df['avg_size_per_cluster'].mean(),
        'median': df['avg_size_per_cluster'].median(),
        'std': df['avg_size_per_cluster'].std(),
        'min': df['avg_size_per_cluster'].min(),
        'max': df['avg_size_per_cluster'].max(),
        'q25': df['avg_size_per_cluster'].quantile(0.25),
        'q75': df['avg_size_per_cluster'].quantile(0.75)
    }
    
    # Identify outliers using IQR method
    print("\n🔍 Identifying outliers...")
    
    def find_outliers(series, column_name):
        Q1 = series.quantile(0.25)
        Q3 = series.quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        
        outliers = df[(series < lower_bound) | (series > upper_bound)]
        return outliers, lower_bound, upper_bound
    
    cluster_outliers, cluster_lower, cluster_upper = find_outliers(df['row_count'], 'row_count')
    sequence_outliers, seq_lower, seq_upper = find_outliers(df['total_size'], 'total_size')
    avg_size_outliers, avg_lower, avg_upper = find_outliers(df['avg_size_per_cluster'], 'avg_size_per_cluster')
    
    # Create summary report
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    summary_file = os.path.join(output_dir, 'size_distribution_summary.txt')
    
    with open(summary_file, 'w') as f:
        f.write("18S EukCensus Size Distribution Analysis\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis date: {timestamp}\n\n")
        
        f.write("CLUSTER COUNT STATISTICS:\n")
        f.write(f"  Mean: {cluster_stats['mean']:.1f}\n")
        f.write(f"  Median: {cluster_stats['median']:.1f}\n")
        f.write(f"  Standard Deviation: {cluster_stats['std']:.1f}\n")
        f.write(f"  Min: {cluster_stats['min']:,}\n")
        f.write(f"  Max: {cluster_stats['max']:,}\n")
        f.write(f"  25th Percentile: {cluster_stats['q25']:.1f}\n")
        f.write(f"  75th Percentile: {cluster_stats['q75']:.1f}\n")
        f.write(f"  Outliers (IQR method): {len(cluster_outliers)} combinations\n\n")
        
        f.write("SEQUENCE COUNT STATISTICS:\n")
        f.write(f"  Mean: {sequence_stats['mean']:.1f}\n")
        f.write(f"  Median: {sequence_stats['median']:.1f}\n")
        f.write(f"  Standard Deviation: {sequence_stats['std']:.1f}\n")
        f.write(f"  Min: {sequence_stats['min']:,}\n")
        f.write(f"  Max: {sequence_stats['max']:,}\n")
        f.write(f"  25th Percentile: {sequence_stats['q25']:.1f}\n")
        f.write(f"  75th Percentile: {sequence_stats['q75']:.1f}\n")
        f.write(f"  Outliers (IQR method): {len(sequence_outliers)} combinations\n\n")
        
        f.write("AVERAGE SIZE PER CLUSTER STATISTICS:\n")
        f.write(f"  Mean: {avg_size_stats['mean']:.2f}\n")
        f.write(f"  Median: {avg_size_stats['median']:.2f}\n")
        f.write(f"  Standard Deviation: {avg_size_stats['std']:.2f}\n")
        f.write(f"  Min: {avg_size_stats['min']:.2f}\n")
        f.write(f"  Max: {avg_size_stats['max']:.2f}\n")
        f.write(f"  25th Percentile: {avg_size_stats['q25']:.2f}\n")
        f.write(f"  75th Percentile: {avg_size_stats['q75']:.2f}\n")
        f.write(f"  Outliers (IQR method): {len(avg_size_outliers)} combinations\n\n")
        
        # Top combinations by different metrics
        f.write("TOP 10 BY CLUSTER COUNT:\n")
        top_clusters = df.nlargest(10, 'row_count')
        for i, (idx, row) in enumerate(top_clusters.iterrows(), 1):
            f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['row_count']:,} clusters)\n")
        
        f.write("\nTOP 10 BY SEQUENCE COUNT:\n")
        top_sequences = df.nlargest(10, 'total_size')
        for i, (idx, row) in enumerate(top_sequences.iterrows(), 1):
            f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['total_size']:,} sequences)\n")
        
        f.write("\nTOP 10 BY AVERAGE SIZE PER CLUSTER:\n")
        top_avg_size = df.nlargest(10, 'avg_size_per_cluster')
        for i, (idx, row) in enumerate(top_avg_size.iterrows(), 1):
            f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['avg_size_per_cluster']:.2f} avg size)\n")
        
        f.write("\nSMALLEST 10 BY AVERAGE SIZE PER CLUSTER:\n")
        bottom_avg_size = df.nsmallest(10, 'avg_size_per_cluster')
        for i, (idx, row) in enumerate(bottom_avg_size.iterrows(), 1):
            f.write(f"  {i:2d}. {row['taxonomic_combination']:<60} "
                   f"({row['avg_size_per_cluster']:.2f} avg size)\n")
    
    print(f"💾 Saved summary to: {summary_file}")
    
    # Save outlier files
    outlier_files = [
        (cluster_outliers, 'cluster_count_outliers.csv', 'row_count'),
        (sequence_outliers, 'sequence_count_outliers.csv', 'total_size'),
        (avg_size_outliers, 'avg_size_outliers.csv', 'avg_size_per_cluster')
    ]
    
    for outlier_df, filename, metric in outlier_files:
        if not outlier_df.empty:
            outlier_file = os.path.join(output_dir, filename)
            
            with open(outlier_file, 'w') as f:
                f.write(f"# {metric.replace('_', ' ').title()} Outliers\n")
                f.write(f"# Analysis date: {timestamp}\n")
                f.write(f"# Total outliers: {len(outlier_df):,}\n")
                f.write(f"# Outliers identified using IQR method (1.5 * IQR)\n\n")
                
                # Sort by the metric in descending order
                outlier_df_sorted = outlier_df.sort_values(metric, ascending=False)
                outlier_df_sorted.to_csv(f, index=False)
            
            print(f"💾 Saved {len(outlier_df)} {metric} outliers to: {outlier_file}")
    
    # Analyze size distribution by major taxonomic groups
    print("\n🔍 Analyzing size patterns by major taxonomic groups...")
    
    # Group by major division (first part before .U.)
    df['major_division'] = df['division'].str.split('.U.').str[0]
    
    division_stats = []
    for division in df['major_division'].unique():
        div_data = df[df['major_division'] == division]
        
        if len(div_data) > 0:
            division_stats.append({
                'division': division,
                'combinations': len(div_data),
                'total_clusters': div_data['row_count'].sum(),
                'total_sequences': div_data['total_size'].sum(),
                'avg_clusters_per_combination': div_data['row_count'].mean(),
                'avg_sequences_per_combination': div_data['total_size'].mean(),
                'avg_size_per_cluster': div_data['avg_size_per_cluster'].mean(),
                'median_clusters': div_data['row_count'].median(),
                'median_sequences': div_data['total_size'].median(),
                'median_avg_size': div_data['avg_size_per_cluster'].median()
            })
    
    division_stats_df = pd.DataFrame(division_stats)
    division_stats_df = division_stats_df.sort_values('total_clusters', ascending=False)
    
    division_file = os.path.join(output_dir, 'division_size_statistics.csv')
    with open(division_file, 'w') as f:
        f.write(f"# Division Size Statistics\n")
        f.write(f"# Analysis date: {timestamp}\n")
        f.write(f"# Size distribution statistics by major taxonomic division\n\n")
        division_stats_df.to_csv(f, index=False)
    
    print(f"💾 Saved division statistics to: {division_file}")
    
    return {
        'cluster_stats': cluster_stats,
        'sequence_stats': sequence_stats,
        'avg_size_stats': avg_size_stats,
        'outlier_counts': {
            'cluster_outliers': len(cluster_outliers),
            'sequence_outliers': len(sequence_outliers),
            'avg_size_outliers': len(avg_size_outliers)
        }
    }

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
    
    print("🚀 Starting size distribution analysis...")
    print("=" * 60)
    
    # Analyze size distributions
    stats = analyze_size_distributions(combinations_file, output_dir)
    
    print("\n" + "=" * 60)
    print("✅ Analysis completed successfully!")
    
    # Show quick summary
    print(f"\n📈 Quick Summary:")
    print(f"  Cluster count - Mean: {stats['cluster_stats']['mean']:.1f}, "
          f"Median: {stats['cluster_stats']['median']:.1f}, "
          f"Max: {stats['cluster_stats']['max']:,}")
    print(f"  Sequence count - Mean: {stats['sequence_stats']['mean']:.1f}, "
          f"Median: {stats['sequence_stats']['median']:.1f}, "
          f"Max: {stats['sequence_stats']['max']:,}")
    print(f"  Avg size per cluster - Mean: {stats['avg_size_stats']['mean']:.2f}, "
          f"Median: {stats['avg_size_stats']['median']:.2f}")
    print(f"  Outliers found: {stats['outlier_counts']['cluster_outliers']} cluster, "
          f"{stats['outlier_counts']['sequence_outliers']} sequence, "
          f"{stats['outlier_counts']['avg_size_outliers']} avg size")

if __name__ == "__main__":
    main()
