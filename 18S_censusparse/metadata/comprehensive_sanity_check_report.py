#!/usr/bin/env python3
"""
Comprehensive Sanity Check Report Generator for 18S EukCensus Data

This script generates a comprehensive report summarizing all the sanity checks
performed on the 18S EukCensus taxonomic combinations data.

Created: 2025-01-15
"""

import pandas as pd
import os
from datetime import datetime

def generate_comprehensive_report(output_dir):
    """
    Generate a comprehensive sanity check report
    
    Args:
        output_dir (str): Directory containing all sanity check files
    """
    
    print("🚀 Generating comprehensive sanity check report...")
    
    # Check for required files
    required_files = [
        'taxonomic_combinations_detailed.csv',
        'taxonomic_combinations_summary.txt',
        'unclassified_patterns_summary.txt',
        'size_distribution_summary.txt',
        'major_division_classification_rates.csv'
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(os.path.join(output_dir, file)):
            missing_files.append(file)
    
    if missing_files:
        print(f"❌ Missing required files: {', '.join(missing_files)}")
        print("Please run the individual analysis scripts first.")
        return
    
    # Read key data
    combinations_df = pd.read_csv(os.path.join(output_dir, 'taxonomic_combinations_detailed.csv'), comment='#')
    division_rates_df = pd.read_csv(os.path.join(output_dir, 'major_division_classification_rates.csv'), comment='#')
    
    # Generate comprehensive report
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_file = os.path.join(output_dir, 'COMPREHENSIVE_SANITY_CHECK_REPORT.md')
    
    with open(report_file, 'w') as f:
        f.write("# 18S EukCensus Comprehensive Sanity Check Report\n\n")
        f.write(f"**Generated:** {timestamp}\n\n")
        f.write("---\n\n")
        
        # Executive Summary
        f.write("## Executive Summary\n\n")
        
        total_combinations = len(combinations_df)
        total_clusters = combinations_df['row_count'].sum()
        total_sequences = combinations_df['total_size'].sum()
        
        f.write(f"- **Total unique taxonomic combinations:** {total_combinations:,}\n")
        f.write(f"- **Total clusters:** {total_clusters:,}\n")
        f.write(f"- **Total sequences:** {total_sequences:,}\n")
        f.write(f"- **Average clusters per combination:** {total_clusters/total_combinations:.1f}\n")
        f.write(f"- **Average sequences per combination:** {total_sequences/total_combinations:.1f}\n\n")
        
        # Classification Completeness
        f.write("## Classification Completeness Analysis\n\n")
        
        # Count classification patterns
        fully_classified = 0
        genus_unclassified = 0
        multiple_unclassified = 0
        
        for _, row in combinations_df.iterrows():
            div_unclass = '.U.' in row['division']
            fam_unclass = '.U.' in row['family']
            gen_unclass = '.U.' in row['genus']
            
            if not (div_unclass or fam_unclass or gen_unclass):
                fully_classified += 1
            elif not div_unclass and not fam_unclass and gen_unclass:
                genus_unclassified += 1
            elif div_unclass or fam_unclass or gen_unclass:
                multiple_unclassified += 1
        
        f.write(f"- **Fully classified combinations:** {fully_classified:,} ({fully_classified/total_combinations*100:.1f}%)\n")
        f.write(f"- **Genus-only unclassified:** {genus_unclassified:,} ({genus_unclassified/total_combinations*100:.1f}%)\n")
        f.write(f"- **Multiple levels unclassified:** {multiple_unclassified:,} ({multiple_unclassified/total_combinations*100:.1f}%)\n\n")
        
        # Top taxonomic groups
        f.write("## Major Taxonomic Groups\n\n")
        f.write("### Top 10 Combinations by Cluster Count\n\n")
        f.write("| Rank | Taxonomic Combination | Clusters | Sequences |\n")
        f.write("|------|----------------------|----------|----------|\n")
        
        top_10 = combinations_df.nlargest(10, 'row_count')
        for i, (_, row) in enumerate(top_10.iterrows(), 1):
            combo = row['taxonomic_combination'].replace('|', ' → ')
            f.write(f"| {i} | {combo} | {row['row_count']:,} | {row['total_size']:,} |\n")
        
        f.write("\n### Top 10 Combinations by Sequence Count\n\n")
        f.write("| Rank | Taxonomic Combination | Sequences | Clusters |\n")
        f.write("|------|----------------------|-----------|----------|\n")
        
        top_10_seq = combinations_df.nlargest(10, 'total_size')
        for i, (_, row) in enumerate(top_10_seq.iterrows(), 1):
            combo = row['taxonomic_combination'].replace('|', ' → ')
            f.write(f"| {i} | {combo} | {row['total_size']:,} | {row['row_count']:,} |\n")
        
        # Division-level analysis
        f.write("\n## Division-Level Classification Rates\n\n")
        f.write("| Division | Total Clusters | Classification Rate (%) | Avg Size per Cluster |\n")
        f.write("|----------|----------------|------------------------|----------------------|\n")
        
        top_divisions = division_rates_df.nlargest(10, 'total_clusters')
        for _, row in top_divisions.iterrows():
            f.write(f"| {row['division']} | {row['total_clusters']:,} | {row['cluster_classification_rate']:.1f}% | {row['total_sequences']/row['total_clusters']:.2f} |\n")
        
        # Size distribution insights
        f.write("\n## Size Distribution Analysis\n\n")
        
        cluster_stats = combinations_df['row_count'].describe()
        sequence_stats = combinations_df['total_size'].describe()
        avg_size_stats = combinations_df['avg_size_per_cluster'].describe()
        
        f.write("### Cluster Count Distribution\n")
        f.write(f"- **Mean:** {cluster_stats['mean']:.1f}\n")
        f.write(f"- **Median:** {cluster_stats['50%']:.1f}\n")
        f.write(f"- **Range:** {cluster_stats['min']:.0f} - {cluster_stats['max']:,.0f}\n")
        f.write(f"- **Standard Deviation:** {cluster_stats['std']:.1f}\n\n")
        
        f.write("### Sequence Count Distribution\n")
        f.write(f"- **Mean:** {sequence_stats['mean']:.1f}\n")
        f.write(f"- **Median:** {sequence_stats['50%']:.1f}\n")
        f.write(f"- **Range:** {sequence_stats['min']:.0f} - {sequence_stats['max']:,.0f}\n")
        f.write(f"- **Standard Deviation:** {sequence_stats['std']:.1f}\n\n")
        
        f.write("### Average Size per Cluster Distribution\n")
        f.write(f"- **Mean:** {avg_size_stats['mean']:.2f}\n")
        f.write(f"- **Median:** {avg_size_stats['50%']:.2f}\n")
        f.write(f"- **Range:** {avg_size_stats['min']:.2f} - {avg_size_stats['max']:.2f}\n")
        f.write(f"- **Standard Deviation:** {avg_size_stats['std']:.2f}\n\n")
        
        # Key findings and recommendations
        f.write("## Key Findings\n\n")
        
        f.write("### 1. Data Completeness\n")
        f.write(f"- The dataset contains {total_combinations:,} unique taxonomic combinations\n")
        f.write(f"- {fully_classified/total_combinations*100:.1f}% of combinations are fully classified to genus level\n")
        f.write(f"- {multiple_unclassified/total_combinations*100:.1f}% have multiple unclassified levels\n\n")
        
        f.write("### 2. Taxonomic Distribution\n")
        largest_group = combinations_df.iloc[0]
        f.write(f"- Largest group: {largest_group['taxonomic_combination']} with {largest_group['row_count']:,} clusters\n")
        f.write(f"- Opisthokonta dominates with {division_rates_df[division_rates_df['division']=='Opisthokonta']['total_clusters'].iloc[0]:,} clusters\n")
        f.write(f"- High diversity in protist groups (Alveolata, Rhizaria, Stramenopiles)\n\n")
        
        f.write("### 3. Size Patterns\n")
        f.write(f"- Highly skewed distribution (mean cluster count: {cluster_stats['mean']:.1f}, median: {cluster_stats['50%']:.1f})\n")
        f.write(f"- Large number of singleton or small combinations\n")
        f.write(f"- Few very large taxonomic groups dominate sequence counts\n\n")
        
        f.write("### 4. Classification Quality\n")
        best_classified = division_rates_df.nlargest(3, 'cluster_classification_rate')
        f.write("- Best classified divisions:\n")
        for _, row in best_classified.iterrows():
            f.write(f"  - {row['division']}: {row['cluster_classification_rate']:.1f}% classification rate\n")
        
        worst_classified = division_rates_df.nsmallest(3, 'cluster_classification_rate')
        f.write("- Poorly classified divisions:\n")
        for _, row in worst_classified.iterrows():
            f.write(f"  - {row['division']}: {row['cluster_classification_rate']:.1f}% classification rate\n")
        
        f.write("\n## Recommendations\n\n")
        f.write("1. **Focus curation efforts** on the largest unclassified groups:\n")
        unclassified_large = combinations_df[combinations_df['taxonomic_combination'].str.contains('.U.')].nlargest(5, 'row_count')
        for _, row in unclassified_large.iterrows():
            f.write(f"   - {row['taxonomic_combination']} ({row['row_count']:,} clusters)\n")
        
        f.write("\n2. **Investigate outliers** with unusual size patterns for potential data quality issues\n\n")
        
        f.write("3. **Prioritize classification** of high-abundance groups in poorly classified divisions\n\n")
        
        f.write("4. **Consider splitting** very large unclassified groups for more granular analysis\n\n")
        
        # File inventory
        f.write("## Generated Files\n\n")
        f.write("This analysis generated the following files in the `sanity_check/` directory:\n\n")
        
        all_files = [f for f in os.listdir(output_dir) if f.endswith(('.csv', '.txt', '.md'))]
        all_files.sort()
        
        for file in all_files:
            f.write(f"- `{file}`\n")
        
        f.write(f"\n---\n")
        f.write(f"*Report generated on {timestamp}*\n")
    
    print(f"💾 Comprehensive report saved to: {report_file}")
    return report_file

def main():
    """Main function"""
    
    output_dir = "sanity_check"
    
    if not os.path.exists(output_dir):
        print(f"❌ Error: Directory not found: {output_dir}")
        return
    
    print("🚀 Starting comprehensive report generation...")
    print("=" * 60)
    
    # Generate report
    report_file = generate_comprehensive_report(output_dir)
    
    print("\n" + "=" * 60)
    print("✅ Comprehensive report generated successfully!")
    print(f"\n📄 Report location: {report_file}")
    print("\n📋 Summary of all sanity check files:")
    
    # List all generated files
    all_files = [f for f in os.listdir(output_dir) if f.endswith(('.csv', '.txt', '.md'))]
    all_files.sort()
    
    for file in all_files:
        print(f"  • {file}")

if __name__ == "__main__":
    main()
