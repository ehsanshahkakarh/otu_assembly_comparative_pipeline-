#!/usr/bin/env python3
"""
Unclassified lineage parser for 18S EukCensus data
Analyzes .U. (unclassified) entries at division, family, and genus levels
Shows full lineage information to understand classification differences
"""

import pandas as pd
import os
from collections import defaultdict
from tqdm import tqdm

def analyze_unclassified_entries(input_file, output_dir):
    """
    Analyze unclassified entries and their lineage information
    
    Args:
        input_file (str): Path to input TSV file
        output_dir (str): Directory to save output CSV files
    """
    
    print(f"📊 Reading TSV file: {input_file}")
    
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"✅ Loaded {len(df):,} records")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each taxonomic level
    taxonomic_levels = ['division', 'family', 'genus']
    
    for level in taxonomic_levels:
        print(f"\n🔍 Processing unclassified {level} entries...")
        
        # Filter for unclassified entries (containing .U.)
        unclassified_mask = df[level].str.contains('.U.', na=False)
        unclassified_df = df[unclassified_mask].copy()
        
        if len(unclassified_df) == 0:
            print(f"  No unclassified entries found for {level}")
            continue
            
        print(f"  Found {len(unclassified_df):,} unclassified {level} entries")
        
        # Group by taxonomic name and collect lineage information
        lineage_groups = defaultdict(list)
        
        for idx, row in tqdm(unclassified_df.iterrows(), total=len(unclassified_df), 
                           desc=f"Analyzing {level}"):
            taxonomic_name = row[level]
            if pd.notna(taxonomic_name):
                # Create a lineage signature from division, family, genus
                lineage_sig = f"{row['division']}|{row['family']}|{row['genus']}"
                
                # Get member info
                members = row['members'] if pd.notna(row['members']) else '[]'
                size = row['size'] if pd.notna(row['size']) else 0
                
                lineage_groups[taxonomic_name].append({
                    'lineage_signature': lineage_sig,
                    'division': row['division'],
                    'family': row['family'], 
                    'genus': row['genus'],
                    'members': members,
                    'size': size,
                    'centroid': row['centroid']
                })
        
        # Create summary DataFrame
        summary_data = []
        
        for taxonomic_name, entries in lineage_groups.items():
            # Group by unique lineage signatures
            lineage_counts = defaultdict(lambda: {'count': 0, 'total_size': 0, 'examples': []})
            
            for entry in entries:
                sig = entry['lineage_signature']
                lineage_counts[sig]['count'] += 1
                lineage_counts[sig]['total_size'] += entry['size']
                if len(lineage_counts[sig]['examples']) < 3:  # Keep up to 3 examples
                    lineage_counts[sig]['examples'].append(entry['centroid'])
            
            # Add summary for each unique lineage pattern
            for lineage_sig, stats in lineage_counts.items():
                div, fam, gen = lineage_sig.split('|')
                summary_data.append({
                    'unclassified_name': taxonomic_name,
                    'count': stats['count'],
                    'sequence_count': stats['total_size'],
                    'division': div,
                    'family': fam,
                    'genus': gen,
                    'lineage_signature': lineage_sig,
                    'example_centroids': '; '.join(stats['examples'][:3])
                })
        
        # Convert to DataFrame and sort alphabetically by division, then by count
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values(['division', 'count'], ascending=[True, False])
        
        # Add metadata row at the top
        total_unclassified = len(unclassified_df)
        total_all = len(df)
        metadata_row = pd.DataFrame([{
            'unclassified_name': f'# Unclassified {level} entries: {total_unclassified:,} of {total_all:,} total',
            'count': '',
            'sequence_count': '',
            'division': '',
            'family': '',
            'genus': '',
            'lineage_signature': '',
            'example_centroids': ''
        }])
        
        # Combine metadata and data
        final_df = pd.concat([metadata_row, summary_df], ignore_index=True)
        
        # Save to CSV
        output_file = os.path.join(output_dir, f'unclassified_{level}_lineages.csv')
        final_df.to_csv(output_file, index=False)
        
        print(f"💾 Saved analysis to: {output_file}")
        
        # Show summary statistics
        unique_names = summary_df['unclassified_name'].nunique()
        unique_lineages = summary_df['lineage_signature'].nunique()
        total_entries = summary_df['count'].sum()
        total_sequences = summary_df['sequence_count'].sum()
        
        print(f"📈 Summary for {level}:")
        print(f"  • {unique_names} unique unclassified {level} names")
        print(f"  • {unique_lineages} unique lineage patterns")
        print(f"  • {total_entries:,} total entries")
        print(f"  • {total_sequences:,} total sequences")
        
        # Show top unclassified entries
        print(f"📋 Top unclassified {level} entries:")
        top_entries = summary_df.groupby('unclassified_name').agg({
            'count': 'sum',
            'sequence_count': 'sum'
        }).sort_values('count', ascending=False).head(5)
        
        for name, row in top_entries.iterrows():
            print(f"  • {name}: {row['count']:,} entries, {row['sequence_count']:,} sequences")

def create_unclassified_summary(input_file, output_dir):
    """
    Create a summary of unclassified entries across all taxonomic levels

    Args:
        input_file (str): Path to input TSV file
        output_dir (str): Directory to save output CSV files
    """

    print(f"\n🔍 Creating unclassified summary...")

    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Calculate statistics for each taxonomic level
    taxonomic_levels = ['division', 'family', 'genus']
    summary_data = []

    total_rows = len(df)
    total_sequences = df['size'].sum()

    for level in taxonomic_levels:
        # Count unclassified entries
        unclassified_mask = df[level].str.contains('.U.', na=False)
        unclassified_count = unclassified_mask.sum()
        unclassified_sequences = df[unclassified_mask]['size'].sum()

        # Calculate percentages
        unclassified_pct = (unclassified_count / total_rows) * 100
        unclassified_seq_pct = (unclassified_sequences / total_sequences) * 100

        # Count classified entries
        classified_count = total_rows - unclassified_count
        classified_sequences = total_sequences - unclassified_sequences
        classified_pct = (classified_count / total_rows) * 100
        classified_seq_pct = (classified_sequences / total_sequences) * 100

        summary_data.extend([
            {
                'taxonomic_level': level,
                'classification_status': 'classified',
                'entry_count': classified_count,
                'entry_percentage': round(classified_pct, 1),
                'sequence_count': classified_sequences,
                'sequence_percentage': round(classified_seq_pct, 1)
            },
            {
                'taxonomic_level': level,
                'classification_status': 'unclassified (.U.)',
                'entry_count': unclassified_count,
                'entry_percentage': round(unclassified_pct, 1),
                'sequence_count': unclassified_sequences,
                'sequence_percentage': round(unclassified_seq_pct, 1)
            }
        ])

    # Create DataFrame
    summary_df = pd.DataFrame(summary_data)

    # Add metadata row at the top
    metadata_row = pd.DataFrame([{
        'taxonomic_level': f'# Total metadata: {total_rows:,} entries, {total_sequences:,} sequences',
        'classification_status': '',
        'entry_count': '',
        'entry_percentage': '',
        'sequence_count': '',
        'sequence_percentage': ''
    }])

    # Combine metadata and data
    final_df = pd.concat([metadata_row, summary_df], ignore_index=True)

    # Save to CSV
    output_file = os.path.join(output_dir, 'unclassified_summary_stats.csv')
    final_df.to_csv(output_file, index=False)

    print(f"💾 Saved summary statistics to: {output_file}")

    # Print summary to terminal
    print(f"\n📊 Unclassified Data Summary:")
    print(f"{'Level':<10} {'Status':<15} {'Entries':<10} {'Entry %':<8} {'Sequences':<12} {'Seq %':<6}")
    print("-" * 70)

    for level in taxonomic_levels:
        level_data = summary_df[summary_df['taxonomic_level'] == level]
        for _, row in level_data.iterrows():
            status = row['classification_status']
            entries = f"{row['entry_count']:,}"
            entry_pct = f"{row['entry_percentage']}%"
            sequences = f"{row['sequence_count']:,}"
            seq_pct = f"{row['sequence_percentage']}%"

            print(f"{level:<10} {status:<15} {entries:<10} {entry_pct:<8} {sequences:<12} {seq_pct:<6}")

def main():
    """Main function"""

    # File paths - script is in metadata folder
    input_file = "eukcensus_18S.clusters.97.tsv"
    output_dir = "sanity_check"

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file not found: {input_file}")
        return

    print("🚀 Starting unclassified lineage analysis...")
    print("=" * 60)

    # Analyze unclassified entries
    analyze_unclassified_entries(input_file, output_dir)

    # Create summary statistics
    create_unclassified_summary(input_file, output_dir)

    print("\n" + "=" * 60)
    print("✅ Analysis completed successfully!")

    # Show output files
    print(f"\n📁 Output files created in '{output_dir}':")
    for level in ['division', 'family', 'genus']:
        output_file = os.path.join(output_dir, f'unclassified_{level}_lineages.csv')
        if os.path.exists(output_file):
            print(f"  • {output_file}")

    summary_file = os.path.join(output_dir, 'unclassified_summary_stats.csv')
    if os.path.exists(summary_file):
        print(f"  • {summary_file}")

if __name__ == "__main__":
    main()
