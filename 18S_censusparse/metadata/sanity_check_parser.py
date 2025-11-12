#!/usr/bin/env python3
"""
Sanity check parser for 18S EukCensus data
Parses division, family, and genus columns from TSV file and creates separate CSV files with counts
"""

import pandas as pd
import os
from collections import Counter
from tqdm import tqdm

def parse_taxonomic_levels(input_file, output_dir):
    """
    Parse taxonomic levels from TSV file and create separate CSV files with counts

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
        print(f"\n🔍 Processing {level} level...")

        # Count occurrences and sum sequence counts for each taxonomic name
        level_counts = Counter()
        level_sequence_counts = Counter()

        for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"Counting {level}"):
            taxonomic_name = row[level]
            sequence_size = row['size'] if pd.notna(row['size']) else 0

            if pd.notna(taxonomic_name) and taxonomic_name.strip():  # Skip empty/NaN values
                clean_name = taxonomic_name.strip()
                level_counts[clean_name] += 1
                level_sequence_counts[clean_name] += sequence_size

        # Convert to DataFrame and sort by count (descending)
        count_df = pd.DataFrame([
            {
                'name': name,
                'count': level_counts[name],
                'sequence_count': level_sequence_counts[name]
            }
            for name in sorted(level_counts.keys(), key=lambda x: level_counts[x], reverse=True)
        ])

        # Add metadata row at the top
        metadata_row = pd.DataFrame([{
            'name': f'# Total rows in original metadata: {len(df):,}',
            'count': '',
            'sequence_count': ''
        }])

        # Combine metadata and data
        final_df = pd.concat([metadata_row, count_df], ignore_index=True)

        # Save to CSV
        output_file = os.path.join(output_dir, f'eukcensus_18S_{level}_counts.csv')
        final_df.to_csv(output_file, index=False)

        print(f"💾 Saved {len(count_df)} unique {level} entries to: {output_file}")

        # Show top 10 for sanity check
        print(f"📈 Top 10 {level} by count:")
        for i, name in enumerate(sorted(level_counts.keys(), key=lambda x: level_counts[x], reverse=True)[:10], 1):
            count = level_counts[name]
            seq_count = level_sequence_counts[name]
            print(f"  {i:2d}. {name:<30} ({count:,} occurrences, {seq_count:,} sequences)")

def main():
    """Main function"""
    
    # File paths - script is in metadata folder
    input_file = "eukcensus_18S.clusters.97.tsv"
    output_dir = "sanity_check"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file not found: {input_file}")
        return
    
    print("🚀 Starting 18S EukCensus taxonomic parsing...")
    print("=" * 60)
    
    # Parse taxonomic levels
    parse_taxonomic_levels(input_file, output_dir)
    
    print("\n" + "=" * 60)
    print("✅ Parsing completed successfully!")
    
    # Show output files
    print(f"\n📁 Output files created in '{output_dir}':")
    for level in ['division', 'family', 'genus']:
        output_file = os.path.join(output_dir, f'eukcensus_18S_{level}_counts.csv')
        if os.path.exists(output_file):
            print(f"  • {output_file}")

if __name__ == "__main__":
    main()
