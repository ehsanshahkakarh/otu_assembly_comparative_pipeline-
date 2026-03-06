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
        
        # Count occurrences of each taxonomic name
        level_counts = Counter()
        
        for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"Counting {level}"):
            taxonomic_name = row[level]
            if pd.notna(taxonomic_name) and taxonomic_name.strip():  # Skip empty/NaN values
                level_counts[taxonomic_name.strip()] += 1
        
        # Convert to DataFrame and sort by count (descending)
        count_df = pd.DataFrame([
            {'name': name, 'count': count} 
            for name, count in level_counts.most_common()
        ])
        
        # Save to CSV
        output_file = os.path.join(output_dir, f'eukcensus_18S_{level}_counts.csv')
        count_df.to_csv(output_file, index=False)
        
        print(f"💾 Saved {len(count_df)} unique {level} entries to: {output_file}")
        
        # Show top 10 for sanity check
        print(f"📈 Top 10 {level} by count:")
        for i, (name, count) in enumerate(level_counts.most_common(10), 1):
            print(f"  {i:2d}. {name:<30} ({count:,} occurrences)")

def main():
    """Main function"""

    # File paths - assuming script is in metadata folder
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
