#!/usr/bin/env python3
"""
Sort Merged Files by Novelty Factor
====================================

Reorders merged census files by novelty factor (descending).
High novelty = priority targets for genome sequencing.

Usage:
    python sort_by_novelty.py <input_file.csv>
    
    # Or sort all files in 18S_output/
    python sort_by_novelty.py --all
"""

import sys
import pandas as pd
from pathlib import Path
import argparse


def sort_by_novelty(input_file: Path, output_file: Path = None) -> Path:
    """
    Sort merged file by novelty factor (descending).
    
    Args:
        input_file: Path to merged CSV file
        output_file: Optional output path (default: adds _sorted_by_novelty suffix)
    
    Returns:
        Path to sorted output file
    """
    print(f"Loading: {input_file.name}")
    
    # Load data
    df = pd.read_csv(input_file, dtype={'census_taxid': str})
    
    # Separate finite and infinite novelty values
    df_finite = df[df['novelty_factor'] != float('inf')].copy()
    df_inf = df[df['novelty_factor'] == float('inf')].copy()
    
    # Sort finite values by novelty (descending - highest first)
    df_finite_sorted = df_finite.sort_values('novelty_factor', ascending=False)
    
    # Concatenate: finite (sorted) + infinite
    df_sorted = pd.concat([df_finite_sorted, df_inf])
    
    # Determine output file
    if output_file is None:
        stem = input_file.stem
        output_file = input_file.parent / f"{stem}_sorted_by_novelty.csv"
    
    # Save
    df_sorted.to_csv(output_file, index=False)
    
    print(f"✅ Saved: {output_file.name}")
    print(f"   Total taxa: {len(df_sorted)}")
    print(f"   Matched taxa: {(df_sorted['match_status'] == 'matched').sum()}")
    
    # Show top 5
    print()
    print("   TOP 5 HIGHEST NOVELTY:")
    level_col = [col for col in df.columns if col in ['division', 'family', 'genus']][0]
    
    for i, (_, row) in enumerate(df_sorted.head(5).iterrows(), 1):
        novelty = row['novelty_factor']
        if novelty == float('inf'):
            novelty_str = "∞"
        else:
            novelty_str = f"{novelty:.1f}"
        
        print(f"   {i}. {row[level_col]:25s} | Novelty: {novelty_str:>8s} | OTUs: {row['census_otu_count']:,}")
    
    print()
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Sort merged files by novelty factor')
    parser.add_argument('input_file', nargs='?', help='Input CSV file to sort')
    parser.add_argument('--all', action='store_true', help='Sort all merged files in 18S_output/')
    parser.add_argument('-o', '--output', help='Output file path')
    
    args = parser.parse_args()
    
    if args.all:
        # Sort all merged files in 18S_output/
        output_dir = Path(__file__).parent / '18S_output'
        
        if not output_dir.exists():
            print(f"❌ Error: {output_dir} not found")
            sys.exit(1)
        
        # Find all merged files (not already sorted)
        merged_files = sorted(output_dir.glob('18s_*_merged_*.csv'))
        merged_files = [f for f in merged_files if 'sorted' not in f.name]
        
        if not merged_files:
            print("❌ No merged files found in 18S_output/")
            sys.exit(1)
        
        print("=" * 80)
        print("SORTING ALL MERGED FILES BY NOVELTY")
        print("=" * 80)
        print()
        
        for input_file in merged_files:
            sort_by_novelty(input_file)
        
        print("=" * 80)
        print(f"✅ Sorted {len(merged_files)} files")
        print("=" * 80)
        
    elif args.input_file:
        # Sort single file
        input_file = Path(args.input_file)
        
        if not input_file.exists():
            print(f"❌ Error: File not found: {input_file}")
            sys.exit(1)
        
        output_file = Path(args.output) if args.output else None
        
        print("=" * 80)
        print("SORTING BY NOVELTY FACTOR")
        print("=" * 80)
        print()
        
        sort_by_novelty(input_file, output_file)
        
        print("=" * 80)
    
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

