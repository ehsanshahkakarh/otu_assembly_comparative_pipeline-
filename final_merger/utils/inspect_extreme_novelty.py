#!/usr/bin/env python3
"""
Inspect taxa with extreme novelty factors (>1000).

This script helps understand which taxa have very high novelty factors
and whether this is expected (underrepresented in genome databases).

Usage:
    python inspect_extreme_novelty.py --gene 18S
    python inspect_extreme_novelty.py --gene 16S
    python inspect_extreme_novelty.py --gene both --threshold 100
"""

import pandas as pd
import argparse
from pathlib import Path


def inspect_extreme_novelty(gene: str, threshold: float = 1000):
    """Find and display taxa with extreme novelty factors."""
    
    base_dir = Path(__file__).resolve().parent.parent
    output_dir = base_dir / 'outputs'
    
    levels = ['division', 'family', 'genus']
    
    print(f"\n{'='*80}")
    print(f"Taxa with Novelty Factor > {threshold} ({gene})")
    print(f"{'='*80}\n")
    
    all_extreme = []
    
    for level in levels:
        file_path = output_dir / f'{gene.lower()}_ncbi_merged_{level}.csv'
        
        if not file_path.exists():
            print(f"⚠️  File not found: {file_path}")
            continue
        
        df = pd.read_csv(file_path)
        
        # Filter for extreme novelty factors (excluding inf)
        extreme = df[
            (df['novelty_factor'] > threshold) & 
            (df['novelty_factor'] != float('inf'))
        ].copy()
        
        if len(extreme) > 0:
            extreme['level'] = level
            all_extreme.append(extreme)
            
            print(f"\n{level.upper()} Level: {len(extreme)} taxa")
            print("-" * 80)
            
            # Sort by novelty factor descending
            extreme_sorted = extreme.sort_values('novelty_factor', ascending=False)
            
            for idx, row in extreme_sorted.iterrows():
                taxon = row[df.columns[0]]  # First column is taxon name
                nf = row['novelty_factor']
                census = row['census_otu_count']
                ncbi = row['ncbi_species_count']
                domain = row.get('domain', 'Unknown')
                
                print(f"\n  {taxon}")
                print(f"    Domain: {domain}")
                print(f"    Novelty Factor: {nf:.1f}")
                print(f"    Census OTUs: {census:,}")
                print(f"    NCBI Species: {ncbi}")
                print(f"    → {census:,} census OTUs but only {ncbi} NCBI genome(s)")
    
    # Summary
    if all_extreme:
        combined = pd.concat(all_extreme, ignore_index=True)
        
        print(f"\n\n{'='*80}")
        print(f"SUMMARY")
        print(f"{'='*80}\n")
        print(f"Total taxa with NF > {threshold}: {len(combined)}")
        print(f"\nBy level:")
        for level in levels:
            count = len(combined[combined['level'] == level])
            if count > 0:
                print(f"  {level}: {count}")
        
        print(f"\nTop 5 highest novelty factors:")
        top5 = combined.nlargest(5, 'novelty_factor')
        for idx, row in top5.iterrows():
            taxon = row[combined.columns[0]]
            nf = row['novelty_factor']
            level = row['level']
            print(f"  {taxon} ({level}): {nf:.1f}")
        
        print(f"\n{'='*80}")
        print("INTERPRETATION")
        print(f"{'='*80}\n")
        print("High novelty factors indicate taxa that are:")
        print("  • Well-represented in environmental sequencing (high census OTUs)")
        print("  • Poorly represented in genome databases (low NCBI genomes)")
        print("\nThis is EXPECTED for:")
        print("  • Unculturable organisms")
        print("  • Understudied lineages")
        print("  • Environmental microbes")
        print("  • Novel or rare taxa")
        print("\nThese are HIGH-PRIORITY targets for genome sequencing!")
        print(f"{'='*80}\n")
    else:
        print(f"\n✅ No taxa found with novelty factor > {threshold}")


def main():
    parser = argparse.ArgumentParser(
        description='Inspect taxa with extreme novelty factors',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--gene', required=True, choices=['18S', '16S', 'both'],
                       help='Gene type to inspect')
    parser.add_argument('--threshold', type=float, default=1000,
                       help='Novelty factor threshold (default: 1000)')
    
    args = parser.parse_args()
    
    if args.gene == 'both':
        genes = ['18S', '16S']
    else:
        genes = [args.gene]
    
    for gene in genes:
        inspect_extreme_novelty(gene, args.threshold)


if __name__ == "__main__":
    main()

