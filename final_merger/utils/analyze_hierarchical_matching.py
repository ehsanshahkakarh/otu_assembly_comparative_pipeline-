#!/usr/bin/env python3
"""
Hierarchical Matching Analysis
===============================

Analyzes the results of the new hierarchical matching strategy.
"""

import pandas as pd
from pathlib import Path

def analyze_matching(census_type: str):
    """Analyze hierarchical matching results for a census type."""
    
    print("=" * 80)
    print(f"{census_type} HIERARCHICAL MATCHING ANALYSIS")
    print("=" * 80)
    
    levels = ['division', 'family', 'genus']
    
    for level in levels:
        output_file = Path(f"outputs/{census_type.lower()}_ncbi_merged_{level}.csv")
        
        if not output_file.exists():
            print(f"\n⚠️  File not found: {output_file}")
            continue
        
        df = pd.read_csv(output_file)
        matched = df[df['match_status'] == 'matched']
        
        print(f"\n{'=' * 80}")
        print(f"{level.upper()} LEVEL")
        print(f"{'=' * 80}")
        
        print(f"\nTotal taxa: {len(df)}")
        print(f"Matched taxa: {len(matched)} ({len(matched)/len(df)*100:.1f}%)")
        print(f"Census-only taxa: {len(df) - len(matched)}")
        
        # Count primary match types
        direct_name = len(matched[matched['direct_name_matches'] > 0])
        direct_taxid = len(matched[(matched['direct_name_matches'] == 0) & (matched['direct_taxid_matches'] > 0)])
        hierarchical_taxid = len(matched[(matched['direct_name_matches'] == 0) & (matched['direct_taxid_matches'] == 0) & (matched['hierarchical_taxid_matches'] > 0)])
        lineage_name = len(matched[(matched['direct_name_matches'] == 0) & (matched['direct_taxid_matches'] == 0) & (matched['hierarchical_taxid_matches'] == 0) & (matched['lineage_name_matches'] > 0)])
        
        print(f"\nPrimary match method:")
        print(f"  1. Direct name matches:        {direct_name:5d} taxa ({direct_name/len(matched)*100:5.1f}%)")
        print(f"  2. Direct taxid matches:       {direct_taxid:5d} taxa ({direct_taxid/len(matched)*100:5.1f}%)")
        print(f"  3. Hierarchical taxid matches: {hierarchical_taxid:5d} taxa ({hierarchical_taxid/len(matched)*100:5.1f}%)")
        print(f"  4. Lineage name matches:       {lineage_name:5d} taxa ({lineage_name/len(matched)*100:5.1f}%)")
        
        print(f"\nTotal match counts (can overlap):")
        print(f"  Direct name matches:        {matched['direct_name_matches'].sum():,}")
        print(f"  Direct taxid matches:       {matched['direct_taxid_matches'].sum():,}")
        print(f"  Hierarchical taxid matches: {matched['hierarchical_taxid_matches'].sum():,}")
        print(f"  Lineage name matches:       {matched['lineage_name_matches'].sum():,}")
        
        # Show examples of hierarchical taxid matches
        hierarchical_examples = matched[(matched['direct_name_matches'] == 0) & (matched['direct_taxid_matches'] == 0) & (matched['hierarchical_taxid_matches'] > 0)]
        
        if len(hierarchical_examples) > 0:
            print(f"\n🔍 Examples of hierarchical taxid matches:")
            for idx, (_, row) in enumerate(hierarchical_examples.head(5).iterrows(), 1):
                taxon_col = level
                print(f"  {idx}. {row[taxon_col]:30s} (taxid: {row['census_taxid']}) → {row['hierarchical_taxid_matches']} NCBI matches")
        
        # Show examples of lineage name only matches
        lineage_only_examples = matched[(matched['direct_name_matches'] == 0) & (matched['direct_taxid_matches'] == 0) & (matched['hierarchical_taxid_matches'] == 0) & (matched['lineage_name_matches'] > 0)]
        
        if len(lineage_only_examples) > 0:
            print(f"\n📝 Examples of lineage name only matches:")
            for idx, (_, row) in enumerate(lineage_only_examples.head(5).iterrows(), 1):
                taxon_col = level
                print(f"  {idx}. {row[taxon_col]:30s} (taxid: {row['census_taxid']}) → {row['lineage_name_matches']} NCBI matches")


def main():
    """Main analysis function."""
    
    # Analyze 18S
    analyze_matching('18S')
    
    print("\n\n")
    
    # Analyze 16S
    analyze_matching('16S')
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    main()

