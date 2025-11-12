#!/usr/bin/env python3
"""
Novelty Factor Family Mapper for 18S NCBI Results

This script maps high novelty genera to their families using the taxonomic combinations file.
It identifies genera with the highest novelty factors and shows their family assignments.
"""

import pandas as pd
import sys
from pathlib import Path

def map_novelty_genera_to_families():
    """Map high novelty genera to their families"""
    
    # File paths
    script_dir = Path(__file__).resolve().parent
    genus_file = script_dir / "18s_ncbi_merged_clean_genus.csv"
    combinations_file = script_dir.parent.parent.parent / "18S_censusparse" / "metadata" / "sanity_check" / "taxonomic_combinations_detailed.csv"
    
    print("🔍 High Novelty Factor Genera and Their Families")
    print("=" * 60)
    
    # Load genus data
    print(f"📊 Loading genus data from: {genus_file}")
    genus_df = pd.read_csv(genus_file)
    
    # Load taxonomic combinations
    print(f"📋 Loading taxonomic combinations from: {combinations_file}")
    # Skip comment lines and load the CSV
    combinations_df = pd.read_csv(combinations_file, comment='#')
    
    # Filter out infinite novelty factors and sort by novelty factor (descending)
    finite_novelty_df = genus_df[genus_df['novelty_factor'] != float('inf')]
    top_novelty_genera = finite_novelty_df.nlargest(20, 'novelty_factor')

    # Also get the infinite novelty genera count
    infinite_count = len(genus_df[genus_df['novelty_factor'] == float('inf')])
    print(f"📈 Found {infinite_count} genera with infinite novelty (no NCBI species matches)")
    print(f"📊 Analyzing top 20 genera with finite novelty factors...")
    
    print(f"\n🎯 Top 20 Genera by Novelty Factor:")
    print(f"{'Rank':<4} {'Genus':<20} {'Novelty':<10} {'OTUs':<6} {'NCBI_Spp':<8} {'Family':<25} {'Division':<20}")
    print("-" * 100)
    
    results = []
    
    for rank, (_, genus_row) in enumerate(top_novelty_genera.iterrows(), 1):
        genus_name = genus_row['genus']
        novelty_factor = genus_row['novelty_factor']
        otu_count = genus_row['census_otu_count']
        ncbi_species_count = genus_row['ncbi_species_count']
        
        # Find family for this genus in combinations file
        genus_matches = combinations_df[combinations_df['genus'] == genus_name]
        
        if not genus_matches.empty:
            # Get the first match (should be the most abundant combination)
            family = genus_matches.iloc[0]['family']
            division = genus_matches.iloc[0]['division']
            row_count = genus_matches.iloc[0]['row_count']
        else:
            # Try to find partial matches (genus might have .U.genus suffix)
            partial_matches = combinations_df[combinations_df['genus'].str.contains(genus_name, na=False)]
            if not partial_matches.empty:
                family = partial_matches.iloc[0]['family']
                division = partial_matches.iloc[0]['division']
                row_count = partial_matches.iloc[0]['row_count']
            else:
                family = "NOT_FOUND"
                division = "NOT_FOUND"
                row_count = 0
        
        print(f"{rank:<4} {genus_name:<20} {novelty_factor:<10.1f} {otu_count:<6} {ncbi_species_count:<8} {family:<25} {division:<20}")

        results.append({
            'rank': rank,
            'genus': genus_name,
            'novelty_factor': novelty_factor,
            'census_otu_count': otu_count,
            'ncbi_species_count': ncbi_species_count,
            'family': family,
            'division': division,
            'combination_row_count': row_count
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    output_file = script_dir / "high_novelty_genera_families.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n💾 Results saved to: {output_file}")
    
    # Summary by family
    print(f"\n📊 Summary by Family (Top Novelty Genera):")
    family_summary = results_df[results_df['family'] != 'NOT_FOUND'].groupby('family').agg({
        'genus': 'count',
        'novelty_factor': 'mean',
        'census_otu_count': 'sum'
    }).round(1).sort_values('novelty_factor', ascending=False)
    
    print(f"{'Family':<25} {'Genera':<6} {'Avg Novelty':<12} {'Total OTUs':<10}")
    print("-" * 55)
    for family, row in family_summary.iterrows():
        print(f"{family:<25} {row['genus']:<6} {row['novelty_factor']:<12} {row['census_otu_count']:<10}")

if __name__ == "__main__":
    map_novelty_genera_to_families()
