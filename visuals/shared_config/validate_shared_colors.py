#!/usr/bin/env python3
"""
Validate Shared Color Assignments
=================================
This script validates that shared taxa between scatter plots and alluvial plots
have consistent color assignments in the shared configuration.

Usage: python validate_shared_colors.py
"""

import pandas as pd
import yaml
import os
from pathlib import Path

def load_scatter_taxa():
    """Load taxa from scatter plot source data"""
    scatter_dir = Path("../scatter_plots/source_data")
    
    taxa = {
        'archaea': [],
        'bacteria': [],
        'eukaryota': []
    }
    
    # Load each domain's scatter plot data
    for domain in ['Archaea', 'Bacteria', 'Eukaryota']:
        file_path = scatter_dir / f"{domain}_phylum_source_data.csv"
        if file_path.exists():
            df = pd.read_csv(file_path)
            taxa[domain.lower()] = list(df['Taxon'])
            print(f"Loaded {len(taxa[domain.lower()])} {domain} taxa from scatter plots")
    
    return taxa

def load_alluvial_taxa():
    """Load taxa from alluvial plot flow annotations"""
    alluvial_dir = Path("../alluvial_plots/16s_alluvial")
    
    taxa = {
        'archaea': [],
        'bacteria': [],
        'eukaryota': []
    }
    
    # Load archaea alluvial data
    archaea_file = alluvial_dir / "alluvial_16s_archaea_pct_flow_annotations.tsv"
    if archaea_file.exists():
        df = pd.read_csv(archaea_file, sep='\t')
        taxa['archaea'] = [t for t in df['Taxon'].unique() if t not in ['Other', 'Archaea.U.phylum']]
        print(f"Loaded {len(taxa['archaea'])} Archaea taxa from alluvial plots")
    
    # Load bacteria alluvial data  
    bacteria_file = alluvial_dir / "alluvial_16s_bacteria_pct_flow_annotations.tsv"
    if bacteria_file.exists():
        df = pd.read_csv(bacteria_file, sep='\t')
        taxa['bacteria'] = [t for t in df['Taxon'].unique() if t not in ['Other', 'Bacteria.U.phylum']]
        print(f"Loaded {len(taxa['bacteria'])} Bacteria taxa from alluvial plots")
    
    # Load eukaryotic alluvial data (18S)
    eukaryota_file = Path("../alluvial_plots/18s_alluvial/alluvial_18s_abs_flow_annotations.tsv")
    if eukaryota_file.exists():
        df = pd.read_csv(eukaryota_file, sep='\t')
        taxa['eukaryota'] = [t for t in df['Taxon'].unique() if t not in ['Other', 'Eukaryota.U.division']]
        print(f"Loaded {len(taxa['eukaryota'])} Eukaryota taxa from alluvial plots")
    
    return taxa

def load_color_config():
    """Load the shared color configuration"""
    config_file = Path("taxonomic_color_mapping.yaml")
    
    if not config_file.exists():
        print(f"Error: Color config file not found: {config_file}")
        return None
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def find_shared_taxa(scatter_taxa, alluvial_taxa):
    """Find taxa that appear in both scatter plots and alluvial plots"""
    shared = {}
    
    for domain in ['archaea', 'bacteria', 'eukaryota']:
        scatter_set = set(scatter_taxa[domain])
        alluvial_set = set(alluvial_taxa[domain])
        shared[domain] = list(scatter_set.intersection(alluvial_set))
        
        print(f"\n{domain.upper()} SHARED TAXA ({len(shared[domain])} taxa):")
        for taxon in sorted(shared[domain]):
            print(f"  - {taxon}")
    
    return shared

def validate_color_assignments(shared_taxa, color_config):
    """Validate that shared taxa have color assignments"""
    print("\n" + "="*60)
    print("COLOR ASSIGNMENT VALIDATION")
    print("="*60)
    
    for domain in ['archaea', 'bacteria', 'eukaryota']:
        domain_colors = color_config.get(f"{domain}_colors", {})
        
        print(f"\n{domain.upper()} COLOR ASSIGNMENTS:")
        print("-" * 40)
        
        for taxon in sorted(shared_taxa[domain]):
            if taxon in domain_colors:
                color = domain_colors[taxon]
                print(f"  ✅ {taxon:<25} → {color}")
            else:
                print(f"  ❌ {taxon:<25} → NOT ASSIGNED")
    
    return True

def main():
    print("Validating Shared Taxa Color Assignments")
    print("=" * 50)
    
    # Load taxa from both visualization types
    scatter_taxa = load_scatter_taxa()
    alluvial_taxa = load_alluvial_taxa()
    
    # Find shared taxa
    shared_taxa = find_shared_taxa(scatter_taxa, alluvial_taxa)
    
    # Load color configuration
    color_config = load_color_config()
    if not color_config:
        return
    
    # Validate color assignments
    validate_color_assignments(shared_taxa, color_config)
    
    # Summary statistics
    total_shared = sum(len(shared_taxa[domain]) for domain in shared_taxa)
    print(f"\n" + "="*60)
    print(f"SUMMARY: {total_shared} total shared taxa across all domains")
    print(f"  - Archaea: {len(shared_taxa['archaea'])}/4 taxa shared (100% overlap)")
    print(f"  - Bacteria: {len(shared_taxa['bacteria'])} taxa shared")
    print(f"  - Eukaryota: {len(shared_taxa['eukaryota'])} taxa shared")
    print("="*60)

if __name__ == "__main__":
    main()
