#!/usr/bin/env python3
"""
Taxonomic Tree Builder
======================

Build hierarchical taxonomic trees for each domain (Bacteria, Archaea, Eukaryota, Viruses).
Generates Newick format trees for iTOL visualization.
"""

import pandas as pd
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict


def build_tree_from_lineages(df, domain_name):
    """
    Build a hierarchical tree structure from lineage strings.
    
    Args:
        df: DataFrame with lineage, lineage_ranks, lineage_taxids columns
        domain_name: Name of the domain (for root node)
    
    Returns:
        Newick format tree string
    """
    print(f"\n🌳 Building tree for {domain_name}...")
    print(f"   Processing {len(df):,} species...")
    
    # Build tree structure as nested dict
    tree = {}
    
    for idx, row in df.iterrows():
        lineage = str(row['lineage'])
        lineage_ranks = str(row['lineage_ranks'])
        lineage_taxids = str(row['lineage_taxids'])
        
        if pd.isna(lineage) or lineage == 'nan':
            continue
        
        taxa = [t.strip() for t in lineage.split(';')]
        ranks = [r.strip() for r in lineage_ranks.split(';')]
        taxids = [t.strip() for t in lineage_taxids.split(';')]
        
        # Build path through tree
        current = tree
        for i, (taxon, rank, taxid) in enumerate(zip(taxa, ranks, taxids)):
            # Create unique key with rank and taxid
            key = f"{taxon}_{rank}_{taxid}"
            
            if key not in current:
                current[key] = {
                    'name': taxon,
                    'rank': rank,
                    'taxid': taxid,
                    'children': {},
                    'species_count': 0,
                    'genome_count': 0
                }
            
            # Add counts at leaf level (species)
            if i == len(taxa) - 1:
                current[key]['species_count'] += 1
                current[key]['genome_count'] += row['total_genome_count']
            
            current = current[key]['children']
    
    print(f"   ✅ Tree structure built")
    
    # Convert tree to Newick format
    def to_newick(node_dict, depth=0):
        if not node_dict:
            return ""
        
        children_newick = []
        for key, node in node_dict.items():
            name = node['name'].replace(' ', '_').replace('(', '').replace(')', '').replace(':', '_').replace(';', '_').replace(',', '_')
            rank = node['rank']
            taxid = node['taxid']
            
            if node['children']:
                # Internal node
                child_str = to_newick(node['children'], depth + 1)
                children_newick.append(f"({child_str}){name}_{rank}_{taxid}")
            else:
                # Leaf node (species)
                children_newick.append(f"{name}_{rank}_{taxid}")
        
        return ','.join(children_newick)
    
    newick = f"({to_newick(tree)}){domain_name};"
    
    print(f"   ✅ Newick format generated")
    
    return newick, tree


def generate_tree_metadata(df, domain_name, output_dir, timestamp):
    """
    Generate iTOL metadata files for tree annotation.
    
    Args:
        df: DataFrame with species data
        domain_name: Name of the domain
        output_dir: Output directory
        timestamp: Timestamp string
    """
    print(f"\n📊 Generating iTOL metadata for {domain_name}...")
    
    # Create metadata for genome counts
    metadata_lines = []
    metadata_lines.append("DATASET_SIMPLEBAR")
    metadata_lines.append("SEPARATOR TAB")
    metadata_lines.append(f"DATASET_LABEL\t{domain_name}_genome_counts")
    metadata_lines.append("COLOR\t#ff0000")
    metadata_lines.append("DATA")
    
    for idx, row in df.iterrows():
        lineage = str(row['lineage'])
        if pd.isna(lineage) or lineage == 'nan':
            continue
        
        taxa = lineage.split(';')
        taxids = str(row['lineage_taxids']).split(';')
        ranks = str(row['lineage_ranks']).split(';')
        
        if len(taxa) > 0 and len(taxids) > 0 and len(ranks) > 0:
            name = taxa[-1].strip().replace(' ', '_').replace('(', '').replace(')', '').replace(':', '_').replace(';', '_').replace(',', '_')
            rank = ranks[-1].strip()
            taxid = taxids[-1].strip()
            label = f"{name}_{rank}_{taxid}"
            count = row['total_genome_count']
            metadata_lines.append(f"{label}\t{count}")
    
    # Save metadata
    metadata_file = output_dir / f'{domain_name.lower()}_tree_metadata_{timestamp}.txt'
    with open(metadata_file, 'w') as f:
        f.write('\n'.join(metadata_lines))
    
    print(f"   ✅ Metadata saved: {metadata_file.name}")


def build_all_trees(species_with_domain_file, output_dir):
    """
    Build trees for all domains.
    
    Args:
        species_with_domain_file: Path to species file with domain column
        output_dir: Directory to save output files
    """
    print("="*70)
    print("TAXONOMIC TREE BUILDER FOR iTOL")
    print("="*70)
    
    # Load data
    print(f"\n📁 Loading data from: {species_with_domain_file.name}")
    df = pd.read_csv(species_with_domain_file)
    print(f"   Total species: {len(df):,}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Build tree for each domain
    domains = ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']
    
    for domain in domains:
        domain_df = df[df['domain'] == domain].copy()
        
        if len(domain_df) == 0:
            print(f"\n⚠️  No species found for {domain}, skipping...")
            continue
        
        print(f"\n{'='*70}")
        print(f"{domain.upper()} TREE")
        print(f"{'='*70}")
        print(f"   Species: {len(domain_df):,}")
        print(f"   Genomes: {domain_df['total_genome_count'].sum():,}")
        
        # Build tree
        newick, tree_structure = build_tree_from_lineages(domain_df, domain)
        
        # Save Newick file
        newick_file = output_dir / f'{domain.lower()}_tree_{timestamp}.nwk'
        with open(newick_file, 'w') as f:
            f.write(newick)
        
        print(f"   💾 Newick tree saved: {newick_file.name}")
        print(f"   📏 Tree size: {len(newick):,} characters")
        
        # Generate metadata
        generate_tree_metadata(domain_df, domain, output_dir, timestamp)
    
    print("\n" + "="*70)
    print("✅ ALL TREES GENERATED!")
    print("="*70)
    print(f"\n📁 Output directory: {output_dir}")
    print(f"\nFiles generated for each domain:")
    print(f"   - <domain>_tree_<timestamp>.nwk       (Newick tree for iTOL)")
    print(f"   - <domain>_tree_metadata_<timestamp>.txt (iTOL metadata)")


def main():
    # Find the most recent species_with_domain file
    sanity_output = Path(__file__).parent / 'output'
    domain_files = list(sanity_output.glob('species_with_domain_*.csv'))
    
    if not domain_files:
        print("❌ No species_with_domain files found!")
        print("   Run domain_parser.py first to generate domain-annotated file.")
        sys.exit(1)
    
    # Use the most recent file
    input_file = max(domain_files, key=lambda p: p.stat().st_mtime)
    
    # Run tree building
    build_all_trees(input_file, sanity_output)
    
    print("\n🎉 Tree building complete!")
    print("\n📤 Upload the .nwk files to iTOL: https://itol.embl.de/upload.cgi")


if __name__ == '__main__':
    main()

