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
    Collapses to genus level (stops at genus, doesn't include species).

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
    genus_stats = {}  # Track species count and genome count per genus

    for idx, row in df.iterrows():
        lineage = str(row['lineage'])
        lineage_ranks = str(row['lineage_ranks'])
        lineage_taxids = str(row['lineage_taxids'])

        if pd.isna(lineage) or lineage == 'nan':
            continue

        taxa = [t.strip() for t in lineage.split(';')]
        ranks = [r.strip() for r in lineage_ranks.split(';')]
        taxids = [t.strip() for t in lineage_taxids.split(';')]

        # Build path through tree, but STOP at genus level
        current = tree
        genus_key = None

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

            # Add counts to this node
            current[key]['species_count'] += 1
            current[key]['genome_count'] += row['total_genome_count']

            # STOP at genus level - don't add species as children
            if rank == 'genus':
                genus_key = key
                # Track genus statistics
                if genus_key not in genus_stats:
                    genus_stats[genus_key] = {
                        'name': taxon,
                        'taxid': taxid,
                        'species_count': 0,
                        'genome_count': 0
                    }
                genus_stats[genus_key]['species_count'] += 1
                genus_stats[genus_key]['genome_count'] += row['total_genome_count']
                break

            current = current[key]['children']

    print(f"   ✅ Tree structure built (collapsed to genus level)")
    print(f"   📊 Found {len(genus_stats):,} unique genera")

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
                # Leaf node (genus)
                children_newick.append(f"{name}_{rank}_{taxid}")

        return ','.join(children_newick)

    newick = f"({to_newick(tree)}){domain_name};"

    print(f"   ✅ Newick format generated")

    return newick, genus_stats


def generate_tree_metadata(genus_stats, domain_name, output_dir, timestamp):
    """
    Generate iTOL metadata files for tree annotation.
    Creates bar graph showing species count per genus.

    Args:
        genus_stats: Dictionary with genus statistics
        domain_name: Name of the domain
        output_dir: Output directory
        timestamp: Timestamp string
    """
    print(f"\n📊 Generating iTOL metadata for {domain_name}...")

    # Create metadata for species counts per genus
    metadata_lines = []
    metadata_lines.append("DATASET_SIMPLEBAR")
    metadata_lines.append("SEPARATOR TAB")
    metadata_lines.append(f"DATASET_LABEL\t{domain_name}_species_per_genus")
    metadata_lines.append("COLOR\t#4a90e2")
    metadata_lines.append("LEGEND_TITLE\tSpecies Count")
    metadata_lines.append("LEGEND_SHAPES\t1")
    metadata_lines.append("LEGEND_COLORS\t#4a90e2")
    metadata_lines.append("LEGEND_LABELS\tNumber of species")
    metadata_lines.append("DATA")

    # Add data for each genus
    for genus_key, stats in genus_stats.items():
        species_count = stats['species_count']
        metadata_lines.append(f"{genus_key}\t{species_count}")

    # Save metadata
    metadata_file = output_dir / f'{domain_name.lower()}_species_count_{timestamp}.txt'
    with open(metadata_file, 'w') as f:
        f.write('\n'.join(metadata_lines))

    print(f"   ✅ Species count metadata saved: {metadata_file.name}")
    print(f"   📊 {len(genus_stats):,} genera with species counts")


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
        newick, genus_stats = build_tree_from_lineages(domain_df, domain)

        # Save Newick file
        newick_file = output_dir / f'{domain.lower()}_tree_{timestamp}.nwk'
        with open(newick_file, 'w') as f:
            f.write(newick)

        print(f"   💾 Newick tree saved: {newick_file.name}")
        print(f"   📏 Tree size: {len(newick):,} characters")

        # Generate metadata (species count per genus bar graph)
        generate_tree_metadata(genus_stats, domain, output_dir, timestamp)
    
    print("\n" + "="*70)
    print("✅ ALL TREES GENERATED!")
    print("="*70)
    print(f"\n📁 Output directory: {output_dir}")
    print(f"\nFiles generated for each domain:")
    print(f"   - <domain>_tree_<timestamp>.nwk              (Newick tree for iTOL)")
    print(f"   - <domain>_species_count_<timestamp>.txt     (Species count bar graph)")
    print(f"\n💡 Trees are collapsed to genus level")
    print(f"   Bar graphs show number of species per genus")


def main():
    # Find the most recent species_with_domain file
    sanity_output = Path(__file__).parent / 'output'

    # Look in both output/ and output/csv/ directories
    domain_files = list(sanity_output.glob('species_with_domain_*.csv'))
    domain_files.extend(list(sanity_output.glob('csv/species_with_domain_*.csv')))

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

