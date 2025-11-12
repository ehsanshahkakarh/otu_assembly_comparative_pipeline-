#!/usr/bin/env python3
"""
Enhanced hierarchical analysis of 18S NCBI merged data.
Creates a comprehensive taxonomic hierarchy with novelty factors and counts.

Shows:
- Phyla with novelty factors, OTU counts, size counts, NCBI counts
- Families under each phylum (ordered by novelty factor)
- Genera under each family (ordered by novelty factor)
- Handles mismatches gracefully

Created: 2025-10-13
"""

import pandas as pd
import os
import subprocess
import tempfile
from datetime import datetime
# from collections import defaultdict  # Not used in current implementation

def load_merged_data():
    """Load the 18S NCBI merged data files."""
    base_dir = ".."  # Go up one level from txt_sanity to csv_results
    
    # Load the three taxonomic levels
    phylum_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_phylum.csv"))
    family_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_family.csv"))
    genus_df = pd.read_csv(os.path.join(base_dir, "18s_ncbi_merged_clean_genus.csv"))
    
    # Filter for matched entries only
    phylum_matched = phylum_df[phylum_df['match_status'] == 'matched'].copy()
    family_matched = family_df[family_df['match_status'] == 'matched'].copy()
    genus_matched = genus_df[genus_df['match_status'] == 'matched'].copy()
    
    return phylum_matched, family_matched, genus_matched

def calculate_novelty_factor(row):
    """Calculate novelty factor (OTU count / NCBI species count)."""
    if row['ncbi_species_count'] > 0:
        return round(row['census_otu_count'] / row['ncbi_species_count'], 2)
    else:
        return float('inf')  # Infinite novelty for taxa with no NCBI species

def get_lineages_for_taxa(taxa_list):
    """Get lineages for a list of taxa using taxonkit."""
    if not taxa_list:
        return {}

    name_to_lineage = {}

    try:
        # Create temporary file with taxa names
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            temp_filename = temp_file.name
            for taxon in taxa_list:
                temp_file.write(f"{taxon}\n")

        # First get taxids
        cmd1 = ['taxonkit', 'name2taxid', temp_filename]
        result1 = subprocess.run(cmd1, capture_output=True, text=True, check=True)

        # Parse name2taxid results to get taxids
        name_to_taxid = {}
        for line in result1.stdout.strip().split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()
                    if taxid and taxid != '0':
                        name_to_taxid[name] = taxid

        # Get lineages for valid taxids
        if name_to_taxid:
            # Create temporary file with taxids
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as taxid_file:
                taxid_filename = taxid_file.name
                for taxid in name_to_taxid.values():
                    taxid_file.write(f"{taxid}\n")

            # Get lineages
            cmd2 = ['taxonkit', 'lineage', '-R', taxid_filename]
            result2 = subprocess.run(cmd2, capture_output=True, text=True, check=True)

            # Parse lineage results
            taxid_to_lineage = {}
            for line in result2.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        taxid = parts[0].strip()
                        lineage = parts[1].strip() if len(parts) > 1 else ""
                        if lineage:
                            taxid_to_lineage[taxid] = lineage

            # Map back to names
            for name, taxid in name_to_taxid.items():
                if taxid in taxid_to_lineage:
                    name_to_lineage[name] = taxid_to_lineage[taxid]

            # Clean up taxid file
            os.unlink(taxid_filename)

        # Clean up name file
        os.unlink(temp_filename)

    except Exception as e:
        print(f"Error getting lineages: {e}")

    return name_to_lineage

def find_parent_in_lineage(lineage, parent_list):
    """Find which parent taxon appears in the lineage."""
    if not lineage:
        return None
    
    lineage_lower = lineage.lower()
    for parent in parent_list:
        if parent.lower() in lineage_lower:
            return parent
    return None

def prepare_taxa_data(df, level_name):
    """Prepare taxa data with novelty factors and sorting."""
    taxa_data = []
    
    for _, row in df.iterrows():
        novelty = calculate_novelty_factor(row)
        taxa_data.append({
            'name': row[level_name],
            'census_otu_count': row['census_otu_count'],
            'census_size_count': row.get('census_size_count', 0),  # Handle missing column
            'ncbi_species_count': row['ncbi_species_count'],
            'novelty_factor': novelty
        })
    
    # Sort by novelty factor (descending), then by OTU count (descending)
    taxa_data.sort(key=lambda x: (-x['novelty_factor'] if x['novelty_factor'] != float('inf') else -999999, 
                                  -x['census_otu_count']))
    
    return taxa_data

def create_hierarchical_structure(phylum_data, family_data, genus_data):
    """Create the hierarchical structure with proper parent-child relationships."""
    
    print("Getting lineages for taxonomic mapping...")
    
    # Get all unique taxa names
    all_phyla = [p['name'] for p in phylum_data]
    all_families = [f['name'] for f in family_data]
    all_genera = [g['name'] for g in genus_data]
    
    # Get lineages
    family_lineages = get_lineages_for_taxa(all_families)
    genus_lineages = get_lineages_for_taxa(all_genera)
    
    print(f"Got lineages for {len(family_lineages)}/{len(all_families)} families")
    print(f"Got lineages for {len(genus_lineages)}/{len(all_genera)} genera")
    
    # Create hierarchical structure
    hierarchy = {}
    unmatched_families = []
    unmatched_genera = []
    
    # Initialize phyla
    for phylum in phylum_data:
        hierarchy[phylum['name']] = {
            'data': phylum,
            'families': {},
            'orphan_genera': []  # Genera that match phylum but no family
        }
    
    # Map families to phyla
    for family in family_data:
        family_name = family['name']
        lineage = family_lineages.get(family_name, "")
        parent_phylum = find_parent_in_lineage(lineage, all_phyla)
        
        if parent_phylum and parent_phylum in hierarchy:
            hierarchy[parent_phylum]['families'][family_name] = {
                'data': family,
                'genera': []
            }
        else:
            unmatched_families.append(family)
    
    # Map genera to families and phyla
    for genus in genus_data:
        genus_name = genus['name']
        lineage = genus_lineages.get(genus_name, "")
        
        # Try to find parent family first
        parent_family = find_parent_in_lineage(lineage, all_families)
        parent_phylum = find_parent_in_lineage(lineage, all_phyla)
        
        placed = False
        
        # If genus matches a family, place it under that family
        if parent_family:
            for phylum_name, phylum_info in hierarchy.items():
                if parent_family in phylum_info['families']:
                    phylum_info['families'][parent_family]['genera'].append(genus)
                    placed = True
                    break
        
        # If genus matches phylum but not family, place as orphan under phylum
        if not placed and parent_phylum and parent_phylum in hierarchy:
            hierarchy[parent_phylum]['orphan_genera'].append(genus)
            placed = True
        
        # If genus doesn't match anything, add to unmatched
        if not placed:
            unmatched_genera.append(genus)
    
    return hierarchy, unmatched_families, unmatched_genera

def format_novelty(novelty):
    """Format novelty factor for display."""
    if novelty == float('inf'):
        return "∞"
    else:
        return f"{novelty:.2f}"

def write_hierarchical_output(hierarchy, unmatched_families, unmatched_genera, output_file):
    """Write the hierarchical structure to output file."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(output_file, 'w') as f:
        f.write("18S NCBI Enhanced Hierarchical Analysis\n")
        f.write(f"Generated: {timestamp}\n")
        f.write("Organized by taxonomic hierarchy with novelty factors\n")
        f.write("Novelty Factor = Census_OTU_Count / NCBI_Species_Count\n")
        f.write("=" * 100 + "\n\n")

        # Sort phyla by novelty factor
        sorted_phyla = sorted(hierarchy.items(),
                            key=lambda x: (-x[1]['data']['novelty_factor'] if x[1]['data']['novelty_factor'] != float('inf') else -999999,
                                         -x[1]['data']['census_otu_count']))

        for phylum_name, phylum_info in sorted_phyla:
            phylum_data = phylum_info['data']

            # Write phylum header
            f.write(f"PHYLUM: {phylum_name}\n")
            f.write("=" * 80 + "\n")
            f.write(f"Novelty Factor: {format_novelty(phylum_data['novelty_factor'])} | ")
            f.write(f"OTU Count: {phylum_data['census_otu_count']} | ")
            f.write(f"Size Count: {phylum_data['census_size_count']} | ")
            f.write(f"NCBI Species: {phylum_data['ncbi_species_count']}\n")
            f.write("-" * 80 + "\n\n")

            # Count families and genera
            family_count = len(phylum_info['families'])
            orphan_genus_count = len(phylum_info['orphan_genera'])
            total_genus_count = sum(len(fam_info['genera']) for fam_info in phylum_info['families'].values()) + orphan_genus_count

            f.write(f"  Families: {family_count} | Genera: {total_genus_count} (Orphan genera: {orphan_genus_count})\n\n")

            # Write families (sorted by novelty factor)
            if phylum_info['families']:
                sorted_families = sorted(phylum_info['families'].items(),
                                       key=lambda x: (-x[1]['data']['novelty_factor'] if x[1]['data']['novelty_factor'] != float('inf') else -999999,
                                                    -x[1]['data']['census_otu_count']))

                for family_name, family_info in sorted_families:
                    family_data = family_info['data']
                    genus_count = len(family_info['genera'])

                    f.write(f"    FAMILY: {family_name}\n")
                    f.write(f"    Novelty: {format_novelty(family_data['novelty_factor'])} | ")
                    f.write(f"OTU: {family_data['census_otu_count']} | ")
                    f.write(f"Size: {family_data['census_size_count']} | ")
                    f.write(f"NCBI: {family_data['ncbi_species_count']} | ")
                    f.write(f"Genera: {genus_count}\n")

                    # Write genera under this family (sorted by novelty factor)
                    if family_info['genera']:
                        sorted_genera = sorted(family_info['genera'],
                                             key=lambda x: (-x['novelty_factor'] if x['novelty_factor'] != float('inf') else -999999,
                                                           -x['census_otu_count']))

                        for genus in sorted_genera:
                            f.write(f"      → GENUS: {genus['name']} | ")
                            f.write(f"Novelty: {format_novelty(genus['novelty_factor'])} | ")
                            f.write(f"OTU: {genus['census_otu_count']} | ")
                            f.write(f"Size: {genus['census_size_count']} | ")
                            f.write(f"NCBI: {genus['ncbi_species_count']}\n")
                    else:
                        f.write(f"      (No genera found under this family)\n")

                    f.write("\n")

            # Write orphan genera (genera that match phylum but no family)
            if phylum_info['orphan_genera']:
                f.write(f"    ORPHAN GENERA (match phylum but no family):\n")
                f.write(f"    " + "-" * 60 + "\n")

                sorted_orphans = sorted(phylum_info['orphan_genera'],
                                      key=lambda x: (-x['novelty_factor'] if x['novelty_factor'] != float('inf') else -999999,
                                                    -x['census_otu_count']))

                for genus in sorted_orphans:
                    f.write(f"      → GENUS: {genus['name']} | ")
                    f.write(f"Novelty: {format_novelty(genus['novelty_factor'])} | ")
                    f.write(f"OTU: {genus['census_otu_count']} | ")
                    f.write(f"Size: {genus['census_size_count']} | ")
                    f.write(f"NCBI: {genus['ncbi_species_count']}\n")
                f.write("\n")

            f.write("=" * 100 + "\n\n")

        # Write unmatched families
        if unmatched_families:
            f.write("UNMATCHED FAMILIES (no phylum found in lineage)\n")
            f.write("=" * 80 + "\n")
            f.write(f"Total unmatched families: {len(unmatched_families)}\n\n")

            sorted_unmatched_fam = sorted(unmatched_families,
                                        key=lambda x: (-x['novelty_factor'] if x['novelty_factor'] != float('inf') else -999999,
                                                     -x['census_otu_count']))

            for family in sorted_unmatched_fam:
                f.write(f"  FAMILY: {family['name']} | ")
                f.write(f"Novelty: {format_novelty(family['novelty_factor'])} | ")
                f.write(f"OTU: {family['census_otu_count']} | ")
                f.write(f"Size: {family['census_size_count']} | ")
                f.write(f"NCBI: {family['ncbi_species_count']}\n")

            f.write("\n" + "=" * 100 + "\n\n")

        # Write unmatched genera
        if unmatched_genera:
            f.write("UNMATCHED GENERA (no phylum or family found in lineage)\n")
            f.write("=" * 80 + "\n")
            f.write(f"Total unmatched genera: {len(unmatched_genera)}\n\n")

            sorted_unmatched_gen = sorted(unmatched_genera,
                                        key=lambda x: (-x['novelty_factor'] if x['novelty_factor'] != float('inf') else -999999,
                                                     -x['census_otu_count']))

            for genus in sorted_unmatched_gen:
                f.write(f"  GENUS: {genus['name']} | ")
                f.write(f"Novelty: {format_novelty(genus['novelty_factor'])} | ")
                f.write(f"OTU: {genus['census_otu_count']} | ")
                f.write(f"Size: {genus['census_size_count']} | ")
                f.write(f"NCBI: {genus['ncbi_species_count']}\n")

            f.write("\n" + "=" * 100 + "\n")

def main():
    """Main function to run the enhanced hierarchical analysis."""
    print("Enhanced 18S NCBI Hierarchical Analysis")
    print("=" * 60)

    # Load data
    print("Loading merged data...")
    phylum_df, family_df, genus_df = load_merged_data()

    print(f"Loaded:")
    print(f"  - Matched Phyla: {len(phylum_df)}")
    print(f"  - Matched Families: {len(family_df)}")
    print(f"  - Matched Genera: {len(genus_df)}")
    print()

    # Prepare data with novelty factors
    print("Calculating novelty factors...")
    phylum_data = prepare_taxa_data(phylum_df, 'phylum')
    family_data = prepare_taxa_data(family_df, 'family')
    genus_data = prepare_taxa_data(genus_df, 'genus')

    # Create hierarchical structure
    print("Building hierarchical structure...")
    hierarchy, unmatched_families, unmatched_genera = create_hierarchical_structure(
        phylum_data, family_data, genus_data)

    # Write output
    output_file = "18s_ncbi_enhanced_hierarchical_analysis.txt"
    print(f"Writing results to {output_file}...")
    write_hierarchical_output(hierarchy, unmatched_families, unmatched_genera, output_file)

    # Print summary
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print(f"Results saved to: {output_file}")
    print()
    print("Summary:")
    print(f"  - Phyla processed: {len(hierarchy)}")
    print(f"  - Families mapped to phyla: {sum(len(p['families']) for p in hierarchy.values())}")
    print(f"  - Genera mapped to families: {sum(sum(len(f['genera']) for f in p['families'].values()) for p in hierarchy.values())}")
    print(f"  - Orphan genera (phylum only): {sum(len(p['orphan_genera']) for p in hierarchy.values())}")
    print(f"  - Unmatched families: {len(unmatched_families)}")
    print(f"  - Unmatched genera: {len(unmatched_genera)}")

if __name__ == "__main__":
    main()
