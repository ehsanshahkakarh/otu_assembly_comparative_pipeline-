#!/usr/bin/env python3
"""
Extract NCBI Overrepresented Taxonomic Diversity Using TaxIDs
============================================================
Extract all genera from NCBI data for overrepresented taxonomic groups,
using TaxIDs for precise matching instead of name searching.

Usage:
    python extract_ncbi_taxid_diversity.py --gene 16S --category overrepresented
    python extract_ncbi_taxid_diversity.py --gene 18S --category overrepresented
"""

import pandas as pd
import json
import os
import argparse
from collections import defaultdict, Counter
from pathlib import Path

class NCBITaxIDExtractor:
    def __init__(self):
        # Define NCBI CSV file paths
        self.ncbi_files = {
            "family": "../../ncbi_parse/csv_ncbi/ncbi_family_counts.csv",
            "genus": "../../ncbi_parse/csv_ncbi/ncbi_genus_counts.csv"
        }
        
        self.json_base_path = "source_data/taxa_extraction"
        
        # Domain mappings for gene types
        self.gene_domains = {
            "16S": ["bacteria", "archaea"],
            "18S": ["eukaryota"]
        }
        
    def load_taxa_from_json_with_taxids(self, category, domain, rank):
        """Load taxa with TaxIDs from JSON file."""
        json_file = f"{self.json_base_path}/{category}/{domain}/{rank}_{category}_taxa_with_taxids.json"
        
        if not os.path.exists(json_file):
            print(f"❌ TaxID JSON file not found: {json_file}")
            return []
            
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        taxa_list = data.get('taxa_with_taxids', [])
        print(f"📋 Loaded {len(taxa_list)} {category} {domain} {rank} taxa with TaxIDs from JSON")
        return taxa_list
    
    def load_ncbi_data(self, rank):
        """Load NCBI CSV data."""
        csv_file = self.ncbi_files[rank]
        
        if not os.path.exists(csv_file):
            print(f"❌ NCBI CSV file not found: {csv_file}")
            return None
            
        df = pd.read_csv(csv_file)
        print(f"📊 Loaded NCBI {rank} data: {len(df):,} total {rank} entries")
        return df
    
    def extract_genera_by_taxid(self, taxon_data, genus_df, gene, category, domain):
        """Extract all genera belonging to specified TaxIDs from NCBI data."""
        print(f"\n🔍 EXTRACTING NCBI DIVERSITY FOR {category.upper()} {domain.upper()} TAXONOMIC GROUPS (TaxID-based)")
        print("=" * 90)
        
        all_taxon_data = {}
        total_genomes_found = 0
        total_species_found = 0
        
        for taxon_info in taxon_data:
            taxon_name = taxon_info['name']
            taxon_taxid = taxon_info['taxid']
            
            print(f"\n🏷️  TAXONOMIC GROUP: {taxon_name} (TaxID: {taxon_taxid})")
            print("-" * 70)
            
            # Search for genera that have this TaxID in their lineage_taxids
            matching_genera = genus_df[
                genus_df['lineage_taxids'].astype(str).str.contains(str(taxon_taxid), na=False)
            ].copy()
            
            if matching_genera.empty:
                print(f"   ❌ No genera found with TaxID '{taxon_taxid}' in lineage")
                continue
            
            # Sort genera by genome count
            matching_genera = matching_genera.sort_values('genus_genome_count', ascending=False)
            
            # Calculate total statistics for this taxonomic group
            total_taxon_genomes = matching_genera['genus_genome_count'].sum()
            total_taxon_species = matching_genera['genus_species_count'].sum()
            
            total_genomes_found += total_taxon_genomes
            total_species_found += total_taxon_species
            
            print(f"   📊 Total genomes: {total_taxon_genomes:,} (across {len(matching_genera):,} genera)")
            print(f"   🌳 TAXONOMIC DIVERSITY:")
            print(f"      Genera: {len(matching_genera)} | Species: {total_taxon_species:,}")
            
            print(f"\n   🔬 TOP GENERA (by genome count):")
            genus_data = []
            
            for i, (_, genus_row) in enumerate(matching_genera.iterrows(), 1):
                genus_name = genus_row['genus']
                genus_genomes = genus_row['genus_genome_count']
                genus_species = genus_row['genus_species_count']
                genus_taxid = genus_row['taxid']
                genus_lineage = genus_row['lineage']
                
                genome_pct = (genus_genomes / total_taxon_genomes) * 100 if total_taxon_genomes > 0 else 0
                species_pct = (genus_species / total_taxon_species) * 100 if total_taxon_species > 0 else 0
                
                genus_data.append({
                    'genus': genus_name,
                    'genomes': int(genus_genomes),
                    'species': int(genus_species),
                    'taxid': int(genus_taxid),
                    'lineage': genus_lineage,
                    'genome_pct': genome_pct,
                    'species_pct': species_pct
                })
                
                if i <= 10:  # Show top 10
                    print(f"      {genus_name}: {genus_genomes:,} genomes ({genome_pct:.1f}%)")
            
            if len(genus_data) > 10:
                remaining = len(genus_data) - 10
                remaining_genomes = sum(g['genomes'] for g in genus_data[10:])
                print(f"      ... and {remaining} more genera ({remaining_genomes:,} genomes)")
            
            # Show detailed breakdown for top 3 genera
            print(f"\n   🔍 DETAILED BREAKDOWN (Top 3 genera):")
            for i, genus in enumerate(genus_data[:3], 1):
                print(f"\n      🧬 GENUS: {genus['genus']} ({genus['genomes']:,} genomes)")
                print(f"         TaxID: {genus['taxid']}")
                # Show truncated lineage for context
                lineage_parts = genus['lineage'].split(';')
                if len(lineage_parts) > 6:
                    context_lineage = ';'.join(lineage_parts[-4:])
                    print(f"         Context: ...;{context_lineage}")
                else:
                    print(f"         Lineage: {genus['lineage']}")
                print(f"         Genomes: {genus['genomes']:,} | Species: {genus['species']:,}")
            
            # Concentration analysis
            if len(genus_data) > 0:
                top_genus = genus_data[0]
                concentration = top_genus['genome_pct']
                
                if len(genus_data) > 1:
                    top_3_genomes = sum(g['genomes'] for g in genus_data[:3])
                    top_3_pct = (top_3_genomes / total_taxon_genomes) * 100 if total_taxon_genomes > 0 else 0
                    
                    print(f"\n   🎯 CONCENTRATION ANALYSIS:")
                    print(f"      Top genus: {top_genus['genus']} ({concentration:.1f}% of taxonomic group)")
                    print(f"      Top 3 genera: {top_3_pct:.1f}% of taxonomic group genomes")
                    print(f"      Genus diversity: {len(genus_data)} genera total")
                    
                    if concentration > 80:
                        concentration_level = "🔥 HIGHLY CONCENTRATED"
                        print(f"      🔥 HIGHLY CONCENTRATED: {top_genus['genus']} dominates this taxonomic group!")
                    elif top_3_pct > 80:
                        concentration_level = "📊 MODERATELY CONCENTRATED"
                        print(f"      📊 MODERATELY CONCENTRATED: Top 3 genera dominate")
                    else:
                        concentration_level = "🌳 DIVERSE"
                        print(f"      🌳 DIVERSE: Genomes spread across many genera")
                else:
                    concentration_level = "🔥 HIGHLY CONCENTRATED"
                    print(f"\n   🎯 SINGLE GENUS GROUP: {top_genus['genus']} (100% of taxonomic group)")
            else:
                concentration_level = "Unknown"
            
            # Store data for summary
            all_taxon_data[taxon_name] = {
                'taxid': taxon_taxid,
                'total_genomes': int(total_taxon_genomes),
                'total_genera': len(matching_genera),
                'genera_count': len(genus_data),
                'species_count': int(total_taxon_species),
                'genus_data': genus_data,
                'top_genus': genus_data[0] if genus_data else None,
                'concentration': genus_data[0]['genome_pct'] if genus_data else 0,
                'concentration_level': concentration_level
            }
            
            print("=" * 70)
        
        # Overall summary
        if all_taxon_data:
            print(f"\n🎯 OVERALL SUMMARY FOR {gene} {category.upper()} {domain.upper()} TAXONOMIC GROUPS:")
            print("=" * 90)
            print(f"Taxonomic groups analyzed: {len(all_taxon_data)}")
            print(f"Total genomes found: {total_genomes_found:,}")
            print(f"Average genomes per group: {total_genomes_found/len(all_taxon_data):,.1f}")
            
            # Sort taxonomic groups by genome count
            sorted_taxa = sorted(all_taxon_data.items(), 
                               key=lambda x: x[1]['total_genomes'], reverse=True)
            
            print(f"\n📊 TAXONOMIC GROUP RANKING (by genome count):")
            for i, (taxon, data) in enumerate(sorted_taxa, 1):
                top_genus = data['top_genus']
                concentration = data['concentration']
                
                print(f"   [{i:2d}] {taxon}: {data['total_genomes']:,} genomes")
                print(f"        Diversity: {data['genera_count']} genera, {data['species_count']} species")
                if top_genus:
                    print(f"        Top genus: {top_genus['genus']} ({concentration:.1f}% genomes)")
                print(f"        {data['concentration_level']}")
        
        return all_taxon_data

    def run_extraction(self, gene, category, output_file):
        """Run the NCBI TaxID-based diversity extraction for specified parameters."""
        print(f"🚀 STARTING NCBI TAXID-BASED OVERREPRESENTED DIVERSITY EXTRACTION")
        print(f"Gene: {gene} | Category: {category}")
        print(f"Output: {output_file}")

        if category != "overrepresented":
            print(f"❌ This script only works with category='overrepresented', got '{category}'")
            return

        # Load NCBI data
        genus_df = self.load_ncbi_data("genus")

        if genus_df is None:
            return

        # Determine domains to search
        domains = self.gene_domains.get(gene, [])

        all_results = {}

        # Extract for each domain
        for domain in domains:
            # Load taxonomic group data with TaxIDs from JSON
            taxon_data = self.load_taxa_from_json_with_taxids(category, domain, "family")
            if not taxon_data:
                continue

            # Filter NCBI data by domain
            if domain == "bacteria":
                domain_genus_df = genus_df[genus_df['domain'] == 'Bacteria'].copy()
            elif domain == "archaea":
                domain_genus_df = genus_df[genus_df['domain'] == 'Archaea'].copy()
            elif domain == "eukaryota":
                domain_genus_df = genus_df[genus_df['domain'] == 'Eukaryota'].copy()
            else:
                print(f"❌ Unknown domain: {domain}")
                continue

            print(f"\n📊 Domain {domain.upper()}: {len(domain_genus_df):,} genera")

            # Extract diversity for these taxonomic groups using TaxIDs
            domain_results = self.extract_genera_by_taxid(
                taxon_data, domain_genus_df, gene, category, domain
            )

            if domain_results:
                all_results[domain] = domain_results

        return all_results

def main():
    parser = argparse.ArgumentParser(description='Extract NCBI overrepresented taxonomic diversity using TaxIDs')
    parser.add_argument('--gene', choices=['16S', '18S'], required=True,
                       help='Gene type (16S for prokaryotes, 18S for eukaryotes)')
    parser.add_argument('--category', choices=['overrepresented'], required=True,
                       help='Taxa category (only overrepresented supported)')
    parser.add_argument('--output',
                       help='Output file path (default: family_diversity_results/{gene}_{category}_family_diversity.txt)')

    args = parser.parse_args()

    # Set default output path
    if args.output is None:
        output_dir = "family_diversity_results"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        args.output = f"{output_dir}/{args.gene}_{args.category}_family_diversity.txt"

    extractor = NCBITaxIDExtractor()
    results = extractor.run_extraction(args.gene, args.category, args.output)

    return results

if __name__ == "__main__":
    main()
