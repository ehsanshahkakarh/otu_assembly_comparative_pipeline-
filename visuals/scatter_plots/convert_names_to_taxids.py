#!/usr/bin/env python3
"""
Convert Taxonomic Names to TaxIDs in JSON Files
==============================================
Look up TaxIDs for taxonomic names in the overrepresented JSON files
and create new JSON files with both names and TaxIDs.
"""

import pandas as pd
import json
import os
from pathlib import Path

class TaxonomicNameToTaxIDConverter:
    def __init__(self):
        # Load NCBI data
        self.family_df = None
        self.genus_df = None
        self.load_ncbi_data()
        
    def load_ncbi_data(self):
        """Load NCBI CSV data."""
        family_file = "../../ncbi_parse/csv_ncbi/ncbi_family_counts.csv"
        genus_file = "../../ncbi_parse/csv_ncbi/ncbi_genus_counts.csv"
        
        if os.path.exists(family_file):
            self.family_df = pd.read_csv(family_file)
            print(f"📊 Loaded NCBI family data: {len(self.family_df):,} entries")
        else:
            print(f"❌ Family file not found: {family_file}")
            
        if os.path.exists(genus_file):
            self.genus_df = pd.read_csv(genus_file)
            print(f"📊 Loaded NCBI genus data: {len(self.genus_df):,} entries")
        else:
            print(f"❌ Genus file not found: {genus_file}")
    
    def find_taxid_for_name(self, taxon_name, domain=None):
        """Find TaxID for a taxonomic name by searching lineages."""
        results = []
        
        # Search in family data first
        if self.family_df is not None:
            family_matches = self.family_df[
                (self.family_df['family'] == taxon_name) |
                (self.family_df['lineage'].str.contains(taxon_name, na=False, case=False))
            ].copy()
            
            if domain:
                family_matches = family_matches[family_matches['domain'] == domain]
            
            for _, row in family_matches.iterrows():
                results.append({
                    'name': taxon_name,
                    'taxid': int(row['taxid']),
                    'matched_name': row['family'],
                    'domain': row['domain'],
                    'lineage': row['lineage'],
                    'source': 'family',
                    'genomes': int(row['family_genome_count']),
                    'species': int(row['family_species_count'])
                })
        
        # Search in genus data for higher-order taxa
        if self.genus_df is not None:
            genus_matches = self.genus_df[
                self.genus_df['lineage'].str.contains(taxon_name, na=False, case=False)
            ].copy()
            
            if domain:
                genus_matches = genus_matches[genus_matches['domain'] == domain]
            
            # Group by the taxonomic level that matches our search term
            # and find the most representative TaxID
            if not genus_matches.empty:
                # Get the TaxID that appears most frequently or has most genomes
                best_match = genus_matches.loc[genus_matches['genus_genome_count'].idxmax()]
                
                # Extract the specific TaxID for this taxonomic level from lineage_taxids
                lineage_parts = best_match['lineage'].split(';')
                lineage_taxids = str(best_match['lineage_taxids']).split(';')
                
                # Find which part of the lineage matches our search term
                matching_taxid = None
                for i, part in enumerate(lineage_parts):
                    if taxon_name.lower() in part.lower():
                        if i < len(lineage_taxids):
                            try:
                                matching_taxid = int(lineage_taxids[i])
                                break
                            except (ValueError, IndexError):
                                continue
                
                if matching_taxid is None:
                    matching_taxid = int(best_match['taxid'])  # Fallback to genus taxid
                
                # Calculate total genomes for this taxonomic group
                total_genomes = genus_matches['genus_genome_count'].sum()
                total_species = genus_matches['genus_species_count'].sum()
                
                results.append({
                    'name': taxon_name,
                    'taxid': int(matching_taxid),
                    'matched_name': f"{taxon_name} (from lineage)",
                    'domain': best_match['domain'],
                    'lineage': best_match['lineage'],
                    'source': 'genus_lineage',
                    'genomes': int(total_genomes),
                    'species': int(total_species),
                    'genera_count': int(len(genus_matches))
                })
        
        return results
    
    def convert_json_file(self, json_path, domain):
        """Convert a single JSON file to include TaxIDs."""
        print(f"\n🔍 Processing: {json_path}")
        
        # Load the JSON file
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        taxa_list = data.get('taxa', [])
        print(f"📋 Found {len(taxa_list)} taxa to convert")
        
        # Convert each taxon name to TaxID
        converted_taxa = []
        not_found = []
        
        for taxon_name in taxa_list:
            results = self.find_taxid_for_name(taxon_name, domain)
            
            if results:
                # Take the best result (highest genome count)
                best_result = max(results, key=lambda x: x.get('genomes', 0))
                converted_taxa.append(best_result)
                print(f"   ✅ {taxon_name} → TaxID: {best_result['taxid']} ({best_result['genomes']:,} genomes)")
            else:
                not_found.append(taxon_name)
                print(f"   ❌ {taxon_name} → No TaxID found")
        
        # Create new JSON structure
        new_data = {
            'taxa_with_taxids': converted_taxa,
            'original_taxa_names': taxa_list,
            'not_found': not_found,
            'conversion_stats': {
                'total_taxa': len(taxa_list),
                'converted': len(converted_taxa),
                'not_found': len(not_found),
                'domain': domain
            }
        }
        
        # Save the new JSON file
        new_json_path = json_path.replace('.json', '_with_taxids.json')
        with open(new_json_path, 'w') as f:
            json.dump(new_data, f, indent=2)
        
        print(f"💾 Saved: {new_json_path}")
        return new_json_path, converted_taxa, not_found
    
    def convert_all_overrepresented_files(self):
        """Convert all overrepresented JSON files."""
        print("🚀 CONVERTING OVERREPRESENTED TAXA NAMES TO TAXIDS")
        print("="*60)
        
        base_path = "source_data/taxa_extraction/overrepresented"
        
        # Define the files to convert
        files_to_convert = [
            ("bacteria/family_overrepresented_taxa.json", "Bacteria"),
            ("archaea/family_overrepresented_taxa.json", "Archaea"),
            ("eukaryota/family_overrepresented_taxa.json", "Eukaryota")
        ]
        
        all_results = {}
        
        for file_path, domain in files_to_convert:
            full_path = os.path.join(base_path, file_path)
            
            if os.path.exists(full_path):
                new_file, converted, not_found = self.convert_json_file(full_path, domain)
                all_results[domain] = {
                    'original_file': full_path,
                    'new_file': new_file,
                    'converted': len(converted),
                    'not_found': len(not_found),
                    'not_found_names': not_found
                }
            else:
                print(f"❌ File not found: {full_path}")
        
        # Summary
        print(f"\n🎯 CONVERSION SUMMARY:")
        print("="*40)
        
        total_converted = 0
        total_not_found = 0
        
        for domain, results in all_results.items():
            print(f"{domain}:")
            print(f"   ✅ Converted: {results['converted']}")
            print(f"   ❌ Not found: {results['not_found']}")
            if results['not_found_names']:
                print(f"   Missing: {', '.join(results['not_found_names'])}")
            
            total_converted += results['converted']
            total_not_found += results['not_found']
        
        print(f"\nOverall:")
        print(f"   ✅ Total converted: {total_converted}")
        print(f"   ❌ Total not found: {total_not_found}")
        
        return all_results

def main():
    """Main function."""
    converter = TaxonomicNameToTaxIDConverter()
    
    if converter.family_df is None or converter.genus_df is None:
        print("❌ Could not load NCBI data files")
        return
    
    results = converter.convert_all_overrepresented_files()
    
    print(f"\n📁 New JSON files with TaxIDs have been created!")
    print("Next step: Update the extraction script to use TaxIDs instead of name matching.")

if __name__ == "__main__":
    main()
