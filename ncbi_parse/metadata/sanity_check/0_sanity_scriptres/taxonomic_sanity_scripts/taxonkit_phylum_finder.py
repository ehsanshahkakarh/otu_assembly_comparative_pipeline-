#!/usr/bin/env python3
"""
Script to extract taxids from unmapped entries, run them through taxonkit,
and identify entries that have phyla in their taxonomic lineage.
"""

import pandas as pd
import subprocess
import os
import sys
from tqdm import tqdm
import tempfile
from pathlib import Path

def extract_taxids_from_unmapped(unmapped_file):
    """Extract unique taxids from the unmapped entries file."""
    print(f"📊 Reading unmapped entries from: {unmapped_file}")
    
    try:
        # Read the unmapped entries file
        df = pd.read_csv(unmapped_file, comment='#')
        
        # Extract unique taxids (both taxid and species_taxid columns)
        taxids = set()
        
        if 'taxid' in df.columns:
            taxids.update(df['taxid'].dropna().astype(int))
        
        if 'species_taxid' in df.columns:
            taxids.update(df['species_taxid'].dropna().astype(int))
        
        taxids = sorted(list(taxids))
        print(f"✅ Extracted {len(taxids)} unique taxids")
        
        return taxids, df
        
    except Exception as e:
        print(f"❌ Error reading unmapped file: {e}")
        sys.exit(1)

def run_taxonkit_lineage(taxids, output_dir):
    """Run taxonkit lineage command on the taxids."""
    print(f"🔍 Running taxonkit lineage on {len(taxids)} taxids...")
    
    # Create temporary file with taxids
    temp_taxids_file = os.path.join(output_dir, "temp_taxids.txt")
    temp_lineage_file = os.path.join(output_dir, "temp_lineages.txt")
    
    try:
        # Write taxids to temporary file
        with open(temp_taxids_file, 'w') as f:
            for taxid in taxids:
                f.write(f"{taxid}\n")
        
        print(f"📝 Created temporary taxid file: {temp_taxids_file}")
        
        # Run taxonkit lineage command
        cmd = f"cat {temp_taxids_file} | taxonkit lineage -R > {temp_lineage_file}"
        print(f"🚀 Running command: {cmd}")
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ Error running taxonkit: {result.stderr}")
            return None
        
        print(f"✅ Taxonkit completed successfully")
        print(f"📄 Lineage results saved to: {temp_lineage_file}")
        
        return temp_lineage_file
        
    except Exception as e:
        print(f"❌ Error running taxonkit: {e}")
        return None

def parse_taxonkit_results(lineage_file):
    """Parse taxonkit lineage results and identify entries with phyla."""
    print(f"📖 Parsing taxonkit results from: {lineage_file}")
    
    results = []
    phylum_found_count = 0
    
    try:
        with open(lineage_file, 'r') as f:
            for line in tqdm(f, desc="Processing lineage results"):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 3:
                    taxid = parts[0]
                    lineage = parts[1] if len(parts) > 1 else ""
                    ranks = parts[2] if len(parts) > 2 else ""
                    
                    # Check if 'phylum' is in the ranks
                    has_phylum = 'phylum' in ranks.lower()
                    
                    # Extract phylum name if present
                    phylum_name = ""
                    if has_phylum and lineage and ranks:
                        lineage_parts = lineage.split(';')
                        rank_parts = ranks.split(';')
                        
                        for i, rank in enumerate(rank_parts):
                            if rank.lower().strip() == 'phylum' and i < len(lineage_parts):
                                phylum_name = lineage_parts[i].strip()
                                break
                    
                    results.append({
                        'taxid': int(taxid),
                        'lineage': lineage,
                        'ranks': ranks,
                        'has_phylum': has_phylum,
                        'phylum_name': phylum_name
                    })
                    
                    if has_phylum:
                        phylum_found_count += 1
        
        print(f"✅ Processed {len(results)} lineage results")
        print(f"🎯 Found {phylum_found_count} entries with phylum information")
        
        return results
        
    except Exception as e:
        print(f"❌ Error parsing taxonkit results: {e}")
        return []

def merge_with_unmapped_data(taxonkit_results, unmapped_df):
    """Merge taxonkit results with original unmapped data."""
    print("🔗 Merging taxonkit results with unmapped entry data...")
    
    # Convert taxonkit results to DataFrame
    taxonkit_df = pd.DataFrame(taxonkit_results)
    
    # Merge with unmapped data on taxid
    merged_df = unmapped_df.merge(
        taxonkit_df, 
        left_on='taxid', 
        right_on='taxid', 
        how='left'
    )
    
    # Also try merging on species_taxid for entries where taxid didn't match
    if 'species_taxid' in unmapped_df.columns:
        # For entries that didn't get lineage info from taxid, try species_taxid
        no_lineage_mask = merged_df['lineage'].isna()
        species_merge = unmapped_df[no_lineage_mask].merge(
            taxonkit_df,
            left_on='species_taxid',
            right_on='taxid',
            how='left',
            suffixes=('', '_species')
        )
        
        # Update the merged dataframe with species_taxid results
        for idx in species_merge.index:
            if not pd.isna(species_merge.loc[idx, 'lineage_species']):
                merged_df.loc[idx, 'lineage'] = species_merge.loc[idx, 'lineage_species']
                merged_df.loc[idx, 'ranks'] = species_merge.loc[idx, 'ranks_species']
                merged_df.loc[idx, 'has_phylum'] = species_merge.loc[idx, 'has_phylum_species']
                merged_df.loc[idx, 'phylum_name'] = species_merge.loc[idx, 'phylum_name_species']
    
    print(f"✅ Merged data for {len(merged_df)} entries")
    
    return merged_df

def save_phylum_found_entries(merged_df, output_dir):
    """Save entries that have phylum information."""
    # Filter entries with phylum information
    phylum_entries = merged_df[merged_df['has_phylum'] == True].copy()
    
    if len(phylum_entries) == 0:
        print("⚠️  No entries with phylum information found")
        return
    
    # Sort by phylum name and organism name
    phylum_entries = phylum_entries.sort_values(['phylum_name', 'organism_name'])
    
    output_file = os.path.join(output_dir, "phylum_found_entries.csv")
    
    # Select relevant columns for output
    output_columns = [
        'assembly_accession', 'taxid', 'species_taxid', 'organism_name', 
        'infraspecific_name', 'isolate', 'assembly_level', 'genome_rep', 
        'asm_name', 'phylum_name', 'lineage', 'ranks'
    ]
    
    # Only include columns that exist
    available_columns = [col for col in output_columns if col in phylum_entries.columns]
    output_df = phylum_entries[available_columns]
    
    # Save with header statistics
    with open(output_file, 'w') as f:
        f.write("# UNMAPPED ENTRIES WITH PHYLUM INFORMATION FOUND VIA TAXONKIT\n")
        f.write(f"# Total entries with phylum found: {len(phylum_entries)}\n")
        f.write(f"# Unique phylums found: {phylum_entries['phylum_name'].nunique()}\n")
        f.write("#\n")
        f.write("# TOP 10 PHYLUMS BY ENTRY COUNT:\n")
        
        phylum_counts = phylum_entries['phylum_name'].value_counts().head(10)
        for i, (phylum, count) in enumerate(phylum_counts.items(), 1):
            f.write(f"# {i:2d}. {phylum:<25} {count:>6} entries\n")
        
        f.write("#\n")
        f.write("# DATA COLUMNS: " + ", ".join(available_columns) + "\n")
        f.write("#\n")
        
        # Write CSV data
        output_df.to_csv(f, index=False)
    
    print(f"💾 Saved {len(phylum_entries)} entries with phylum info to: {output_file}")
    
    # Show summary
    print(f"\n📊 PHYLUM SUMMARY:")
    print(f"   Total entries with phylum: {len(phylum_entries)}")
    print(f"   Unique phylums found: {phylum_entries['phylum_name'].nunique()}")
    print(f"\n🏆 Top 5 phylums:")
    for i, (phylum, count) in enumerate(phylum_counts.head().items(), 1):
        print(f"   {i}. {phylum}: {count} entries")

def cleanup_temp_files(output_dir):
    """Clean up temporary files."""
    temp_files = [
        os.path.join(output_dir, "temp_taxids.txt"),
        os.path.join(output_dir, "temp_lineages.txt")
    ]
    
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"🗑️  Cleaned up: {temp_file}")

def main():
    # Setup paths using flexible path handling
    script_dir = Path(__file__).resolve().parent

    # File paths (parent directory is 0_sanity_scriptres)
    unmapped_file = script_dir.parent / "phylum_sanity" / "phylum_unmapped_entries.csv"
    output_dir = script_dir.parent / "phylum_sanity"
    
    print("🚀 Starting taxonkit phylum finder...")
    print(f"📁 Input file: {unmapped_file}")
    print(f"📁 Output directory: {output_dir}")
    
    # Check if input file exists
    if not os.path.exists(unmapped_file):
        print(f"❌ Input file not found: {unmapped_file}")
        sys.exit(1)
    
    # Extract taxids from unmapped entries
    taxids, unmapped_df = extract_taxids_from_unmapped(unmapped_file)
    
    # Run taxonkit lineage
    lineage_file = run_taxonkit_lineage(taxids, output_dir)
    if not lineage_file:
        sys.exit(1)
    
    # Parse taxonkit results
    taxonkit_results = parse_taxonkit_results(lineage_file)
    if not taxonkit_results:
        sys.exit(1)
    
    # Merge with unmapped data
    merged_df = merge_with_unmapped_data(taxonkit_results, unmapped_df)
    
    # Save entries with phylum information
    save_phylum_found_entries(merged_df, output_dir)
    
    # Cleanup temporary files
    cleanup_temp_files(output_dir)
    
    print("\n✅ Taxonkit phylum finder completed successfully!")

if __name__ == "__main__":
    main()
