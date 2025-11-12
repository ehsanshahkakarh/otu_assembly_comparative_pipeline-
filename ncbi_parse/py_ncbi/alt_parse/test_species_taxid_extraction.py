#!/usr/bin/env python3
"""
Test Species Taxid Extraction

This script tests the species_taxid extraction logic to ensure it works correctly
with the NCBI assembly data.

Usage:
    python test_species_taxid_extraction.py --test-extraction
    python test_species_taxid_extraction.py --compare-approaches
"""

import subprocess
import tempfile
import os
import pandas as pd
import argparse
from pathlib import Path

def test_taxonkit_species_extraction():
    """Test the taxonkit species_taxid extraction with sample data."""
    print("🧪 Testing taxonkit species_taxid extraction")
    print("=" * 60)
    
    # Sample taxids from your terminal output
    sample_taxids = [
        "58343",    # Streptomyces canus
        "68286",    # Streptomyces zaomyceticus  
        "1892",     # Streptomyces anulatus
        "285562",   # Streptomyces coelicoflavus
        "2903613",  # Streptomyces sp. NBC_00029
        "2903614",  # Streptomyces sp. NBC_00035
    ]
    
    print(f"Testing with {len(sample_taxids)} sample taxids:")
    for taxid in sample_taxids:
        print(f"   {taxid}")
    
    # Create temporary file with taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in sample_taxids:
            temp_file.write(f"{taxid}\n")
    
    try:
        # Test taxonkit reformat command
        print(f"\n🔍 Running taxonkit reformat...")
        result = subprocess.run(
            ["taxonkit", "reformat", "--taxid-field", "1", "--format", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode == 0 and result.stdout.strip():
            print(f"✅ Taxonkit command successful")
            
            # Parse the output
            taxid_to_species_taxid = {}
            lines = result.stdout.strip().split('\n')
            
            print(f"\n📊 Parsing results:")
            for i, line in enumerate(lines):
                parts = line.split('\t')
                if len(parts) >= 8:
                    taxid = parts[0]
                    kingdom = parts[0] if len(parts) > 0 else ""
                    phylum = parts[1] if len(parts) > 1 else ""
                    class_name = parts[2] if len(parts) > 2 else ""
                    order = parts[3] if len(parts) > 3 else ""
                    family = parts[4] if len(parts) > 4 else ""
                    genus = parts[5] if len(parts) > 5 else ""
                    species = parts[6] if len(parts) > 6 else ""
                    taxid_full = parts[7] if len(parts) > 7 else ""
                    
                    print(f"\n   Taxid {taxid}:")
                    print(f"      Species: {species}")
                    print(f"      Full taxid: {taxid_full}")
                    
                    # Extract species taxid from species info
                    if "[taxid:" in species:
                        species_taxid = species.split("[taxid:")[1].split("]")[0]
                        taxid_to_species_taxid[taxid] = species_taxid
                        print(f"      Species taxid: {species_taxid}")
                    else:
                        print(f"      ⚠️  No species taxid found in: {species}")
            
            print(f"\n✅ Successfully extracted {len(taxid_to_species_taxid)} species taxids:")
            for taxid, species_taxid in taxid_to_species_taxid.items():
                print(f"   {taxid} → {species_taxid}")
            
            return taxid_to_species_taxid
            
        else:
            print(f"❌ Taxonkit command failed")
            print(f"Return code: {result.returncode}")
            print(f"Error: {result.stderr}")
            return {}
            
    except Exception as e:
        print(f"❌ Error running taxonkit: {e}")
        return {}
    finally:
        try:
            os.unlink(temp_filename)
        except:
            pass

def compare_subsetting_approaches():
    """Compare different species subsetting approaches."""
    print("\n🔍 Comparing Species Subsetting Approaches")
    print("=" * 60)
    
    # Sample organism names from your data
    sample_data = [
        {"taxid": "58343", "organism_name": "Streptomyces canus strain=NBC_00002"},
        {"taxid": "58343", "organism_name": "Streptomyces canus strain=NBC_00087"},
        {"taxid": "68286", "organism_name": "Streptomyces zaomyceticus strain=NBC_00088"},
        {"taxid": "1892", "organism_name": "Streptomyces anulatus strain=NBC_00107"},
        {"taxid": "285562", "organism_name": "Streptomyces coelicoflavus strain=NBC_00099"},
        {"taxid": "2903613", "organism_name": "Streptomyces sp. NBC_00029 strain=NBC_00029"},
        {"taxid": "2903614", "organism_name": "Streptomyces sp. NBC_00035 strain=NBC_00035"},
        {"taxid": "2903615", "organism_name": "Streptomyces sp. NBC_00040 strain=NBC_00040"},
    ]
    
    df = pd.DataFrame(sample_data)
    print(f"📋 Sample data ({len(df)} entries):")
    for _, row in df.iterrows():
        print(f"   {row['taxid']}: {row['organism_name']}")
    
    # Approach 1: Using organism_name directly (old script behavior)
    print(f"\n1️⃣  Approach 1: Using organism_name directly (strain-level)")
    organism_name_groups = df.groupby('organism_name').size()
    print(f"   Unique organism names: {len(organism_name_groups)}")
    print(f"   Result: Keeps {len(organism_name_groups)} entries (one per strain)")
    
    # Approach 2: Extracting species from organism_name
    print(f"\n2️⃣  Approach 2: Extracting species from organism_name")
    df['extracted_species'] = df['organism_name'].str.split().str[:2].str.join(' ')
    species_groups = df.groupby('extracted_species').size()
    print(f"   Unique extracted species: {len(species_groups)}")
    print(f"   Species groups:")
    for species, count in species_groups.items():
        print(f"      {species}: {count} strains")
    print(f"   Result: Keeps {len(species_groups)} entries (one per species)")
    
    # Approach 3: Using species_taxid (if we had it)
    print(f"\n3️⃣  Approach 3: Using species_taxid (most accurate)")
    
    # Simulate species_taxid mapping based on known relationships
    # In reality, this would come from taxonkit
    species_taxid_mapping = {
        "58343": "58343",      # Streptomyces canus (both strains same species)
        "68286": "68286",      # Streptomyces zaomyceticus
        "1892": "1892",        # Streptomyces anulatus
        "285562": "285562",    # Streptomyces coelicoflavus
        "2903613": "1883",     # Streptomyces sp. (likely maps to genus Streptomyces)
        "2903614": "1883",     # Streptomyces sp. (likely maps to genus Streptomyces)
        "2903615": "1883",     # Streptomyces sp. (likely maps to genus Streptomyces)
    }
    
    df['species_taxid'] = df['taxid'].map(species_taxid_mapping)
    species_taxid_groups = df.groupby('species_taxid').size()
    print(f"   Unique species taxids: {len(species_taxid_groups)}")
    print(f"   Species taxid groups:")
    for species_taxid, count in species_taxid_groups.items():
        print(f"      {species_taxid}: {count} strains")
    print(f"   Result: Keeps {len(species_taxid_groups)} entries (one per true species)")
    
    # Summary
    print(f"\n📊 Summary:")
    print(f"   Organism name approach: {len(organism_name_groups)} entries (strain-level)")
    print(f"   Extracted species approach: {len(species_groups)} entries (species-level)")
    print(f"   Species taxid approach: {len(species_taxid_groups)} entries (true species-level)")
    
    print(f"\n💡 Recommendation:")
    print(f"   Use species_taxid for the most biologically accurate species-level subsetting")
    print(f"   This avoids issues with strain naming conventions and provides taxonomic accuracy")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Test species taxid extraction and subsetting approaches")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--test-extraction", action="store_true", help="Test taxonkit species_taxid extraction")
    group.add_argument("--compare-approaches", action="store_true", help="Compare different subsetting approaches")
    group.add_argument("--test-both", action="store_true", help="Run both tests")
    
    args = parser.parse_args()
    
    if args.test_extraction or args.test_both:
        species_mapping = test_taxonkit_species_extraction()
    
    if args.compare_approaches or args.test_both:
        compare_subsetting_approaches()
    
    if args.test_both:
        print(f"\n🎉 Both tests completed!")
        print(f"The updated genus scripts will now use species_taxid for true species-level subsetting")

if __name__ == "__main__":
    main()
