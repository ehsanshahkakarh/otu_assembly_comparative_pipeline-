#!/usr/bin/env python3
"""
Simple test script to debug lineage generation
"""

import subprocess
import os
import pandas as pd

# Set up environment
env = os.environ.copy()
env["TAXONKIT_DB"] = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"

def test_name_to_taxid():
    """Test name to taxid conversion"""
    print("🔍 Testing name to taxid conversion...")
    
    test_names = ["Opisthokonta", "Pythium", "Agaricomycetes", "Centraloplasthelida", "Insecta"]
    
    result = subprocess.run(
        ["taxonkit", "name2taxid"],
        input="\n".join(test_names),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env
    )
    
    print(f"Return code: {result.returncode}")
    print("Results:")
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        for i, line in enumerate(lines):
            parts = line.strip().split('\t')
            name = test_names[i] if i < len(test_names) else "Unknown"
            taxid = parts[1] if len(parts) >= 2 and parts[1] != "0" else "NOT_FOUND"
            print(f"  {name} → {taxid}")
    
    if result.stderr:
        print(f"Errors: {result.stderr}")
    
    return result

def test_lineage_generation():
    """Test lineage generation"""
    print("\n🌳 Testing lineage generation...")
    
    test_taxids = ["33154", "4797", "155619"]
    
    result = subprocess.run(
        ["taxonkit", "lineage", "--show-name"],
        input="\n".join(test_taxids),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env
    )
    
    print(f"Return code: {result.returncode}")
    print("Results:")
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                taxid = parts[0]
                lineage = parts[2]  # This should be the lineage
                print(f"  {taxid} → {lineage}")
            else:
                print(f"  Unexpected format: {line}")
    
    if result.stderr:
        print(f"Errors: {result.stderr}")
    
    return result

def test_simple_processing():
    """Test simple data processing"""
    print("\n📊 Testing simple data processing...")
    
    # Create test data
    test_data = {
        'taxon_name': ['Opisthokonta', 'Pythium', 'Agaricomycetes', 'Centraloplasthelida'],
        'member_size': [100, 50, 75, 25],
        'occurrence_count': [10, 5, 8, 3]
    }
    
    df = pd.DataFrame(test_data)
    print("Test DataFrame:")
    print(df)
    
    # Test name to taxid mapping
    name_to_taxid = {}
    
    result = subprocess.run(
        ["taxonkit", "name2taxid"],
        input="\n".join(df['taxon_name'].tolist()),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env
    )
    
    if result.returncode == 0 and result.stdout:
        lines = result.stdout.strip().split('\n')
        for i, line in enumerate(lines):
            if i < len(df):
                parts = line.strip().split('\t')
                name = df.iloc[i]['taxon_name']
                if len(parts) >= 2 and parts[1] != "0":
                    name_to_taxid[name] = parts[1]
                else:
                    name_to_taxid[name] = 'NA'
    
    # Add taxids to dataframe
    df['taxid'] = df['taxon_name'].map(name_to_taxid).fillna('NA')
    print("\nWith taxids:")
    print(df)
    
    # Test lineage generation for valid taxids
    valid_taxids = df[df['taxid'] != 'NA']['taxid'].tolist()
    
    if valid_taxids:
        result = subprocess.run(
            ["taxonkit", "lineage", "--show-name"],
            input="\n".join(valid_taxids),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        taxid_to_lineage = {}
        if result.returncode == 0 and result.stdout:
            lines = result.stdout.strip().split('\n')
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    taxid = parts[0]
                    lineage = parts[2]
                    taxid_to_lineage[taxid] = lineage
        
        # Add lineages to dataframe
        df['lineage'] = df['taxid'].map(taxid_to_lineage).fillna('')
        df['match_type'] = df['taxid'].apply(lambda x: 'full' if x != 'NA' else 'none')
        
        print("\nFinal result:")
        print(df[['taxon_name', 'taxid', 'match_type', 'lineage', 'member_size', 'occurrence_count']])

if __name__ == "__main__":
    print("🧪 Testing taxonkit functionality...")
    
    test_name_to_taxid()
    test_lineage_generation()
    test_simple_processing()
    
    print("\n✅ Test complete!")
