#!/usr/bin/env python3
"""
Extract Names with Taxonomy

This script extracts the Name_to_Use column along with taxonomic information
from the Eukprot_included_datasets.txt file and attempts to map them to NCBI taxids.

Usage:
    python extract_names_with_taxonomy.py
"""

import pandas as pd
import os
import subprocess
import logging
from pathlib import Path
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("extract_names_with_taxonomy.log"),
        logging.StreamHandler()
    ]
)

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
EUKPROT_DATA = SCRIPT_DIR / "Eukprot_included_datasets.txt"
OUTPUT_DIR = SCRIPT_DIR.parent / "csv_eukprot"
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")

# Ensure output directory exists
OUTPUT_DIR.mkdir(exist_ok=True)

def load_eukprot_data():
    """Load the EukProt dataset."""
    try:
        df = pd.read_csv(EUKPROT_DATA, sep='\t')
        logging.info(f"✅ Loaded EukProt data with {len(df)} rows")
        return df
    except Exception as e:
        logging.error(f"❌ Error loading EukProt data: {e}")
        sys.exit(1)

def extract_taxonomy_info(df):
    """
    Extract Name_to_Use along with taxonomic information.
    Returns a DataFrame with the extracted information.
    """
    # Create a new DataFrame to store the results
    result = pd.DataFrame()
    
    # Copy the Name_to_Use column
    result['organism_name'] = df['Name_to_Use']
    
    # Extract genus and species (if available)
    result['genus'] = df['Genus_UniEuk']
    result['species_epithet'] = df['Epithet_UniEuk']
    
    # Extract taxonomic ranks
    result['supergroup'] = df['Supergroup_UniEuk']
    result['taxogroup1'] = df['Taxogroup1_UniEuk']
    result['phylum'] = df['Taxogroup2_UniEuk']  # Using Taxogroup2 as phylum
    
    # Extract strain information
    result['strain'] = df['Strain']
    
    # Extract EukProt ID
    result['eukprot_id'] = df['EukProt_ID']
    
    # Extract full taxonomy string
    result['full_taxonomy'] = df['Taxonomy_UniEuk']
    
    # Clean up the data
    for col in result.columns:
        if result[col].dtype == 'object':
            result[col] = result[col].str.replace("'", "").str.strip()
    
    logging.info(f"✅ Extracted taxonomy information for {len(result)} organisms")
    return result

def try_get_taxid(name, env=None):
    """
    Try to get NCBI taxid for a name using taxonkit.
    Returns the taxid if found, None otherwise.
    """
    if not name or pd.isna(name):
        return None
        
    if env is None:
        env = os.environ.copy()
        env["TAXONKIT_DB"] = str(TAXDUMP_DIR)
    
    try:
        # Run taxonkit name2taxid
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=name,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1].strip() and parts[1] != "0":
                return parts[1].strip()
    except Exception:
        pass
    
    return None

def add_taxids(df):
    """
    Add NCBI taxids to the DataFrame.
    Tries to match at different taxonomic levels.
    """
    # Set up environment for taxonkit
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)
    
    # Check if taxonkit and database are available
    try:
        result = subprocess.run(
            ["taxonkit", "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        if result.returncode != 0:
            logging.warning("⚠️ taxonkit not available or database not found. Skipping taxid mapping.")
            return df
        logging.info(f"✅ Using taxonkit {result.stdout.strip()}")
    except Exception:
        logging.warning("⚠️ taxonkit not available. Skipping taxid mapping.")
        return df
    
    # Add taxid columns
    df['organism_taxid'] = None
    df['genus_taxid'] = None
    df['phylum_taxid'] = None
    
    # Try to get taxids for each row
    for idx, row in df.iterrows():
        # Try organism name first
        if not pd.isna(row['organism_name']):
            organism_taxid = try_get_taxid(row['organism_name'], env)
            if organism_taxid:
                df.at[idx, 'organism_taxid'] = organism_taxid
        
        # Try genus
        if not pd.isna(row['genus']):
            genus_taxid = try_get_taxid(row['genus'], env)
            if genus_taxid:
                df.at[idx, 'genus_taxid'] = genus_taxid
        
        # Try phylum
        if not pd.isna(row['phylum']):
            phylum_taxid = try_get_taxid(row['phylum'], env)
            if phylum_taxid:
                df.at[idx, 'phylum_taxid'] = phylum_taxid
    
    # Count how many taxids were found
    organism_count = df['organism_taxid'].notna().sum()
    genus_count = df['genus_taxid'].notna().sum()
    phylum_count = df['phylum_taxid'].notna().sum()
    
    logging.info(f"✅ Found taxids for {organism_count} organisms, {genus_count} genera, and {phylum_count} phyla")
    
    return df

def main():
    """Main function to extract names with taxonomy."""
    try:
        # Load EukProt data
        df = load_eukprot_data()
        
        # Extract taxonomy information
        taxonomy_df = extract_taxonomy_info(df)
        
        # Try to add taxids
        taxonomy_df = add_taxids(taxonomy_df)
        
        # Save to CSV
        output_file = OUTPUT_DIR / "eukprot_names_with_taxonomy.csv"
        taxonomy_df.to_csv(output_file, index=False)
        logging.info(f"✅ Saved taxonomy information to {output_file}")
        
        # Also save a simple list of names
        names_file = OUTPUT_DIR / "eukprot_names_to_use.txt"
        with open(names_file, 'w') as f:
            for name in taxonomy_df['organism_name'].dropna():
                f.write(f"{name}\n")
        logging.info(f"✅ Saved list of names to {names_file}")
        
    except Exception as e:
        logging.error(f"❌ Error in main function: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
