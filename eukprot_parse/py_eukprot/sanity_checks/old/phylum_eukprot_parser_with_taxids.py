#!/usr/bin/env python3
"""
EukProt Phylum Parser with NCBI Taxids

This script processes the EukProt dataset to extract phylum information
and adds NCBI taxids to the output files using the taxon mapping.

Usage:
    python phylum_eukprot_parser_with_taxids.py

The script will:
1. Load the EukProt dataset
2. Extract phylum information
3. Add NCBI taxids from the mapping file
4. Save the results to CSV files
"""

import pandas as pd
import os
import logging
from pathlib import Path
import sys
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("phylum_eukprot_parser_with_taxids.log"),
        logging.StreamHandler()
    ]
)

# Get absolute path to the current script
script_dir = Path(__file__).resolve().parent

# Define paths
csv_output_dir = script_dir.parent / "csv_eukprot"
input_file = script_dir / "Eukprot_included_datasets.txt"
output_file = csv_output_dir / "eukprot_phylum_counts.csv"
detailed_output_file = csv_output_dir / "eukprot_phylum_with_accessions.csv"
taxid_map_file = csv_output_dir / "taxon_maps" / "eukprot_phylum_taxid_map.csv"

# Ensure output directory exists
csv_output_dir.mkdir(exist_ok=True)
(csv_output_dir / "taxon_maps").mkdir(exist_ok=True)

def load_taxid_map():
    """Load the phylum to taxid mapping file if it exists."""
    try:
        if os.path.exists(taxid_map_file):
            taxid_map = pd.read_csv(taxid_map_file)
            logging.info(f"✅ Loaded taxid map with {len(taxid_map)} phyla")
            return taxid_map
        else:
            logging.warning(f"⚠️ Taxid map file not found at {taxid_map_file}")
            logging.warning("⚠️ Run eukprot_taxon_mapper.py first to generate the mapping")
            return None
    except Exception as e:
        logging.error(f"❌ Error loading taxid map: {e}")
        return None

def process_eukprot_data():
    """Process the EukProt dataset file to extract phylum information with taxids."""
    try:
        # Load EukProt dataset file
        df = pd.read_csv(input_file, sep='\t')
        logging.info(f"✅ Successfully loaded EukProt file with {len(df)} rows")

        # Extract phylum directly from Taxogroup2_UniEuk
        df["phylum"] = df["Taxogroup2_UniEuk"].astype(str).str.strip().str.replace("'", "")
        
        # Add domain column (all are Eukaryota)
        df["domain"] = "Eukaryota"
        
        # Create accession column from EukProt_ID
        df["accession_clean"] = df["EukProt_ID"]

        # Load taxid mapping if available
        taxid_map = load_taxid_map()
        
        # Add taxids if mapping is available
        if taxid_map is not None:
            # Create a dictionary for faster lookups
            taxid_dict = dict(zip(taxid_map["taxon_name"], taxid_map["taxid"]))
            
            # Add taxid column to the dataframe
            df["taxid"] = df["phylum"].map(taxid_dict)
            
            # Count how many phyla have taxids
            phyla_with_taxids = df["taxid"].notna().sum()
            logging.info(f"✅ Added taxids to {phyla_with_taxids} out of {df['phylum'].notna().sum()} phyla")
            
            # Add NCBI rank information if available
            if "ncbi_rank" in taxid_map.columns:
                ncbi_rank_dict = dict(zip(taxid_map["taxon_name"], taxid_map["ncbi_rank"]))
                df["ncbi_rank"] = df["phylum"].map(ncbi_rank_dict)
        
        # Save detailed file with all accessions
        detailed_df = df.dropna(subset=["phylum"])
        
        # Select columns for detailed output
        if taxid_map is not None:
            detailed_df[["accession_clean", "phylum", "domain", "taxid", "ncbi_rank"]].to_csv(detailed_output_file, index=False)
        else:
            detailed_df[["accession_clean", "phylum", "domain"]].to_csv(detailed_output_file, index=False)
            
        logging.info(f"✅ Saved detailed file with {len(detailed_df)} entries to {detailed_output_file}")
        
        # Drop rows with null phylum values
        df = df.dropna(subset=["phylum"])
        logging.info(f"✅ After dropping null phylum values: {len(df)} rows")

        # Group by phylum and domain, then count
        if taxid_map is not None:
            counts = df.groupby(['phylum', 'domain', 'taxid', 'ncbi_rank']).size().reset_index(name='eukprot_genome_count')
        else:
            counts = df.groupby(['phylum', 'domain']).size().reset_index(name='eukprot_genome_count')
            
        logging.info(f"✅ Found {len(counts)} unique phylum combinations")

        # Add accession lists to the counts dataframe
        if taxid_map is not None:
            accession_lists = df.groupby(['phylum', 'domain', 'taxid', 'ncbi_rank'])["accession_clean"].apply(list).reset_index()
            counts = pd.merge(counts, accession_lists, on=['phylum', 'domain', 'taxid', 'ncbi_rank'])
        else:
            accession_lists = df.groupby(['phylum', 'domain'])["accession_clean"].apply(list).reset_index()
            counts = pd.merge(counts, accession_lists, on=['phylum', 'domain'])

        # Sort by count (descending)
        counts = counts.sort_values('eukprot_genome_count', ascending=False)

        # Save
        counts.to_csv(output_file, index=False)
        logging.info(f"✅ Saved phylum counts to {output_file}")
        
        # Print preview
        logging.info("\n🔝 Top 10 phyla:")
        if taxid_map is not None:
            logging.info(counts[["phylum", "domain", "taxid", "eukprot_genome_count"]].head(10))
        else:
            logging.info(counts[["phylum", "domain", "eukprot_genome_count"]].head(10))
        
        # Generate verification files for taxonomic validation
        generate_verification_files(counts)
        
    except Exception as e:
        logging.error(f"❌ Error: {e}")

def generate_verification_files(counts_df):
    """Generate files for taxonomic verification."""
    try:
        # Create directory for verification files
        verify_dir = csv_output_dir / "sanity_check" / "taxon_names"
        verify_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a file with just the phylum names and taxids for verification
        verify_file = verify_dir / "eukaryotic_phyla.csv"
        
        if "taxid" in counts_df.columns:
            # If we have taxids, include them
            phyla_for_verify = counts_df[["phylum", "taxid"]].drop_duplicates()
            phyla_for_verify.columns = ["taxon_name", "taxid"]
        else:
            # Otherwise just use the names
            phyla_for_verify = pd.DataFrame({"taxon_name": counts_df["phylum"].unique()})
        
        phyla_for_verify.to_csv(verify_file, index=False)
        logging.info(f"✅ Generated verification file with {len(phyla_for_verify)} phyla at {verify_file}")
        
    except Exception as e:
        logging.error(f"❌ Error generating verification files: {e}")

if __name__ == "__main__":
    process_eukprot_data()
