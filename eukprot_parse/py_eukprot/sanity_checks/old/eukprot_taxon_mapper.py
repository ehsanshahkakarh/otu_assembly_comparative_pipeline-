#!/usr/bin/env python3
"""
EukProt Taxon Mapper

This script creates a mapping between EukProt taxonomic names and NCBI taxids.
It uses taxonkit to perform the mapping and handles various taxonomic ranks.

Usage:
    python eukprot_taxon_mapper.py

The script will:
1. Load the EukProt taxonomic data
2. Map taxonomic names to NCBI taxids
3. Save the mapping to CSV files for later use
"""

import pandas as pd
import os
import subprocess
import logging
import re
import difflib
from pathlib import Path
import sys
from tqdm import tqdm
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("eukprot_taxon_mapper.log"),
        logging.StreamHandler()
    ]
)

# Constants
TAXONKIT_DB = "/clusterfs/jgi/groups/science/homes/ehsankakar/.taxonkit/"
EUKPROT_DATA = "Eukprot_included_datasets.txt"
OUTPUT_DIR = Path("../csv_eukprot/taxon_maps")
RANKS = ["phylum", "family", "genus"]

def ensure_taxonkit_available():
    """Check if taxonkit is available and the database exists."""
    try:
        # Check if taxonkit is in PATH
        result = subprocess.run(["which", "taxonkit"], 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, 
                               text=True)
        
        if result.returncode != 0:
            logging.error("taxonkit not found in PATH. Please install taxonkit.")
            sys.exit(1)
            
        # Check if database files exist
        if not os.path.exists(TAXONKIT_DB):
            logging.error(f"taxonkit database directory not found at {TAXONKIT_DB}")
            sys.exit(1)
            
        nodes_file = os.path.join(TAXONKIT_DB, "nodes.dmp")
        names_file = os.path.join(TAXONKIT_DB, "names.dmp")
        
        if not (os.path.exists(nodes_file) and os.path.exists(names_file)):
            logging.error(f"Required taxonkit database files not found in {TAXONKIT_DB}")
            sys.exit(1)
            
        logging.info("✅ taxonkit and database files are available")
        
    except Exception as e:
        logging.error(f"Error checking taxonkit: {e}")
        sys.exit(1)

def load_eukprot_data():
    """Load the EukProt dataset and extract taxonomic information."""
    try:
        df = pd.read_csv(EUKPROT_DATA, sep='\t')
        logging.info(f"✅ Loaded EukProt data with {len(df)} rows")
        return df
    except Exception as e:
        logging.error(f"Error loading EukProt data: {e}")
        sys.exit(1)

def extract_taxonomic_names(df, rank):
    """Extract taxonomic names for a specific rank from the EukProt dataset."""
    if rank == "phylum":
        # Extract phylum from Taxogroup2_UniEuk
        taxa = df["Taxogroup2_UniEuk"].dropna().unique()
    elif rank == "family":
        # Extract family from Taxonomy_UniEuk (index 4 if present)
        taxa = []
        for taxonomy in df["Taxonomy_UniEuk"].dropna():
            parts = taxonomy.split(';')
            if len(parts) > 4:
                family = parts[4].replace("'", "").strip()
                if family and not family.startswith("'") and not family.startswith('"'):
                    taxa.append(family)
        taxa = list(set(taxa))
    elif rank == "genus":
        # Extract genus from Genus_UniEuk
        taxa = df["Genus_UniEuk"].dropna().unique()
    else:
        logging.error(f"Unsupported rank: {rank}")
        return []
    
    # Clean up taxa names
    taxa = [t.strip().replace("'", "") for t in taxa if isinstance(t, str) and t.strip()]
    logging.info(f"✅ Extracted {len(taxa)} unique {rank} names")
    return sorted(taxa)

def get_taxid_with_variants(taxon_name, env=None):
    """
    Try to get taxid using different variants of the name.
    Returns (taxid, matched_name) tuple.
    """
    if env is None:
        env = os.environ.copy()
        env["TAXONKIT_DB"] = TAXONKIT_DB
    
    # Try different variants of the name
    variants = [
        taxon_name,  # Original
        taxon_name.lower(),  # Lowercase
        taxon_name.capitalize(),  # Capitalized
        re.sub(r'[-_\s]+', ' ', taxon_name),  # Replace hyphens/underscores with spaces
        re.sub(r'[-_\s]+', '', taxon_name),  # Remove hyphens/underscores/spaces
    ]
    
    # Add variants with common suffix transformations for taxonomic names
    if taxon_name.endswith('ae'):
        variants.append(taxon_name[:-2] + 'a')  # e.g., Chromulinaceae -> Chromulina
    if taxon_name.endswith('ales'):
        variants.append(taxon_name[:-4] + 'a')  # e.g., Chromulinales -> Chromulina
    if taxon_name.endswith('ida'):
        variants.append(taxon_name[:-3] + 'a')  # e.g., Chromulinida -> Chromulina
    
    for variant in variants:
        try:
            # Run taxonkit name2taxid
            result = subprocess.run(
                ["taxonkit", "name2taxid"],
                input=variant,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 2 and parts[1].strip() and parts[1] != "0":
                    return parts[1].strip(), variant
        except Exception as e:
            logging.debug(f"Error with variant '{variant}': {e}")
    
    return None, None

def map_taxa_to_taxids(taxa, rank):
    """Map taxonomic names to NCBI taxids using taxonkit."""
    env = os.environ.copy()
    env["TAXONKIT_DB"] = TAXONKIT_DB
    
    results = []
    not_found = []
    
    for taxon in tqdm(taxa, desc=f"Mapping {rank} names to taxids"):
        taxid, matched_variant = get_taxid_with_variants(taxon, env)
        
        if taxid:
            # Verify the rank using taxonkit lineage
            try:
                lineage_result = subprocess.run(
                    ["taxonkit", "lineage", "--show-rank", taxid],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env
                )
                
                if lineage_result.returncode == 0 and lineage_result.stdout.strip():
                    lineage_parts = lineage_result.stdout.strip().split('\t')
                    if len(lineage_parts) >= 3:
                        taxid_rank = lineage_parts[2].strip()
                        lineage = lineage_parts[1].strip()
                        
                        # Check if this is a eukaryotic taxon
                        is_eukaryotic = "Eukaryota" in lineage
                        
                        results.append({
                            "taxon_name": taxon,
                            "matched_variant": matched_variant,
                            "taxid": taxid,
                            "ncbi_rank": taxid_rank,
                            "is_eukaryotic": is_eukaryotic,
                            "lineage": lineage
                        })
                        continue
            except Exception as e:
                logging.debug(f"Error verifying rank for {taxon}: {e}")
        
        # If we get here, the taxon wasn't found or had issues
        not_found.append(taxon)
    
    logging.info(f"✅ Successfully mapped {len(results)} out of {len(taxa)} {rank} names")
    if not_found:
        logging.warning(f"❌ Could not map {len(not_found)} {rank} names")
        
    return results, not_found

def fuzzy_match_remaining(not_found, rank):
    """Try to fuzzy match remaining taxa that weren't found."""
    # First, get a comprehensive list of taxonomic names from NCBI
    env = os.environ.copy()
    env["TAXONKIT_DB"] = TAXONKIT_DB
    
    try:
        # Get all names at the specified rank
        logging.info(f"Retrieving all {rank} names from NCBI taxonomy...")
        
        # Use taxonkit list to get all taxids at the specified rank
        list_result = subprocess.run(
            ["taxonkit", "list", "--rank", rank],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        if list_result.returncode != 0:
            logging.error(f"Error retrieving {rank} names: {list_result.stderr}")
            return []
        
        # Extract taxids
        taxids = [line.strip() for line in list_result.stdout.strip().split('\n')]
        
        # Get names for these taxids
        names_result = subprocess.run(
            ["taxonkit", "name"],
            input='\n'.join(taxids),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        if names_result.returncode != 0:
            logging.error(f"Error retrieving names: {names_result.stderr}")
            return []
        
        # Create a dictionary of name -> taxid
        name_to_taxid = {}
        for line in names_result.stdout.strip().split('\n'):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name_to_taxid[parts[1].lower()] = parts[0]
        
        # Now try fuzzy matching
        fuzzy_matches = []
        for taxon in tqdm(not_found, desc=f"Fuzzy matching {rank} names"):
            best_match = None
            best_score = 0
            
            for ncbi_name in name_to_taxid.keys():
                score = difflib.SequenceMatcher(None, taxon.lower(), ncbi_name).ratio()
                if score > 0.9 and score > best_score:  # 90% similarity threshold
                    best_match = ncbi_name
                    best_score = score
            
            if best_match:
                taxid = name_to_taxid[best_match]
                
                # Verify it's eukaryotic
                lineage_result = subprocess.run(
                    ["taxonkit", "lineage", taxid],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env
                )
                
                if lineage_result.returncode == 0 and "Eukaryota" in lineage_result.stdout:
                    fuzzy_matches.append({
                        "taxon_name": taxon,
                        "matched_variant": best_match,
                        "taxid": taxid,
                        "ncbi_rank": rank,
                        "is_eukaryotic": True,
                        "lineage": lineage_result.stdout.strip().split('\t')[1],
                        "match_score": best_score,
                        "is_fuzzy_match": True
                    })
        
        logging.info(f"✅ Found {len(fuzzy_matches)} fuzzy matches for {rank} names")
        return fuzzy_matches
        
    except Exception as e:
        logging.error(f"Error during fuzzy matching: {e}")
        return []

def main():
    """Main function to run the EukProt taxon mapping process."""
    start_time = time.time()
    
    # Ensure taxonkit is available
    ensure_taxonkit_available()
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load EukProt data
    df = load_eukprot_data()
    
    # Process each taxonomic rank
    for rank in RANKS:
        logging.info(f"Processing {rank} rank...")
        
        # Extract taxonomic names
        taxa = extract_taxonomic_names(df, rank)
        
        # Map to taxids
        mapped_taxa, not_found = map_taxa_to_taxids(taxa, rank)
        
        # Try fuzzy matching for remaining taxa
        fuzzy_matches = fuzzy_match_remaining(not_found, rank)
        
        # Combine results
        all_results = mapped_taxa + fuzzy_matches
        
        # Save results
        if all_results:
            result_df = pd.DataFrame(all_results)
            output_file = OUTPUT_DIR / f"eukprot_{rank}_taxid_map.csv"
            result_df.to_csv(output_file, index=False)
            logging.info(f"✅ Saved {len(result_df)} {rank} mappings to {output_file}")
            
            # Also save not found taxa
            if not_found:
                not_found_df = pd.DataFrame({"taxon_name": not_found})
                not_found_file = OUTPUT_DIR / f"eukprot_{rank}_not_found.csv"
                not_found_df.to_csv(not_found_file, index=False)
                logging.info(f"✅ Saved {len(not_found)} unmapped {rank} names to {not_found_file}")
    
    elapsed_time = time.time() - start_time
    logging.info(f"✅ Completed EukProt taxon mapping in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
