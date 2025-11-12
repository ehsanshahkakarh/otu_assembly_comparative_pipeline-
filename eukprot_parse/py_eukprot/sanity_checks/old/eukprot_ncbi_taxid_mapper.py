#!/usr/bin/env python3
"""
EukProt NCBI Taxid Mapper

This script maps EukProt taxonomic names (genus, family, phylum) to NCBI taxids
using the NCBI taxdump files. It employs fuzzy matching to handle naming differences
while respecting the lineage information in the EukProt metadata.

Usage:
    python eukprot_ncbi_taxid_mapper.py

Output:
    CSV files with taxon name to taxid mappings for each taxonomic rank
"""

import pandas as pd
import os
import re
import sys
import logging
from pathlib import Path
import difflib
from tqdm import tqdm
import subprocess

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("eukprot_ncbi_taxid_mapper.log"),
        logging.StreamHandler()
    ]
)

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
EUKPROT_DATA = SCRIPT_DIR / "Eukprot_included_datasets.txt"
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")
OUTPUT_DIR = SCRIPT_DIR.parent / "csv_eukprot" / "taxon_maps"
VERIFICATION_DIR = SCRIPT_DIR.parent / "csv_eukprot" / "sanity_check" / "taxon_names"

# Ensure output directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
VERIFICATION_DIR.mkdir(parents=True, exist_ok=True)

# Taxonomic ranks to process
RANKS = ["phylum", "family", "genus"]

def load_eukprot_data():
    """Load the EukProt dataset."""
    try:
        df = pd.read_csv(EUKPROT_DATA, sep='\t')
        logging.info(f"✅ Loaded EukProt data with {len(df)} rows")
        return df
    except Exception as e:
        logging.error(f"❌ Error loading EukProt data: {e}")
        sys.exit(1)

def load_ncbi_names():
    """Load NCBI names.dmp file into a DataFrame."""
    try:
        names_file = TAXDUMP_DIR / "names.dmp"
        
        if not names_file.exists():
            logging.error(f"❌ NCBI names.dmp file not found at {names_file}")
            sys.exit(1)
            
        logging.info(f"Loading NCBI names from {names_file}...")
        
        # Read names.dmp with custom delimiter
        names_df = pd.read_csv(
            names_file, 
            sep='|', 
            header=None,
            names=['taxid', 'name', 'unique_name', 'name_class'],
            engine='python'
        )
        
        # Clean up whitespace
        names_df = names_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        
        # Filter to scientific names only
        scientific_names = names_df[names_df['name_class'] == 'scientific name']
        
        logging.info(f"✅ Loaded {len(scientific_names)} scientific names from NCBI taxonomy")
        return scientific_names
    
    except Exception as e:
        logging.error(f"❌ Error loading NCBI names: {e}")
        sys.exit(1)

def load_ncbi_nodes():
    """Load NCBI nodes.dmp file into a DataFrame."""
    try:
        nodes_file = TAXDUMP_DIR / "nodes.dmp"
        
        if not nodes_file.exists():
            logging.error(f"❌ NCBI nodes.dmp file not found at {nodes_file}")
            sys.exit(1)
            
        logging.info(f"Loading NCBI nodes from {nodes_file}...")
        
        # Read nodes.dmp with custom delimiter
        nodes_df = pd.read_csv(
            nodes_file, 
            sep='|', 
            header=None,
            names=['taxid', 'parent_taxid', 'rank', 'embl_code', 'division_id', 
                   'inherited_div_flag', 'genetic_code_id', 'inherited_gc_flag',
                   'mitochondrial_genetic_code_id', 'inherited_mgc_flag',
                   'genbank_hidden_flag', 'hidden_subtree_root_flag', 
                   'comments', 'plastid_genetic_code_id', 'inherited_pgc_flag',
                   'specified_species', 'inherited_specified_species',
                   'genbank_hidden_flag2', 'hidden_subtree_root_flag2'],
            engine='python'
        )
        
        # Clean up whitespace
        nodes_df = nodes_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        
        logging.info(f"✅ Loaded {len(nodes_df)} nodes from NCBI taxonomy")
        return nodes_df
    
    except Exception as e:
        logging.error(f"❌ Error loading NCBI nodes: {e}")
        sys.exit(1)

def extract_eukprot_taxa(df, rank):
    """Extract taxonomic names for a specific rank from the EukProt dataset."""
    if rank == "phylum":
        # Extract phylum from Taxogroup2_UniEuk
        taxa = df["Taxogroup2_UniEuk"].dropna().astype(str)
        taxa = taxa[taxa != "nan"].unique()
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
        taxa = df["Genus_UniEuk"].dropna().astype(str)
        taxa = taxa[taxa != "nan"].unique()
    else:
        logging.error(f"❌ Unsupported rank: {rank}")
        return []
    
    # Clean up taxa names
    taxa = [t.strip().replace("'", "") for t in taxa if isinstance(t, str) and t.strip()]
    logging.info(f"✅ Extracted {len(taxa)} unique {rank} names from EukProt")
    return sorted(taxa)

def get_lineage_from_taxid(taxid, nodes_df, names_df):
    """Get the full lineage for a taxid."""
    lineage = []
    current_taxid = taxid
    
    # Prevent infinite loops
    visited = set()
    
    while current_taxid != "1" and current_taxid not in visited:
        visited.add(current_taxid)
        
        # Get the name for this taxid
        name_row = names_df[names_df['taxid'] == current_taxid]
        if not name_row.empty:
            lineage.append(name_row.iloc[0]['name'])
        
        # Get the parent taxid
        parent_row = nodes_df[nodes_df['taxid'] == current_taxid]
        if parent_row.empty:
            break
            
        current_taxid = parent_row.iloc[0]['parent_taxid']
    
    # Reverse to get root->leaf order
    return lineage[::-1]

def is_eukaryotic(taxid, nodes_df, names_df):
    """Check if a taxid belongs to Eukaryota."""
    lineage = get_lineage_from_taxid(taxid, nodes_df, names_df)
    return "Eukaryota" in lineage

def get_rank(taxid, nodes_df):
    """Get the taxonomic rank for a taxid."""
    row = nodes_df[nodes_df['taxid'] == taxid]
    if row.empty:
        return None
    return row.iloc[0]['rank']

def fuzzy_match_taxa(eukprot_taxa, ncbi_names_df, nodes_df, rank):
    """
    Match EukProt taxa to NCBI taxids using fuzzy matching.
    Returns a DataFrame with the mapping results.
    """
    # Filter NCBI names by rank if possible
    if rank in ["phylum", "family", "genus"]:
        rank_nodes = nodes_df[nodes_df['rank'] == rank]
        rank_taxids = set(rank_nodes['taxid'])
        rank_names = ncbi_names_df[ncbi_names_df['taxid'].isin(rank_taxids)]
    else:
        rank_names = ncbi_names_df
    
    logging.info(f"Found {len(rank_names)} NCBI names at rank '{rank}'")
    
    # Create a dictionary of lowercase name -> original name for faster lookups
    name_map = {name.lower(): (taxid, name) for taxid, name in zip(rank_names['taxid'], rank_names['name'])}
    
    results = []
    for taxon in tqdm(eukprot_taxa, desc=f"Matching {rank} names"):
        # Try exact match first (case insensitive)
        taxon_lower = taxon.lower()
        if taxon_lower in name_map:
            taxid, ncbi_name = name_map[taxon_lower]
            if is_eukaryotic(taxid, nodes_df, ncbi_names_df):
                results.append({
                    "taxon_name": taxon,
                    "taxid": taxid,
                    "ncbi_name": ncbi_name,
                    "match_type": "exact",
                    "match_score": 1.0,
                    "ncbi_rank": get_rank(taxid, nodes_df)
                })
                continue
        
        # Try various transformations
        variants = [
            taxon_lower,
            taxon_lower.replace("-", ""),
            taxon_lower.replace("_", ""),
            taxon_lower.replace(" ", ""),
            re.sub(r'[^a-z0-9]', '', taxon_lower)
        ]
        
        # Add variants with common suffix transformations
        if taxon_lower.endswith('ae'):
            variants.append(taxon_lower[:-2] + 'a')
        if taxon_lower.endswith('ales'):
            variants.append(taxon_lower[:-4] + 'a')
        if taxon_lower.endswith('ida'):
            variants.append(taxon_lower[:-3] + 'a')
        
        matched = False
        for variant in variants:
            if variant in name_map:
                taxid, ncbi_name = name_map[variant]
                if is_eukaryotic(taxid, nodes_df, ncbi_names_df):
                    results.append({
                        "taxon_name": taxon,
                        "taxid": taxid,
                        "ncbi_name": ncbi_name,
                        "match_type": "variant",
                        "match_score": 0.9,
                        "ncbi_rank": get_rank(taxid, nodes_df)
                    })
                    matched = True
                    break
        
        if matched:
            continue
            
        # If still no match, try fuzzy matching
        best_match = None
        best_score = 0
        best_taxid = None
        
        # Only check against names that are at least 3 characters long
        if len(taxon_lower) >= 3:
            for ncbi_lower, (taxid, ncbi_name) in name_map.items():
                # Skip very short names to avoid false matches
                if len(ncbi_lower) < 3:
                    continue
                    
                # Calculate similarity score
                score = difflib.SequenceMatcher(None, taxon_lower, ncbi_lower).ratio()
                
                # Use a threshold based on name length
                threshold = 0.85 if len(taxon_lower) > 5 else 0.9
                
                if score > threshold and score > best_score:
                    # Verify it's eukaryotic
                    if is_eukaryotic(taxid, nodes_df, ncbi_names_df):
                        best_match = ncbi_name
                        best_score = score
                        best_taxid = taxid
        
        if best_match:
            results.append({
                "taxon_name": taxon,
                "taxid": best_taxid,
                "ncbi_name": best_match,
                "match_type": "fuzzy",
                "match_score": best_score,
                "ncbi_rank": get_rank(best_taxid, nodes_df)
            })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Add unmapped taxa
    if not results_df.empty:
        mapped_taxa = set(results_df['taxon_name'])
        unmapped = [t for t in eukprot_taxa if t not in mapped_taxa]
        
        for taxon in unmapped:
            results_df = pd.concat([results_df, pd.DataFrame([{
                "taxon_name": taxon,
                "taxid": None,
                "ncbi_name": None,
                "match_type": "unmapped",
                "match_score": 0.0,
                "ncbi_rank": None
            }])], ignore_index=True)
    else:
        # If no matches at all, create DataFrame with all unmapped
        data = []
        for taxon in eukprot_taxa:
            data.append({
                "taxon_name": taxon,
                "taxid": None,
                "ncbi_name": None,
                "match_type": "unmapped",
                "match_score": 0.0,
                "ncbi_rank": None
            })
        results_df = pd.DataFrame(data)
    
    # Log results
    mapped_count = len(results_df[results_df['match_type'] != 'unmapped'])
    logging.info(f"✅ Mapped {mapped_count} out of {len(eukprot_taxa)} {rank} names ({mapped_count/len(eukprot_taxa)*100:.1f}%)")
    
    return results_df

def create_verification_files(mapping_df, rank):
    """Create verification files for the taxonomic names."""
    try:
        # Create the verification file with taxon_name and taxid columns
        verify_file = VERIFICATION_DIR / f"eukaryotic_{rank}s.csv"
        
        # Select only the mapped taxa (with taxids)
        mapped_df = mapping_df[mapping_df['taxid'].notna()].copy()
        
        if not mapped_df.empty:
            # Create the verification file
            verify_df = mapped_df[['taxon_name', 'taxid']].copy()
            verify_df['taxid'] = verify_df['taxid'].astype(str)
            verify_df.to_csv(verify_file, index=False)
            logging.info(f"✅ Created verification file with {len(verify_df)} {rank} names at {verify_file}")
        else:
            # Create an empty file with headers
            pd.DataFrame(columns=['taxon_name', 'taxid']).to_csv(verify_file, index=False)
            logging.warning(f"⚠️ Created empty verification file at {verify_file} (no mapped taxa)")
    
    except Exception as e:
        logging.error(f"❌ Error creating verification file: {e}")

def main():
    """Main function to run the EukProt NCBI taxid mapping process."""
    try:
        # Load EukProt data
        eukprot_df = load_eukprot_data()
        
        # Load NCBI taxonomy data
        ncbi_names_df = load_ncbi_names()
        ncbi_nodes_df = load_ncbi_nodes()
        
        # Process each taxonomic rank
        for rank in RANKS:
            logging.info(f"Processing {rank} rank...")
            
            # Extract taxa from EukProt
            eukprot_taxa = extract_eukprot_taxa(eukprot_df, rank)
            
            if not eukprot_taxa:
                logging.warning(f"⚠️ No {rank} taxa found in EukProt data")
                continue
            
            # Match taxa to NCBI taxids
            mapping_df = fuzzy_match_taxa(eukprot_taxa, ncbi_names_df, ncbi_nodes_df, rank)
            
            # Save mapping results
            output_file = OUTPUT_DIR / f"eukprot_{rank}_taxid_map.csv"
            mapping_df.to_csv(output_file, index=False)
            logging.info(f"✅ Saved {rank} mapping to {output_file}")
            
            # Create verification files
            create_verification_files(mapping_df, rank)
        
        logging.info("✅ EukProt NCBI taxid mapping completed successfully")
        
    except Exception as e:
        logging.error(f"❌ Error in main function: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
