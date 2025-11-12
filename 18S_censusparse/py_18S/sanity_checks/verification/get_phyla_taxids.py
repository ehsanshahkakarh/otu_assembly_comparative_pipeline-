#!/usr/bin/env python3
"""
Get NCBI Taxids for Phyla Names

This script reads phyla names from phyla_names_for_taxonkit.txt,
uses taxonkit to get their NCBI taxids, and outputs the taxids to taxid.txt.

Usage:
    python get_phyla_taxids.py

Input:
    phyla_names_for_taxonkit.txt: A text file with one phylum name per line

Output:
    taxid.txt: A text file with one taxid per line
"""

import os
import subprocess
from tqdm import tqdm

def get_taxid_for_name(taxon_name, env):
    """
    Get NCBI taxid for a taxon name using taxonkit.
    
    Args:
        taxon_name: The taxon name to look up
        env: Environment variables for subprocess
        
    Returns:
        The taxid as a string, or "NA" if not found
    """
    try:
        # Run taxonkit name2taxid
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=taxon_name,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0":
                return parts[1]
    except Exception as e:
        print(f"Error getting taxid for {taxon_name}: {e}")
    
    return "NA"

def get_phyla_taxids():
    """
    Get taxids for phyla names and save to a text file.
    """
    # Input file
    input_file = "phyla_names_for_taxonkit.txt"
    
    # Output file
    output_file = "taxid.txt"
    
    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        return
    
    try:
        # Read phyla names from input file
        with open(input_file, 'r') as f:
            phyla_names = [line.strip() for line in f if line.strip()]
        
        print(f"Read {len(phyla_names)} phyla names from {input_file}")
        
        # Get taxids for each phylum
        taxids = []
        found_count = 0
        not_found_count = 0
        
        print("Getting taxids for phyla names...")
        for name in tqdm(phyla_names):
            taxid = get_taxid_for_name(name, env)
            taxids.append(taxid)
            
            if taxid != "NA":
                found_count += 1
            else:
                not_found_count += 1
        
        # Write taxids to output file
        with open(output_file, 'w') as f:
            for taxid in taxids:
                if taxid != "NA":  # Only write valid taxids
                    f.write(f"{taxid}\n")
        
        print(f"Successfully wrote {found_count} taxids to {output_file}")
        print(f"Could not find taxids for {not_found_count} phyla names")
        
        # Print a summary of the results
        print("\nSummary:")
        print(f"Total phyla names: {len(phyla_names)}")
        print(f"Taxids found: {found_count} ({found_count/len(phyla_names)*100:.1f}%)")
        print(f"Taxids not found: {not_found_count} ({not_found_count/len(phyla_names)*100:.1f}%)")
        
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    try:
        from tqdm import tqdm
    except ImportError:
        # Define a simple replacement if tqdm is not available
        def tqdm(iterable, **kwargs):
            print("Note: tqdm not available, using simple progress reporting")
            total = len(iterable)
            for i, item in enumerate(iterable):
                if i % 10 == 0 or i == total - 1:
                    print(f"Progress: {i+1}/{total} ({(i+1)/total*100:.1f}%)")
                yield item
    
    get_phyla_taxids()
