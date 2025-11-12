#!/usr/bin/env python3
"""
Get Lineage and Rank Information for Taxids

This script reads taxids from taxid.txt and uses taxonkit
to get the lineage and rank information for each taxid.

Usage:
    python get_lineage_ranks.py

Input:
    taxid.txt: A text file with one taxid per line

Output:
    lineage_ranks.txt: A tab-separated file with taxid, lineage, and rank columns
"""

import os
import subprocess
import tempfile

def get_lineage_ranks():
    """
    Get lineage and rank information for taxids using taxonkit.
    """
    # Input and output files
    input_file = "taxid.txt"
    output_file = "lineage_ranks.txt"
    
    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        return
    
    try:
        # Read taxids from input file
        with open(input_file, 'r') as f:
            taxids = [line.strip() for line in f if line.strip()]
        
        print(f"Read {len(taxids)} taxids from {input_file}")
        
        # Create a temporary file with all taxids
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_path = temp_file.name
            for taxid in taxids:
                temp_file.write(f"{taxid}\n")
        
        try:
            # Run taxonkit lineage with rank information
            print("Running taxonkit to get lineage and rank information...")
            result = subprocess.run(
                ["taxonkit", "lineage", "-r", temp_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )
            
            if result.returncode != 0:
                print(f"Error running taxonkit: {result.stderr}")
                return
            
            # Parse the output and write to file
            with open(output_file, 'w') as f:
                # Write header
                f.write("taxid\tlineage\trank\n")
                
                # Process each line of output
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        taxid = parts[0]
                        lineage = parts[1]
                        rank = parts[2]
                        f.write(f"{taxid}\t{lineage}\t{rank}\n")
                    elif len(parts) == 2:
                        # If no rank information, use empty string
                        taxid = parts[0]
                        lineage = parts[1]
                        rank = ""
                        f.write(f"{taxid}\t{lineage}\t{rank}\n")
            
            print(f"Successfully wrote lineage and rank information to {output_file}")
            
            # Show a preview of the results
            print("\nPreview of results:")
            with open(output_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines[:6]):  # Show header + first 5 entries
                    print(line.strip())
                if len(lines) > 6:
                    print(f"... and {len(lines) - 6} more entries")
        
        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as e:
                print(f"Warning: Could not remove temporary file: {e}")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    get_lineage_ranks()
