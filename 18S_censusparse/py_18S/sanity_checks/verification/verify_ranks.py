#!/usr/bin/env python3
"""
Verify Taxonomic Ranks

This script verifies taxonomic rank information for taxids using taxonkit.
It reads taxids from a file and checks their taxonomic ranks and lineages.

Usage:
    python verify_ranks.py [input_file]

Input:
    input_file: A text file with one taxid per line (default: taxid.txt)

Output:
    rank_verification.txt: A tab-separated file with taxid, lineage, and rank verification
"""

import os
import subprocess
import tempfile
import sys

def verify_ranks(input_file="taxid.txt"):
    """
    Verify taxonomic ranks for taxids using taxonkit.
    
    Args:
        input_file: Path to file containing taxids (one per line)
    """
    # Output file
    output_file = "rank_verification.txt"
    
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
        
        if not taxids:
            print("No valid taxids found in input file")
            return
        
        # Create a temporary file with all taxids
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_path = temp_file.name
            for taxid in taxids:
                temp_file.write(f"{taxid}\n")
        
        try:
            # Run taxonkit lineage with rank information
            print("Running taxonkit to verify ranks and lineages...")
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
            verified_count = 0
            failed_count = 0
            
            with open(output_file, 'w') as f:
                # Write header
                f.write("taxid\tlineage\trank\tstatus\n")
                
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
                        status = "verified"
                        verified_count += 1
                        f.write(f"{taxid}\t{lineage}\t{rank}\t{status}\n")
                    elif len(parts) == 2:
                        taxid = parts[0]
                        lineage = parts[1]
                        rank = "unknown"
                        status = "partial"
                        verified_count += 1
                        f.write(f"{taxid}\t{lineage}\t{rank}\t{status}\n")
                    else:
                        # Handle failed lookups
                        if parts:
                            taxid = parts[0]
                            lineage = "failed"
                            rank = "failed"
                            status = "failed"
                            failed_count += 1
                            f.write(f"{taxid}\t{lineage}\t{rank}\t{status}\n")
            
            # Add any taxids that weren't found in the output
            output_taxids = set()
            with open(output_file, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    if line.strip():
                        output_taxids.add(line.split('\t')[0])
            
            missing_taxids = set(taxids) - output_taxids
            if missing_taxids:
                with open(output_file, 'a') as f:
                    for taxid in missing_taxids:
                        f.write(f"{taxid}\tnot_found\tnot_found\tfailed\n")
                        failed_count += 1
            
            print(f"Successfully wrote rank verification to {output_file}")
            print(f"Verified: {verified_count}, Failed: {failed_count}")
            
            # Show a summary of the results
            print("\nSummary of rank verification:")
            rank_counts = {}
            status_counts = {}
            
            with open(output_file, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            rank = parts[2]
                            status = parts[3]
                            rank_counts[rank] = rank_counts.get(rank, 0) + 1
                            status_counts[status] = status_counts.get(status, 0) + 1
            
            print("\nRank distribution:")
            for rank, count in sorted(rank_counts.items()):
                print(f"  {rank}: {count}")
            
            print("\nStatus distribution:")
            for status, count in sorted(status_counts.items()):
                print(f"  {status}: {count}")
            
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

def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "taxid.txt"
    
    verify_ranks(input_file)

if __name__ == "__main__":
    main()
