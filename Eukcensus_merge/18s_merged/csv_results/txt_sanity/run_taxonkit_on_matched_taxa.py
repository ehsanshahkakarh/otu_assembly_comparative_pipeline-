#!/usr/bin/env python3
"""
Run taxonkit commands on matched taxa from 18S NCBI files.
Gets taxids and lineages for all phyla, families, and genera.

Created: 2025-01-12
"""

import subprocess
import tempfile
import os
from datetime import datetime

def read_taxa_from_file(filename):
    """Read taxa names from a .txt file, skipping header lines."""
    taxa = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
        # Skip header lines (until we find the separator line)
        start_reading = False
        for line in lines:
            line = line.strip()
            if line.startswith('----'):
                start_reading = True
                continue
            
            if start_reading and line and not line.startswith('='):
                # Extract taxon name (first column)
                parts = line.split('\t')
                if len(parts) >= 3:  # Ensure we have all columns
                    taxon = parts[0].strip()
                    if taxon:
                        taxa.append(taxon)
    
    return taxa

def run_taxonkit_name2taxid(taxa_list):
    """Run taxonkit name2taxid on a list of taxa."""
    if not taxa_list:
        return {}
    
    name_to_taxid = {}
    
    try:
        # Create temporary file with taxa names
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            temp_filename = temp_file.name
            for taxon in taxa_list:
                temp_file.write(f"{taxon}\n")
        
        # Run taxonkit name2taxid
        cmd = ['taxonkit', 'name2taxid', temp_filename]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse results
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()
                    if taxid and taxid != 'NA':
                        name_to_taxid[name] = taxid
        
        # Clean up
        os.unlink(temp_filename)
        
    except Exception as e:
        print(f"Error running taxonkit name2taxid: {e}")
    
    return name_to_taxid

def run_taxonkit_lineage(taxids):
    """Run taxonkit lineage on a list of taxids."""
    if not taxids:
        return {}
    
    taxid_to_lineage = {}
    
    try:
        # Create temporary file with taxids
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            temp_filename = temp_file.name
            for taxid in taxids:
                temp_file.write(f"{taxid}\n")
        
        # Run taxonkit lineage
        cmd = ['taxonkit', 'lineage', temp_filename]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse results
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip() if len(parts) > 1 else ""
                    taxid_to_lineage[taxid] = lineage
        
        # Clean up
        os.unlink(temp_filename)
        
    except Exception as e:
        print(f"Error running taxonkit lineage: {e}")
    
    return taxid_to_lineage

def process_taxa_file(filename, output_suffix):
    """Process a single taxa file and generate taxonkit results."""
    print(f"Processing {filename}...")
    
    # Read taxa from file
    taxa = read_taxa_from_file(filename)
    print(f"  Found {len(taxa)} taxa")
    
    if not taxa:
        print(f"  No taxa found in {filename}")
        return
    
    # Get taxids
    print("  Getting taxids...")
    name_to_taxid = run_taxonkit_name2taxid(taxa)
    print(f"  Found taxids for {len(name_to_taxid)}/{len(taxa)} taxa")
    
    # Get lineages
    print("  Getting lineages...")
    taxids = list(name_to_taxid.values())
    taxid_to_lineage = run_taxonkit_lineage(taxids)
    print(f"  Found lineages for {len(taxid_to_lineage)}/{len(taxids)} taxids")
    
    # Write results
    output_filename = f"18s_ncbi_taxonkit_{output_suffix}.txt"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    with open(output_filename, 'w') as f:
        f.write(f"18S NCBI Taxonkit Results - {output_suffix.title()}\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"Total taxa processed: {len(taxa)}\n")
        f.write(f"Taxa with taxids: {len(name_to_taxid)}\n")
        f.write(f"Taxa with lineages: {len(taxid_to_lineage)}\n")
        f.write("=" * 80 + "\n\n")
        f.write("Taxon_Name\tTaxID\tLineage\n")
        f.write("-" * 80 + "\n")
        
        for taxon in taxa:
            taxid = name_to_taxid.get(taxon, "NA")
            lineage = taxid_to_lineage.get(taxid, "NA") if taxid != "NA" else "NA"
            f.write(f"{taxon}\t{taxid}\t{lineage}\n")
    
    print(f"  Results saved to: {output_filename}")
    return output_filename

def main():
    """Main function to process all taxa files."""
    print("Running taxonkit on matched 18S NCBI taxa...")
    print("=" * 60)
    
    # Define input files and their output suffixes
    files_to_process = [
        ("18s_ncbi_matched_phyla.txt", "phyla"),
        ("18s_ncbi_matched_families.txt", "families"),
        ("18s_ncbi_matched_genera.txt", "genera")
    ]
    
    output_files = []
    
    for input_file, suffix in files_to_process:
        if os.path.exists(input_file):
            output_file = process_taxa_file(input_file, suffix)
            if output_file:
                output_files.append(output_file)
            print()
        else:
            print(f"Warning: {input_file} not found, skipping...")
            print()
    
    print("=" * 60)
    print("Taxonkit processing complete!")
    print(f"Generated {len(output_files)} output files:")
    for output_file in output_files:
        print(f"  - {output_file}")

if __name__ == "__main__":
    main()
