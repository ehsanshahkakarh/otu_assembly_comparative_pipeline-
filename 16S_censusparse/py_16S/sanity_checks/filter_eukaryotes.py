#!/usr/bin/env python3
"""
Filter Eukaryotic Entries from EukCensus Genus CSV

This script filters the eukcensus_by_genus.csv file to remove any rows that contain
"Bacteria" or "Archaea" in them. It also checks if the taxon names are Eukaryotic
using taxonkit and removes any non-Eukaryotic entries.

Output file:
- eukcensus_by_genus_eukaryotes_only.csv
"""

import pandas as pd
import os
import subprocess
import tempfile
import csv
import multiprocessing as mp
from functools import partial
import math
from tqdm import tqdm

def extract_organism_name(taxon_name):
    """
    Extract the organism name part from a taxon name, removing organelle information.

    Args:
        taxon_name: The taxon name to process

    Returns:
        The organism name without organelle information
    """
    # Handle names with organelle information
    if "." in taxon_name:
        # For names like "Genus_species.Mitochondria", get the part before the dot
        parts = taxon_name.split(".")
        return parts[0]

    return taxon_name

def get_taxid_for_name(name, env):
    """
    Get NCBI taxid for a name using taxonkit name2taxid.

    Args:
        name: Taxon name to look up
        env: Environment variables for subprocess

    Returns:
        Taxid if found, None otherwise
    """
    # Extract organism name without organelle information
    organism_name = extract_organism_name(name)

    # Replace underscores with spaces for better matching
    clean_name = organism_name.replace("_", " ")

    try:
        # Run taxonkit name2taxid
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=clean_name,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0":
                return parts[1]

        # If full name doesn't match, try genus
        if "_" in organism_name or " " in clean_name:
            genus = organism_name.split("_")[0] if "_" in organism_name else clean_name.split(" ")[0]
            genus_result = subprocess.run(
                ["taxonkit", "name2taxid"],
                input=genus,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )

            if genus_result.returncode == 0 and genus_result.stdout.strip():
                genus_parts = genus_result.stdout.strip().split('\t')
                if len(genus_parts) >= 2 and genus_parts[1] != "0":
                    return genus_parts[1]

    except Exception:
        pass

    return None

def process_taxon_batch(taxon_batch, env):
    """
    Process a batch of taxon names to check if they are eukaryotic.

    Args:
        taxon_batch: List of taxon names to check
        env: Environment variables for subprocess

    Returns:
        Dictionary mapping taxon names to boolean (True if Eukaryotic, False otherwise)
    """
    # Create a temporary file with the taxon names
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        taxon_list_path = temp_file.name
        for name in taxon_batch:
            # Clean the name by removing organelle information and underscores
            organism_name = extract_organism_name(name)
            clean_name = organism_name.replace("_", " ")
            temp_file.write(f"{clean_name}\n")

    taxon_is_eukaryotic = {}

    try:
        # Run taxonkit name2taxid for the batch
        name2taxid_path = taxon_list_path + "_name2taxid.txt"
        with open(name2taxid_path, 'w') as output_file:
            subprocess.run(
                ["taxonkit", "name2taxid", taxon_list_path],
                stdout=output_file,
                stderr=subprocess.PIPE,
                env=env
            )

        # Parse the name2taxid results
        taxid_map = {}
        with open(name2taxid_path, 'r') as f:
            for i, line in enumerate(f):
                if i < len(taxon_batch):  # Ensure we don't go out of bounds
                    parts = line.strip().split('\t')
                    if len(parts) >= 2 and parts[1] != "0":
                        taxid_map[taxon_batch[i]] = parts[1]

        # For names that didn't get a taxid, try genus fallback
        genus_fallback_names = []
        for name in taxon_batch:
            if name not in taxid_map:
                organism_name = extract_organism_name(name)
                if "_" in organism_name:
                    genus = organism_name.split("_")[0]
                    genus_fallback_names.append((name, genus))

        # Process genus fallbacks if any
        if genus_fallback_names:
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as genus_file:
                genus_list_path = genus_file.name
                for _, genus in genus_fallback_names:
                    genus_file.write(f"{genus}\n")

            genus_taxid_path = genus_list_path + "_genus_taxid.txt"
            with open(genus_taxid_path, 'w') as output_file:
                subprocess.run(
                    ["taxonkit", "name2taxid", genus_list_path],
                    stdout=output_file,
                    stderr=subprocess.PIPE,
                    env=env
                )

            # Parse the genus taxid results
            with open(genus_taxid_path, 'r') as f:
                for i, line in enumerate(f):
                    if i < len(genus_fallback_names):
                        parts = line.strip().split('\t')
                        if len(parts) >= 2 and parts[1] != "0":
                            name = genus_fallback_names[i][0]
                            taxid_map[name] = parts[1]

            # Clean up genus files
            os.unlink(genus_list_path)
            if os.path.exists(genus_taxid_path):
                os.unlink(genus_taxid_path)

        # Get lineages for all taxids
        if taxid_map:
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as taxid_file:
                taxid_list_path = taxid_file.name
                for taxid in taxid_map.values():
                    taxid_file.write(f"{taxid}\n")

            lineage_path = taxid_list_path + "_lineage.txt"
            with open(lineage_path, 'w') as output_file:
                subprocess.run(
                    ["taxonkit", "lineage", taxid_list_path],
                    stdout=output_file,
                    stderr=subprocess.PIPE,
                    env=env
                )

            # Parse the lineage results
            taxid_to_lineage = {}
            with open(lineage_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        taxid = parts[0]
                        lineage = parts[1]
                        taxid_to_lineage[taxid] = "Eukaryota" in lineage

            # Map results back to taxon names
            for name, taxid in taxid_map.items():
                if taxid in taxid_to_lineage:
                    taxon_is_eukaryotic[name] = taxid_to_lineage[taxid]
                else:
                    # If lineage can't be determined, assume it's eukaryotic
                    taxon_is_eukaryotic[name] = True

            # Clean up taxid files
            os.unlink(taxid_list_path)
            if os.path.exists(lineage_path):
                os.unlink(lineage_path)

        # For any names that still don't have a result, assume they're eukaryotic
        for name in taxon_batch:
            if name not in taxon_is_eukaryotic:
                taxon_is_eukaryotic[name] = True

    except Exception as e:
        # If there's an error, assume all names in the batch are eukaryotic
        for name in taxon_batch:
            taxon_is_eukaryotic[name] = True

    finally:
        # Clean up temporary files
        try:
            os.unlink(taxon_list_path)
            if os.path.exists(name2taxid_path):
                os.unlink(name2taxid_path)
        except Exception:
            pass

    return taxon_is_eukaryotic

def check_if_eukaryotic(taxon_names):
    """
    Check if taxon names belong to Eukaryota using taxonkit.
    Uses batch processing and multiprocessing for better performance.

    For names with organelle information (e.g., "Genus_species.Mitochondria"),
    checks if the organism part (before the period) is eukaryotic.

    Args:
        taxon_names: List of taxon names to check

    Returns:
        Dictionary mapping taxon names to boolean (True if Eukaryotic, False otherwise)
    """
    if not taxon_names:
        return {}

    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir

    # Determine the number of workers and batch size
    num_workers = min(os.cpu_count() or 4, 8)  # Use at most 8 workers
    batch_size = max(50, math.ceil(len(taxon_names) / (num_workers * 2)))  # Ensure enough batches

    # Split the taxon names into batches
    batches = [taxon_names[i:i + batch_size] for i in range(0, len(taxon_names), batch_size)]
    print(f"Processing {len(taxon_names)} taxon names in {len(batches)} batches using {num_workers} workers")

    # Process batches in parallel
    results = {}

    # Also check organism parts for names with organelle information
    organism_parts = {}
    for name in taxon_names:
        if "." in name:
            organism_part = extract_organism_name(name)
            organism_parts[name] = organism_part

    with mp.Pool(processes=num_workers) as pool:
        # Create a partial function with the environment
        process_batch = partial(process_taxon_batch, env=env)

        # Process batches in parallel with progress bar
        for batch_result in tqdm(pool.imap(process_batch, batches), total=len(batches), desc="Processing taxon batches"):
            results.update(batch_result)

    # For names with organelle information, also check if the organism part is eukaryotic
    if organism_parts:
        print(f"Checking {len(organism_parts)} organism parts for names with organelle information...")
        organism_names = list(set(organism_parts.values()))
        organism_batches = [organism_names[i:i + batch_size] for i in range(0, len(organism_names), batch_size)]

        organism_results = {}
        with mp.Pool(processes=num_workers) as pool:
            for batch_result in tqdm(pool.imap(process_batch, organism_batches),
                                    total=len(organism_batches),
                                    desc="Processing organism parts"):
                organism_results.update(batch_result)

        # If the organism part is eukaryotic, mark the full name as eukaryotic
        for full_name, organism_part in organism_parts.items():
            if organism_part in organism_results and organism_results[organism_part]:
                results[full_name] = True
                print(f"✅ Keeping {full_name} because organism part {organism_part} is eukaryotic")

    print(f"Completed processing {len(results)} taxon names")
    return results

def main():
    # Input and output file paths
    input_file = "eukcensus_by_genus.csv"
    output_file = "eukcensus_by_genus_eukaryotes_only.csv"

    print(f"Reading input file: {input_file}")

    # Read the CSV file
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    print(f"Successfully read {len(df)} rows from the input file")

    # Filter out rows with "Bacteria" or "Archaea" in the taxon_name
    print("Filtering out rows with 'Bacteria' or 'Archaea' in the taxon_name...")

    # Create a list to store indices of rows to remove
    indices_to_remove = []

    # Check each row with progress bar
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Checking for bacterial/archaeal terms"):
        taxon_name = row['taxon_name']

        # Special case: If the name has organelle information (contains a period)
        # We need to check if the organism part (before the period) is eukaryotic
        if "." in taxon_name:
            organism_part = extract_organism_name(taxon_name)
            # Skip this check - we'll handle it in the eukaryotic check phase
            continue

        # Check if the taxon name contains bacterial or archaeal terms
        if any(term in taxon_name.lower() for term in ['bacill', 'cocc', 'bacterium', 'bacteria', 'archaea', 'vibrio']):
            indices_to_remove.append(idx)
            # Don't print each removal to reduce verbosity

    # Remove the identified rows
    filtered_df = df.drop(indices_to_remove)

    print(f"Removed {len(df) - len(filtered_df)} rows with 'Bacteria' or 'Archaea' in the name")
    print(f"Remaining rows: {len(filtered_df)}")

    # Check if taxon names are Eukaryotic
    print("Checking if taxon names are Eukaryotic...")
    taxon_names = filtered_df['taxon_name'].unique().tolist()
    taxon_is_eukaryotic = check_if_eukaryotic(taxon_names)

    print(f"Checked {len(taxon_is_eukaryotic)} unique taxon names")

    # Filter out non-Eukaryotic entries based on taxon name
    print("Filtering out non-Eukaryotic entries...")
    eukaryotic_rows = []
    non_eukaryotic_entries = []

    for _, row in tqdm(filtered_df.iterrows(), total=len(filtered_df), desc="Filtering entries"):
        taxon_name = row['taxon_name']

        # Special case for organelle entries: check if the organism part is eukaryotic
        if "." in taxon_name:
            organism_part = extract_organism_name(taxon_name)
            # If the organism part is in our dictionary and is eukaryotic, keep the row
            if organism_part in taxon_is_eukaryotic and taxon_is_eukaryotic[organism_part]:
                print(f"✅ Keeping organelle entry: {taxon_name} (organism part: {organism_part} is eukaryotic)")
                eukaryotic_rows.append(row)
                continue

        # If taxon name is not in our dictionary, keep the row (benefit of the doubt)
        if taxon_name not in taxon_is_eukaryotic:
            eukaryotic_rows.append(row)
        # If taxon name is Eukaryotic, keep the row
        elif taxon_is_eukaryotic[taxon_name]:
            eukaryotic_rows.append(row)
        # Otherwise, it's non-Eukaryotic, so skip it
        else:
            non_eukaryotic_entries.append((row['taxid'], taxon_name))

    # Create a new DataFrame with only Eukaryotic rows
    eukaryotic_df = pd.DataFrame(eukaryotic_rows)

    print(f"Removed {len(filtered_df) - len(eukaryotic_df)} non-Eukaryotic entries based on taxid")
    print(f"Final number of rows: {len(eukaryotic_df)}")

    # Save the filtered data to the output file
    eukaryotic_df.to_csv(output_file, index=False)

    print(f"Saved filtered data to {output_file}")

    # Print some statistics
    print("\nStatistics:")
    print(f"Original rows: {len(df)}")
    print(f"Rows after filtering out 'Bacteria' and 'Archaea': {len(filtered_df)}")
    print(f"Final Eukaryotic rows: {len(eukaryotic_df)}")

    # Print some examples of non-Eukaryotic entries that were removed
    if non_eukaryotic_entries:
        print("\nExamples of non-Eukaryotic entries that were removed:")
        for i, (taxid, taxon_name) in enumerate(non_eukaryotic_entries[:10]):
            print(f"{i+1}. Taxid: {taxid}, Taxon: {taxon_name}")

        if len(non_eukaryotic_entries) > 10:
            print(f"... and {len(non_eukaryotic_entries) - 10} more")

if __name__ == "__main__":
    main()
