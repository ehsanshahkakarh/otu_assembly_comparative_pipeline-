#!/usr/bin/env python3
"""
Map Organism Names to NCBI Taxids (Parallel Version)

This script takes the organism names from the Name_to_Use column in the
Eukprot_included_datasets.txt file, finds their corresponding NCBI taxids
using the taxonkit name2taxid function, and outputs a new CSV file with
the original names, their taxids, and match type.

For names that don't have a matching taxid, "FAILED" is added in the taxid column.

The script implements a multi-stage fallback strategy:
1. First tries the full name match
2. If that fails, tries extracting and using just the genus name
3. For remaining failures, attempts using the Previous_Names column if available
4. Finally, tries one more approach using the first word as genus

This version uses parallel processing to speed up the taxid lookup process.

Usage:
    python map_names_to_taxids_parallel.py

Output:
    A CSV file named 'eukprot_names_with_taxids.csv' with columns:
    'Name_to_use', 'taxid', and 'match_type'.

    The match_type column indicates how the match was made:
    - 'full': Full name match
    - 'genus': Genus-only match
    - 'previous_name': Match using previous name
    - 'none': No match found
"""

import pandas as pd
import os
import subprocess
import logging
from pathlib import Path
import sys
import multiprocessing
from functools import partial
from concurrent.futures import ProcessPoolExecutor
import time
import re

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("tqdm not available, progress bars will be disabled")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("map_names_to_taxids_parallel.log"),
        logging.StreamHandler()
    ]
)

# NCBI taxdump directory
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")

def extract_genus_from_name(name):
    """
    Extract the genus part from a name.

    This function tries multiple approaches:
    1. First, look for 'sp' pattern (e.g., "Acanthamoeba sp")
    2. If that fails, try to get the first word as the genus

    Examples:
    - "Acanthamoeba sp" -> "Acanthamoeba"
    - "Acanthamoeba sp." -> "Acanthamoeba"
    - "Acanthamoeba sp ATCC" -> "Acanthamoeba"
    - "Acanthamoeba castellanii" -> "Acanthamoeba"
    """
    if not name or pd.isna(name):
        return None

    # First approach: Look for " sp" with a space before it
    # This will match " sp", " sp.", " sp ", etc.
    sp_match = re.search(r'(.+?)\s+sp\b.*', name, re.IGNORECASE)
    if sp_match:
        genus = sp_match.group(1).strip()
        return genus

    # Second approach: Try to get the first word as the genus
    # This assumes binomial nomenclature (Genus species)
    parts = name.split()
    if len(parts) >= 1:
        # Check if the first word looks like a genus (capitalized)
        first_word = parts[0].strip()
        if first_word and first_word[0].isupper():
            return first_word

    return None

def get_taxid_for_name(name, env=None, try_genus=True):
    """
    Get NCBI taxid for a name using taxonkit name2taxid.
    If try_genus is True and the name contains 'sp', will also try to get the taxid for just the genus.
    Returns the taxid if found, "FAILED" otherwise.
    """
    if not name or pd.isna(name):
        return "FAILED"

    if env is None:
        env = os.environ.copy()
        env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    # First try the full name
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
    except Exception as e:
        logging.debug(f"Error getting taxid for full name '{name}': {e}")

    # If we get here, the full name lookup failed
    # If try_genus is True, try extracting and looking up just the genus
    if try_genus:
        genus = extract_genus_from_name(name)
        if genus:
            logging.debug(f"Trying genus extraction for '{name}' -> '{genus}'")
            try:
                # Run taxonkit name2taxid for the genus
                result = subprocess.run(
                    ["taxonkit", "name2taxid"],
                    input=genus,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env
                )

                if result.returncode == 0 and result.stdout.strip():
                    parts = result.stdout.strip().split('\t')
                    if len(parts) >= 2 and parts[1].strip() and parts[1] != "0":
                        logging.debug(f"Found taxid for genus '{genus}': {parts[1].strip()}")
                        return parts[1].strip()
            except Exception as e:
                logging.debug(f"Error getting taxid for genus '{genus}': {e}")

    return "FAILED"

def process_batch(names_batch, previous_names_map=None):
    """
    Process a batch of names.
    Returns a list of (name, taxid, match_type) tuples.

    Args:
        names_batch: List of names to process
        previous_names_map: Dictionary mapping current names to previous names
    """
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    results = []
    for name in names_batch:
        # First try the full name
        taxid = get_taxid_for_name(name, env, try_genus=False)
        match_type = "full"

        # If that fails, try to extract and use the genus
        if taxid == "FAILED":
            genus = extract_genus_from_name(name)
            if genus:
                genus_taxid = get_taxid_for_name(genus, env, try_genus=False)
                if genus_taxid != "FAILED":
                    taxid = genus_taxid
                    match_type = "genus"

        # If still failed and we have a previous name, try that
        if taxid == "FAILED" and previous_names_map and name in previous_names_map and previous_names_map[name]:
            previous_names_str = previous_names_map[name]

            # Split by comma to handle multiple previous names
            previous_names = [pn.strip() for pn in previous_names_str.split(',')]
            logging.debug(f"Trying {len(previous_names)} previous name(s) for '{name}': {previous_names}")

            # Try each previous name in order
            for previous_name in previous_names:
                if not previous_name:
                    continue

                logging.debug(f"Trying previous name: '{previous_name}'")

                # Try the full previous name
                prev_taxid = get_taxid_for_name(previous_name, env, try_genus=False)
                if prev_taxid != "FAILED":
                    taxid = prev_taxid
                    match_type = "previous_name"
                    logging.debug(f"Found taxid using previous name '{previous_name}': {taxid}")
                    break  # Stop once we find a match
                else:
                    # Try the genus of the previous name
                    prev_genus = extract_genus_from_name(previous_name)
                    if prev_genus:
                        prev_genus_taxid = get_taxid_for_name(prev_genus, env, try_genus=False)
                        if prev_genus_taxid != "FAILED":
                            taxid = prev_genus_taxid
                            match_type = "previous_name_genus"
                            logging.debug(f"Found taxid using previous name genus '{prev_genus}': {taxid}")
                            break  # Stop once we find a match

        # If still failed, match_type remains "none"
        if taxid == "FAILED":
            match_type = "none"

        results.append((name, taxid, match_type))
    return results

def main():
    """Map organism names to NCBI taxids using parallel processing."""
    start_time = time.time()

    # Get the current directory
    script_dir = Path(__file__).resolve().parent

    # Input file path - either use the CSV from extract_name_to_use.py or directly from the source
    csv_file = script_dir / "eukprot_names_to_use.csv"
    txt_file = script_dir / "Eukprot_included_datasets.txt"

    # Output file path (save in the current directory)
    output_file = script_dir / "eukprot_names_with_taxids.csv"

    # Log the output file path
    logging.info(f"Output will be saved to: {output_file}")

    try:
        # Check if taxonkit is available
        try:
            result = subprocess.run(
                ["which", "taxonkit"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.returncode != 0:
                logging.error("❌ taxonkit not found in PATH. Please install taxonkit.")
                sys.exit(1)
            logging.info("✅ taxonkit is available")
        except Exception as e:
            logging.error(f"❌ Error checking taxonkit: {e}")
            sys.exit(1)

        # Load the data
        previous_names_map = {}  # Initialize empty map

        if csv_file.exists():
            logging.info(f"Loading data from {csv_file}...")
            df = pd.read_csv(csv_file)
            names_column = "Name_to_use"

            # Check if Previous_Names column exists and create a mapping
            if "Previous_Names" in df.columns:
                logging.info("Found Previous_Names column, creating mapping for alternative names")
                for idx, row in df.iterrows():
                    if pd.notna(row["Previous_Names"]) and row["Previous_Names"].strip():
                        previous_names_map[row[names_column]] = row["Previous_Names"].strip()

                logging.info(f"Created mapping with {len(previous_names_map)} alternative names")
                # Preview some mappings
                preview_count = min(5, len(previous_names_map))
                if preview_count > 0:
                    logging.info("Preview of name mappings:")
                    for i, (current, previous) in enumerate(list(previous_names_map.items())[:preview_count]):
                        logging.info(f"  {i+1}. {current} → {previous}")
        elif txt_file.exists():
            logging.info(f"Loading data from {txt_file}...")
            df = pd.read_csv(txt_file, sep='\t')
            # Extract the Name_to_Use column and replace underscores with spaces
            names = df["Name_to_Use"].dropna()
            names = names.str.replace('_', ' ')

            # Check if Previous_Names column exists in the original file
            if "Previous_Names" in df.columns:
                logging.info("Found Previous_Names column in the original file")
                # Create a mapping from current names to previous names
                for idx, row in df.iterrows():
                    if pd.notna(row["Previous_Names"]) and row["Previous_Names"].strip():
                        current_name = row["Name_to_Use"].replace('_', ' ')
                        previous_name = row["Previous_Names"].replace('_', ' ')
                        previous_names_map[current_name] = previous_name

                logging.info(f"Created mapping with {len(previous_names_map)} alternative names")

            # Create a new DataFrame with just the Name_to_use column
            df = pd.DataFrame({"Name_to_use": names})
            names_column = "Name_to_use"
        else:
            logging.error("❌ No input file found.")
            sys.exit(1)

        # Create a new DataFrame for the results
        results_df = pd.DataFrame()
        results_df[names_column] = df[names_column]
        results_df["taxid"] = "FAILED"  # Initialize all with FAILED

        # Get list of names to process
        names_list = df[names_column].tolist()

        # Determine the number of processes to use (use 80% of available cores)
        num_processes = max(1, int(multiprocessing.cpu_count() * 0.8))
        logging.info(f"Using {num_processes} processes for parallel processing")

        # Determine batch size (aim for ~100 batches total)
        batch_size = max(1, len(names_list) // (num_processes * 10))
        logging.info(f"Processing in batches of {batch_size} names")

        # Split names into batches
        batches = [names_list[i:i + batch_size] for i in range(0, len(names_list), batch_size)]
        logging.info(f"Split into {len(batches)} batches")

        # Process batches in parallel
        all_results = []

        logging.info(f"Processing {len(names_list)} organism names in parallel...")
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            # Submit all batches to the executor, passing the previous_names_map
            futures = [executor.submit(process_batch, batch, previous_names_map) for batch in batches]

            # Process results as they complete
            if TQDM_AVAILABLE:
                for future in tqdm(futures, total=len(futures), desc="Processing batches"):
                    batch_results = future.result()
                    all_results.extend(batch_results)
            else:
                for i, future in enumerate(futures):
                    batch_results = future.result()
                    all_results.extend(batch_results)
                    if (i + 1) % 10 == 0 or i + 1 == len(futures):
                        logging.info(f"Processed {i + 1}/{len(futures)} batches ({(i + 1)/len(futures)*100:.1f}%)")

        # Update the results DataFrame
        success_count = 0
        genus_match_count = 0
        previous_name_match_count = 0
        previous_name_genus_match_count = 0
        fail_count = 0

        # Add a column to track the match type
        results_df["match_type"] = "none"

        # First pass: Process all results
        for name, taxid, match_type in all_results:
            # Find the index of this name in the DataFrame
            idx = df.index[df[names_column] == name].tolist()
            if idx:
                results_df.at[idx[0], "taxid"] = taxid
                results_df.at[idx[0], "match_type"] = match_type

                if taxid != "FAILED":
                    if match_type == "genus":
                        genus_match_count += 1
                    elif match_type == "previous_name":
                        previous_name_match_count += 1
                    elif match_type == "previous_name_genus":
                        previous_name_genus_match_count += 1
                    success_count += 1
                else:
                    fail_count += 1

        # Second pass: Try one more time for failed entries by extracting the first word as genus
        failed_entries = results_df[results_df["taxid"] == "FAILED"]
        if not failed_entries.empty:
            logging.info(f"Trying one more approach for {len(failed_entries)} failed entries...")

            env = os.environ.copy()
            env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

            for idx, row in failed_entries.iterrows():
                name = row[names_column]
                parts = name.split()

                # Try the first word if it looks like a genus (capitalized)
                if len(parts) >= 1 and parts[0][0].isupper():
                    genus = parts[0]
                    genus_taxid = get_taxid_for_name(genus, env, try_genus=False)

                    if genus_taxid != "FAILED":
                        results_df.at[idx, "taxid"] = genus_taxid
                        results_df.at[idx, "match_type"] = "genus"
                        genus_match_count += 1
                        success_count += 1
                        fail_count -= 1

            logging.info(f"After second pass: {success_count} successes, {fail_count} failures")

        # Save to CSV file
        results_df.to_csv(output_file, index=False)

        # Calculate elapsed time
        elapsed_time = time.time() - start_time

        logging.info(f"✅ Successfully processed {len(df)} organism names in {elapsed_time:.2f} seconds")
        logging.info(f"✅ Found taxids for {success_count} names ({success_count/len(df)*100:.1f}%)")

        # Calculate full matches (total success minus all other match types)
        full_match_count = success_count - genus_match_count - previous_name_match_count - previous_name_genus_match_count

        logging.info(f"  - Full name matches: {full_match_count} ({full_match_count/len(df)*100:.1f}%)")
        logging.info(f"  - Genus-only matches: {genus_match_count} ({genus_match_count/len(df)*100:.1f}%)")

        if previous_name_match_count > 0:
            logging.info(f"  - Previous name matches: {previous_name_match_count} ({previous_name_match_count/len(df)*100:.1f}%)")

        if previous_name_genus_match_count > 0:
            logging.info(f"  - Previous name genus matches: {previous_name_genus_match_count} ({previous_name_genus_match_count/len(df)*100:.1f}%)")

        logging.info(f"❌ Failed to find taxids for {fail_count} names ({fail_count/len(df)*100:.1f}%)")
        logging.info(f"✅ Results saved to {output_file}")

        # Print the first 10 results as a preview
        logging.info("\nPreview of results:")
        preview = results_df.head(10)
        for _, row in preview.iterrows():
            if row["taxid"] == "FAILED":
                status = "❌"
                match_info = ""
            else:
                status = "✅"
                match_type = row["match_type"]
                if match_type == "genus":
                    match_info = " (genus match)"
                elif match_type == "previous_name":
                    match_info = " (previous name match)"
                elif match_type == "previous_name_genus":
                    match_info = " (previous name genus match)"
                else:
                    match_info = " (full match)"
            logging.info(f"{status} {row[names_column]}: {row['taxid']}{match_info}")

        if len(results_df) > 10:
            logging.info(f"... and {len(results_df) - 10} more")

    except Exception as e:
        logging.error(f"❌ Error: {e}")
        import traceback
        logging.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
