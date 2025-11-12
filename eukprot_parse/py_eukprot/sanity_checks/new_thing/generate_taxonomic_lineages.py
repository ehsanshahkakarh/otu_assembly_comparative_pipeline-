#!/usr/bin/env python3
"""
Generate Taxonomic Lineages from Taxids

This script takes the taxids from the eukprot_names_with_taxids.csv file
and uses taxonkit to generate complete taxonomic lineages for each taxid.
It outputs a new CSV file with the original data plus the taxonomic lineage
information.

Usage:
    python generate_taxonomic_lineages.py

Output:
    A CSV file named 'eukprot_with_lineages.csv' containing the original data
    plus taxonomic lineage information.
"""

import pandas as pd
import os
import subprocess
import logging
from pathlib import Path
import sys
import time
import tempfile
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("generate_taxonomic_lineages.log"),
        logging.StreamHandler()
    ]
)

# NCBI taxdump directory
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")

def get_lineage_for_taxid(taxid, env=None):
    """
    Get the taxonomic lineage for a taxid using taxonkit lineage.
    Returns a dictionary with lineage information.
    """
    if taxid == "FAILED" or not taxid:
        return {
            "lineage": "",
            "lineage_taxids": "",
            "kingdom": "",
            "phylum": "",
            "class": "",
            "order": "",
            "family": "",
            "genus": "",
            "species": ""
        }

    if env is None:
        env = os.environ.copy()
        env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    try:
        # First, just get the basic lineage to make sure it works
        basic_cmd = ["taxonkit", "lineage", taxid]
        logging.debug(f"Running command: {' '.join(basic_cmd)}")

        basic_result = subprocess.run(
            basic_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if basic_result.returncode != 0 or not basic_result.stdout.strip():
            logging.debug(f"Basic lineage command failed for taxid {taxid}: {basic_result.stderr}")
            return {
                "lineage": f"Error: {basic_result.stderr}",
                "lineage_taxids": "",
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "species": ""
            }

        # If basic command works, try the full command
        full_cmd = ["taxonkit", "lineage", "--show-name", "--show-rank", "--show-lineage-taxids", taxid]
        logging.debug(f"Running command: {' '.join(full_cmd)}")

        result = subprocess.run(
            full_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')

            # Debug output to see what we're getting
            logging.debug(f"Lineage output for taxid {taxid}: {result.stdout}")

            # Make sure we have enough parts
            if len(parts) < 4:
                logging.debug(f"Not enough parts in lineage output for taxid {taxid}: {parts}")
                lineage = parts[1] if len(parts) > 1 else ""
                lineage_taxids = ""
            else:
                lineage = parts[1]
                lineage_taxids = parts[3]

            # Try to get formatted ranks
            try:
                format_cmd = ["taxonkit", "reformat", "--format", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}"]
                logging.debug(f"Running command: {' '.join(format_cmd)}")

                format_result = subprocess.run(
                    format_cmd,
                    input=result.stdout,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env
                )

                if format_result.returncode == 0 and format_result.stdout.strip():
                    logging.debug(f"Reformat output for taxid {taxid}: {format_result.stdout}")
                    format_parts = format_result.stdout.strip().split('\t')

                    if len(format_parts) >= 8:
                        return {
                            "lineage": lineage,
                            "lineage_taxids": lineage_taxids,
                            "kingdom": format_parts[1],
                            "phylum": format_parts[2],
                            "class": format_parts[3],
                            "order": format_parts[4],
                            "family": format_parts[5],
                            "genus": format_parts[6],
                            "species": format_parts[7]
                        }
                    else:
                        logging.debug(f"Not enough parts in reformat output for taxid {taxid}: {format_parts}")
                else:
                    logging.debug(f"Reformat command failed for taxid {taxid}: {format_result.stderr}")
            except Exception as e:
                logging.debug(f"Error in reformat for taxid {taxid}: {e}")

            # If we get here, the reformat failed but we still have the basic lineage
            return {
                "lineage": lineage,
                "lineage_taxids": lineage_taxids,
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "species": ""
            }
    except Exception as e:
        logging.debug(f"Error getting lineage for taxid {taxid}: {e}")
        import traceback
        logging.debug(traceback.format_exc())

    # If we get here, something went wrong
    return {
        "lineage": f"Error processing taxid: {taxid}",
        "lineage_taxids": "",
        "kingdom": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": ""
    }

def process_taxids_batch(taxids_batch):
    """
    Process a batch of taxids to get their lineages.
    Returns a list of (taxid, lineage_info) tuples.
    """
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    # First, try a simple test to make sure taxonkit is working
    test_taxid = "9606"  # Human
    try:
        test_result = subprocess.run(
            ["taxonkit", "lineage", test_taxid],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        if test_result.returncode == 0 and test_result.stdout.strip():
            logging.debug(f"Taxonkit test successful: {test_result.stdout.strip()}")
        else:
            logging.error(f"Taxonkit test failed: {test_result.stderr}")
            logging.error("This may indicate a problem with the taxonkit installation or database")
    except Exception as e:
        logging.error(f"Error running taxonkit test: {e}")

    results = []
    for taxid in taxids_batch:
        try:
            lineage_info = get_lineage_for_taxid(taxid, env)
            results.append((taxid, lineage_info))
        except Exception as e:
            logging.error(f"Error processing taxid {taxid}: {e}")
            results.append((taxid, {
                "lineage": f"Error: {str(e)}",
                "lineage_taxids": "",
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "species": ""
            }))
    return results

def process_taxids_with_file(taxids):
    """
    Process taxids by writing them to a file and using taxonkit batch processing.
    This is more efficient than processing one by one.
    """
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    # Create a temporary file with the taxids
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_path = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        # Run taxonkit lineage on the file
        lineage_cmd = ["taxonkit", "lineage", "--show-name", "--show-rank", "--show-lineage-taxids", temp_path]
        logging.info(f"Running command: {' '.join(lineage_cmd)}")

        lineage_result = subprocess.run(
            lineage_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if lineage_result.returncode != 0:
            logging.error(f"Error running taxonkit lineage: {lineage_result.stderr}")
            return {}

        # Save the lineage output to a temporary file
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as lineage_file:
            lineage_path = lineage_file.name
            lineage_file.write(lineage_result.stdout)

        # Run taxonkit reformat on the lineage output with proper rank extraction
        # Use the --fill-miss-rank option to ensure consistent output format
        reformat_cmd = ["taxonkit", "reformat", "--format", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}",
                        "--fill-miss-rank", "--miss-rank-repl", "", lineage_path]
        logging.info(f"Running command: {' '.join(reformat_cmd)}")

        reformat_result = subprocess.run(
            reformat_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if reformat_result.returncode != 0:
            logging.error(f"Error running taxonkit reformat: {reformat_result.stderr}")

        # Also run taxonkit lineage with --show-rank to get rank information
        rank_cmd = ["taxonkit", "lineage", "--show-name", "--show-rank", lineage_path]
        logging.info(f"Running command: {' '.join(rank_cmd)}")

        rank_result = subprocess.run(
            rank_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if rank_result.returncode != 0:
            logging.error(f"Error running taxonkit with rank info: {rank_result.stderr}")

        # Parse the results and initialize the taxid_to_lineage dictionary
        taxid_to_lineage = {}

        # First parse the lineage output for basic lineage information
        lineage_lines = lineage_result.stdout.strip().split('\n')
        for line in lineage_lines:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                taxid = parts[0]
                lineage = parts[1]
                rank = parts[2] if len(parts) > 2 else ""
                lineage_taxids = parts[3] if len(parts) > 3 else ""

                # Replace semicolons with spaces in the lineage
                lineage_with_spaces = lineage.replace(';', ' ')
                lineage_taxids_with_spaces = lineage_taxids.replace(';', ' ')

                # Initialize the entry with basic lineage info
                taxid_to_lineage[taxid] = {
                    "lineage": lineage_with_spaces,
                    "rank": rank,
                    "lineage_taxids": lineage_taxids_with_spaces,
                    "kingdom": "",
                    "phylum": "",
                    "class": "",
                    "order": "",
                    "family": "",
                    "genus": "",
                    "species": ""
                }

        # Extract rank information from the reformat output
        if reformat_result.returncode == 0:
            reformat_lines = reformat_result.stdout.strip().split('\n')
            for line in reformat_lines:
                parts = line.strip().split('\t')
                if len(parts) >= 8:  # taxid + 7 ranks
                    taxid = parts[0]
                    if taxid in taxid_to_lineage:
                        # Store the full lineage parts for each rank
                        kingdom_full = parts[1].strip() if parts[1].strip() else ""
                        phylum_full = parts[2].strip() if parts[2].strip() else ""
                        class_full = parts[3].strip() if parts[3].strip() else ""
                        order_full = parts[4].strip() if parts[4].strip() else ""
                        family_full = parts[5].strip() if parts[5].strip() else ""
                        genus_full = parts[6].strip() if parts[6].strip() else ""
                        species_full = parts[7].strip() if parts[7].strip() else ""

                        # Also extract just the last part of each rank (the actual rank name)
                        kingdom = kingdom_full.split()[-1] if kingdom_full else ""
                        phylum = phylum_full.split()[-1] if phylum_full else ""
                        class_ = class_full.split()[-1] if class_full else ""
                        order = order_full.split()[-1] if order_full else ""
                        family = family_full.split()[-1] if family_full else ""
                        genus = genus_full.split()[-1] if genus_full else ""
                        species = species_full.split()[-1] if species_full else ""

                        taxid_to_lineage[taxid].update({
                            "kingdom": kingdom,
                            "phylum": phylum,
                            "class": class_,
                            "order": order,
                            "family": family,
                            "genus": genus,
                            "species": species,
                            "kingdom_full": kingdom_full,
                            "phylum_full": phylum_full,
                            "class_full": class_full,
                            "order_full": order_full,
                            "family_full": family_full,
                            "genus_full": genus_full,
                            "species_full": species_full
                        })

        # If reformat failed, try to extract ranks from the lineage information
        elif rank_result.returncode == 0:
            rank_lines = rank_result.stdout.strip().split('\n')
            for line in rank_lines:
                parts = line.strip().split('\t')
                if len(parts) >= 3:  # taxid + lineage + rank
                    taxid = parts[0]
                    lineage = parts[1]
                    rank = parts[2]

                    if taxid in taxid_to_lineage:
                        # Extract ranks from the lineage
                        lineage_parts = lineage.split(';')

                        # Iteratively extract ranks in decreasing order
                        # Start with the assumption that Eukaryota is the domain
                        if "Eukaryota" in lineage_parts:
                            # Find the index of Eukaryota
                            eukaryota_idx = lineage_parts.index("Eukaryota")

                            # Kingdom is typically the first major group after Eukaryota
                            kingdom_idx = eukaryota_idx + 1
                            if kingdom_idx < len(lineage_parts):
                                # Extract just the name, not the full path
                                taxid_to_lineage[taxid]["kingdom"] = lineage_parts[kingdom_idx].strip()

                                # Phylum is typically after kingdom
                                phylum_idx = kingdom_idx + 1
                                if phylum_idx < len(lineage_parts):
                                    taxid_to_lineage[taxid]["phylum"] = lineage_parts[phylum_idx].strip()

                                    # Class is typically after phylum
                                    class_idx = phylum_idx + 1
                                    if class_idx < len(lineage_parts):
                                        taxid_to_lineage[taxid]["class"] = lineage_parts[class_idx].strip()

                                        # Order is typically after class
                                        order_idx = class_idx + 1
                                        if order_idx < len(lineage_parts):
                                            taxid_to_lineage[taxid]["order"] = lineage_parts[order_idx].strip()

                                            # Family is typically after order
                                            family_idx = order_idx + 1
                                            if family_idx < len(lineage_parts):
                                                taxid_to_lineage[taxid]["family"] = lineage_parts[family_idx].strip()

                                                # Genus is typically after family
                                                genus_idx = family_idx + 1
                                                if genus_idx < len(lineage_parts):
                                                    taxid_to_lineage[taxid]["genus"] = lineage_parts[genus_idx].strip()

                                                    # Species is typically after genus
                                                    species_idx = genus_idx + 1
                                                    if species_idx < len(lineage_parts):
                                                        taxid_to_lineage[taxid]["species"] = lineage_parts[species_idx].strip()

        return taxid_to_lineage

    finally:
        # Clean up temporary files
        try:
            os.unlink(temp_path)
            if 'lineage_path' in locals():
                os.unlink(lineage_path)
        except Exception as e:
            logging.debug(f"Error cleaning up temporary files: {e}")

def main():
    """Generate taxonomic lineages for taxids."""
    start_time = time.time()

    # Get the current directory
    script_dir = Path(__file__).resolve().parent

    # Input file path
    input_file = script_dir / "eukprot_names_with_taxids.csv"

    # Output file path
    output_file = script_dir / "eukprot_with_lineages.csv"

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
        if not input_file.exists():
            logging.error(f"❌ Input file not found: {input_file}")
            sys.exit(1)

        logging.info(f"Loading data from {input_file}...")
        df = pd.read_csv(input_file)

        # Get column names
        if len(df.columns) < 2:
            logging.error(f"❌ Input file does not have the expected format. Found columns: {df.columns}")
            sys.exit(1)

        name_column = df.columns[0]  # First column should be the name
        taxid_column = df.columns[1]  # Second column should be the taxid

        # Filter out failed taxids
        valid_df = df[df[taxid_column] != "FAILED"].copy()
        failed_df = df[df[taxid_column] == "FAILED"].copy()

        logging.info(f"Found {len(valid_df)} entries with valid taxids")
        logging.info(f"Found {len(failed_df)} entries with failed taxids")

        # Get unique taxids to process
        unique_taxids = valid_df[taxid_column].unique()
        logging.info(f"Found {len(unique_taxids)} unique taxids to process")

        # Process all taxids at once using file-based approach
        logging.info(f"Processing {len(unique_taxids)} unique taxids using file-based approach")

        # Convert taxids to strings
        unique_taxids_str = [str(taxid) for taxid in unique_taxids]

        # Process taxids
        taxid_to_lineage = process_taxids_with_file(unique_taxids_str)

        if not taxid_to_lineage:
            logging.error("Failed to get lineage information. Trying individual processing as fallback.")

            # Fallback to individual processing
            taxid_to_lineage = {}
            batch_size = 10
            batches = [unique_taxids_str[i:i + batch_size] for i in range(0, len(unique_taxids_str), batch_size)]

            for batch in tqdm(batches, desc="Processing taxids (fallback)"):
                batch_results = process_taxids_batch(batch)
                for taxid, lineage_info in batch_results:
                    taxid_to_lineage[taxid] = lineage_info

        # Add lineage information to the dataframe
        for column in ["kingdom", "phylum", "class", "order", "family", "genus", "species",
                      "kingdom_full", "phylum_full", "class_full", "order_full", "family_full", "genus_full", "species_full",
                      "lineage", "lineage_taxids"]:
            valid_df[column] = ""

        # Count how many taxids have lineage info
        taxids_with_lineage = sum(1 for taxid in taxid_to_lineage if taxid_to_lineage[taxid].get("lineage"))
        logging.info(f"Got lineage information for {taxids_with_lineage} out of {len(unique_taxids)} taxids")

        # Add the lineage info to the dataframe
        for i, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Adding lineage info to dataframe"):
            taxid = str(row[taxid_column])  # Convert to string to ensure matching
            if taxid in taxid_to_lineage:
                lineage_info = taxid_to_lineage[taxid]
                for column, value in lineage_info.items():
                    if column != "rank":  # Skip the rank column
                        valid_df.at[i, column] = value

        # Add empty lineage columns to failed entries
        for column in ["kingdom", "phylum", "class", "order", "family", "genus", "species", "lineage", "lineage_taxids"]:
            failed_df[column] = ""

        # Combine valid and failed entries
        result_df = pd.concat([valid_df, failed_df], ignore_index=True)

        # Now create a new column that contains the full lineage up to the match_type
        result_df["full_lineage_to_match"] = ""

        # Process each row to create the appropriate lineage
        for i, row in result_df.iterrows():
            match_type = row["match_type"] if "match_type" in result_df.columns else ""

            if match_type == "genus":
                # Include lineage up to genus
                lineage_parts = []
                if row["kingdom_full"]:
                    lineage_parts.append(row["kingdom_full"])
                if row["phylum_full"]:
                    lineage_parts.append(row["phylum_full"])
                if row["class_full"]:
                    lineage_parts.append(row["class_full"])
                if row["order_full"]:
                    lineage_parts.append(row["order_full"])
                if row["family_full"]:
                    lineage_parts.append(row["family_full"])
                if row["genus_full"]:
                    lineage_parts.append(row["genus_full"])

                result_df.at[i, "full_lineage_to_match"] = " ".join(lineage_parts)
            elif match_type == "full":
                # Include full lineage up to species
                lineage_parts = []
                if row["kingdom_full"]:
                    lineage_parts.append(row["kingdom_full"])
                if row["phylum_full"]:
                    lineage_parts.append(row["phylum_full"])
                if row["class_full"]:
                    lineage_parts.append(row["class_full"])
                if row["order_full"]:
                    lineage_parts.append(row["order_full"])
                if row["family_full"]:
                    lineage_parts.append(row["family_full"])
                if row["genus_full"]:
                    lineage_parts.append(row["genus_full"])
                if row["species_full"]:
                    lineage_parts.append(row["species_full"])

                result_df.at[i, "full_lineage_to_match"] = " ".join(lineage_parts)
            else:
                # For failed matches or unknown match types, use the original lineage
                result_df.at[i, "full_lineage_to_match"] = row["lineage"]

        # Reorder columns to match the requested format
        column_order = [name_column, taxid_column, "match_type",
                        "kingdom", "phylum", "class", "order", "family", "genus", "species",
                        "full_lineage_to_match", "lineage", "lineage_taxids"]

        # Make sure all columns exist
        for col in column_order:
            if col not in result_df.columns:
                result_df[col] = ""

        # Reorder columns
        result_df = result_df[column_order]

        # Remove the individual full rank columns as they're no longer needed
        for col in ["kingdom_full", "phylum_full", "class_full", "order_full", "family_full", "genus_full", "species_full"]:
            if col in result_df.columns:
                result_df = result_df.drop(columns=[col])

        # Replace any remaining semicolons with spaces in all columns
        for col in result_df.columns:
            if result_df[col].dtype == 'object':  # Only process string columns
                result_df[col] = result_df[col].astype(str).apply(lambda x: x.replace(';', ' '))

        # Save to CSV file
        result_df.to_csv(output_file, index=False)

        # Calculate elapsed time
        elapsed_time = time.time() - start_time

        logging.info(f"✅ Successfully processed {len(unique_taxids)} unique taxids in {elapsed_time:.2f} seconds")
        logging.info(f"✅ Results saved to {output_file}")

        # Print a preview of the results
        logging.info("\nPreview of results:")
        preview = result_df.head(5)
        for _, row in preview.iterrows():
            name = row[name_column]
            taxid = row[taxid_column]
            if taxid != "FAILED":
                lineage = row["lineage"]
                logging.info(f"✅ {name} (taxid: {taxid}): {lineage}")
            else:
                logging.info(f"❌ {name}: No lineage (FAILED taxid)")

        if len(result_df) > 5:
            logging.info(f"... and {len(result_df) - 5} more")

    except Exception as e:
        logging.error(f"❌ Error: {e}")
        import traceback
        logging.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
