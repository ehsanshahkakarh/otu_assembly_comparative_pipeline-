#!/usr/bin/env python3
"""
Taxonkit utilities module for 18S Census Parser

Handles NCBI taxonkit integration with 4-tier fallback system for taxid resolution.
"""

import os
import subprocess
import tempfile
import logging
import math
import concurrent.futures
from typing import Dict, List, Tuple

from .taxon_cleaner import (
    clean_taxon_name,
    strip_trailing_numbers,
    extract_genus,
    extract_taxa_from_hyphenated
)


def process_taxon_batch(taxon_batch: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    Process a batch of taxon names to get their taxids using 4-tier fallback system.

    Tier 1: Direct lookup with cleaned names
    Tier 2: Genus fallback for species-level entries
    Tier 3: Number stripping fallback
    Tier 4: Hyphenated pattern extraction fallback

    Args:
        taxon_batch: List of taxon names to process

    Returns:
        Dictionary mapping taxon names to (taxid, method) tuples
    """
    # Set up environment - taxonkit uses its own built-in NCBI database
    env = os.environ.copy()

    results = {}

    # Clean the names
    cleaned_names = [clean_taxon_name(name) for name in taxon_batch]

    # TIER 1: Run taxonkit name2taxid for the batch (direct lookup)
    try:
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input="\n".join(cleaned_names),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0:
            # Check stderr for multi-match warnings
            multi_match_names = set()
            if result.stderr:
                for line in result.stderr.split('\n'):
                    if 'multiple TaxIds found for' in line:
                        # Extract the name from warning like: "multiple TaxIds found for 'Craniata'"
                        import re
                        match = re.search(r"'([^']+)'", line)
                        if match:
                            multi_match_names.add(match.group(1))

            # Parse the output
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[1] != "0" and parts[1].strip():
                    # Use the name from the output, not the index, to avoid misalignment
                    output_name = parts[0]
                    taxid_raw = parts[1]

                    # Clean the taxid - remove ALL whitespace including newlines
                    # This handles cases where taxonkit returns multi-line results
                    taxid = taxid_raw.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
                    taxid = ' '.join(taxid.split()).strip()

                    # Skip if the cleaned taxid contains spaces (indicates multiple taxids were concatenated)
                    if ' ' in taxid:
                        logging.warning(f"Skipping multi-match result for '{output_name}' (taxid contains multiple values: '{taxid_raw}')")
                        continue

                    # Skip if this name had multiple matches (detected in stderr)
                    if output_name in multi_match_names:
                        logging.warning(f"Skipping '{output_name}' - multiple taxids found by taxonkit")
                        continue

                    # Find the corresponding original name (handle cleaning)
                    original_name = None
                    for orig_name in taxon_batch:
                        if clean_taxon_name(orig_name) == output_name:
                            original_name = orig_name
                            break

                    if original_name:
                        # Only store if we haven't already found a match (first match wins)
                        if original_name not in results:
                            results[original_name] = (taxid, "direct")
    except Exception:
        pass

    # TIER 2: For names that didn't get a match, try genus fallback
    for name in taxon_batch:
        if name not in results:
            genus = extract_genus(name)
            if genus:
                try:
                    genus_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=genus.replace("_", " "),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if genus_result.returncode == 0 and genus_result.stdout.strip():
                        genus_parts = genus_result.stdout.strip().split('\t')
                        if len(genus_parts) >= 2 and genus_parts[1] != "0" and genus_parts[1].strip():
                            # Strip the taxid to remove any newlines or extra whitespace
                            results[name] = (genus_parts[1].strip(), "genus_fallback")
                except Exception:
                    pass

    # TIER 3: For names that still didn't get a match, try stripping numbers as a fallback
    for name in taxon_batch:
        if name not in results:
            # Try stripping numbers from the original name
            stripped_name = strip_trailing_numbers(name.replace("_", " "))
            if stripped_name != name.replace("_", " "):  # Only try if we actually stripped something
                try:
                    stripped_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=stripped_name,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if stripped_result.returncode == 0 and stripped_result.stdout.strip():
                        stripped_parts = stripped_result.stdout.strip().split('\t')
                        if len(stripped_parts) >= 2 and stripped_parts[1] != "0" and stripped_parts[1].strip():
                            # Strip the taxid to remove any newlines or extra whitespace
                            results[name] = (stripped_parts[1].strip(), "number_stripped")
                except Exception:
                    pass

    # TIER 4: For names that still didn't get a match, try extracting from hyphenated patterns as final fallback
    for name in taxon_batch:
        if name not in results:
            # Try extracting taxa from hyphenated names
            extracted_name = extract_taxa_from_hyphenated(name)
            if extracted_name:  # Only try if we actually extracted something
                try:
                    extracted_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=extracted_name.replace("_", " "),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if extracted_result.returncode == 0 and extracted_result.stdout.strip():
                        extracted_parts = extracted_result.stdout.strip().split('\t')
                        if len(extracted_parts) >= 2 and extracted_parts[1] != "0" and extracted_parts[1].strip():
                            # Strip the taxid to remove any newlines or extra whitespace
                            results[name] = (extracted_parts[1].strip(), "hyphenated_extracted")
                except Exception:
                    pass

    return results


def get_taxids_for_names(taxon_names: List[str], rank_name: str = "unknown") -> Tuple[Dict[str, str], Dict[str, dict]]:
    """
    Get NCBI taxids for a list of taxon names using taxonkit.
    Handles underscore removal and genus fallback.
    Uses parallel processing for better performance.

    Args:
        taxon_names: List of taxon names
        rank_name: Name of the taxonomic rank being processed (for logging)

    Returns:
        Tuple of (results_dict, failed_names_dict) where:
        - results_dict: Dictionary mapping taxon names to their taxids
        - failed_names_dict: Dictionary of failed names with failure details
    """
    if not taxon_names:
        return {}, {}

    logging.info(f"Getting taxids for {len(taxon_names)} unique names...")

    # Determine the number of workers and batch size
    num_workers = min(os.cpu_count() or 4, 8)  # Use at most 8 workers
    batch_size = max(1, math.ceil(len(taxon_names) / (num_workers * 4)))  # Ensure enough batches

    # Split the taxon names into batches
    batches = [taxon_names[i:i + batch_size] for i in range(0, len(taxon_names), batch_size)]

    # Process batches in parallel
    results = {}
    failed_names = {}
    direct_match_count = 0
    genus_fallback_count = 0
    number_stripped_count = 0
    hyphenated_extracted_count = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all batches
        future_to_batch = {executor.submit(process_taxon_batch, batch): i for i, batch in enumerate(batches)}

        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_batch):
            batch_results = future.result()

            # Count match types
            for name, (taxid, method) in batch_results.items():
                results[name] = taxid
                if method == "direct":
                    direct_match_count += 1
                elif method == "genus_fallback":
                    genus_fallback_count += 1
                elif method == "number_stripped":
                    number_stripped_count += 1
                elif method == "hyphenated_extracted":
                    hyphenated_extracted_count += 1

    # Collect failed names and add "NA" for names that didn't get a match
    for name in taxon_names:
        if name not in results:
            results[name] = "NA"
            failed_names[name] = {
                'type': 'NO_TAXID_FOUND',
                'details': f'No taxid found for {rank_name} name after direct and genus fallback attempts',
                'taxid': 'NA',
                'rank': rank_name
            }

    # Log statistics
    total = len(taxon_names)
    matched = direct_match_count + genus_fallback_count + number_stripped_count + hyphenated_extracted_count

    logging.info(f"Taxid matching complete: {matched}/{total} matched ({matched/total*100:.1f}%)")

    return results, failed_names


def get_lineages_for_taxids(taxids: List[str], env: dict) -> Dict[str, Tuple[str, str, str]]:
    """
    Get lineages for a list of taxids using taxonkit with temporary file approach.
    Returns lineage, ranks, and taxids in a single pass.

    Args:
        taxids: List of taxids
        env: Environment variables for subprocess

    Returns:
        Dictionary mapping taxids to a tuple of (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}

    lineage_data = {}

    # Create temporary file for taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        # Write taxids to temporary file
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        # Run taxonkit lineage command with -R flag for ranks and -t flag for taxids
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')

            for line in lines:
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 4:  # Expecting taxid, lineage, lineage_taxids, lineage_ranks
                    taxid = parts[0]
                    lineage = parts[1]  # The lineage is in the second column
                    lineage_taxids = parts[2]  # The taxids are in the third column
                    lineage_ranks = parts[3]  # The ranks are in the fourth column

                    # Clean up the data
                    if lineage and lineage != taxid:
                        # Clean and format the lineage components
                        clean_lineage = lineage.strip()
                        clean_ranks = lineage_ranks.strip()
                        clean_taxids = lineage_taxids.strip()

                        # Store the cleaned data
                        lineage_data[taxid] = (
                            clean_lineage,
                            clean_ranks,
                            clean_taxids
                        )

    except Exception as e:
        logging.error(f"Error getting lineages: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except Exception as e:
            logging.warning(f"Error cleaning up temporary file: {e}")

    return lineage_data

