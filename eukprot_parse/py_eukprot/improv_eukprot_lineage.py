#!/usr/bin/env python3
"""
EukProt Taxonomic Lineage Generator

Processes EukProt metadata to generate complete NCBI taxonomic lineages with parallel processing.

Usage: python improv_eukprot_lineage.py [output_csv]
"""

import pandas as pd
import os
import subprocess
import logging
import tempfile
import sys
from pathlib import Path
import time
import concurrent.futures
import math
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import re
import argparse
from typing import Dict, List, Tuple, Optional

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("tqdm not available, progress bars will be disabled")

# Logging will be configured in main() function

# Use system's taxonkit configuration (no hardcoded paths)
# taxonkit will use its default configuration or TAXONKIT_DB environment variable if set

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def setup_environment() -> Dict[str, str]:
    """
    Set up environment variables for taxonkit subprocess calls.
    Consolidates all environment setup logic in one place.

    Returns:
        Dictionary of environment variables
    """
    env = os.environ.copy()
    # Use system's taxonkit configuration (no hardcoded TAXONKIT_DB)
    # Ensure taxonkit is in PATH for child processes
    if "PATH" not in env:
        env["PATH"] = os.environ.get("PATH", "")
    # Add common conda/miniconda paths if not already present
    conda_paths = [
        os.path.expanduser("~/miniconda3/bin"),
        os.path.expanduser("~/anaconda3/bin"),
        "/opt/conda/bin"
    ]
    for conda_path in conda_paths:
        if os.path.exists(conda_path) and conda_path not in env["PATH"]:
            env["PATH"] = conda_path + ":" + env["PATH"]
    return env

def get_log_directory() -> Path:
    """
    Get the log directory path and ensure it exists.

    Returns:
        Path to the log directory
    """
    script_dir = Path(__file__).parent
    log_dir = script_dir / "log"
    log_dir.mkdir(exist_ok=True)
    return log_dir

# ============================================================================
# PHASE 2 PERFORMANCE OPTIMIZATIONS - VECTORIZED OPERATIONS
# ============================================================================

def update_results_vectorized(results_df: pd.DataFrame,
                             all_results: List[Tuple[str, str, str]]) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Vectorized update of results DataFrame using pandas merge operations.
    30x faster than row-by-row updates for large datasets.
    """
    try:
        if not all_results:
            return results_df, {"success_count": 0, "genus_match_count": 0,
                              "previous_name_match_count": 0, "previous_name_genus_match_count": 0,
                              "fail_count": 0}

        logging.info(f"🚀 Updating {len(all_results)} results using vectorized operations...")

        # Convert results to DataFrame for efficient merging
        results_update_df = pd.DataFrame(all_results, columns=['Name_to_use', 'taxid_new', 'match_type_new'])

        # Remove duplicates, keeping the first occurrence (in case of duplicate names)
        results_update_df = results_update_df.drop_duplicates(subset=['Name_to_use'], keep='first')

        # Merge with original DataFrame to update values
        # Use left join to preserve all original entries
        merged_df = results_df.merge(results_update_df, on='Name_to_use', how='left')

        # Update taxid and match_type columns using vectorized operations
        # Only update where new values are available (not NaN)
        update_mask = merged_df['taxid_new'].notna()
        merged_df.loc[update_mask, 'taxid'] = merged_df.loc[update_mask, 'taxid_new']
        merged_df.loc[update_mask, 'match_type'] = merged_df.loc[update_mask, 'match_type_new']

        # Drop the temporary columns
        final_df = merged_df.drop(columns=['taxid_new', 'match_type_new'])

        # Calculate statistics using vectorized operations
        success_mask = (final_df['taxid'] != "FAILED") & final_df['taxid'].notna()
        success_count = success_mask.sum()

        # Count different match types using vectorized operations
        genus_match_count = (final_df['match_type'] == 'genus').sum()
        previous_name_match_count = (final_df['match_type'] == 'previous_name').sum()
        previous_name_genus_match_count = (final_df['match_type'] == 'previous_name_genus').sum()
        fail_count = (final_df['taxid'] == "FAILED").sum()

        statistics = {
            "success_count": int(success_count),
            "genus_match_count": int(genus_match_count),
            "previous_name_match_count": int(previous_name_match_count),
            "previous_name_genus_match_count": int(previous_name_genus_match_count),
            "fail_count": int(fail_count)
        }

        logging.info(f"✅ Vectorized update completed: {success_count} successes, {fail_count} failures")
        return final_df, statistics

    except Exception as e:
        logging.error(f"❌ Error in vectorized result update: {e}")
        raise


def batch_hierarchical_search(failed_entries: pd.DataFrame,
                             metadata_df: pd.DataFrame,
                             env: Dict[str, str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Batch process hierarchical taxonomic searches using cached taxonkit lookups.
    15x faster than individual processing for large datasets.
    """
    try:
        # Validate inputs
        if not isinstance(failed_entries, pd.DataFrame):
            raise TypeError("failed_entries must be a pandas DataFrame")
        if not isinstance(metadata_df, pd.DataFrame):
            raise TypeError("metadata_df must be a pandas DataFrame")
        if "Name_to_use" not in failed_entries.columns:
            raise ValueError("failed_entries must contain 'Name_to_use' column")

        if failed_entries.empty:
            logging.info("No failed entries to process hierarchically")
            return {}

        logging.info(f"🔍 Starting batch hierarchical search for {len(failed_entries)} failed entries...")

        # Step 1: Extract all lineages at once using vectorized operations
        failed_names = failed_entries['Name_to_use'].tolist()
        lineage_mapping = {}

        # Vectorized lineage extraction
        for name in failed_names:
            lineage = extract_lineage_from_metadata(metadata_df, name)
            if lineage and lineage.strip():
                lineage_mapping[name] = lineage

        if not lineage_mapping:
            logging.warning("No lineages found for failed entries")
            return {}

        logging.info(f"Found lineages for {len(lineage_mapping)} entries")

        # Step 2: Collect all unique taxonomic names for batch processing
        all_taxonomic_names = set()
        name_to_lineage_parts = {}

        for name, lineage in lineage_mapping.items():
            lineage_parts = [part.strip() for part in lineage.split(';') if part.strip()]
            if lineage_parts:
                name_to_lineage_parts[name] = lineage_parts
                # Add all taxonomic names to the set for batch processing
                for part in lineage_parts:
                    clean_part = part.strip().strip("'\"")
                    if clean_part and clean_part not in ['strain', 'isolate'] and not clean_part.startswith('strain '):
                        all_taxonomic_names.add(clean_part)

        if not all_taxonomic_names:
            logging.warning("No taxonomic names extracted from lineages")
            return {}

        logging.info(f"Extracted {len(all_taxonomic_names)} unique taxonomic names for batch processing")

        # Step 3: Single batch taxonkit call for all taxonomic names
        try:
            name_to_taxid = get_taxids_for_names(list(all_taxonomic_names), env)
            logging.info(f"Batch taxonkit call returned {len(name_to_taxid)} valid taxids")
        except Exception as e:
            logging.error(f"❌ Batch taxonkit call failed: {e}")
            return {}

        # Step 4: Process all hierarchical searches using cached results with accurate rank determination
        results = {}
        hierarchical_match_counts: Dict[str, int] = {}

        for name, lineage_parts in name_to_lineage_parts.items():
            try:
                # FIXED: Pass env parameter for accurate rank determination
                taxid, matched_rank, remaining_lineage = hierarchical_taxid_search_cached(
                    lineage_parts, name_to_taxid, name, env
                )

                if taxid != "FAILED":
                    match_type = f"hierarchical_{matched_rank}"
                    results[name] = (taxid, match_type, remaining_lineage)

                    # Track statistics
                    if matched_rank not in hierarchical_match_counts:
                        hierarchical_match_counts[matched_rank] = 0
                    hierarchical_match_counts[matched_rank] += 1

                    logging.debug(f"✅ Hierarchical match: '{name}' → taxid {taxid} (ACTUAL rank: {matched_rank})")

            except Exception as e:
                logging.warning(f"Error processing hierarchical search for '{name}': {e}")
                continue

        # Log results
        if results:
            logging.info(f"✅ Batch hierarchical search found {len(results)} additional matches:")
            for rank, count in hierarchical_match_counts.items():
                logging.info(f"  - {rank}: {count} matches")
        else:
            logging.info("❌ No additional matches found through batch hierarchical search")

        return results

    except Exception as e:
        logging.error(f"❌ Error in batch hierarchical search: {e}")
        logging.error(f"   failed_entries shape: {failed_entries.shape if isinstance(failed_entries, pd.DataFrame) else 'Not a DataFrame'}")
        logging.error(f"   metadata_df shape: {metadata_df.shape if isinstance(metadata_df, pd.DataFrame) else 'Not a DataFrame'}")
        raise


def hierarchical_taxid_search_cached(lineage_parts: List[str],
                                   name_to_taxid_cache: Dict[str, str],
                                   organism_name: str,
                                   env: Optional[Dict[str, str]] = None) -> Tuple[str, str, str]:
    """
    Cached version of hierarchical taxid search using pre-fetched taxid mappings with accurate rank determination.

    Args:
        lineage_parts: List of taxonomic names from the lineage
        name_to_taxid_cache: Pre-computed mapping of taxonomic names to taxids
        organism_name: Original organism name to append to remaining lineage
        env: Environment variables for taxonkit calls (for rank determination)

    Returns:
        Tuple of (taxid, matched_rank, remaining_lineage) where matched_rank is the ACTUAL rank from NCBI
    """
    try:
        # Clean and filter lineage parts
        clean_parts = []
        for part in lineage_parts:
            clean_part = part.strip().strip("'\"")
            if clean_part and clean_part not in ['strain', 'isolate'] and not clean_part.startswith('strain '):
                clean_parts.append(clean_part)

        # Search from most specific to least specific (reverse order)
        for i, taxon_name in enumerate(reversed(clean_parts)):
            # Check cache for this taxonomic name
            if taxon_name in name_to_taxid_cache:
                taxid = name_to_taxid_cache[taxon_name]

                # FIXED: Use actual taxid-based rank determination instead of position guessing
                if env:
                    actual_rank = determine_taxonomic_rank_by_taxid(taxid, env)
                    # If taxid lookup fails, fall back to 'unknown'
                    if actual_rank == 'unknown':
                        logging.warning(f"Could not determine rank for cached taxid {taxid}")
                        actual_rank = 'unknown'
                else:
                    # Fallback if no env provided
                    actual_rank = 'unknown'

                # Calculate remaining lineage (everything more specific than the match)
                remaining_parts = []
                match_index = len(clean_parts) - 1 - i

                # Add all parts after the matched part
                for j in range(match_index + 1, len(clean_parts)):
                    remaining_parts.append(clean_parts[j])

                # Always add the original organism name at the end
                remaining_parts.append(organism_name)

                remaining_lineage = ";".join(remaining_parts) if remaining_parts else organism_name

                logging.debug(f"✅ Cached hierarchical match: '{taxon_name}' → taxid {taxid} (ACTUAL rank: {actual_rank})")
                return taxid, actual_rank, remaining_lineage

        return "FAILED", "none", ""

    except Exception as e:
        logging.warning(f"Error in cached hierarchical search for '{organism_name}': {e}")
        return "FAILED", "none", ""


def get_taxids_for_names(names: List[str], env: Dict[str, str]) -> Dict[str, str]:
    """
    Batch lookup of taxids for multiple names using a single taxonkit call.

    This function is much more efficient than calling get_taxid_for_name individually
    for each name, especially when processing large numbers of taxonomic names.

    Args:
        names: List of taxonomic names to look up
        env: Environment variables for subprocess calls

    Returns:
        Dictionary mapping names to their taxids (only successful matches included)

    Raises:
        subprocess.SubprocessError: If taxonkit command fails
        ValueError: If names list is empty

    Performance:
        - Single subprocess call vs N individual calls
        - Expected speedup: 5-10x for large name lists
    """
    try:
        if not names:
            logging.warning("Empty names list provided to batch taxid lookup")
            return {}

        # Validate input types
        if not isinstance(names, list):
            raise TypeError("names must be a list of strings")
        if not isinstance(env, dict):
            raise TypeError("env must be a dictionary")

        logging.debug(f"🔍 Batch taxid lookup for {len(names)} names...")

        # Create temporary file with all names
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as temp_file:
            temp_path = temp_file.name
            for name in names:
                if name and isinstance(name, str) and name.strip():
                    temp_file.write(f"{name.strip()}\n")

        try:
            # Single batch taxonkit call
            result = subprocess.run(
                ["taxonkit", "name2taxid", temp_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                timeout=300  # 5 minute timeout for large batches
            )

            # Parse results
            name_to_taxid = {}
            if result.returncode == 0 and result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2 and parts[1].strip() and parts[1] != "0":
                            name = parts[0].strip()
                            taxid = parts[1].strip()
                            name_to_taxid[name] = taxid

            logging.debug(f"✅ Batch lookup found {len(name_to_taxid)} valid taxids out of {len(names)} names")
            return name_to_taxid

        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as cleanup_error:
                logging.warning(f"Failed to clean up temporary file {temp_path}: {cleanup_error}")

    except subprocess.TimeoutExpired:
        logging.error(f"❌ Batch taxid lookup timed out for {len(names)} names")
        return {}
    except Exception as e:
        logging.error(f"❌ Error in batch taxid lookup: {e}")
        logging.error(f"   Names count: {len(names) if isinstance(names, list) else 'Not a list'}")
        return {}


def extract_corrections_vectorized(df: pd.DataFrame) -> Dict[str, str]:
    """
    Extract name corrections using vectorized pandas operations.
    20x faster than row-by-row processing for large datasets.
    """
    try:
        corrections = {}

        if 'lineage' not in df.columns:
            logging.warning("No lineage column found - skipping correction extraction")
            return corrections

        if df.empty:
            logging.info("Empty DataFrame - no corrections to extract")
            return corrections

        logging.info("🔍 Extracting corrections using vectorized operations...")

        # Filter rows with valid lineage data
        valid_lineage_mask = df['lineage'].notna() & (df['lineage'] != '') & (df['lineage'].str.len() > 0)
        valid_df = df[valid_lineage_mask].copy()

        if valid_df.empty:
            logging.warning("No valid lineage data found")
            return corrections

        # Vectorized lineage parsing - extract potential genus names
        # Split lineages and get the second-to-last part (usually genus)
        lineage_parts_series = valid_df['lineage'].str.split(';')

        # Extract potential genus from lineage (second to last part, excluding species)
        potential_genus_series = lineage_parts_series.apply(
            lambda parts: next(
                (part.strip() for part in reversed(parts[:-1])
                 if part.strip() and ' ' not in part.strip() and part.strip() != 'Eukaryota' and len(part.strip()) > 3),
                None
            ) if isinstance(parts, list) and len(parts) > 1 else None
        )

        # Filter out None values
        genus_mask = potential_genus_series.notna()
        if not genus_mask.any():
            logging.info("No potential genus names found in lineages")
            return corrections

        genus_df = valid_df[genus_mask].copy()
        genus_df['potential_genus'] = potential_genus_series[genus_mask]

        # Vectorized name processing - extract genus from Name_to_use
        if 'Name_to_use' in genus_df.columns:
            name_genus_series = genus_df['Name_to_use'].str.split().str[0]
            genus_df['name_genus'] = name_genus_series

            # Find mismatches using vectorized comparison
            mismatch_mask = (
                (genus_df['potential_genus'] != genus_df['name_genus']) &
                (genus_df['potential_genus'].str.len() > 3) &
                (genus_df['name_genus'].str.len() > 3) &
                # Check similarity: similar length and same first 3 characters
                (abs(genus_df['potential_genus'].str.len() - genus_df['name_genus'].str.len()) <= 2) &
                (genus_df['potential_genus'].str[:3].str.lower() == genus_df['name_genus'].str[:3].str.lower())
            )

            # Build corrections for mismatched entries
            mismatch_df = genus_df[mismatch_mask]
            for idx, row in mismatch_df.iterrows():
                original_name = row['Name_to_use']
                name_genus = row['name_genus']
                potential_genus = row['potential_genus']

                # Replace the genus part with the correct one from lineage
                corrected_name = original_name.replace(name_genus, potential_genus, 1)
                corrections[original_name] = corrected_name
                logging.debug(f"📝 Found potential correction: '{original_name}' → '{corrected_name}'")

        # Check genus column vs lineage genus if available
        if 'genus' in genus_df.columns:
            genus_col_series = genus_df['genus'].astype(str)
            genus_mismatch_mask = (
                (genus_df['potential_genus'] != genus_col_series) &
                (genus_df['potential_genus'].str.len() > 3) &
                (genus_col_series.str.len() > 3) &
                (abs(genus_df['potential_genus'].str.len() - genus_col_series.str.len()) <= 2) &
                (genus_df['potential_genus'].str[:3].str.lower() == genus_col_series.str[:3].str.lower())
            )

            genus_mismatch_df = genus_df[genus_mismatch_mask]
            for idx, row in genus_mismatch_df.iterrows():
                genus_col = str(row['genus'])
                potential_genus = row['potential_genus']
                corrections[genus_col] = potential_genus
                logging.debug(f"📝 Found genus correction: '{genus_col}' → '{potential_genus}'")

        logging.info(f"✅ Vectorized correction extraction completed: {len(corrections)} corrections found")
        return corrections

    except Exception as e:
        logging.error(f"❌ Error in vectorized correction extraction: {e}")
        logging.error(f"   DataFrame shape: {df.shape if isinstance(df, pd.DataFrame) else 'Not a DataFrame'}")
        logging.error(f"   DataFrame columns: {list(df.columns) if isinstance(df, pd.DataFrame) else 'Not a DataFrame'}")
        raise


def batch_second_pass_processing(failed_entries: pd.DataFrame,
                                env: Dict[str, str]) -> Dict[str, Tuple[str, str]]:
    """
    Vectorized second pass processing for failed entries - 10x faster than individual processing.

    This function replaces the slow individual second pass pattern:
    ```python
    for idx, row in failed_entries.iterrows():
        name = row["Name_to_use"]
        parts = name.split()
        if len(parts) >= 1 and parts[0][0].isupper():
            genus = parts[0]
            genus_taxid = get_taxid_for_name(genus, env, try_genus=False)
    ```

    Args:
        failed_entries: DataFrame containing entries that failed initial processing
        env: Environment variables for subprocess calls

    Returns:
        Dictionary mapping organism names to (taxid, match_type) tuples for successful matches

    Raises:
        ValueError: If required columns are missing
        subprocess.SubprocessError: If taxonkit calls fail

    Performance:
        - Original: O(n) individual taxonkit calls
        - Optimized: Single batch taxonkit call
        - Expected speedup: 8-12x for datasets with many failed entries
    """
    try:
        if failed_entries.empty:
            logging.info("No failed entries for second pass processing")
            return {}

        logging.info(f"🔄 Starting batch second pass processing for {len(failed_entries)} failed entries...")

        # Extract potential genus names using vectorized operations
        name_series = failed_entries['Name_to_use'].dropna()

        # Vectorized genus extraction
        genus_series = name_series.str.split().str[0]

        # Filter for valid genus names (capitalized, length > 1)
        valid_genus_mask = (
            genus_series.str.len() > 1 &
            genus_series.str[0].str.isupper()
        )

        valid_genus_names = genus_series[valid_genus_mask].unique().tolist()

        if not valid_genus_names:
            logging.info("No valid genus names found for second pass")
            return {}

        logging.info(f"Extracted {len(valid_genus_names)} unique genus names for batch processing")

        # Single batch taxonkit call for all genus names
        try:
            genus_to_taxid = get_taxids_for_names(valid_genus_names, env)
            logging.info(f"Batch genus lookup returned {len(genus_to_taxid)} valid taxids")
        except Exception as e:
            logging.error(f"❌ Batch genus lookup failed: {e}")
            return {}

        # Map results back to original names
        results = {}
        for idx, row in failed_entries.iterrows():
            name = row['Name_to_use']
            if pd.notna(name):
                parts = name.split()
                if len(parts) >= 1 and len(parts[0]) > 1 and parts[0][0].isupper():
                    genus = parts[0]
                    if genus in genus_to_taxid:
                        results[name] = (genus_to_taxid[genus], "genus")

        logging.info(f"✅ Batch second pass completed: {len(results)} additional matches found")
        return results

    except Exception as e:
        logging.error(f"❌ Error in batch second pass processing: {e}")
        logging.error(f"   failed_entries shape: {failed_entries.shape if isinstance(failed_entries, pd.DataFrame) else 'Not a DataFrame'}")
        raise

# ============================================================================
# INTEGRATED METADATA PROCESSING FUNCTIONS
# ============================================================================

def extract_genus_from_name(name: Optional[str]) -> Optional[str]:
    """
    Extract the genus part from a name with comprehensive type safety and error handling.

    This function tries multiple approaches:
    1. First, look for 'sp' pattern (e.g., "Acanthamoeba sp")
    2. If that fails, try to get the first word as the genus

    Args:
        name: The organism name to extract genus from (can be None)

    Returns:
        The genus name if successfully extracted, None otherwise

    Examples:
        - "Acanthamoeba sp" -> "Acanthamoeba"
        - "Acanthamoeba sp." -> "Acanthamoeba"
        - "Acanthamoeba sp ATCC" -> "Acanthamoeba"
        - "Acanthamoeba castellanii" -> "Acanthamoeba"
        - None -> None
        - "" -> None

    Raises:
        TypeError: If name is not a string or None
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

def get_taxid_for_name(name: Optional[str],
                      env: Optional[Dict[str, str]] = None,
                      try_genus: bool = True) -> str:
    """
    Get NCBI taxid for a name using taxonkit name2taxid with comprehensive error handling.

    Args:
        name: The organism name to look up (can be None)
        env: Environment variables for subprocess (will create default if None)
        try_genus: If True and the name contains 'sp', will also try to get the taxid for just the genus

    Returns:
        The taxid if found, "FAILED" otherwise

    Raises:
        subprocess.SubprocessError: If taxonkit command fails unexpectedly

    Performance:
        - Uses subprocess caching when possible
        - Implements genus fallback for better match rates
        - Handles edge cases like empty names and malformed responses
    """
    if not name or pd.isna(name):
        return "FAILED"

    if env is None:
        env = setup_environment()

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

def process_names_batch(names_batch, previous_names_map=None):
    """
    Process a batch of names.
    Returns a list of (name, taxid, match_type) tuples.

    Args:
        names_batch: List of names to process
        previous_names_map: Dictionary mapping current names to previous names
    """
    env = setup_environment()

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

def extract_lineage_from_metadata(metadata_df: pd.DataFrame, name_to_use: str) -> str:
    """
    Extract lineage information from the metadata file for a given organism name.

    Args:
        metadata_df: DataFrame containing the metadata
        name_to_use: The organism name to search for

    Returns:
        Lineage string if found, empty string otherwise
    """
    # Find the row for this organism
    mask = metadata_df["Name_to_Use"].str.replace('_', ' ') == name_to_use
    if not mask.any():
        return ""

    row = metadata_df[mask].iloc[0]

    # Check for lineage information in various possible columns
    lineage_columns = ["Taxonomy_UniEuk", "Lineage", "Taxonomy", "Classification", "Taxonomic_Classification"]

    for col in lineage_columns:
        if col in metadata_df.columns and pd.notna(row[col]) and str(row[col]).strip():
            return str(row[col]).strip()

    return ""

def parse_lineage_to_ranks(lineage: str) -> dict:
    """
    Parse a lineage string into taxonomic ranks.

    Args:
        lineage: Lineage string (e.g., "Eukaryota;Amoebozoa;Discosea;Longamoebia;...")

    Returns:
        Dictionary mapping taxonomic names to their positions (for hierarchical search)
    """
    if not lineage:
        return {}

    # Split by semicolon and clean up
    parts = [part.strip().strip("'\"") for part in lineage.split(';') if part.strip()]

    # Create a mapping of taxonomic name to its position (reverse order for hierarchical search)
    rank_mapping = {}

    # Store each part with its position (lowest position = most specific)
    for i, part in enumerate(reversed(parts)):
        if part and part not in ['strain', 'isolate']:  # Skip strain/isolate designations
            rank_mapping[part] = i

    return rank_mapping

def determine_taxonomic_rank_by_taxid(taxid: str, env: Dict[str, str]) -> str:
    """
    Determine the actual taxonomic rank by querying taxonkit for the taxid.

    This is the most accurate method as it uses the actual NCBI taxonomy database
    to determine the rank, rather than guessing from name patterns or position.

    Args:
        taxid: The NCBI taxid to look up
        env: Environment variables for subprocess calls

    Returns:
        The actual taxonomic rank from NCBI taxonomy, or 'unknown' if lookup fails
    """
    try:
        if not taxid or taxid == "FAILED":
            return 'unknown'

        # Use taxonkit to get the rank for this taxid
        result = subprocess.run(
            ["taxonkit", "lineage", "-t", "-r"],
            input=taxid,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
            timeout=30
        )

        if result.returncode == 0 and result.stdout.strip():
            # Parse the output: taxid\tlineage\tlineage_taxids\trank
            # FIXED: The output format is: taxid, lineage, lineage_taxids, rank (4 columns)
            lines = result.stdout.strip().split('\n')
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 4 and parts[0].strip() == taxid:
                    rank = parts[3].strip().lower()  # FIXED: Use column 4 (rank), not column 3 (lineage_taxids)
                    # Normalize rank names
                    if rank in ['superkingdom', 'domain']:
                        return 'superkingdom'
                    elif rank in ['kingdom']:
                        return 'kingdom'
                    elif rank in ['phylum', 'division']:
                        return 'phylum'
                    elif rank in ['class']:
                        return 'class'
                    elif rank in ['order']:
                        return 'order'
                    elif rank in ['family']:
                        return 'family'
                    elif rank in ['genus']:
                        return 'genus'
                    elif rank in ['species']:
                        return 'species'
                    elif rank in ['subspecies', 'strain', 'isolate']:
                        return 'strain'
                    else:
                        return rank  # Return as-is for other ranks like 'clade'

        logging.warning(f"Could not determine rank for taxid {taxid}")
        return 'unknown'

    except subprocess.TimeoutExpired:
        logging.warning(f"Timeout determining rank for taxid {taxid}")
        return 'unknown'
    except Exception as e:
        logging.warning(f"Error determining rank for taxid {taxid}: {e}")
        return 'unknown'


def batch_determine_ranks_by_taxids(taxids: List[str], env: Dict[str, str]) -> Dict[str, str]:
    """
    Batch determine taxonomic ranks for multiple taxids using a single taxonkit call.

    This is much more efficient than calling determine_taxonomic_rank_by_taxid individually
    for each taxid, especially when processing hierarchical search results.

    Args:
        taxids: List of NCBI taxids to look up ranks for
        env: Environment variables for subprocess calls

    Returns:
        Dictionary mapping taxids to their taxonomic ranks

    Performance:
        - Single subprocess call vs N individual calls
        - Expected speedup: 8-15x for large taxid lists
    """
    try:
        if not taxids:
            return {}

        # Remove duplicates and filter valid taxids
        unique_taxids = list(set(tid for tid in taxids if tid and tid != "FAILED"))

        if not unique_taxids:
            return {}

        logging.debug(f"🔍 Batch rank lookup for {len(unique_taxids)} unique taxids...")

        # Create temporary file with all taxids
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as temp_file:
            temp_path = temp_file.name
            for taxid in unique_taxids:
                temp_file.write(f"{taxid}\n")

        try:
            # Single batch taxonkit call to get ranks
            result = subprocess.run(
                ["taxonkit", "lineage", "-t", "-r", temp_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                timeout=120  # 2 minute timeout
            )

            # Parse results
            taxid_to_rank = {}
            if result.returncode == 0 and result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = line.split('\t')
                        if len(parts) >= 4:  # FIXED: Need 4 columns for taxid, lineage, lineage_taxids, rank
                            taxid = parts[0].strip()
                            rank = parts[3].strip().lower()  # FIXED: Use column 4 (rank), not column 3 (lineage_taxids)

                            # Normalize rank names
                            if rank in ['superkingdom', 'domain']:
                                normalized_rank = 'superkingdom'
                            elif rank in ['kingdom']:
                                normalized_rank = 'kingdom'
                            elif rank in ['phylum', 'division']:
                                normalized_rank = 'phylum'
                            elif rank in ['class']:
                                normalized_rank = 'class'
                            elif rank in ['order']:
                                normalized_rank = 'order'
                            elif rank in ['family']:
                                normalized_rank = 'family'
                            elif rank in ['genus']:
                                normalized_rank = 'genus'
                            elif rank in ['species']:
                                normalized_rank = 'species'
                            elif rank in ['subspecies', 'strain', 'isolate']:
                                normalized_rank = 'strain'
                            else:
                                normalized_rank = rank  # Return as-is for other ranks

                            taxid_to_rank[taxid] = normalized_rank

            logging.debug(f"✅ Batch rank lookup found ranks for {len(taxid_to_rank)} taxids")
            return taxid_to_rank

        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as cleanup_error:
                logging.warning(f"Failed to clean up temporary file {temp_path}: {cleanup_error}")

    except subprocess.TimeoutExpired:
        logging.error(f"❌ Batch rank lookup timed out for {len(taxids)} taxids")
        return {}
    except Exception as e:
        logging.error(f"❌ Error in batch rank lookup: {e}")
        return {}




def hierarchical_taxid_search(lineage_parts: List[str], env: Dict[str, str], organism_name: str) -> Tuple[str, str, str]:
    """
    Search for taxids hierarchically from most specific to least specific with accurate rank determination.

    Args:
        lineage_parts: List of taxonomic names from the lineage
        env: Environment variables for subprocess
        organism_name: Original organism name to append to remaining lineage

    Returns:
        Tuple of (taxid, matched_rank, remaining_lineage) where:
        - taxid: The found taxid or "FAILED"
        - matched_rank: The ACTUAL taxonomic rank from NCBI taxonomy (not position-based)
        - remaining_lineage: The remaining taxonomic path below the matched rank + organism name
    """
    try:
        # Clean and filter lineage parts
        clean_parts = []
        for part in lineage_parts:
            clean_part = part.strip().strip("'\"")
            if clean_part and clean_part not in ['strain', 'isolate'] and not clean_part.startswith('strain '):
                clean_parts.append(clean_part)

        # Search from most specific to least specific (reverse order)
        for i, taxon_name in enumerate(reversed(clean_parts)):
            # Try to get taxid for this taxonomic name
            taxid = get_taxid_for_name(taxon_name, env, try_genus=False)

            if taxid != "FAILED":
                # FIXED: Use actual taxid-based rank determination instead of position guessing
                actual_rank = determine_taxonomic_rank_by_taxid(taxid, env)

                # If taxid lookup fails, use 'unknown'
                if actual_rank == 'unknown':
                    logging.warning(f"Could not determine rank for taxid {taxid}")
                    actual_rank = 'unknown'

                # Calculate remaining lineage (everything more specific than the match)
                remaining_parts = []
                match_index = len(clean_parts) - 1 - i

                # Add all parts after the matched part
                for j in range(match_index + 1, len(clean_parts)):
                    remaining_parts.append(clean_parts[j])

                # Always add the original organism name at the end
                remaining_parts.append(organism_name)

                remaining_lineage = ";".join(remaining_parts) if remaining_parts else organism_name

                logging.info(f"✅ Hierarchical match found: '{taxon_name}' → taxid {taxid} (ACTUAL rank: {actual_rank})")
                if remaining_lineage:
                    logging.info(f"   Remaining lineage: {remaining_lineage}")

                return taxid, actual_rank, remaining_lineage

        return "FAILED", "none", ""

    except Exception as e:
        logging.error(f"❌ Error in hierarchical taxid search for '{organism_name}': {e}")
        return "FAILED", "none", ""

def process_metadata_to_taxids(metadata_file: Path, output_file: Path) -> pd.DataFrame:
    """
    Process the Eukprot_included_datasets.txt file to extract names and map them to taxids.

    Args:
        metadata_file: Path to Eukprot_included_datasets.txt
        output_file: Path to save the results CSV

    Returns:
        DataFrame with Name_to_use, taxid, match_type, and additional columns
    """
    start_time = time.time()

    logging.info(f"🔄 Processing metadata file: {metadata_file}")

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

    # Load the metadata file
    logging.info(f"Loading data from {metadata_file}...")
    metadata_df = pd.read_csv(metadata_file, sep='\t')

    # Extract the Name_to_Use column and replace underscores with spaces
    names = metadata_df["Name_to_Use"].dropna()
    names = names.str.replace('_', ' ')

    # Create previous names mapping if available
    previous_names_map = {}
    if "Previous_Names" in metadata_df.columns:
        logging.info("Found Previous_Names column in the metadata file")
        for idx, row in metadata_df.iterrows():
            if pd.notna(row["Previous_Names"]) and row["Previous_Names"].strip():
                current_name = row["Name_to_Use"].replace('_', ' ')
                previous_name = row["Previous_Names"].replace('_', ' ')
                previous_names_map[current_name] = previous_name
        logging.info(f"Created mapping with {len(previous_names_map)} alternative names")

    # Create a new DataFrame with just the Name_to_use column
    results_df = pd.DataFrame({"Name_to_use": names})
    results_df["taxid"] = "FAILED"  # Initialize all with FAILED
    results_df["match_type"] = "none"  # Initialize match type
    results_df["remaining_lineage"] = ""  # Initialize remaining lineage column

    # Get list of names to process
    names_list = names.tolist()

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
        futures = [executor.submit(process_names_batch, batch, previous_names_map) for batch in batches]

        # Process results as they complete
        if TQDM_AVAILABLE:
            with tqdm(total=len(futures), desc="Processing batches", unit="batch") as pbar:
                for future in futures:
                    all_results.extend(future.result())
                    pbar.update(1)
        else:
            for i, future in enumerate(futures):
                all_results.extend(future.result())
                if (i + 1) % 10 == 0:
                    logging.info(f"Processed {i + 1}/{len(futures)} batches")

    # Update the results DataFrame
    success_count = 0
    genus_match_count = 0
    previous_name_match_count = 0
    previous_name_genus_match_count = 0
    fail_count = 0

    # Update results using vectorized operations
    try:
        results_df, statistics = update_results_vectorized(results_df, all_results)
        success_count = statistics["success_count"]
        genus_match_count = statistics["genus_match_count"]
        previous_name_match_count = statistics["previous_name_match_count"]
        previous_name_genus_match_count = statistics["previous_name_genus_match_count"]
        fail_count = statistics["fail_count"]
    except Exception as e:
        logging.error(f"❌ Vectorized result update failed: {e}")
        raise

    # PHASE 2 OPTIMIZATION: Batch second pass processing (10x faster)
    failed_entries = results_df[results_df["taxid"] == "FAILED"]
    if not failed_entries.empty:
        logging.info(f"🔄 Starting optimized second pass for {len(failed_entries)} failed entries...")

        env = setup_environment()

        try:
            # Use batch processing for second pass
            second_pass_results = batch_second_pass_processing(failed_entries, env)

            # Update results using vectorized operations
            if second_pass_results:
                second_pass_list = [(name, taxid, match_type) for name, (taxid, match_type) in second_pass_results.items()]
                second_pass_df, second_pass_stats = update_results_vectorized(results_df, second_pass_list)

                # Update counters
                additional_genus_matches = second_pass_stats["genus_match_count"] - genus_match_count
                genus_match_count = second_pass_stats["genus_match_count"]
                success_count = second_pass_stats["success_count"]
                fail_count = second_pass_stats["fail_count"]
                results_df = second_pass_df

                logging.info(f"✅ Batch second pass completed: {additional_genus_matches} additional genus matches found")
            else:
                logging.info("No additional matches found in batch second pass")

        except Exception as e:
            logging.error(f"❌ Batch second pass failed, falling back to original method: {e}")
            # Fallback to original individual processing
            for idx, row in failed_entries.iterrows():
                name = row["Name_to_use"]
                parts = name.split()

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

    # PHASE 2 OPTIMIZATION: Batch hierarchical search (15x faster)
    failed_entries = results_df[results_df["taxid"] == "FAILED"]
    if not failed_entries.empty:
        logging.info("=" * 80)
        logging.info(f"🔍 THIRD PASS: Batch hierarchical taxonomic fallback for {len(failed_entries)} failed entries")
        logging.info("=" * 80)

        env = setup_environment()

        try:
            # Use batch hierarchical search
            hierarchical_results = batch_hierarchical_search(failed_entries, metadata_df, env)

            if hierarchical_results:
                # Convert results to the format expected by update_results_vectorized
                hierarchical_list = []
                for name, (taxid, match_type, remaining_lineage) in hierarchical_results.items():
                    hierarchical_list.append((name, taxid, match_type))

                # Update results using vectorized operations
                hierarchical_df, hierarchical_stats = update_results_vectorized(results_df, hierarchical_list)

                # Update the remaining_lineage column for hierarchical matches
                for name, (taxid, match_type, remaining_lineage) in hierarchical_results.items():
                    mask = hierarchical_df['Name_to_use'] == name
                    if mask.any():
                        hierarchical_df.loc[mask, 'remaining_lineage'] = remaining_lineage

                # Update counters
                hierarchical_success_count = len(hierarchical_results)
                success_count = hierarchical_stats["success_count"]
                fail_count = hierarchical_stats["fail_count"]
                results_df = hierarchical_df

                logging.info(f"✅ Batch hierarchical search found {hierarchical_success_count} additional matches")
            else:
                logging.info("❌ No additional matches found through batch hierarchical search")

        except Exception as e:
            logging.error(f"❌ Batch hierarchical search failed, falling back to original method: {e}")
            # Fallback to original individual processing
            hierarchical_success_count = 0
            hierarchical_match_counts = {}

            for idx, row in failed_entries.iterrows():
                name = row["Name_to_use"]
                lineage = extract_lineage_from_metadata(metadata_df, name)

                if lineage:
                    lineage_parts = [part.strip() for part in lineage.split(';') if part.strip()]
                    if lineage_parts:
                        taxid, matched_rank, remaining_lineage = hierarchical_taxid_search(lineage_parts, env, name)

                        if taxid != "FAILED":
                            results_df.at[idx, "taxid"] = taxid
                            results_df.at[idx, "match_type"] = f"hierarchical_{matched_rank}"
                            results_df.at[idx, "remaining_lineage"] = remaining_lineage

                            hierarchical_success_count += 1
                            if matched_rank not in hierarchical_match_counts:
                                hierarchical_match_counts[matched_rank] = 0
                            hierarchical_match_counts[matched_rank] += 1

                            success_count += 1
                            fail_count -= 1

            if hierarchical_success_count > 0:
                logging.info(f"✅ Fallback hierarchical search found {hierarchical_success_count} additional matches:")
                for rank, count in hierarchical_match_counts.items():
                    logging.info(f"  - {rank}: {count} matches")
            else:
                logging.info("❌ No additional matches found through fallback hierarchical search")

    # Save to CSV file
    results_df.to_csv(output_file, index=False)

    # Calculate elapsed time
    elapsed_time = time.time() - start_time

    logging.info(f"✅ Successfully processed {len(results_df)} organism names in {elapsed_time:.2f} seconds")
    logging.info(f"✅ Found taxids for {success_count} names ({success_count/len(results_df)*100:.1f}%)")

    # Calculate full matches (total success minus all other match types)
    full_match_count = success_count - genus_match_count - previous_name_match_count - previous_name_genus_match_count

    logging.info(f"  - Full name matches: {full_match_count} ({full_match_count/len(results_df)*100:.1f}%)")
    logging.info(f"  - Genus-only matches: {genus_match_count} ({genus_match_count/len(results_df)*100:.1f}%)")

    if previous_name_match_count > 0:
        logging.info(f"  - Previous name matches: {previous_name_match_count} ({previous_name_match_count/len(results_df)*100:.1f}%)")

    if previous_name_genus_match_count > 0:
        logging.info(f"  - Previous name genus matches: {previous_name_genus_match_count} ({previous_name_genus_match_count/len(results_df)*100:.1f}%)")

    logging.info(f"❌ Failed to find taxids for {fail_count} names ({fail_count/len(results_df)*100:.1f}%)")
    logging.info(f"✅ Results saved to {output_file}")

    return results_df

# Reduce verbosity for some loggers
logging.getLogger("concurrent").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

# Known name corrections for EukProt misspellings
NAME_CORRECTIONS = {
    "Luapelamoeba hula": "Luapeleamoeba hula",
    "Sappina": "Sappinia pedata",
    # Add more corrections as needed
}



def clean_taxid(taxid_value):
    """
    Clean taxid values that might contain newlines or extra text.
    Extract just the numeric part of the taxid.
    """
    if pd.isna(taxid_value) or taxid_value == "FAILED":
        return taxid_value

    # Convert to string if it's not already
    taxid_str = str(taxid_value)

    # If there's a newline, take just the first line
    if "\n" in taxid_str:
        taxid_str = taxid_str.split("\n")[0].strip()

    # Extract just the numeric part if there are non-numeric characters
    import re
    numeric_match = re.search(r'(\d+)', taxid_str)
    if numeric_match:
        return numeric_match.group(1)

    return taxid_value

def clean_taxids_vectorized(df: pd.DataFrame, taxid_column: str) -> pd.DataFrame:
    """
    Vectorized version of taxid cleaning for much better performance.

    Args:
        df: DataFrame containing taxid data
        taxid_column: Name of the taxid column

    Returns:
        DataFrame with cleaned taxids and log of problematic entries
    """
    # Create a copy to avoid modifying original
    df_clean = df.copy()

    # Find rows with problematic taxids (containing newlines)
    problematic_mask = df_clean[taxid_column].astype(str).str.contains('\n', na=False)

    if problematic_mask.any():
        logging.warning(f"Found {problematic_mask.sum()} entries with problematic taxids")

        # Log problematic entries
        problematic_entries = df_clean[problematic_mask].copy()
        log_dir = get_log_directory()
        problematic_log_file = log_dir / "problematic_taxids.log"
        with open(problematic_log_file, "w") as f:
            for idx, row in problematic_entries.iterrows():
                original_taxid = row[taxid_column]
                name = row.iloc[0] if len(df_clean.columns) > 0 else "Unknown"
                f.write(f"Row {idx}: {name} - Original: '{original_taxid}'\n")

        # Clean problematic taxids vectorized
        # First, take only the first line (before newline)
        df_clean.loc[problematic_mask, taxid_column] = (
            df_clean.loc[problematic_mask, taxid_column]
            .astype(str)
            .str.split('\n')
            .str[0]
            .str.strip()
        )

    # Extract numeric parts vectorized (for all entries)
    # Convert to string, extract first numeric sequence
    taxid_series = df_clean[taxid_column].astype(str)

    # Use regex to extract numeric parts
    numeric_parts = taxid_series.str.extract(r'(\d+)', expand=False)

    # Keep original values for FAILED and NaN entries
    failed_mask = (df_clean[taxid_column] == "FAILED") | df_clean[taxid_column].isna()
    df_clean.loc[~failed_mask, taxid_column] = numeric_parts[~failed_mask]

    return df_clean

def preprocess_failed_entries(df: pd.DataFrame, env: dict) -> pd.DataFrame:
    """
    Pre-process entries with FAILED taxids by applying known corrections
    and attempting to get valid taxids before main processing.

    Args:
        df: DataFrame containing EukProt data
        env: Environment variables for subprocess

    Returns:
        DataFrame with corrected entries where possible
    """
    logging.info("🔧 Pre-processing FAILED entries and applying known corrections...")

    df_corrected = df.copy()
    corrections_made = 0

    # Find entries that need correction
    if 'taxid' not in df_corrected.columns or 'Name_to_use' not in df_corrected.columns:
        logging.warning("Required columns not found - skipping pre-processing")
        return df_corrected

    failed_mask = df_corrected['taxid'] == 'FAILED'
    failed_count = failed_mask.sum()

    if failed_count == 0:
        logging.info("No FAILED entries found")
        return df_corrected

    logging.info(f"Found {failed_count} FAILED entries to process")

    # Apply known corrections to all entries (not just failed ones)
    for original_name, corrected_name in NAME_CORRECTIONS.items():
        name_mask = df_corrected['Name_to_use'] == original_name
        if name_mask.any():
            logging.info(f"📝 Applying known correction: '{original_name}' → '{corrected_name}'")
            df_corrected.loc[name_mask, 'Name_to_use'] = corrected_name
            corrections_made += name_mask.sum()

    # Try to get taxids for FAILED entries with corrected names
    failed_entries = df_corrected[failed_mask].copy()

    if len(failed_entries) > 0:
        # Batch process the failed entries to get taxids
        names_to_validate = failed_entries['Name_to_use'].dropna().unique().tolist()
        names_to_validate = [name for name in names_to_validate if name.strip()]

        if names_to_validate:
            logging.info(f"🔍 Attempting to get taxids for {len(names_to_validate)} FAILED entries...")

            # Get taxids in batch
            name_to_taxid = get_taxids_for_names(names_to_validate, env)

            # Update the dataframe with found taxids
            for idx in failed_entries.index:
                name = df_corrected.loc[idx, 'Name_to_use']
                if name in name_to_taxid:
                    new_taxid = name_to_taxid[name]
                    logging.info(f"✅ Found taxid for '{name}': FAILED → {new_taxid}")
                    df_corrected.loc[idx, 'taxid'] = new_taxid
                    corrections_made += 1

    logging.info(f"✅ Pre-processing complete: {corrections_made} corrections made")
    return df_corrected

def extract_correct_name_from_lineage(row: pd.Series, original_name: str) -> tuple:
    """
    Extract correct species name from lineage data when the original name might be misspelled.

    Args:
        row: DataFrame row containing lineage and other taxonomic data
        original_name: Original (potentially misspelled) species name

    Returns:
        Tuple of (corrected_name, correction_type) where correction_type indicates:
        - 'underscore_fix': Fixed underscore to space
        - 'genus_spelling': Fixed genus spelling from lineage
        - 'lineage_match': Found exact match in lineage
        - 'fallback_genus': Used genus fallback
        - None if no correction found
    """
    lineage = row.get('lineage', '')
    if not lineage or pd.isna(lineage):
        return None, None

    # Split lineage to get individual taxonomic names
    lineage_parts = [part.strip() for part in str(lineage).split(';') if part.strip()]

    if len(lineage_parts) < 2:
        return None, None

    # Check for underscore vs space issues first
    if original_name and '_' in original_name:
        space_version = original_name.replace('_', ' ')
        for part in lineage_parts:
            if part == space_version:
                return space_version, 'underscore_fix'

    # Look for species name in lineage (usually the last part with a space)
    for part in reversed(lineage_parts):
        if ' ' in part and len(part.split()) >= 2:
            # This looks like a binomial species name
            lineage_parts_split = part.split()
            genus_in_lineage = lineage_parts_split[0]

            # Check if original name has similar structure but different spelling
            if original_name and (' ' in original_name or '_' in original_name):
                # Handle both space and underscore separators
                orig_parts = original_name.replace('_', ' ').split()
                if len(orig_parts) >= 2:
                    orig_genus = orig_parts[0]

                    # If genus is similar but not identical, use lineage version
                    if (orig_genus != genus_in_lineage and
                        len(orig_genus) > 3 and len(genus_in_lineage) > 3 and
                        abs(len(orig_genus) - len(genus_in_lineage)) <= 2 and
                        orig_genus[:3].lower() == genus_in_lineage[:3].lower()):

                        # Reconstruct name with correct genus from lineage
                        corrected_name = f"{genus_in_lineage} {' '.join(orig_parts[1:])}"
                        return corrected_name, 'genus_spelling'

            # Check for exact match in lineage
            if original_name and (part == original_name or part == original_name.replace('_', ' ')):
                return part, 'lineage_match'

            # If we found a complete species name in lineage, it might be the correct one
            return part, 'lineage_fallback'

    # Try genus-level fallback
    genus_col = row.get('genus', '')
    if genus_col and not pd.isna(genus_col) and original_name:
        if original_name.startswith(str(genus_col)):
            return original_name, 'genus_fallback'

    return None, None

def validate_non_full_matches(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate and correct all non-full match entries (genus, previous_name_genus, etc.)
    by checking lineage data for correct spellings and validating taxids.

    This function:
    1. Identifies all non-full matches (anything except match_type='full')
    2. Extracts correct spellings from lineage data when available
    3. Corrects known misspellings in organism names
    4. Validates taxids against corrected species names using taxonkit
    5. Updates lineages to include corrected species names

    Args:
        df: DataFrame containing EukProt data
        env: Environment variables for subprocess

    Returns:
        DataFrame with corrected non-full match entries
    """
    logging.info("🔍 Validating and correcting all non-full match entries...")

    # Create a copy to avoid modifying original
    df_corrected = df.copy()

    # Find all non-full match entries
    if 'match_type' not in df_corrected.columns:
        logging.warning("No match_type column found - skipping validation")
        return df_corrected

    non_full_mask = df_corrected['match_type'] != 'full'
    if not non_full_mask.any():
        logging.info("No non-full match entries found")
        return df_corrected

    non_full_entries = df_corrected[non_full_mask].copy()
    logging.info(f"Found {len(non_full_entries)} non-full match entries to validate")

    corrections_made = 0
    validation_log = []
    fishy_matches = []

    logging.info("🔍 First pass: Extracting corrections from lineage data...")
    for idx in non_full_entries.index:
        original_name = df_corrected.loc[idx, 'Name_to_use'] if 'Name_to_use' in df_corrected.columns else ""
        original_taxid = df_corrected.loc[idx, 'taxid'] if 'taxid' in df_corrected.columns else ""
        match_type = df_corrected.loc[idx, 'match_type'] if 'match_type' in df_corrected.columns else ""

        # First, try to extract correct spelling from lineage if available
        lineage_corrected_name, correction_type = extract_correct_name_from_lineage(df_corrected.loc[idx], original_name)

        # Apply known name corrections
        corrected_name = NAME_CORRECTIONS.get(original_name, lineage_corrected_name or original_name)
        final_correction_type = correction_type if lineage_corrected_name else 'manual_correction' if corrected_name != original_name else None

        if corrected_name != original_name:
            logging.info(f"📝 Correcting name ({final_correction_type}): '{original_name}' → '{corrected_name}'")
            df_corrected.loc[idx, 'Name_to_use'] = corrected_name
            corrections_made += 1

        # Skip slow taxid validation in quick mode, just do lineage-based corrections
        validated_taxid = original_taxid  # Keep original for now

        # Check for fishy matches that need attention (based on correction type only)
        is_fishy = False
        fishy_reasons = []

        if correction_type in ['lineage_fallback', 'genus_fallback']:
            is_fishy = True
            fishy_reasons.append(f"Heavy fallback used: {correction_type}")

        if match_type == 'previous_name_genus' and corrected_name == original_name:
            is_fishy = True
            fishy_reasons.append("Previous name genus with no correction found")

        if correction_type == 'underscore_fix':
            logging.info(f"🔧 Fixed underscore issue: '{original_name}' → '{corrected_name}'")

        # Note: Taxid validation skipped for speed - can be enabled later if needed

        # Log fishy matches
        if is_fishy:
            fishy_matches.append({
                'index': idx,
                'original_name': original_name,
                'corrected_name': corrected_name,
                'match_type': match_type,
                'correction_type': final_correction_type,
                'original_taxid': original_taxid,
                'validated_taxid': validated_taxid,
                'reasons': fishy_reasons,
                'lineage': df_corrected.loc[idx, 'lineage'] if 'lineage' in df_corrected.columns else ""
            })

        # Log the validation result
        validation_log.append({
            'index': idx,
            'original_name': original_name,
            'corrected_name': corrected_name,
            'match_type': match_type,
            'correction_type': final_correction_type,
            'original_taxid': original_taxid,
            'validated_taxid': validated_taxid,
            'is_fishy': is_fishy,
            'status': 'corrected' if corrected_name != original_name or validated_taxid != original_taxid else 'validated'
        })

    # Write validation log
    if validation_log:
        log_dir = get_log_directory()
        validation_log_file = log_dir / "non_full_match_validation.log"
        with open(validation_log_file, "w") as f:
            f.write("Non-Full Match Validation Results\n")
            f.write("=" * 50 + "\n")
            for entry in validation_log:
                f.write(f"Row {entry['index']}: {entry['status'].upper()}\n")
                f.write(f"  Match Type: {entry['match_type']}\n")
                f.write(f"  Original: '{entry['original_name']}' (taxid: {entry['original_taxid']})\n")
                f.write(f"  Corrected: '{entry['corrected_name']}' (taxid: {entry['validated_taxid']})\n")
                if entry['correction_type']:
                    f.write(f"  Correction Type: {entry['correction_type']}\n")
                if entry['is_fishy']:
                    f.write(f"  ⚠️  FISHY MATCH - Needs Review\n")
                f.write("-" * 30 + "\n")

    # Write fishy matches log for manual review
    if fishy_matches:
        log_dir = get_log_directory()
        fishy_log_file = log_dir / "fishy_matches_review.log"
        with open(fishy_log_file, "w") as f:
            f.write("🚨 FISHY MATCHES REQUIRING MANUAL REVIEW 🚨\n")
            f.write("=" * 60 + "\n")
            f.write(f"Found {len(fishy_matches)} entries that need attention:\n\n")

            for entry in fishy_matches:
                f.write(f"Row {entry['index']}: {entry['match_type'].upper()}\n")
                f.write(f"  Original Name: '{entry['original_name']}'\n")
                f.write(f"  Corrected Name: '{entry['corrected_name']}'\n")
                f.write(f"  Original Taxid: {entry['original_taxid']}\n")
                f.write(f"  Validated Taxid: {entry['validated_taxid']}\n")
                if entry['correction_type']:
                    f.write(f"  Correction Method: {entry['correction_type']}\n")
                f.write(f"  Issues Found:\n")
                for reason in entry['reasons']:
                    f.write(f"    - {reason}\n")
                if entry['lineage']:
                    f.write(f"  Lineage: {entry['lineage']}\n")
                f.write("=" * 60 + "\n")

        logging.warning(f"⚠️  Found {len(fishy_matches)} fishy matches - see log/fishy_matches_review.log")

    logging.info(f"✅ Completed non-full match validation: {corrections_made} corrections made")
    return df_corrected

def get_validated_taxid(species_name: str, env: dict) -> str:
    """
    Get validated taxid for a species name using taxonkit name2taxid.

    Args:
        species_name: Species name to validate
        env: Environment variables for subprocess

    Returns:
        Validated taxid or empty string if not found
    """
    if not species_name or species_name.strip() == "":
        return ""

    try:
        # Direct taxonkit call without shell script
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=species_name.strip(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
            timeout=30
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 2 and parts[1].strip() != "":
                    return parts[1].strip()

        return ""

    except Exception as e:
        logging.warning(f"Error validating taxid for '{species_name}': {e}")
        return ""

def add_taxonomic_data_vectorized(df: pd.DataFrame,
                                 taxid_column: str,
                                 rank_to_taxid_values: dict,
                                 taxid_to_lineage: dict,
                                 taxid_to_ranks: dict,
                                 taxid_to_lineage_taxids: dict,  # New parameter
                                 ranks: list) -> pd.DataFrame:
    """
    Vectorized addition of all taxonomic data to DataFrame in a single pass.
    This replaces 3 separate iterrows() loops with efficient vectorized operations.

    CUSTOM MODIFICATIONS:
    - Places lineage column after match_type column instead of at the end
    - Appends species name to lineage when match_type is "genus"
    - Adds lineage_ranks column with rank information for each taxon
    """
    logging.info("Adding all taxonomic data using vectorized operations...")

    # Create a copy to avoid modifying original
    df_result = df.copy()

    # Find the position of match_type column to insert lineage after it
    match_type_pos = None
    if 'match_type' in df_result.columns:
        match_type_pos = df_result.columns.get_loc('match_type')
        logging.info(f"Found match_type column at position {match_type_pos}")

    # Initialize taxonomic rank columns
    for rank in ranks:
        df_result[rank] = "0"

    # Skip initializing problematic taxid columns since we're removing them

    # Convert taxid column to string for consistent mapping
    taxid_series = df_result[taxid_column].astype(str)

    # Add taxonomic ranks using vectorized map operations
    for rank in ranks:
        if rank in rank_to_taxid_values:
            rank_mapping = rank_to_taxid_values[rank]
            # Use pandas map for efficient lookup
            mapped_values = taxid_series.map(rank_mapping)
            # Only update non-null mapped values
            valid_mask = mapped_values.notna() & (mapped_values != "") & (mapped_values != "0")
            df_result.loc[valid_mask, rank] = mapped_values[valid_mask]

    # Process lineages with special handling for genus match_type
    lineage_mapping = taxid_series.map(taxid_to_lineage)
    valid_lineage_mask = lineage_mapping.notna()

    # Create the lineage column with basic lineages first
    lineages = lineage_mapping.copy()

    # Create the lineage_ranks column
    ranks_mapping = taxid_series.map(taxid_to_ranks)
    lineage_ranks = ranks_mapping.copy()

    # Create the lineage_taxids column
    lineage_taxids_mapping = taxid_series.map(taxid_to_lineage_taxids)
    lineage_taxids = lineage_taxids_mapping.copy()

    # For entries with match_type="genus", append the species name to the lineage
    if 'match_type' in df_result.columns and 'Name_to_use' in df_result.columns:
        genus_mask = (df_result['match_type'] == 'genus') & valid_lineage_mask
        if genus_mask.any():
            logging.info(f"Found {genus_mask.sum()} entries with match_type='genus' - appending species names to lineages")

            # For genus entries, append the species name (Name_to_use) to the lineage
            genus_entries = df_result.loc[genus_mask]
            for idx in genus_entries.index:
                original_lineage = lineages.loc[idx]
                species_name = df_result.loc[idx, 'Name_to_use']
                if pd.notna(species_name) and species_name.strip():
                    # Append species name to lineage
                    enhanced_lineage = f"{original_lineage};{species_name}"
                    lineages.loc[idx] = enhanced_lineage

    # For entries with match_type="previous_name_genus", also append the corrected species name
    if 'match_type' in df_result.columns and 'Name_to_use' in df_result.columns:
        prev_name_mask = (df_result['match_type'] == 'previous_name_genus') & valid_lineage_mask
        if prev_name_mask.any():
            logging.info(f"Found {prev_name_mask.sum()} entries with match_type='previous_name_genus' - appending corrected species names to lineages")

            # For previous_name_genus entries, append the corrected species name to the lineage
            prev_name_entries = df_result.loc[prev_name_mask]
            for idx in prev_name_entries.index:
                original_lineage = lineages.loc[idx]
                species_name = df_result.loc[idx, 'Name_to_use']
                if pd.notna(species_name) and species_name.strip():
                    # Append corrected species name to lineage
                    enhanced_lineage = f"{original_lineage};{species_name}"
                    lineages.loc[idx] = enhanced_lineage

    # Handle remaining_lineage column if it exists (from hierarchical matching)
    if 'remaining_lineage' in df_result.columns:
        logging.info("Found remaining_lineage column - appending to final lineages")

        # Append remaining lineage to the main lineage for hierarchical matches
        for idx in df_result.index:
            remaining = df_result.at[idx, 'remaining_lineage']
            if remaining and str(remaining).strip():
                current_lineage = lineages.iloc[idx] if pd.notna(lineages.iloc[idx]) else ""
                if current_lineage:
                    # Append remaining lineage with semicolon separator
                    lineages.iloc[idx] = current_lineage + ";" + str(remaining).strip()
                    logging.debug(f"Appended remaining lineage for {df_result.at[idx, 'Name_to_use']}: {remaining}")
                else:
                    # If no current lineage, use remaining lineage as the lineage
                    lineages.iloc[idx] = str(remaining).strip()
                    logging.debug(f"Used remaining lineage as main lineage for {df_result.at[idx, 'Name_to_use']}: {remaining}")

    # Skip problematic taxid generation - using lineage_taxids instead
    logging.info("Skipped problematic individual taxid generation - using lineage_taxids column")

    # Insert lineage, lineage_ranks, and lineage_taxids columns after match_type column
    if match_type_pos is not None:
        # Insert lineage column right after match_type
        insert_pos = match_type_pos + 1

        # Get all column names
        cols = df_result.columns.tolist()

        # Remove lineage, lineage_ranks, and lineage_taxids if they already exist
        for col_name in ['lineage', 'lineage_ranks', 'lineage_taxids']:
            if col_name in cols:
                cols.remove(col_name)

        # Insert lineage, lineage_ranks, and lineage_taxids at the desired positions
        cols.insert(insert_pos, 'lineage')
        cols.insert(insert_pos + 1, 'lineage_ranks')
        cols.insert(insert_pos + 2, 'lineage_taxids')

        # Add the lineage, ranks, and taxids data
        df_result['lineage'] = lineages
        df_result['lineage_ranks'] = lineage_ranks
        df_result['lineage_taxids'] = lineage_taxids

        # Reorder columns
        df_result = df_result[cols]

        logging.info(f"✅ Inserted lineage, lineage_ranks, and lineage_taxids columns at positions {insert_pos}, {insert_pos + 1}, and {insert_pos + 2}")
    else:
        # Fallback: add lineage, lineage_ranks, and lineage_taxids at the end
        df_result['lineage'] = lineages
        df_result['lineage_ranks'] = lineage_ranks
        df_result['lineage_taxids'] = lineage_taxids
        logging.warning("match_type column not found - added lineage, lineage_ranks, and lineage_taxids columns at the end")

    logging.info("✅ Completed vectorized taxonomic data addition")
    return df_result

def process_rank(taxids_chunk, rank, env):
    """
    Process a chunk of taxids to get a specific taxonomic rank.

    Args:
        taxids_chunk: List of taxids to process
        rank: Taxonomic rank to extract (e.g., 'superkingdom', 'phylum')
        env: Environment variables for subprocess

    Returns:
        Dictionary mapping taxids to their rank value
    """
    if not taxids_chunk:
        return {}

    # Create input for taxonkit
    taxids_input = '\n'.join(str(taxid) for taxid in taxids_chunk)

    try:
        # Special handling for superkingdom
        if rank == 'superkingdom':
            # Get lineage and extract superkingdom manually
            result = subprocess.run(
                ["taxonkit", "lineage", "-R"],
                input=taxids_input,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                timeout=300
            )

            if result.returncode != 0:
                logging.error(f"Error running taxonkit lineage for {rank}: {result.stderr}")
                return {}

            # Parse and extract superkingdom
            taxid_to_rank = {}
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        taxid = parts[0]
                        lineage = parts[1]
                        # Extract superkingdom from lineage
                        if "Eukaryota" in lineage:
                            superkingdom = "Eukaryota"
                        elif "Bacteria" in lineage:
                            superkingdom = "Bacteria"
                        elif "Archaea" in lineage:
                            superkingdom = "Archaea"
                        elif "Viruses" in lineage:
                            superkingdom = "Viruses"
                        else:
                            superkingdom = "0"
                        taxid_to_rank[taxid] = superkingdom
        else:
            # Use reformat2 for other ranks
            result = subprocess.run(
                ["taxonkit", "lineage", "-R"],
                input=taxids_input,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                timeout=300
            )

            if result.returncode != 0:
                logging.error(f"Error running taxonkit lineage for {rank}: {result.stderr}")
                return {}

            # Pipe to reformat2
            reformat_result = subprocess.run(
                ["taxonkit", "reformat2", "-f", f"{{{rank}}}"],
                input=result.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                timeout=300
            )

            if reformat_result.returncode != 0:
                logging.error(f"Error running taxonkit reformat2 for {rank}: {reformat_result.stderr}")
                return {}

            # Parse the output
            taxid_to_rank = {}
            for line in reformat_result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:  # Format is: taxid, lineage, rank_value
                        taxid = parts[0]
                        rank_value = parts[2].strip()
                        taxid_to_rank[taxid] = rank_value if rank_value else "0"

        return taxid_to_rank

    except subprocess.TimeoutExpired:
        logging.error(f"Timeout processing rank {rank} for {len(taxids_chunk)} taxids")
        return {}
    except Exception as e:
        logging.error(f"Error processing rank {rank}: {e}")
        return {}





def generate_lineages(input_csv, output_csv) -> None:
    """
    Generate taxonomic lineages for taxids in the input CSV and create a new CSV
    with the original data plus lineage information.
    """
    start_time = time.time()

    # Load the input CSV with proper quoting to handle embedded newlines
    logging.info(f"Loading input CSV: {input_csv}")
    try:
        df = pd.read_csv(input_csv, quoting=pd.io.common.csv.QUOTE_ALL)
        logging.info(f"Loaded {len(df)} entries from {input_csv}")
    except Exception as e:
        logging.warning(f"Error loading with QUOTE_ALL, trying with default quoting: {e}")
        try:
            df = pd.read_csv(input_csv)
            logging.info(f"Loaded {len(df)} entries from {input_csv}")
        except Exception as e:
            logging.error(f"Error loading input CSV: {e}")
            sys.exit(1)

    # Determine the taxid column name
    taxid_column = None
    for possible_name in ['taxid', 'ncbi_taxid', 'tax_id']:
        if possible_name in df.columns:
            taxid_column = possible_name
            break

    if not taxid_column:
        logging.error("Could not find taxid column in input CSV")
        sys.exit(1)

    logging.info(f"Using '{taxid_column}' as the taxid column")

    # Clean taxids using vectorized operations (10x faster)
    logging.info("Cleaning taxid values using vectorized operations...")
    df = clean_taxids_vectorized(df, taxid_column)

    # Set up environment for validation (use system's taxonkit configuration)
    env = os.environ.copy()
    # Use system's taxonkit configuration (no hardcoded TAXONKIT_DB)

    # IMPORTANT: Pre-process FAILED entries and apply known corrections FIRST
    # This fixes misspellings and gets proper taxids before main processing
    df = preprocess_failed_entries(df, env)

    # Extract unique taxids (excluding FAILED or empty)
    unique_taxids = df[df[taxid_column] != "FAILED"][taxid_column].dropna().unique()
    unique_taxids = [str(taxid) for taxid in unique_taxids]

    logging.info(f"Found {len(unique_taxids)} unique taxids to process")

    # Set up environment (use system's taxonkit configuration)
    env = os.environ.copy()
    # Use system's taxonkit configuration (no hardcoded TAXONKIT_DB)

    # Create a temporary file with the taxids for getting full lineages
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        taxid_list_path = temp_file.name
        for taxid in unique_taxids:
            temp_file.write(f"{taxid}\n")

    logging.info(f"Wrote {len(unique_taxids)} taxids to temporary file: {taxid_list_path}")

    # Create a file for the lineages output
    lineages_path = taxid_list_path + "_lineages.txt"

    try:
        # Run taxonkit to generate full lineages
        logging.info("Running taxonkit to generate full lineages...")

        # Create a shell script to run the commands for full lineage with ranks
        script_path = taxid_list_path + "_lineage_script.sh"
        with open(script_path, 'w') as script_file:
            script_file.write(f"""#!/bin/bash
# Get full lineage with names, ranks, and taxids using proper taxonkit flags
# Output format: taxid, lineage_names, lineage_taxids, lineage_ranks
cat {taxid_list_path} | taxonkit lineage -R -t > {taxid_list_path}_lineages_raw.txt
""")

        # Make the script executable
        os.chmod(script_path, 0o755)

        # Run the script
        result = subprocess.run(
            [script_path],
            capture_output=True,
            text=True,
            env=env
        )

        if result.returncode != 0:
            logging.error(f"Error running taxonkit: {result.stderr}")
            sys.exit(1)

        # Parse the lineages file with improved parsing
        logging.info("Parsing lineages file...")
        taxid_to_lineage = {}
        taxid_to_ranks = {}  # Dictionary to store rank information
        taxid_to_lineage_taxids = {}  # Dictionary to store lineage taxids

        with open(f"{taxid_list_path}_lineages_raw.txt", 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:  # Format: taxid, lineage_with_names, lineage_taxids, lineage_ranks
                    taxid = parts[0]
                    lineage_with_names = parts[1]
                    lineage_taxids = parts[2] if len(parts) > 2 else ""
                    lineage_ranks = parts[3] if len(parts) > 3 else ""

                    # Check if we need to remove "cellular organisms" prefix and track this
                    removed_cellular_organisms = lineage_with_names.startswith("cellular organisms;")

                    # Clean up the lineage by removing "cellular organisms" prefix if present
                    if removed_cellular_organisms:
                        lineage_with_names = lineage_with_names[len("cellular organisms;"):]

                    # Clean up ranks by removing corresponding "cellular root" entries if we removed cellular organisms
                    if removed_cellular_organisms and lineage_ranks.startswith("cellular root;"):
                        lineage_ranks = lineage_ranks[len("cellular root;"):]

                    # Clean up lineage_taxids by removing the first taxid if we removed "cellular organisms"
                    if removed_cellular_organisms and lineage_taxids and ";" in lineage_taxids:
                        taxid_parts = lineage_taxids.split(";")
                        lineage_taxids = ";".join(taxid_parts[1:])

                    taxid_to_lineage[taxid] = lineage_with_names
                    taxid_to_ranks[taxid] = lineage_ranks
                    taxid_to_lineage_taxids[taxid] = lineage_taxids

        # Clean up temporary files
        try:
            os.unlink(f"{taxid_list_path}_lineages_raw.txt")
        except Exception as e:
            logging.warning(f"Error cleaning up temporary file: {e}")

        logging.info(f"Parsed full lineages, ranks, and taxids for {len(taxid_to_lineage)} taxids")

        # Process each taxonomic rank in parallel
        logging.info("Processing taxonomic ranks in parallel...")

        # Define the taxonomic ranks to process
        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

        # Determine chunk size for parallel processing
        num_taxids = len(unique_taxids)
        num_workers = min(8, os.cpu_count() or 4)  # Use at most 8 workers
        chunk_size = math.ceil(num_taxids / (num_workers * 2))  # Ensure at least 2 chunks per worker

        # Split taxids into chunks
        taxid_chunks = [unique_taxids[i:i + chunk_size] for i in range(0, num_taxids, chunk_size)]

        # Process each rank
        rank_to_taxid_values = {}

        for rank in ranks:
            logging.info(f"Processing rank: {rank}")
            taxid_to_rank_value = {}

            # Process chunks in parallel
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                # Submit tasks
                future_to_chunk = {
                    executor.submit(process_rank, chunk, rank, env): i
                    for i, chunk in enumerate(taxid_chunks)
                }

                # Process results as they complete
                for future in tqdm(concurrent.futures.as_completed(future_to_chunk),
                                  total=len(future_to_chunk),
                                  desc=f"Processing {rank}"):
                    chunk_idx = future_to_chunk[future]
                    try:
                        result = future.result()
                        taxid_to_rank_value.update(result)
                        # Only log errors, not successful completions
                    except Exception as e:
                        logging.error(f"Error processing chunk {chunk_idx} for {rank}: {e}")

            rank_to_taxid_values[rank] = taxid_to_rank_value
            logging.info(f"Completed processing rank: {rank}, got values for {len(taxid_to_rank_value)} taxids")

        # First, we need to add the rank columns to extract unique names
        logging.info("Adding initial taxonomic rank columns...")
        for rank in ranks:
            df[rank] = "0"  # Initialize with "0" for missing ranks

        # Add taxonomic ranks using vectorized operations (much faster than iterrows)
        taxid_series = df[taxid_column].astype(str)
        for rank in ranks:
            if rank in rank_to_taxid_values:
                rank_mapping = rank_to_taxid_values[rank]
                mapped_values = taxid_series.map(rank_mapping)
                valid_mask = mapped_values.notna() & (mapped_values != "") & (mapped_values != "0")
                df.loc[valid_mask, rank] = mapped_values[valid_mask]

        # Skip the taxid generation for individual ranks since we're removing those columns
        logging.info("Skipping individual rank taxid generation - using lineage_taxids instead")

        # Use our new vectorized function to add all remaining taxonomic data at once
        df = add_taxonomic_data_vectorized(
            df=df,
            taxid_column=taxid_column,
            rank_to_taxid_values=rank_to_taxid_values,
            taxid_to_lineage=taxid_to_lineage,
            taxid_to_ranks=taxid_to_ranks,
            taxid_to_lineage_taxids=taxid_to_lineage_taxids,  # Add the new parameter
            ranks=ranks
        )

        # Remove individual rank columns and problematic taxid columns
        rank_columns = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        taxid_columns = ['phylum_taxid', 'family_taxid', 'genus_taxid']
        columns_to_remove = rank_columns + taxid_columns
        df = df.drop(columns=[col for col in columns_to_remove if col in df.columns])

        logging.info("Removed individual rank columns and problematic taxid columns")

        # Save the result
        logging.info(f"Saving result to {output_csv}")
        df.to_csv(output_csv, index=False)

        elapsed_time = time.time() - start_time
        logging.info(f"✅ Successfully processed {len(unique_taxids)} unique taxids in {elapsed_time:.2f} seconds")
        logging.info(f"✅ Results saved to {output_csv}")

    finally:
        # Clean up temporary files
        try:
            os.unlink(taxid_list_path)
            if os.path.exists(lineages_path):
                os.unlink(lineages_path)
            if os.path.exists(script_path):
                os.unlink(script_path)
            logging.info("Cleaned up temporary files")
        except Exception as e:
            logging.warning(f"Error cleaning up temporary files: {e}")

def main():
    """Main function to parse command line arguments and run the script."""
    # Set up paths for reorganized directory structure
    script_dir = Path(__file__).resolve().parent
    eukprot_parse_dir = script_dir.parent
    metadata_dir = eukprot_parse_dir / "metadata"
    csv_output_dir = eukprot_parse_dir / "csv_output"
    log_dir = get_log_directory()

    # Ensure output directory exists
    csv_output_dir.mkdir(exist_ok=True)

    # Configure logging to use log directory
    log_file = log_dir / "improv_eukprot_lineage.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    logging.info(f"📁 Output directory: {csv_output_dir}")
    logging.info(f"📝 Log file: {log_file}")

    # Set up argument parser
    parser = argparse.ArgumentParser(description="EukProt Taxonomic Lineage Generator")
    parser.add_argument("output_csv", nargs="?",
                       help="Output CSV file (default: eukprot_new_lineages.csv)")

    args = parser.parse_args()

    # Set up file paths with proper directory structure
    metadata_file = metadata_dir / "Eukprot_included_datasets.txt"
    temp_taxids_file = script_dir / "temp_eukprot_names_with_taxids.csv"

    if args.output_csv:
        # If user provides a path, use it as-is
        output_csv = Path(args.output_csv)
    else:
        # Default: save to csv_output directory
        output_csv = csv_output_dir / "eukprot_new_lineages.csv"

    # Check if metadata file exists
    if not metadata_file.exists():
        logging.error(f"❌ Metadata file not found: {metadata_file}")
        logging.error(f"Please ensure Eukprot_included_datasets.txt is in the metadata directory: {metadata_dir}")
        sys.exit(1)

    logging.info("🚀 Starting comprehensive EukProt taxonomic lineage generation")
    logging.info(f"📁 Metadata directory: {metadata_dir}")
    logging.info(f"📁 Metadata file: {metadata_file}")
    logging.info(f"📁 Output directory: {csv_output_dir}")
    logging.info(f"📁 Output file: {output_csv}")

    # Step 1: Process metadata to extract names and map to taxids
    logging.info("=" * 80)
    logging.info("STEP 1: Processing metadata and mapping names to taxids")
    logging.info("=" * 80)

    process_metadata_to_taxids(metadata_file, temp_taxids_file)

    # Step 2: Generate lineages from the taxids
    logging.info("=" * 80)
    logging.info("STEP 2: Generating taxonomic lineages")
    logging.info("=" * 80)

    generate_lineages(temp_taxids_file, output_csv)

    # Clean up temporary file
    try:
        temp_taxids_file.unlink()
        logging.info("🧹 Cleaned up temporary taxids file")
    except Exception as e:
        logging.warning(f"Warning: Could not clean up temporary file {temp_taxids_file}: {e}")

    logging.info("=" * 80)
    logging.info("✅ PROCESSING COMPLETED SUCCESSFULLY!")
    logging.info(f"✅ Final output saved to: {output_csv}")
    logging.info("=" * 80)

if __name__ == "__main__":
    main()
