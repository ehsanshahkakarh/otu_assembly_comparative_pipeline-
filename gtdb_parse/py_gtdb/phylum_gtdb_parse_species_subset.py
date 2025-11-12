#!/usr/bin/env python3
"""
GTDB Phylum Parser (Species Subset Version)

This script processes GTDB taxonomy files to extract phylum-level information.
It subsets data to one genome per species, then counts the number of species per phylum
and creates output files with both summary statistics and detailed accession information.

Input:
- GTDB bacterial taxonomy file (00bac120_taxonomy.tsv)
- GTDB archaeal taxonomy file (00ar53_taxonomy.tsv)

Output:
- Summary CSV with phylum species counts (gtdb_phylum_species_counts.csv)

Usage:
    python phylum_gtdb_parse.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
"""

import pandas as pd
import argparse
import sys
import os
import logging
import time
import re
import subprocess
import tempfile
from typing import Dict, List, Tuple
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from contextlib import contextmanager
from tqdm import tqdm

# Constants
BAC_TAXONOMY_FILENAME = "00bac120_taxonomy.tsv"
ARC_TAXONOMY_FILENAME = "00ar53_taxonomy.tsv"
SUMMARY_OUTPUT_FILENAME = "gtdb_phylum_species_counts.csv"
# Removed ALL_RANKS_OUTPUT_FILENAME - not needed
TAXID_MAP_FILENAME = "taxid.map"

# Regex patterns for extracting taxonomic ranks
DOMAIN_REGEX = r"d__([A-Za-z0-9_\-]+)"
PHYLUM_REGEX = r"p__([A-Za-z0-9_\-]+)"
CLASS_REGEX = r"c__([A-Za-z0-9_\-]+)"
ORDER_REGEX = r"o__([A-Za-z0-9_\-]+)"
FAMILY_REGEX = r"f__([A-Za-z0-9_\-]+)"
GENUS_REGEX = r"g__([A-Za-z0-9_\-]+)"
SPECIES_REGEX = r"s__([^;]+)"

# Removed accession cleaning - not needed for species subset analysis

# Processing parameters
CHUNK_SIZE = 100000  # Number of rows to process at once

@dataclass
class PerformanceMetrics:
    """Class to store performance metrics for different operations."""
    operation: str
    start_time: float
    end_time: float
    rows_processed: int = 0

    @property
    def duration(self) -> float:
        """Calculate duration in seconds."""
        return self.end_time - self.start_time

    def __str__(self) -> str:
        """String representation of the metrics."""
        base_str = f"{self.operation}: {self.duration:.2f}s"
        if self.rows_processed > 0:
            base_str += f" ({self.rows_processed:,} rows, {self.rows_processed/self.duration:.0f} rows/s)"
        return base_str

@contextmanager
def timer(operation: str, metrics_list: List[PerformanceMetrics]):
    """Context manager for timing operations."""
    start_time = time.time()
    metrics = PerformanceMetrics(operation=operation, start_time=start_time, end_time=start_time)
    try:
        yield metrics
    finally:
        metrics.end_time = time.time()
        metrics_list.append(metrics)

def setup_logging(log_level: str, output_dir: Path) -> None:
    """
    Set up logging configuration with clean terminal output and detailed file logging.

    Args:
        log_level (str): Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        output_dir (Path): Directory to store log files
    """
    # Create logs directory if it doesn't exist
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    # Create log filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"gtdb_phylum_parse_{timestamp}.log"

    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(message)s')  # Clean terminal output

    # Create handlers
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(file_formatter)
    file_handler.setLevel(getattr(logging, log_level))

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(console_formatter)
    console_handler.setLevel(getattr(logging, log_level))

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)  # Capture all levels
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)

    # Log initialization with file path only to console
    print(f"📝 Log file: {log_file}")

def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Parse GTDB taxonomy files for phylum-level information")
    parser.add_argument("--input-dir", help="Directory containing input files (default: script directory)")
    parser.add_argument("--output-dir", help="Directory for output files (default: ../csv_gtdb)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--log-level",
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       default='INFO',
                       help="Set the logging level")
    return parser.parse_args()

def setup_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path]:
    """
    Set up input and output file paths.

    Args:
        args (argparse.Namespace): Command line arguments

    Returns:
        Tuple[Path, Path, Path, Path]: Paths for bacterial taxonomy file, archaeal taxonomy file,
                                      summary output file, and taxid map file
    """
    # Determine current directory and set up paths
    script_dir = Path(__file__).resolve().parent

    # Set input directory (default to ../metadata if not specified)
    input_dir = Path(args.input_dir) if args.input_dir else script_dir.parent / "metadata"

    # Set output directory (default to ../csv_gtdb if not specified)
    output_dir = Path(args.output_dir) if args.output_dir else script_dir.parent / "csv_gtdb"

    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)

    # Define file paths
    bac_file = input_dir / BAC_TAXONOMY_FILENAME
    arc_file = input_dir / ARC_TAXONOMY_FILENAME
    output_file = output_dir / SUMMARY_OUTPUT_FILENAME

    # Define taxid map file path (in the script directory)
    taxid_map_file = script_dir / TAXID_MAP_FILENAME

    # Check if taxid map file exists
    if not taxid_map_file.exists():
        print(f"⚠️ Warning: Taxid map file not found at {taxid_map_file}")
        print("  Taxid information will not be included in the output.")
        taxid_map_file = None

    return bac_file, arc_file, output_file, taxid_map_file

def load_taxid_map(file_path: Path) -> Dict[str, str]:
    """
    Load taxid.map file into a dictionary mapping accession to taxid.

    Args:
        file_path (Path): Path to the taxid.map file

    Returns:
        Dict[str, str]: Dictionary mapping accession to taxid
    """
    if file_path is None or not file_path.exists():
        logging.warning("Taxid map file not provided or does not exist")
        return {}

    try:
        taxid_map = {}

        # Process the file
        with tqdm(desc=f"Loading taxid map", unit="entries") as pbar:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split('\t')
                    if len(parts) != 2:
                        continue

                    accession, taxid = parts
                    # Use accession as-is for species subset analysis
                    taxid_map[accession.strip()] = taxid
                    pbar.update(1)

        return taxid_map

    except Exception as e:
        logging.error(f"Error loading taxid map: {str(e)}")
        return {}

def load_taxonomy_file(file_path: Path, use_chunks: bool = True) -> pd.DataFrame:
    """
    Load a GTDB taxonomy file, optionally in chunks for large files.

    Args:
        file_path (Path): Path to the taxonomy file
        use_chunks (bool): Whether to use chunked processing for large files

    Returns:
        pd.DataFrame: DataFrame containing the taxonomy data

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    try:
        logging.debug(f"Loading taxonomy file: {file_path}")

        if not use_chunks:
            # Load the entire file at once
            return pd.read_csv(file_path, sep="\t", header=None, names=["accession", "taxonomy"])

        # Process in chunks for large files
        chunks = []
        total_rows = 0

        # First count the lines to set up the progress bar
        with open(file_path, 'r') as f:
            for count, _ in enumerate(f):
                pass
        total_lines = count + 1

        # Now process in chunks with a progress bar
        with tqdm(total=total_lines, desc=f"Loading {file_path.name}", unit="rows") as pbar:
            for chunk in pd.read_csv(file_path, sep="\t", header=None,
                                    names=["accession", "taxonomy"],
                                    chunksize=CHUNK_SIZE):
                chunks.append(chunk)
                total_rows += len(chunk)
                pbar.update(len(chunk))

        logging.debug(f"Loaded {total_rows} rows from {file_path}")
        return pd.concat(chunks, ignore_index=True)

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        raise
    except pd.errors.EmptyDataError:
        logging.error(f"File is empty: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading {file_path}: {str(e)}")
        raise

def process_taxonomy_data(df: pd.DataFrame, taxid_map: Dict[str, str] = None) -> pd.DataFrame:
    """
    Process taxonomy data to extract all taxonomic ranks and add taxids if available.

    Args:
        df (pd.DataFrame): DataFrame with taxonomy strings
        taxid_map (Dict[str, str], optional): Dictionary mapping accessions to taxids

    Returns:
        pd.DataFrame: DataFrame with extracted taxonomic information for all ranks
    """
    logging.debug("Processing taxonomy data")

    # Use original accession for species subset analysis (no cleaning needed)

    # Extract all taxonomic ranks using regex, fill missing values with "Unclassified"
    print("Extracting taxonomic ranks...")
    ranks = [
        ("domain", DOMAIN_REGEX),
        ("phylum", PHYLUM_REGEX),
        ("class", CLASS_REGEX),
        ("order", ORDER_REGEX),
        ("family", FAMILY_REGEX),
        ("genus", GENUS_REGEX),
        ("species", SPECIES_REGEX)
    ]

    for rank_name, regex_pattern in tqdm(ranks, desc="Extracting ranks", unit="rank"):
        df[rank_name] = df["taxonomy"].str.extract(regex_pattern).fillna("Unclassified")
        # Strip whitespace from extracted taxonomic names
        df[rank_name] = df[rank_name].str.strip()

    # Remove rows with missing phylum (for backward compatibility)
    df = df.dropna(subset=["phylum"])
    logging.debug(f"Processed {len(df)} valid entries with phylum information")

    # Assert that there are no missing accessions
    assert df[["accession"]].isna().sum().sum() == 0, "Missing accessions detected!"

    # Add taxids if available
    if taxid_map and len(taxid_map) > 0:
        print("Adding taxids from taxid map...")

        # Create a new column for taxids, initialized with NaN
        df["taxid"] = pd.NA

        # Map accessions to taxids
        matched_count = 0
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Mapping taxids", unit="entry"):
            accession = row["accession"]
            if accession in taxid_map:
                df.at[idx, "taxid"] = taxid_map[accession]
                matched_count += 1

        # Report matching statistics
        match_percentage = (matched_count / len(df)) * 100
        print(f"✅ Matched {matched_count} of {len(df)} entries with taxids ({match_percentage:.2f}%)")

    return df

def standardize_phylum_name(phylum: str) -> str:
    """
    Standardize phylum names by removing suffixes like _A, _B, etc.

    Args:
        phylum (str): Original phylum name

    Returns:
        str: Standardized phylum name
    """
    # Check if the phylum name ends with _[A-Z]
    if re.search(r'_[A-Z]$', phylum):
        # Remove the suffix
        return phylum.rsplit('_', 1)[0]
    return phylum

def merge_related_phyla(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge counts for related phyla (e.g., Bacillota and Bacillota_A).

    Args:
        counts_df (pd.DataFrame): DataFrame with phylum counts

    Returns:
        pd.DataFrame: DataFrame with merged phylum counts
    """
    logging.info("Merging related phyla...")

    # Create a copy to avoid modifying the original
    df = counts_df.copy()

    # Create a new column with standardized phylum names
    df['phylum_base'] = df['phylum'].apply(standardize_phylum_name)

    # Find phyla that need to be merged (those where phylum_base != phylum)
    phyla_to_merge = df[df['phylum_base'] != df['phylum']]['phylum_base'].unique()

    # Log the phyla that will be merged (only in debug mode)
    if len(phyla_to_merge) > 0:
        logging.debug(f"Found {len(phyla_to_merge)} phyla to merge")
        for phylum in phyla_to_merge:
            variants = df[df['phylum_base'] == phylum]['phylum'].unique()
            logging.debug(f"  {phylum}: {', '.join(variants)}")

    if len(phyla_to_merge) == 0:
        logging.info("No related phyla found to merge")
        return counts_df

    # For each base phylum that has variants
    merged_rows = []
    for base_phylum in tqdm(phyla_to_merge, desc="Merging phyla", unit="phylum"):
        # Get all rows for this base phylum and its variants
        related_rows = df[df['phylum_base'] == base_phylum]

        # Group by domain and sum the genome counts
        for domain, domain_group in related_rows.groupby('domain'):
            # Sum the genome counts
            total_count = domain_group['gtdb_genome_count'].sum()

            # Create a new row with the merged data (no accession tracking for species subset)
            merged_row = {
                'phylum': base_phylum,
                'domain': domain,
                'gtdb_genome_count': total_count,
                'phylum_base': base_phylum  # Keep this for filtering later
            }

            # Combine taxids if available
            if 'taxids' in domain_group.columns:
                all_taxids = []
                for taxid_list in domain_group['taxids']:
                    if isinstance(taxid_list, list):
                        all_taxids.extend(taxid_list)
                merged_row['taxids'] = all_taxids
            merged_rows.append(merged_row)

            # Log the merge
            variants = domain_group['phylum'].tolist()
            logging.info(f"Merged {len(variants)} variants of {base_phylum} ({domain}): {', '.join(variants)}")

    # Create a DataFrame from the merged rows
    merged_df = pd.DataFrame(merged_rows)

    # Remove the original rows for the merged phyla (both variants and base phyla)
    df = df[~(df['phylum_base'].isin(phyla_to_merge))]

    # Add the merged rows
    df = pd.concat([df, merged_df], ignore_index=True)

    # Remove the temporary column and sort by genome count
    df = df.drop(columns=['phylum_base'])
    df = df.sort_values('gtdb_genome_count', ascending=False)

    logging.info(f"Merged {len(phyla_to_merge)} phylum groups")
    return df

def validate_taxonomy(df: pd.DataFrame) -> bool:
    """
    Validate taxonomic hierarchy consistency.

    Args:
        df (pd.DataFrame): DataFrame with taxonomic information

    Returns:
        bool: True if taxonomy is consistent, False otherwise
    """
    logging.info("Validating taxonomic hierarchy consistency...")
    inconsistencies = []

    # Define the taxonomic hierarchy
    hierarchy = [
        ("genus", "family"),
        ("family", "order"),
        ("order", "class"),
        ("class", "phylum"),
        ("phylum", "domain")
    ]

    # Check each level of the hierarchy
    for lower_rank, higher_rank in tqdm(hierarchy, desc="Validating taxonomy", unit="level"):
        # Skip ranks with all "Unclassified" values
        if (df[lower_rank] == "Unclassified").all() or (df[higher_rank] == "Unclassified").all():
            continue

        # Group by the lower rank and check if all entries have the same higher rank
        rank_groups = df.groupby(lower_rank)[higher_rank].unique()
        inconsistent_ranks = rank_groups[rank_groups.apply(len) > 1]

        if len(inconsistent_ranks) > 0:
            logging.warning(f"Found {len(inconsistent_ranks)} {lower_rank} entries with inconsistent {higher_rank} assignments")
            for rank, higher_ranks in inconsistent_ranks.items():
                inconsistency = {
                    "lower_rank": lower_rank,
                    "lower_value": rank,
                    "higher_rank": higher_rank,
                    "higher_values": list(higher_ranks)
                }
                inconsistencies.append(inconsistency)
                logging.warning(f"  {rank} ({lower_rank}) has multiple {higher_rank} assignments: {', '.join(higher_ranks)}")

    # Return True if no inconsistencies were found
    is_consistent = len(inconsistencies) == 0
    if is_consistent:
        logging.info("✅ Taxonomic hierarchy is consistent")
    else:
        logging.warning(f"⚠️ Found {len(inconsistencies)} taxonomic inconsistencies")

    return is_consistent

def get_taxids_from_names_gtdb(phylum_names: List[str]) -> Dict[str, str]:
    """
    Get taxids for phylum names using GTDB taxdump via taxonkit name2taxid.

    Args:
        phylum_names: List of phylum names

    Returns:
        Dictionary mapping phylum names to their GTDB taxids
    """
    if not phylum_names:
        return {}

    print(f"Getting GTDB taxids for {len(phylum_names)} phylum names...")

    # Set up environment to use GTDB taxdump
    script_dir = Path(__file__).resolve().parent
    gtdb_taxdump_dir = script_dir.parent / "taxdump_gtdp" / "gtdb-taxdump-R226"
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(gtdb_taxdump_dir)

    name_to_taxid = {}

    # Create temporary file for phylum names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in phylum_names:
            temp_file.write(f"{name}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit name2taxid {temp_filename} (using GTDB taxdump)")

        # Run taxonkit name2taxid command with GTDB database
        result = subprocess.run(
            ["taxonkit", "name2taxid", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} GTDB name2taxid results")

            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()

                    if taxid and taxid != "":
                        name_to_taxid[name] = taxid

                        # Show first few examples
                        if i < 3:
                            print(f"  Example: {name} -> {taxid} (GTDB)")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No results from GTDB taxonkit name2taxid")

    except Exception as e:
        print(f"Error running GTDB taxonkit name2taxid: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(name_to_taxid)} phylum name->GTDB taxid mappings")
    return name_to_taxid

def get_lineages_from_gtdb_taxids(taxids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Get lineages for GTDB taxids using GTDB taxdump via taxonkit lineage -R -t.

    Args:
        taxids: List of GTDB taxids

    Returns:
        Dictionary mapping taxids to (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}

    print(f"Getting GTDB lineages for {len(taxids)} taxids...")

    # Set up environment to use GTDB taxdump
    script_dir = Path(__file__).resolve().parent
    gtdb_taxdump_dir = script_dir.parent / "taxdump_gtdp" / "gtdb-taxdump-R226"
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(gtdb_taxdump_dir)

    lineage_data = {}

    # Create temporary file for taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit lineage -R -t {temp_filename} (using GTDB taxdump)")

        # Run taxonkit lineage command with GTDB database
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} GTDB lineage results")

            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip()
                    lineage_taxids = parts[2].strip()
                    lineage_ranks = parts[3].strip()

                    # Store the lineage data
                    if lineage and lineage != taxid:
                        lineage_data[taxid] = (
                            lineage,
                            lineage_ranks,
                            lineage_taxids
                        )

                        # Show first few examples
                        if i < 3:
                            print(f"  Example: GTDB Taxid {taxid}:")
                            print(f"    Lineage: {lineage}")
                            print(f"    Ranks: {lineage_ranks}")
                            print(f"    Taxids: {lineage_taxids}")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No lineage results from GTDB taxonkit")

    except Exception as e:
        print(f"Error running GTDB taxonkit lineage: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(lineage_data)} GTDB lineages")
    return lineage_data

def get_taxids_from_names_ncbi_fallback(phylum_names: List[str]) -> Dict[str, str]:
    """
    Fallback function to get taxids for phylum names using NCBI taxonkit.
    Used for names that failed GTDB taxdump lookup.

    Args:
        phylum_names: List of phylum names that failed GTDB lookup

    Returns:
        Dictionary mapping phylum names to their NCBI taxids
    """
    if not phylum_names:
        return {}

    print(f"🔄 Trying NCBI fallback for {len(phylum_names)} phylum names...")

    name_to_taxid = {}

    # Create temporary file for names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in phylum_names:
            temp_file.write(f"{name}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit name2taxid {temp_filename} (using NCBI)")

        # Run taxonkit name2taxid command with standard NCBI database
        result = subprocess.run(
            ["taxonkit", "name2taxid", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} NCBI fallback results")

            for line in lines:
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()

                    # Only store if we got a valid taxid
                    if taxid and taxid != name:
                        name_to_taxid[name] = taxid
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No results from NCBI fallback")

    except Exception as e:
        print(f"Error running NCBI fallback: {str(e)}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"✅ NCBI fallback found {len(name_to_taxid)} additional mappings")
    return name_to_taxid

def get_lineages_from_ncbi_taxids(taxids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Get lineages for NCBI taxids using standard NCBI taxonkit lineage -R -t.

    Args:
        taxids: List of NCBI taxids

    Returns:
        Dictionary mapping taxids to (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}

    print(f"Getting NCBI lineages for {len(taxids)} fallback taxids...")

    lineage_data = {}

    # Create temporary file for taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit lineage -R -t {temp_filename} (using NCBI)")

        # Run taxonkit lineage command with standard NCBI database
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} NCBI lineage results")

            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip()
                    lineage_taxids = parts[2].strip()
                    lineage_ranks = parts[3].strip()

                    # Store the lineage data
                    if lineage and lineage != taxid:
                        lineage_data[taxid] = (
                            lineage,
                            lineage_ranks,
                            lineage_taxids
                        )

                        # Show first few examples
                        if i < 3:
                            print(f"  Example: NCBI Taxid {taxid}:")
                            print(f"    Lineage: {lineage}")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No lineage results from NCBI fallback")

    except Exception as e:
        print(f"Error running NCBI lineage: {str(e)}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"✅ Got NCBI lineages for {len(lineage_data)} fallback taxids")
    return lineage_data

def cleanup_lineage_to_phylum_level(lineage: str, lineage_ranks: str, lineage_taxids: str) -> Tuple[str, str, str]:
    """
    Clean up lineage data to only include taxonomic levels up to phylum level.
    Removes any ranks below phylum (class, order, family, genus, species).

    Args:
        lineage: Full lineage string
        lineage_ranks: Full lineage ranks string
        lineage_taxids: Full lineage taxids string

    Returns:
        Tuple of (cleaned_lineage, cleaned_ranks, cleaned_taxids) with only phylum-level and above
    """
    if not lineage or not lineage_ranks or not lineage_taxids:
        return "", "", ""

    try:
        lineage_parts = lineage.split(';')
        ranks_parts = lineage_ranks.split(';')
        taxids_parts = lineage_taxids.split(';')

        # Find the phylum rank position
        phylum_idx = -1
        for i, rank in enumerate(ranks_parts):
            if rank.lower() == 'phylum':
                phylum_idx = i
                break

        # If we found a phylum rank, truncate everything after it
        if phylum_idx >= 0:
            # Keep everything up to and including phylum level
            cleaned_lineage = ';'.join(lineage_parts[:phylum_idx + 1])
            cleaned_ranks = ';'.join(ranks_parts[:phylum_idx + 1])
            cleaned_taxids = ';'.join(taxids_parts[:phylum_idx + 1])

            return cleaned_lineage, cleaned_ranks, cleaned_taxids

        # If no phylum rank found, return original (this preserves existing behavior)
        return lineage, lineage_ranks, lineage_taxids

    except Exception as e:
        print(f"Warning: Error cleaning lineage: {e}")
        return lineage, lineage_ranks, lineage_taxids

def main() -> int:
    """
    Main function to run the script.

    Returns:
        int: Exit code (0 for success, non-zero for error)
    """
    # Parse command line arguments
    args = parse_arguments()
    metrics_list: List[PerformanceMetrics] = []
    script_start_time = time.time()

    try:
        # Set up paths
        bac_file, arc_file, output_file, taxid_map_file = setup_paths(args)

        # Set up logging
        setup_logging(args.log_level, Path(args.output_dir) if args.output_dir else Path(__file__).resolve().parent.parent / "csv_gtdb")

        # print("Loading data files...")

        # Load taxid map if available
        taxid_map = {}
        if taxid_map_file:
            with timer("Load taxid map", metrics_list) as metrics:
                taxid_map = load_taxid_map(taxid_map_file)
                metrics.rows_processed = len(taxid_map)
            # print(f"✅ Loaded taxid mappings")

        # Load bacterial taxonomy file with chunking
        with timer("Load bacterial taxonomy", metrics_list) as metrics:
            bac_df = load_taxonomy_file(bac_file, use_chunks=True)
            metrics.rows_processed = len(bac_df)

        # Load archaeal taxonomy file with chunking
        with timer("Load archaeal taxonomy", metrics_list) as metrics:
            arc_df = load_taxonomy_file(arc_file, use_chunks=True)
            metrics.rows_processed = len(arc_df)

        # Combine dataframes
        with timer("Combine datasets", metrics_list) as metrics:
            df = pd.concat([bac_df, arc_df], ignore_index=True)
            metrics.rows_processed = len(df)
        # print(f"✅ Loaded and combined taxonomy data")

        # Process taxonomy data and add taxids if available
        # print("Processing taxonomy data...")
        with timer("Process taxonomy", metrics_list) as metrics:
            df = process_taxonomy_data(df, taxid_map)
            metrics.rows_processed = len(df)

        # Validate taxonomic hierarchy
        with timer("Validate taxonomy", metrics_list) as metrics:
            validate_taxonomy(df)  # We don't need to store the result
            metrics.rows_processed = len(df)
        # print(f"✅ Processed taxonomy data")

        # Skip all ranks file generation - not needed

        # Group and count
        # print("Counting and processing phyla...")
        with timer("Count genomes", metrics_list) as metrics:
            # Clean phylum names by removing underscores and suffixes (e.g., Bacillota_A -> Bacillota)
            # print("Cleaning phylum names...")
            df['phylum_clean'] = df['phylum'].str.replace(r'_[A-Z]$', '', regex=True)

            # Step 1: Subset to keep one genome per species (best one)
            # print("Subsetting to keep one genome per species (best one)...")
            initial_count = len(df)

            # Debug: Check Wallbacteria species before deduplication
            wallbacteria_before = df[df['phylum'] == 'Wallbacteria']
            # if not wallbacteria_before.empty:
            #     print(f"🔍 DEBUG: Wallbacteria before deduplication: {len(wallbacteria_before)} genomes")
            #     unique_species = wallbacteria_before['species'].unique()
            #     print(f"🔍 DEBUG: Wallbacteria unique species: {len(unique_species)}")
            #     for i, species in enumerate(unique_species):
            #         print(f"  {i+1}. '{species}'")

            species_subset_df = df.drop_duplicates(subset="species", keep="first")

        # Species subset version - no accession file generation
        # print("✅ Species subset version - no accession file generated")

        # Continue with counting
        with timer("Count genomes continued", metrics_list) as metrics:
            # print(f"✅ Reduced from {initial_count:,} genomes to {len(species_subset_df):,} unique species")

            # Debug: Check Wallbacteria species after deduplication
            wallbacteria_after = species_subset_df[species_subset_df['phylum'] == 'Wallbacteria']
            # if not wallbacteria_after.empty:
            #     print(f"🔍 DEBUG: Wallbacteria after deduplication: {len(wallbacteria_after)} genomes")
            #     unique_species_after = wallbacteria_after['species'].unique()
            #     print(f"🔍 DEBUG: Wallbacteria unique species after: {len(unique_species_after)}")
            #     for i, species in enumerate(unique_species_after):
            #         print(f"  {i+1}. '{species}'")

            # Step 2: Group by cleaned phylum and domain, and count unique species
            # print("Counting species by cleaned phylum and domain...")
            counts = species_subset_df.groupby(["phylum_clean", "domain"]).size().reset_index(name="phylum_species_count")
            counts = counts.rename(columns={'phylum_clean': 'phylum'})
            counts = counts.sort_values("phylum_species_count", ascending=False)

            # Get unique phylum names for taxonkit processing
            unique_phyla = counts['phylum'].unique().tolist()
            # print(f"🔍 Found {len(unique_phyla)} unique cleaned phyla")

            # Show some examples of cleaned names
            original_phyla = df['phylum'].unique()[:5]
            cleaned_phyla = df['phylum_clean'].unique()[:5]
            # print("📝 Examples of name cleaning:")
            # for orig, clean in zip(original_phyla, cleaned_phyla):
            #     if orig != clean:
            #         print(f"  {orig} → {clean}")
            #     else:
            #         print(f"  {orig} (no change)")

            # Get taxids for phylum names using GTDB taxdump
            # print("Getting GTDB taxids for phylum names...")
            phylum_to_taxid = get_taxids_from_names_gtdb(unique_phyla)

            # Add taxid column to counts
            counts['taxid'] = counts['phylum'].map(phylum_to_taxid)

            # Log failed GTDB taxid matches
            failed_matches = counts[counts['taxid'].isna()]
            if not failed_matches.empty:
                # print(f"⚠️  Failed to get GTDB taxids for {len(failed_matches)} phyla")

                # Show percentage statistics
                total_phyla = len(counts)
                failed_count = len(failed_matches)
                success_rate = ((total_phyla - failed_count) / total_phyla) * 100
                # print(f"📊 GTDB success rate: {success_rate:.1f}% ({total_phyla - failed_count}/{total_phyla})")

                # Show top failed matches by species count
                top_failed = failed_matches.nlargest(5, 'phylum_species_count')
                # print("🔍 Top failed matches by species count:")
                # for _, row in top_failed.iterrows():
                #     print(f"  {row['phylum']} ({row['domain']}): {row['phylum_species_count']} species")

                # Try NCBI fallback for failed names
                failed_names = failed_matches['phylum'].unique().tolist()
                # print(f"🔄 Trying NCBI fallback for {len(failed_names)} failed names...")
                ncbi_fallback_taxids = get_taxids_from_names_ncbi_fallback(failed_names)

                # Update counts with fallback results
                fallback_success_count = 0
                for idx, row in failed_matches.iterrows():
                    phylum_name = row['phylum']
                    if phylum_name in ncbi_fallback_taxids:
                        counts.at[idx, 'taxid'] = ncbi_fallback_taxids[phylum_name]
                        fallback_success_count += 1

                # print(f"✅ NCBI fallback found {fallback_success_count} additional taxids")

                # Log final failed matches (after fallback) with enhanced statistics
                final_failed = counts[counts['taxid'].isna()]
                if not final_failed.empty:
                    final_log_file = "gtdb_phylum_final_failed_taxids.log"

                    # Calculate comprehensive statistics
                    total_phyla = len(counts)
                    final_failed_count = len(final_failed)
                    final_success_count = total_phyla - final_failed_count
                    final_success_rate = (final_success_count / total_phyla) * 100
                    final_failure_rate = (final_failed_count / total_phyla) * 100

                    # Calculate species impact
                    total_species = counts['phylum_species_count'].sum()
                    failed_species = final_failed['phylum_species_count'].sum()
                    lost_species_percentage = (failed_species / total_species) * 100

                    with open(final_log_file, 'w') as f:
                        f.write("# GTDB Phylum Names That Failed Both GTDB and NCBI Taxid Lookup\n")
                        f.write("# These names could not be found in either GTDB taxdump or NCBI taxonomy\n")
                        f.write("# Generated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")

                        f.write("=== PHYLUM LEVEL MAPPING STATISTICS ===\n")
                        f.write(f"Total phyla processed: {total_phyla}\n")
                        f.write(f"Successfully mapped: {final_success_count} ({final_success_rate:.1f}%)\n")
                        f.write(f"Failed to map: {final_failed_count} ({final_failure_rate:.1f}%)\n\n")

                        f.write("=== SPECIES IMPACT ANALYSIS ===\n")
                        f.write(f"Total species across all phyla: {total_species:,}\n")
                        f.write(f"Species in failed phyla: {failed_species:,}\n")
                        f.write(f"Percentage of species lost due to failed mapping: {lost_species_percentage:.1f}%\n\n")

                        f.write("=== FAILED PHYLA DETAILS ===\n")
                        f.write("# Format: phylum_name,domain,species_count,percentage_of_total_species\n")

                        for _, row in final_failed.iterrows():
                            species_pct = (row['phylum_species_count'] / total_species) * 100
                            f.write(f"{row['phylum']},{row['domain']},{row['phylum_species_count']},{species_pct:.2f}%\n")

                    print(f"📝 Enhanced failed mapping log saved to: {final_log_file}")
                    print(f"📊 Data loss summary: {final_failed_count}/{total_phyla} phyla ({final_failure_rate:.1f}%), {failed_species:,}/{total_species:,} species ({lost_species_percentage:.1f}%)")

                    # Show final statistics
                    total_phyla = len(counts)
                    final_failed_count = len(final_failed)
                    final_success_rate = ((total_phyla - final_failed_count) / total_phyla) * 100
                    print(f"📊 Final success rate: {final_success_rate:.1f}% ({total_phyla - final_failed_count}/{total_phyla})")
                    print(f"📊 GTDB: {len(phylum_to_taxid)}, NCBI fallback: {fallback_success_count}, Final failed: {final_failed_count}")

                    # Show top final failed matches by species count
                    if len(final_failed) > 0:
                        top_final_failed = final_failed.nlargest(5, 'phylum_species_count')
                        print("🔍 Top final failed matches by species count:")
                        for _, row in top_final_failed.iterrows():
                            print(f"  {row['phylum']} ({row['domain']}): {row['phylum_species_count']} species")

            # Get lineages for all available taxids (both GTDB and NCBI fallback)
            print("Getting lineages for all available taxids...")

            # Collect all taxids that we have (both GTDB and NCBI fallback)
            available_taxids = counts[counts['taxid'].notna()]['taxid'].unique().tolist()
            available_taxids = [str(taxid) for taxid in available_taxids]  # Ensure string format

            # print(f"🔍 Found {len(available_taxids)} unique taxids to process for lineages")

            # First try GTDB lineages for GTDB taxids
            gtdb_taxids = []
            ncbi_taxids = []

            # Separate GTDB and NCBI taxids based on source
            for _, row in counts.iterrows():
                if pd.notna(row['taxid']):
                    phylum_name = row['phylum']
                    taxid = str(row['taxid'])

                    # If this phylum was found in GTDB, use GTDB lineage
                    if phylum_name in phylum_to_taxid:
                        gtdb_taxids.append(taxid)
                    else:
                        # Otherwise it's from NCBI fallback
                        ncbi_taxids.append(taxid)

            # Remove duplicates
            gtdb_taxids = list(set(gtdb_taxids))
            ncbi_taxids = list(set(ncbi_taxids))

            print(f"📊 Taxid breakdown: {len(gtdb_taxids)} GTDB, {len(ncbi_taxids)} NCBI fallback")

            # Get GTDB lineages
            gtdb_lineage_data = {}
            if gtdb_taxids:
                print(f"Getting GTDB lineages for {len(gtdb_taxids)} GTDB taxids...")
                gtdb_lineage_data = get_lineages_from_gtdb_taxids(gtdb_taxids)

            # Get NCBI lineages for fallback taxids
            ncbi_lineage_data = {}
            if ncbi_taxids:
                print(f"Getting NCBI lineages for {len(ncbi_taxids)} NCBI fallback taxids...")
                ncbi_lineage_data = get_lineages_from_ncbi_taxids(ncbi_taxids)

            # Combine lineage data
            all_lineage_data = {**gtdb_lineage_data, **ncbi_lineage_data}
            print(f"✅ Total lineages obtained: {len(all_lineage_data)}")

            # Add lineage columns using taxid lookup with phylum validation
            lineages = []
            lineage_ranks = []
            lineage_taxids = []

            for _, row in counts.iterrows():
                if pd.notna(row['taxid']):
                    taxid = str(row['taxid'])
                    phylum_name = row['phylum']

                    if taxid in all_lineage_data:
                        full_lineage, full_ranks, full_taxids = all_lineage_data[taxid]

                        # Clean up lineage to only include phylum level and above
                        cleaned_lineage, cleaned_ranks, cleaned_taxids = cleanup_lineage_to_phylum_level(
                            full_lineage, full_ranks, full_taxids
                        )

                        lineages.append(cleaned_lineage)
                        lineage_ranks.append(cleaned_ranks)
                        lineage_taxids.append(cleaned_taxids)
                    else:
                        lineages.append("")
                        lineage_ranks.append("")
                        lineage_taxids.append("")
                else:
                    lineages.append("")
                    lineage_ranks.append("")
                    lineage_taxids.append("")

            counts['lineage'] = lineages
            counts['lineage_ranks'] = lineage_ranks
            counts['lineage_taxids'] = lineage_taxids

            # Log lineage success statistics
            lineage_success = len([l for l in lineages if l])
            total_with_taxids = len([t for t in counts['taxid'] if pd.notna(t)])
            if total_with_taxids > 0:
                lineage_success_rate = (lineage_success / total_with_taxids) * 100
                print(f"📊 Lineage success rate: {lineage_success_rate:.1f}% ({lineage_success}/{total_with_taxids})")
            else:
                print("📊 No taxids available for lineage inference")

            # Ensure final counts are sorted by species count
            counts = counts.sort_values("phylum_species_count", ascending=False)

            metrics.rows_processed = len(counts)

        # Save to output directory
        # print("Saving summary file...")
        with timer("Save summary file", metrics_list) as metrics:
            counts.to_csv(output_file, index=False)
            metrics.rows_processed = len(counts)
        # print(f"✅ Saved summary file")

        # Print preview
        print("\nTop 10 phyla by species count:")
        for _, row in counts.head(10).iterrows():
            print(f"  {row['phylum']} ({row['domain']}): {row['phylum_species_count']} species")

        # Print simple performance summary
        total_time = time.time() - script_start_time
        print(f"\nTotal execution time: {total_time:.2f}s")

        return 0

    except FileNotFoundError as e:
        print(f"❌ Error: {str(e)}")
        logging.error(f"File not found: {str(e)}")
        return 1
    except pd.errors.EmptyDataError as e:
        print(f"❌ Error: {str(e)}")
        logging.error(f"Empty data error: {str(e)}")
        return 1
    except AssertionError as e:
        print(f"❌ Data validation error: {str(e)}")
        logging.error(f"Data validation error: {str(e)}")
        return 1
    except ValueError as e:
        print(f"❌ Value error: {str(e)}")
        logging.error(f"Value error: {str(e)}")
        return 1
    except Exception as e:
        print(f"❌ Unexpected error: {str(e)}")
        logging.error(f"Unexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        logging.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())
