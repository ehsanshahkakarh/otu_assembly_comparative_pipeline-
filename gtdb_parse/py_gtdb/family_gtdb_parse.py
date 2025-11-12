#!/usr/bin/env python3
"""
GTDB Family Parser

This script processes GTDB taxonomy files to extract family-level information.
It counts the number of genomes per family and creates output files with
both summary statistics and detailed accession information.

Input:
- GTDB bacterial taxonomy file (00bac120_taxonomy.tsv)
- GTDB archaeal taxonomy file (00ar53_taxonomy.tsv)

Output:
- Summary CSV with family counts (gtdb_family_counts.csv)
- Detailed CSV with accessions (gtdb_family_with_accessions.csv)

Usage:
    python family_gtdb_parse.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
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
SUMMARY_OUTPUT_FILENAME = "gtdb_family_counts.csv"
DETAILED_OUTPUT_FILENAME = "gtdb_family_with_accessions.csv"
# Removed ALL_RANKS_OUTPUT_FILENAME - not needed
TAXID_MAP_FILENAME = "taxid.map"

# Regex patterns for extracting taxonomic ranks
DOMAIN_REGEX = r"d__([A-Za-z0-9_\-]+)"
PHYLUM_REGEX = r"p__([A-Za-z0-9_\-]+)"
CLASS_REGEX = r"c__([A-Za-z0-9_\-]+)"
ORDER_REGEX = r"o__([A-Za-z0-9_\-]+)"
FAMILY_REGEX = r"f__([A-Za-z0-9_\-]+)"
GENUS_REGEX = r"g__([A-Za-z0-9_\-]+)"
SPECIES_REGEX = r"s__([A-Za-z0-9_\-]+)"

# Regex for cleaning accession IDs - remove prefixes like GB_ or RS_
ACCESSION_CLEAN_REGEX = r"^[A-Z]{2}_"

# Processing parameters
CHUNK_SIZE = 100000  # Number of rows to process at once
MIN_GENOME_THRESHOLD = 20  # Minimum number of genomes required for a family to be included

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
    log_file = log_dir / f"gtdb_family_parse_{timestamp}.log"

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
    parser = argparse.ArgumentParser(description="Parse GTDB taxonomy files for family-level information")
    parser.add_argument("--input-dir", help="Directory containing input files (default: script directory)")
    parser.add_argument("--output-dir", help="Directory for output files (default: ../csv_gtdb)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--log-level",
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       default='INFO',
                       help="Set the logging level")
    parser.add_argument("--min-genomes",
                       type=int,
                       default=MIN_GENOME_THRESHOLD,
                       help=f"Minimum number of genomes required for a family to be included (default: {MIN_GENOME_THRESHOLD})")
    return parser.parse_args()

def setup_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path, Path]:
    """
    Set up input and output file paths.

    Args:
        args (argparse.Namespace): Command line arguments

    Returns:
        Tuple[Path, Path, Path, Path, Path]: Paths for bacterial taxonomy file, archaeal taxonomy file,
                                           summary output file, detailed output file, and taxid map file
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
    detailed_output_file = output_dir / DETAILED_OUTPUT_FILENAME

    # Define taxid map file path (in the script directory)
    taxid_map_file = script_dir / TAXID_MAP_FILENAME

    # Check if taxid map file exists
    if not taxid_map_file.exists():
        print(f"⚠️ Warning: Taxid map file not found at {taxid_map_file}")
        print("  Taxid information will not be included in the output.")
        taxid_map_file = None

    return bac_file, arc_file, output_file, detailed_output_file, taxid_map_file

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
                    # Remove the prefix (GB_ or RS_) to match our cleaned accessions
                    accession_clean = re.sub(ACCESSION_CLEAN_REGEX, "", accession).strip()
                    taxid_map[accession_clean] = taxid
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

        # Get total number of lines for progress bar
        with open(file_path, 'r') as f:
            total_lines = sum(1 for _ in f)

        # Read file in chunks with progress bar
        with tqdm(total=total_lines, desc=f"Loading {file_path.name}", unit="lines") as pbar:
            for chunk in pd.read_csv(file_path, sep="\t", header=None, names=["accession", "taxonomy"], chunksize=CHUNK_SIZE):
                chunks.append(chunk)
                total_rows += len(chunk)
                pbar.update(len(chunk))

        # Combine all chunks
        df = pd.concat(chunks, ignore_index=True)
        logging.info(f"Loaded {total_rows:,} rows from {file_path.name}")
        return df

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading {file_path}: {str(e)}")
        raise

def process_taxonomy_data(df: pd.DataFrame, taxid_map: Dict[str, str] = None) -> pd.DataFrame:
    """
    Process taxonomy data to extract family information and other ranks.

    Args:
        df (pd.DataFrame): Input DataFrame with taxonomy data
        taxid_map (Dict[str, str], optional): Dictionary mapping accession to taxid

    Returns:
        pd.DataFrame: Processed DataFrame with family information and other ranks
    """
    # Extract all taxonomic ranks
    df["domain"] = df["taxonomy"].str.extract(DOMAIN_REGEX)
    df["phylum"] = df["taxonomy"].str.extract(PHYLUM_REGEX)
    df["class"] = df["taxonomy"].str.extract(CLASS_REGEX)
    df["order"] = df["taxonomy"].str.extract(ORDER_REGEX)
    df["family"] = df["taxonomy"].str.extract(FAMILY_REGEX)
    df["genus"] = df["taxonomy"].str.extract(GENUS_REGEX)
    df["species"] = df["taxonomy"].str.extract(SPECIES_REGEX)

    # Clean accession by removing prefixes and whitespace
    df["accession_clean"] = df["accession"].str.replace(ACCESSION_CLEAN_REGEX, "", regex=True).str.strip()

    # Add taxid if available
    if taxid_map:
        df["taxid"] = df["accession_clean"].map(taxid_map)

    # Fill missing values with "Unclassified"
    for col in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        df[col] = df[col].fillna("Unclassified")

    return df

def standardize_family_name(family: str) -> str:
    """
    Standardize family name format.

    Args:
        family (str): Family name to standardize

    Returns:
        str: Standardized family name
    """
    if pd.isna(family) or family == "Unclassified":
        return "Unclassified"

    # Remove any trailing numbers or special characters
    family = re.sub(r'[^A-Za-z0-9_\-]+$', '', family)
    
    # Ensure proper capitalization
    family = family.strip()
    if family:
        family = family[0].upper() + family[1:].lower()
    
    return family

def merge_related_families(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge related families based on naming patterns.

    Args:
        counts_df (pd.DataFrame): DataFrame with family counts

    Returns:
        pd.DataFrame: DataFrame with merged families
    """
    # Create a copy to avoid modifying the original
    df = counts_df.copy()

    # Define patterns for related families
    patterns = {
        r'^Enterobacteriaceae.*$': 'Enterobacteriaceae',
        r'^Pseudomonadaceae.*$': 'Pseudomonadaceae',
        r'^Bacillaceae.*$': 'Bacillaceae',
        r'^Staphylococcaceae.*$': 'Staphylococcaceae',
        r'^Streptococcaceae.*$': 'Streptococcaceae',
        r'^Lactobacillaceae.*$': 'Lactobacillaceae',
        r'^Clostridiaceae.*$': 'Clostridiaceae',
        r'^Bacteroidaceae.*$': 'Bacteroidaceae',
        r'^Prevotellaceae.*$': 'Prevotellaceae',
        r'^Ruminococcaceae.*$': 'Ruminococcaceae',
        r'^Lachnospiraceae.*$': 'Lachnospiraceae',
        r'^Bifidobacteriaceae.*$': 'Bifidobacteriaceae',
        r'^Corynebacteriaceae.*$': 'Corynebacteriaceae',
        r'^Mycobacteriaceae.*$': 'Mycobacteriaceae',
        r'^Actinomycetaceae.*$': 'Actinomycetaceae',
        r'^Streptomycetaceae.*$': 'Streptomycetaceae',
        r'^Micrococcaceae.*$': 'Micrococcaceae',
        r'^Propionibacteriaceae.*$': 'Propionibacteriaceae',
        r'^Bifidobacteriaceae.*$': 'Bifidobacteriaceae',
        r'^Corynebacteriaceae.*$': 'Corynebacteriaceae',
        r'^Mycobacteriaceae.*$': 'Mycobacteriaceae',
        r'^Actinomycetaceae.*$': 'Actinomycetaceae',
        r'^Streptomycetaceae.*$': 'Streptomycetaceae',
        r'^Micrococcaceae.*$': 'Micrococcaceae',
        r'^Propionibacteriaceae.*$': 'Propionibacteriaceae'
    }

    # Apply standardization
    df['family'] = df['family'].apply(standardize_family_name)

    # Apply pattern matching and merging
    for pattern, replacement in patterns.items():
        mask = df['family'].str.match(pattern, na=False)
        df.loc[mask, 'family'] = replacement

    # Group by the standardized family names and sum the counts
    df = df.groupby(['family', 'domain']).agg({
        'gtdb_genome_count': 'sum',
        'accession_clean': lambda x: list(set([item for sublist in x for item in sublist])),
        'taxid': lambda x: list(set([item for sublist in x for item in sublist if pd.notna(item)]))
    }).reset_index()

    return df

def validate_taxonomy(df: pd.DataFrame) -> bool:
    """
    Validate taxonomy data for consistency and completeness.

    Args:
        df (pd.DataFrame): DataFrame to validate

    Returns:
        bool: True if validation passes, False otherwise
    """
    try:
        # Check for required columns
        required_columns = ["accession_clean", "domain", "phylum", "class", "order", "family", "genus", "species"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logging.error(f"Missing required columns: {missing_columns}")
            return False

        # Check for missing values in critical columns
        critical_columns = ["accession_clean", "domain", "family"]
        for col in critical_columns:
            null_count = df[col].isna().sum()
            if null_count > 0:
                logging.error(f"Found {null_count} null values in {col}")
                return False

        # Validate domain values
        valid_domains = ["Bacteria", "Archaea"]
        invalid_domains = df[~df["domain"].isin(valid_domains)]["domain"].unique()
        if len(invalid_domains) > 0:
            logging.error(f"Invalid domain values found: {invalid_domains}")
            return False

        # Check for empty strings in critical columns
        for col in critical_columns:
            empty_count = (df[col] == "").sum()
            if empty_count > 0:
                logging.error(f"Found {empty_count} empty strings in {col}")
                return False

        return True

    except Exception as e:
        logging.error(f"Error during validation: {str(e)}")
        return False

def get_taxids_from_names_gtdb(family_names: List[str]) -> Dict[str, str]:
    """
    Get taxids for family names using GTDB taxdump via taxonkit name2taxid.

    Args:
        family_names: List of family names

    Returns:
        Dictionary mapping family names to their GTDB taxids
    """
    if not family_names:
        return {}

    print(f"Getting GTDB taxids for {len(family_names)} family names...")

    # Set up environment to use GTDB taxdump
    script_dir = Path(__file__).resolve().parent
    gtdb_taxdump_dir = script_dir.parent / "taxdump_gtdp" / "gtdb-taxdump-R226"
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(gtdb_taxdump_dir)

    name_to_taxid = {}

    # Create temporary file for family names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in family_names:
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

    print(f"Successfully processed {len(name_to_taxid)} family name->GTDB taxid mappings")
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

def get_taxids_from_names_ncbi_fallback(family_names: List[str]) -> Dict[str, str]:
    """
    Fallback function to get taxids for family names using NCBI taxonkit.
    Used for names that failed GTDB taxdump lookup.

    Args:
        family_names: List of family names that failed GTDB lookup

    Returns:
        Dictionary mapping family names to their NCBI taxids
    """
    if not family_names:
        return {}

    print(f"🔄 Fallback: Getting NCBI taxids for {len(family_names)} failed GTDB names...")

    name_to_taxid = {}

    # Create temporary file for family names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in family_names:
            temp_file.write(f"{name}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit name2taxid {temp_filename} (using NCBI)")

        # Run taxonkit name2taxid command with default NCBI database
        result = subprocess.run(
            ["taxonkit", "name2taxid", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"NCBI Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"NCBI Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} NCBI fallback results")

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
                            print(f"  Fallback success: {name} -> {taxid} (NCBI)")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No results from NCBI fallback taxonkit")

    except Exception as e:
        print(f"Error running NCBI fallback taxonkit: {e}")
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

        print(f"NCBI Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"NCBI Taxonkit stderr: {result.stderr}")

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
        print(f"Error running NCBI fallback lineage: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"✅ Got NCBI lineages for {len(lineage_data)} fallback taxids")
    return lineage_data

def cleanup_lineage_to_family_level(lineage: str, lineage_ranks: str, lineage_taxids: str) -> Tuple[str, str, str]:
    """
    Clean up lineage data to only include taxonomic levels up to family level.
    Removes any ranks below family (genus, species).

    Args:
        lineage: Full lineage string
        lineage_ranks: Full lineage ranks string
        lineage_taxids: Full lineage taxids string

    Returns:
        Tuple of (cleaned_lineage, cleaned_ranks, cleaned_taxids) with only family-level and above
    """
    if not lineage or not lineage_ranks or not lineage_taxids:
        return "", "", ""

    try:
        lineage_parts = lineage.split(';')
        ranks_parts = lineage_ranks.split(';')
        taxids_parts = lineage_taxids.split(';')

        # Find the family rank position
        family_idx = -1
        for i, rank in enumerate(ranks_parts):
            if rank.lower() == 'family':
                family_idx = i
                break

        # If we found a family rank, truncate everything after it
        if family_idx >= 0:
            # Keep everything up to and including family level
            cleaned_lineage = ';'.join(lineage_parts[:family_idx + 1])
            cleaned_ranks = ';'.join(ranks_parts[:family_idx + 1])
            cleaned_taxids = ';'.join(taxids_parts[:family_idx + 1])

            return cleaned_lineage, cleaned_ranks, cleaned_taxids

        # If no family rank found, return original (this preserves existing behavior)
        return lineage, lineage_ranks, lineage_taxids

    except Exception as e:
        logging.warning(f"Error cleaning lineage: {e}")
        return lineage, lineage_ranks, lineage_taxids

def main() -> int:
    """
    Main function to run the script.

    Returns:
        int: Exit code (0 for success, non-zero for error)
    """
    # Parse command line arguments
    args = parse_arguments()

    # Set up logging
    setup_logging(args.log_level, Path(args.output_dir) if args.output_dir else Path(__file__).resolve().parent.parent / "csv_gtdb")

    # Initialize performance metrics list
    metrics_list = []

    try:
        with timer("Total execution", metrics_list):
            # Set up paths
            bac_file, arc_file, output_file, detailed_output_file, taxid_map_file = setup_paths(args)

            # logging.info(f"📂 Input files: {bac_file}, {arc_file}")
            # logging.info(f"📂 Output files: {output_file}, {detailed_output_file}")

            # Load taxid map if available
            taxid_map = load_taxid_map(taxid_map_file)

            # Load bacterial taxonomy
            with timer("Loading bacterial taxonomy", metrics_list) as metrics:
                bac_df = load_taxonomy_file(bac_file)
                metrics.rows_processed = len(bac_df)
                # logging.info(f"✅ Loaded {len(bac_df):,} bacterial entries")

            # Load archaeal taxonomy
            with timer("Loading archaeal taxonomy", metrics_list) as metrics:
                ar_df = load_taxonomy_file(arc_file)
                metrics.rows_processed = len(ar_df)
                # logging.info(f"✅ Loaded {len(ar_df):,} archaeal entries")

            # Process taxonomy data
            with timer("Processing taxonomy data", metrics_list) as metrics:
                bac_processed = process_taxonomy_data(bac_df, taxid_map)
                ar_processed = process_taxonomy_data(ar_df, taxid_map)
                metrics.rows_processed = len(bac_processed) + len(ar_processed)

            # Combine processed data
            df = pd.concat([bac_processed, ar_processed], ignore_index=True)
            # logging.info(f"✅ Combined: {len(df):,} total entries")

            # Validate taxonomy data
            if not validate_taxonomy(df):
                logging.error("❌ Taxonomy validation failed")
                return 1

            # Skip all ranks file generation - not needed

            # Process counts
            with timer("Processing counts", metrics_list) as metrics:
                # Drop entries missing family
                df = df.dropna(subset=["family"])
                # logging.info(f"✅ After dropping null families: {len(df):,} rows")

                # Enable tqdm for pandas operations
                tqdm.pandas(desc="Processing families")

                # Clean family names by removing underscores and suffixes (e.g., Enterobacteriaceae_A -> Enterobacteriaceae)
                # logging.info("Cleaning family names...")
                df['family_clean'] = df['family'].str.replace(r'_[A-Z]$', '', regex=True)

                # Step 1: Count all genomes (no species subsetting)
                # logging.info("Counting all genomes (no species subsetting)...")
                all_genomes_df = df.copy()
                # logging.info(f"✅ Processing {initial_count:,} total genomes (no species subsetting)")

            # Save detailed file with all genomes (no species subsetting)
            with timer("Saving detailed file", metrics_list) as metrics:
                detailed_df = all_genomes_df.dropna(subset=["family"])
                # Sort by domain and family for detailed file
                detailed_df = detailed_df.sort_values(["domain", "family", "taxid"], ascending=[True, True, True])
                detailed_df[["accession_clean", "family", "domain", "taxid"]].to_csv(detailed_output_file, index=False)
                metrics.rows_processed = len(detailed_df)
                logging.info(f"✅ Saved detailed file with {len(detailed_df):,} entries")

            # Continue processing counts
            with timer("Processing counts continued", metrics_list) as metrics:

                # Step 2: Group by cleaned family and domain, and count all genomes
                # logging.info("Counting genomes by cleaned family and domain...")
                counts = all_genomes_df.groupby(["family_clean", "domain"]).size().reset_index(name="family_counts")
                counts = counts.rename(columns={'family_clean': 'family'})
                counts = counts.sort_values("family_counts", ascending=False)

                # Get unique family names for taxonkit processing
                unique_families = counts['family'].unique().tolist()
                # logging.info(f"🔍 Found {len(unique_families)} unique cleaned families")

                # Skip name cleaning examples - not needed for minimal output

                # Get taxids for family names using GTDB taxdump
                logging.info("Getting GTDB taxids for family names...")
                family_to_taxid = get_taxids_from_names_gtdb(unique_families)

                # Add taxid column to counts
                counts['taxid'] = counts['family'].map(family_to_taxid)

                # Try NCBI fallback for failed GTDB matches
                failed_matches = counts[counts['taxid'].isna()]
                if not failed_matches.empty:
                    logging.info(f"⚠️  Failed to get GTDB taxids for {len(failed_matches)} families")

                    # Try NCBI fallback for failed names
                    failed_names = failed_matches['family'].unique().tolist()
                    logging.info(f"🔄 Trying NCBI fallback for {len(failed_names)} failed names...")
                    ncbi_fallback_taxids = get_taxids_from_names_ncbi_fallback(failed_names)

                    # Update counts with fallback results
                    fallback_success_count = 0
                    for idx, row in failed_matches.iterrows():
                        family_name = row['family']
                        if family_name in ncbi_fallback_taxids:
                            counts.at[idx, 'taxid'] = ncbi_fallback_taxids[family_name]
                            fallback_success_count += 1

                    logging.info(f"✅ NCBI fallback found {fallback_success_count} additional taxids")

                    # Log final failed matches (after fallback) with enhanced statistics
                    final_failed = counts[counts['taxid'].isna()]
                    if not final_failed.empty:
                        log_file = "gtdb_family_final_failed_taxids.log"

                        # Calculate comprehensive statistics
                        total_families = len(counts)
                        final_failed_count = len(final_failed)
                        final_success_count = total_families - final_failed_count
                        final_success_rate = (final_success_count / total_families) * 100
                        final_failure_rate = (final_failed_count / total_families) * 100

                        # Calculate genome impact
                        total_genomes = counts['family_counts'].sum()
                        failed_genomes = final_failed['family_counts'].sum()
                        lost_genomes_percentage = (failed_genomes / total_genomes) * 100

                        with open(log_file, 'w') as f:
                            f.write("# GTDB Family Names That Failed Both GTDB and NCBI Taxid Lookup\n")
                            f.write("# These names could not be found in either GTDB taxdump or NCBI taxonomy\n")
                            f.write("# Generated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")

                            f.write("=== FAMILY LEVEL MAPPING STATISTICS ===\n")
                            f.write(f"Total families processed: {total_families}\n")
                            f.write(f"Successfully mapped: {final_success_count} ({final_success_rate:.1f}%)\n")
                            f.write(f"Failed to map: {final_failed_count} ({final_failure_rate:.1f}%)\n\n")

                            f.write("=== GENOME IMPACT ANALYSIS ===\n")
                            f.write(f"Total genomes across all families: {total_genomes:,}\n")
                            f.write(f"Genomes in failed families: {failed_genomes:,}\n")
                            f.write(f"Percentage of genomes lost due to failed mapping: {lost_genomes_percentage:.1f}%\n\n")

                            f.write("=== FAILED FAMILIES DETAILS ===\n")
                            f.write("# Format: family_name,domain,genome_count,percentage_of_total_genomes\n")

                            for _, row in final_failed.iterrows():
                                genome_pct = (row['family_counts'] / total_genomes) * 100
                                f.write(f"{row['family']},{row['domain']},{row['family_counts']},{genome_pct:.2f}%\n")

                        logging.info(f"📝 Enhanced failed mapping log saved to: {log_file}")
                        logging.info(f"📊 Data loss summary: {final_failed_count}/{total_families} families ({final_failure_rate:.1f}%), {failed_genomes:,}/{total_genomes:,} genomes ({lost_genomes_percentage:.1f}%)")

                    # Show final statistics
                    total_families = len(counts)
                    final_failed_count = len(final_failed)
                    final_success_rate = ((total_families - final_failed_count) / total_families) * 100
                    logging.info(f"📊 Final success rate: {final_success_rate:.1f}% ({total_families - final_failed_count}/{total_families})")
                    logging.info(f"📊 GTDB: {len(family_to_taxid)}, NCBI fallback: {fallback_success_count}, Final failed: {final_failed_count}")

                    # Show top final failed matches by genome count
                    if not final_failed.empty:
                        top_failed = final_failed.nlargest(5, 'family_counts')
                        logging.info("🔍 Top final failed matches by genome count:")
                        for _, row in top_failed.iterrows():
                            logging.info(f"  {row['family']} ({row['domain']}): {row['family_counts']} genomes")

                # Get lineages for all available taxids (both GTDB and NCBI fallback)
                logging.info("Getting lineages for all available taxids...")

                # Collect all taxids that we have (both GTDB and NCBI fallback)
                available_taxids = counts[counts['taxid'].notna()]['taxid'].unique().tolist()
                available_taxids = [str(taxid) for taxid in available_taxids]  # Ensure string format

                # logging.info(f"🔍 Found {len(available_taxids)} unique taxids to process for lineages")

                # First try GTDB lineages for GTDB taxids
                gtdb_taxids = []
                ncbi_taxids = []

                # Separate GTDB and NCBI taxids based on source
                for _, row in counts.iterrows():
                    if pd.notna(row['taxid']):
                        family_name = row['family']
                        taxid = str(row['taxid'])

                        # If this family was found in GTDB, use GTDB lineage
                        if family_name in family_to_taxid:
                            gtdb_taxids.append(taxid)
                        else:
                            # Otherwise it's from NCBI fallback
                            ncbi_taxids.append(taxid)

                # Remove duplicates
                gtdb_taxids = list(set(gtdb_taxids))
                ncbi_taxids = list(set(ncbi_taxids))

                logging.info(f"📊 Taxid breakdown: {len(gtdb_taxids)} GTDB, {len(ncbi_taxids)} NCBI fallback")

                # Get GTDB lineages
                gtdb_lineage_data = {}
                if gtdb_taxids:
                    logging.info(f"Getting GTDB lineages for {len(gtdb_taxids)} GTDB taxids...")
                    gtdb_lineage_data = get_lineages_from_gtdb_taxids(gtdb_taxids)

                # Get NCBI lineages for fallback taxids
                ncbi_lineage_data = {}
                if ncbi_taxids:
                    logging.info(f"Getting NCBI lineages for {len(ncbi_taxids)} NCBI fallback taxids...")
                    ncbi_lineage_data = get_lineages_from_ncbi_taxids(ncbi_taxids)

                # Combine lineage data
                all_lineage_data = {**gtdb_lineage_data, **ncbi_lineage_data}
                logging.info(f"✅ Total lineages obtained: {len(all_lineage_data)}")

                # Add lineage columns using taxid lookup with family validation
                lineages = []
                lineage_ranks = []
                lineage_taxids = []

                for _, row in counts.iterrows():
                    if pd.notna(row['taxid']):
                        taxid = str(row['taxid'])
                        family_name = row['family']

                        if taxid in all_lineage_data:
                            full_lineage, full_ranks, full_taxids = all_lineage_data[taxid]

                            # Clean up lineage to only include family level and above
                            cleaned_lineage, cleaned_ranks, cleaned_taxids = cleanup_lineage_to_family_level(
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
                    logging.info(f"📊 Lineage success rate: {lineage_success_rate:.1f}% ({lineage_success}/{total_with_taxids})")
                else:
                    logging.info("📊 No taxids available for lineage inference")

                # Additional grouping by last taxid to merge GTDB-specific family names
                logging.info("Grouping by last lineage taxid to merge related families...")

                # Extract last taxid from lineage_taxids for grouping
                def get_last_taxid(lineage_taxids_str):
                    if pd.isna(lineage_taxids_str) or lineage_taxids_str == "":
                        return None
                    taxids = str(lineage_taxids_str).split(';')
                    return taxids[-1] if taxids else None

                counts['last_taxid'] = counts['lineage_taxids'].apply(get_last_taxid)

                # Group by last_taxid and domain, summing the genome counts
                grouped_counts = []

                for (last_taxid, domain), group in counts.groupby(['last_taxid', 'domain']):
                    if pd.isna(last_taxid):
                        # Keep ungrouped entries as-is
                        for _, row in group.iterrows():
                            grouped_counts.append(row.to_dict())
                    else:
                        # Merge entries with same last_taxid
                        total_genomes = group['family_counts'].sum()

                        # Use the most common family name from the lineage, or first one
                        representative_row = group.iloc[0]

                        # Extract family name from lineage if available
                        lineage = representative_row['lineage']
                        if lineage and ';' in lineage:
                            lineage_parts = lineage.split(';')
                            # Family should be the 4th part (superkingdom;phylum;class;order;family)
                            if len(lineage_parts) >= 5:
                                family_name = lineage_parts[4]  # Fifth part should be family
                            else:
                                family_name = representative_row['family']
                        else:
                            family_name = representative_row['family']

                        grouped_entry = {
                            'family': family_name,
                            'domain': domain,
                            'family_counts': total_genomes,
                            'taxid': last_taxid,
                            'lineage': representative_row['lineage'],
                            'lineage_ranks': representative_row['lineage_ranks'],
                            'lineage_taxids': representative_row['lineage_taxids'],
                            'last_taxid': last_taxid
                        }
                        grouped_counts.append(grouped_entry)

                # Convert back to DataFrame
                counts = pd.DataFrame(grouped_counts)

                # Remove the temporary last_taxid column
                counts = counts.drop('last_taxid', axis=1)

                # Ensure final counts are sorted by genome count
                counts = counts.sort_values("family_counts", ascending=False)

                logging.info(f"✅ After grouping by last taxid: {len(counts)} unique families")

                metrics.rows_processed = len(counts)

            # Save summary
            with timer("Saving summary", metrics_list) as metrics:
                counts.to_csv(output_file, index=False)
                metrics.rows_processed = len(counts)
                # logging.info(f"✅ Saved {len(counts):,} family-domain combinations")

            # Print preview
            logging.info("\n📈 Top 10 families by genome count:")
            for _, row in counts.head(10).iterrows():
                logging.info(f"  {row['family']} ({row['domain']}): {row['family_counts']:,} genomes")

            # Print performance metrics
            logging.info("\n⏱️ Performance Metrics:")
            for metric in metrics_list:
                logging.info(f"  {metric}")

            return 0

    except FileNotFoundError as e:
        logging.error(f"❌ Error: {str(e)}")
        return 1
    except ValueError as e:
        logging.error(f"❌ Error: {str(e)}")
        return 1
    except Exception as e:
        logging.error(f"❌ Unexpected error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
