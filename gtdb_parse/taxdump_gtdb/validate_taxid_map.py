#!/usr/bin/env python3
"""
Validate GTDB taxid.map file against taxonomy dump files.

This script checks:
1. Each taxid in taxid.map exists in nodes.dmp
2. Each taxid has a corresponding name in names.dmp
3. Accession IDs match the expected format
4. Reports any inconsistencies or missing data

Usage:
    python validate_taxid_map.py [--taxid_map PATH] [--nodes PATH] [--names PATH] [--metadata PATH]

Options:
    --taxid_map PATH    Path to taxid.map file [default: gtdb-taxdump-R226/taxid.map]
    --nodes PATH        Path to nodes.dmp file [default: gtdb-taxdump-R226/nodes.dmp]
    --names PATH        Path to names.dmp file [default: gtdb-taxdump-R226/names.dmp]
    --metadata PATH     Path to GTDB metadata file [default: None]
"""

import sys
import re
import argparse
import logging
from pathlib import Path
from typing import Dict, Set, List, Tuple
from datetime import datetime
from tqdm import tqdm

# Logging will be set up in main() after parsing arguments

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Validate GTDB taxid.map file against taxonomy dump files")

    # Get the script's directory
    script_dir = Path(__file__).resolve().parent

    # Set default paths relative to the current directory
    default_taxid_map = str(script_dir / "gtdb-taxdump-R226" / "taxid.map")
    default_nodes = str(script_dir / "gtdb-taxdump-R226" / "nodes.dmp")
    default_names = str(script_dir / "gtdb-taxdump-R226" / "names.dmp")

    parser.add_argument("--taxid_map", type=str, default=default_taxid_map,
                        help=f"Path to taxid.map file [default: {default_taxid_map}]")
    parser.add_argument("--nodes", type=str, default=default_nodes,
                        help=f"Path to nodes.dmp file [default: {default_nodes}]")
    parser.add_argument("--names", type=str, default=default_names,
                        help=f"Path to names.dmp file [default: {default_names}]")
    parser.add_argument("--metadata", type=str, default=None,
                        help="Path to GTDB metadata file (optional)")
    parser.add_argument("--output_dir", type=str, default=str(script_dir),
                        help=f"Directory to save log file [default: {script_dir}]")
    return parser.parse_args()

def load_taxid_map(file_path: str) -> Dict[str, str]:
    """
    Load taxid.map file into a dictionary.

    Args:
        file_path: Path to taxid.map file

    Returns:
        Dictionary mapping accession to taxid
    """
    logging.info(f"Loading taxid.map from {file_path}")
    taxid_map = {}

    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(tqdm(f, desc="Loading taxid.map"), 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    logging.warning(f"Line {line_num} has incorrect format: {line}")
                    continue

                accession, taxid = parts
                taxid_map[accession] = taxid

        logging.info(f"Loaded {len(taxid_map)} entries from taxid.map")
        return taxid_map

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error loading taxid.map: {str(e)}")
        sys.exit(1)

def load_nodes(file_path: str) -> Set[str]:
    """
    Load nodes.dmp file and extract all taxids.

    Args:
        file_path: Path to nodes.dmp file

    Returns:
        Set of all taxids in nodes.dmp
    """
    logging.info(f"Loading nodes.dmp from {file_path}")
    taxids = set()

    try:
        with open(file_path, 'r') as f:
            for line in tqdm(f, desc="Loading nodes.dmp"):
                parts = line.strip().split('|')
                if len(parts) > 1:
                    taxid = parts[0].strip()
                    taxids.add(taxid)

        logging.info(f"Loaded {len(taxids)} taxids from nodes.dmp")
        return taxids

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error loading nodes.dmp: {str(e)}")
        sys.exit(1)

def load_names(file_path: str) -> Set[str]:
    """
    Load names.dmp file and extract all taxids.

    Args:
        file_path: Path to names.dmp file

    Returns:
        Set of all taxids in names.dmp
    """
    logging.info(f"Loading names.dmp from {file_path}")
    taxids = set()

    try:
        with open(file_path, 'r') as f:
            for line in tqdm(f, desc="Loading names.dmp"):
                parts = line.strip().split('|')
                if len(parts) > 1:
                    taxid = parts[0].strip()
                    taxids.add(taxid)

        logging.info(f"Loaded {len(taxids)} taxids from names.dmp")
        return taxids

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error loading names.dmp: {str(e)}")
        sys.exit(1)

def load_metadata_accessions(file_path: str) -> Set[str]:
    """
    Load GTDB metadata file and extract all accessions.

    Args:
        file_path: Path to GTDB metadata file

    Returns:
        Set of all accessions in metadata
    """
    if not file_path:
        logging.info("No metadata file provided, skipping metadata validation")
        return set()

    logging.info(f"Loading metadata from {file_path}")
    accessions = set()

    try:
        with open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            accession_idx = header.index('accession')

            for line in tqdm(f, desc="Loading metadata"):
                parts = line.strip().split('\t')
                if len(parts) > accession_idx:
                    accession = parts[accession_idx]
                    accessions.add(accession)

        logging.info(f"Loaded {len(accessions)} accessions from metadata")
        return accessions

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        return set()
    except ValueError:
        logging.error("Could not find 'accession' column in metadata file")
        return set()
    except Exception as e:
        logging.error(f"Error loading metadata: {str(e)}")
        return set()

def validate_accession_format(accessions: List[str]) -> Tuple[List[str], int]:
    """
    Validate accession format.

    Args:
        accessions: List of accession IDs

    Returns:
        Tuple of (invalid_accessions, count_of_valid)
    """
    logging.info("Validating accession formats")
    valid_pattern = re.compile(r'^(GCA|GCF)_\d{9}\.\d+$')
    invalid_accessions = []
    valid_count = 0

    for accession in tqdm(accessions, desc="Validating accessions"):
        if not valid_pattern.match(accession):
            invalid_accessions.append(accession)
        else:
            valid_count += 1

    return invalid_accessions, valid_count

def validate_taxid_map(taxid_map: Dict[str, str], nodes_taxids: Set[str],
                      names_taxids: Set[str], metadata_accessions: Set[str]) -> None:
    """
    Validate taxid.map against nodes.dmp, names.dmp, and metadata.

    Args:
        taxid_map: Dictionary mapping accession to taxid
        nodes_taxids: Set of all taxids in nodes.dmp
        names_taxids: Set of all taxids in names.dmp
        metadata_accessions: Set of all accessions in metadata
    """
    logging.info("Validating taxid.map against taxonomy files")

    # Check accession format
    accessions = list(taxid_map.keys())
    invalid_accessions, valid_count = validate_accession_format(accessions)

    if invalid_accessions:
        logging.warning(f"Found {len(invalid_accessions)} accessions with invalid format")
        for i, acc in enumerate(invalid_accessions[:10]):
            logging.warning(f"  Invalid accession {i+1}: {acc}")
        if len(invalid_accessions) > 10:
            logging.warning(f"  ... and {len(invalid_accessions) - 10} more")

    logging.info(f"Accession format validation: {valid_count}/{len(accessions)} valid ({valid_count/len(accessions)*100:.2f}%)")

    # Check taxids against nodes.dmp
    missing_in_nodes = []
    for accession, taxid in tqdm(taxid_map.items(), desc="Checking against nodes.dmp"):
        if taxid not in nodes_taxids:
            missing_in_nodes.append((accession, taxid))

    if missing_in_nodes:
        logging.warning(f"Found {len(missing_in_nodes)} taxids missing from nodes.dmp")
        for i, (acc, tid) in enumerate(missing_in_nodes[:10]):
            logging.warning(f"  Missing taxid {i+1}: {acc} -> {tid}")
        if len(missing_in_nodes) > 10:
            logging.warning(f"  ... and {len(missing_in_nodes) - 10} more")

    logging.info(f"Nodes.dmp validation: {len(taxid_map) - len(missing_in_nodes)}/{len(taxid_map)} valid ({(len(taxid_map) - len(missing_in_nodes))/len(taxid_map)*100:.2f}%)")

    # Check taxids against names.dmp
    missing_in_names = []
    for accession, taxid in tqdm(taxid_map.items(), desc="Checking against names.dmp"):
        if taxid not in names_taxids:
            missing_in_names.append((accession, taxid))

    if missing_in_names:
        logging.warning(f"Found {len(missing_in_names)} taxids missing from names.dmp")
        for i, (acc, tid) in enumerate(missing_in_names[:10]):
            logging.warning(f"  Missing taxid {i+1}: {acc} -> {tid}")
        if len(missing_in_names) > 10:
            logging.warning(f"  ... and {len(missing_in_names) - 10} more")

    logging.info(f"Names.dmp validation: {len(taxid_map) - len(missing_in_names)}/{len(taxid_map)} valid ({(len(taxid_map) - len(missing_in_names))/len(taxid_map)*100:.2f}%)")

    # Check accessions against metadata
    if metadata_accessions:
        missing_in_metadata = []
        for accession in tqdm(taxid_map.keys(), desc="Checking against metadata"):
            if accession not in metadata_accessions:
                missing_in_metadata.append(accession)

        if missing_in_metadata:
            logging.warning(f"Found {len(missing_in_metadata)} accessions missing from metadata")
            for i, acc in enumerate(missing_in_metadata[:10]):
                logging.warning(f"  Missing accession {i+1}: {acc}")
            if len(missing_in_metadata) > 10:
                logging.warning(f"  ... and {len(missing_in_metadata) - 10} more")

        logging.info(f"Metadata validation: {len(taxid_map) - len(missing_in_metadata)}/{len(taxid_map)} valid ({(len(taxid_map) - len(missing_in_metadata))/len(taxid_map)*100:.2f}%)")

    # Summary
    logging.info("\nValidation Summary:")
    logging.info(f"Total entries in taxid.map: {len(taxid_map)}")
    logging.info(f"Valid accession formats: {valid_count} ({valid_count/len(accessions)*100:.2f}%)")
    logging.info(f"Taxids found in nodes.dmp: {len(taxid_map) - len(missing_in_nodes)} ({(len(taxid_map) - len(missing_in_nodes))/len(taxid_map)*100:.2f}%)")
    logging.info(f"Taxids found in names.dmp: {len(taxid_map) - len(missing_in_names)} ({(len(taxid_map) - len(missing_in_names))/len(taxid_map)*100:.2f}%)")

    if metadata_accessions:
        logging.info(f"Accessions found in metadata: {len(taxid_map) - len(missing_in_metadata)} ({(len(taxid_map) - len(missing_in_metadata))/len(taxid_map)*100:.2f}%)")

def main():
    """Main function to validate taxid.map."""
    args = parse_arguments()

    # Set up logging with the specified output directory
    log_file = Path(args.output_dir) / f"taxid_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(str(log_file))
        ]
    )

    logging.info(f"Starting validation with the following paths:")
    logging.info(f"  taxid.map: {args.taxid_map}")
    logging.info(f"  nodes.dmp: {args.nodes}")
    logging.info(f"  names.dmp: {args.names}")
    if args.metadata:
        logging.info(f"  metadata: {args.metadata}")
    logging.info(f"  log file: {log_file}")

    # Check if files exist before proceeding
    for file_path, file_name in [
        (args.taxid_map, "taxid.map"),
        (args.nodes, "nodes.dmp"),
        (args.names, "names.dmp")
    ]:
        if not Path(file_path).exists():
            logging.error(f"Error: {file_name} file not found at {file_path}")
            return 1

    if args.metadata and not Path(args.metadata).exists():
        logging.warning(f"Warning: Metadata file not found at {args.metadata}")
        logging.warning("Proceeding without metadata validation")
        args.metadata = None

    # Load files
    taxid_map = load_taxid_map(args.taxid_map)
    nodes_taxids = load_nodes(args.nodes)
    names_taxids = load_names(args.names)
    metadata_accessions = load_metadata_accessions(args.metadata) if args.metadata else set()

    # Validate
    validate_taxid_map(taxid_map, nodes_taxids, names_taxids, metadata_accessions)

    logging.info("Validation complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
