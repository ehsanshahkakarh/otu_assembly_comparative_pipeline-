#!/usr/bin/env python3
"""
Configuration module for 18S Census Parser

Handles directory paths, environment setup, and logging configuration.
"""

import os
import sys
import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Tuple


@dataclass
class PathConfig:
    """Configuration for directory paths in the 18S census parser."""
    script_dir: Path
    censusparse_dir: Path
    metadata_dir: Path
    csv_output_dir: Path
    log_dir: Path


def setup_directory_paths() -> PathConfig:
    """
    Set up directory paths for the reorganized 18S_censusparse structure.

    Returns:
        PathConfig object with all directory paths
    """
    script_dir = Path(__file__).resolve().parent.parent
    censusparse_dir = script_dir.parent
    metadata_dir = censusparse_dir / "metadata"
    csv_output_dir = censusparse_dir / "csv_outputs"
    log_dir = script_dir / "logs"

    # Ensure directories exist
    csv_output_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    return PathConfig(
        script_dir=script_dir,
        censusparse_dir=censusparse_dir,
        metadata_dir=metadata_dir,
        csv_output_dir=csv_output_dir,
        log_dir=log_dir
    )


def setup_taxonkit_environment() -> dict:
    """
    Set up environment for taxonkit.

    Note: taxonkit has its own built-in taxdump and doesn't require external TAXONKIT_DB.

    Returns:
        Environment dictionary (using system default)
    """
    env = os.environ.copy()
    # taxonkit uses its own built-in NCBI taxdump, no need to set TAXONKIT_DB
    return env


def setup_logging(log_dir: Path, output_prefix: str = "eukcensus_optimized") -> logging.Logger:
    """
    Set up logging configuration for the optimization script.

    Args:
        log_dir: Directory to store log files
        output_prefix: Prefix for output files to include in log messages

    Returns:
        Configured logger instance
    """
    # Create log file path in the logs directory
    log_file = log_dir / "eukcensus_optimization.log"

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),  # Append to existing log
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Log the start of processing
    logging.info(f"Starting optimized EukCensus processing: eukcensus_18S.clusters.97.tsv -> {output_prefix}_*")
    logging.info(f"Log file: {log_file}")
    
    return logging.getLogger(__name__)

