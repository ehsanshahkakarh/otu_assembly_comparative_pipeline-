#!/usr/bin/env python3
"""
Data Loader Module
==================

Loads NCBI assembly summary data.
"""

import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger('species_parser')


def load_assembly_data(assembly_file: Path, sample_size: int = None) -> pd.DataFrame:
    """
    Load NCBI assembly summary file.

    Args:
        assembly_file: Path to assembly_summary_genbank.txt
        sample_size: Optional number of records to load for testing

    Returns:
        DataFrame with assembly data
    """
    logger.info(f"Loading assembly data from: {assembly_file}")

    if not assembly_file.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

    # Read the file with tab separator
    # The file has comment lines starting with ##, and the header line starts with #
    # We need to skip lines starting with ## but keep the header line starting with #
    df = pd.read_csv(
        assembly_file,
        sep='\t',
        skiprows=lambda x: x == 0,  # Skip only the first comment line
        low_memory=False
    )

    # Clean column names (remove leading # if present)
    df.columns = df.columns.str.replace('^#', '', regex=True)

    # If sample size specified, take only first N records
    if sample_size:
        df = df.head(sample_size)
        logger.info(f"Loaded sample of {len(df):,} records")
    else:
        logger.info(f"Loaded {len(df):,} total records")

    # Verify required columns exist
    required_cols = ['assembly_accession', 'taxid', 'species_taxid', 'organism_name']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Unique species_taxid values: {df['species_taxid'].nunique():,}")

    return df

