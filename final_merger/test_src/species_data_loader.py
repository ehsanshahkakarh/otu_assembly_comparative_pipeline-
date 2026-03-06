#!/usr/bin/env python3
"""
Species Data Loader Module
===========================

Loads the species_grouped_*.csv file from nev_parse_meth.
"""

import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_species_grouped_data(species_file: Path) -> pd.DataFrame:
    """
    Load species_grouped_*.csv file.
    
    Args:
        species_file: Path to species_grouped_*.csv file
        
    Returns:
        DataFrame with columns:
        - species_taxid
        - species_name
        - total_genome_count
        - isolate_genome_count
        - uncultured_genome_count
        - isolate_genome_percentage
        - species_genome_percentage
        - lineage (semicolon-separated)
        - lineage_ranks (semicolon-separated)
        - lineage_taxids (semicolon-separated)
    
    Raises:
        FileNotFoundError: If species file doesn't exist
        ValueError: If required columns are missing
    """
    logger.info(f"Loading species data from: {species_file.name}")
    
    if not species_file.exists():
        raise FileNotFoundError(f"Species file not found: {species_file}")
    
    df = pd.read_csv(species_file)
    
    # Verify required columns
    required_cols = [
        'species_taxid', 'species_name', 'total_genome_count',
        'isolate_genome_count', 'uncultured_genome_count',
        'lineage', 'lineage_ranks', 'lineage_taxids'
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in species file: {missing_cols}")
    
    # Fill NaN values in string columns
    df['lineage'] = df['lineage'].fillna('')
    df['lineage_ranks'] = df['lineage_ranks'].fillna('')
    df['lineage_taxids'] = df['lineage_taxids'].fillna('')
    df['species_name'] = df['species_name'].fillna('')
    
    logger.info(f"  Loaded {len(df):,} species")
    logger.info(f"  Total genomes: {df['total_genome_count'].sum():,}")
    logger.info(f"  Isolate genomes: {df['isolate_genome_count'].sum():,}")
    logger.info(f"  Uncultured genomes: {df['uncultured_genome_count'].sum():,}")
    
    return df

