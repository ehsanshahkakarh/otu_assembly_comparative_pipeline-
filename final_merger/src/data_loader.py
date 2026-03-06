#!/usr/bin/env python3
"""
Data Loader Module
==================

Loads census and NCBI data files.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def load_census_data(census_file: Path) -> pd.DataFrame:
    """
    Load census data file (18S or 16S).
    
    Args:
        census_file: Path to census CSV file
    
    Returns:
        DataFrame with census data
    
    Raises:
        FileNotFoundError: If census file doesn't exist
        ValueError: If required columns are missing
    """
    logger.info(f"Loading census data from: {census_file.name}")
    
    if not census_file.exists():
        raise FileNotFoundError(f"Census file not found: {census_file}")

    # Read CSV with taxid as string to preserve original format
    df = pd.read_csv(census_file, dtype={'taxid': str})

    # Verify required columns
    required_cols = ['Name_to_use', 'otu_count', 'size_count']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in census file: {missing_cols}")
    
    logger.info(f"  Loaded {len(df):,} census entries")
    logger.info(f"  Total OTUs: {df['otu_count'].sum():,}")
    logger.info(f"  Total size: {df['size_count'].sum():,}")
    
    return df


def load_ncbi_species_data(species_file: Path) -> pd.DataFrame:
    """
    Load NCBI species-level data from species_grouped_*.csv file.

    This loads the raw species-level data from nev_parse_meth output.
    Aggregation by taxonomic level happens during matching.

    Args:
        species_file: Path to species_grouped_*.csv file

    Returns:
        DataFrame with species-level NCBI data containing:
        - species_taxid, species_name
        - total_genome_count, isolate_genome_count, uncultured_genome_count
        - lineage, lineage_ranks, lineage_taxids

    Raises:
        FileNotFoundError: If species file doesn't exist
        ValueError: If required columns are missing
    """
    logger.info(f"Loading NCBI species data from: {species_file.name}")

    if not species_file.exists():
        raise FileNotFoundError(f"NCBI species file not found: {species_file}")

    # Load with low_memory=False to avoid dtype warnings on large file
    df = pd.read_csv(species_file, low_memory=False)

    # Verify required columns
    required_cols = ['species_taxid', 'species_name', 'total_genome_count',
                     'isolate_genome_count', 'uncultured_genome_count',
                     'lineage', 'lineage_ranks', 'lineage_taxids']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in species file: {missing_cols}")

    # Extract domain from lineage
    # Domain is the second element in lineage (after "cellular organisms")
    def extract_domain(lineage_str):
        if pd.isna(lineage_str):
            return 'Unknown'
        parts = str(lineage_str).split(';')
        if len(parts) >= 2:
            return parts[1]  # Second element is domain (Bacteria, Archaea, Eukaryota)
        return 'Unknown'

    df['domain'] = df['lineage'].apply(extract_domain)

    logger.info(f"  Loaded {len(df):,} species")
    logger.info(f"  Total genomes: {df['total_genome_count'].sum():,}")
    logger.info(f"  Isolate genomes: {df['isolate_genome_count'].sum():,}")
    logger.info(f"  Uncultured genomes: {df['uncultured_genome_count'].sum():,}")

    # Show domain breakdown
    domain_counts = df.groupby('domain')['total_genome_count'].sum()
    for domain, count in domain_counts.items():
        logger.info(f"    {domain}: {count:,} genomes")

    return df


def load_isolate_data(isolate_file: Path, level: str) -> Optional[pd.DataFrame]:
    """
    Load NCBI isolate data file (optional).
    
    Args:
        isolate_file: Path to NCBI accessions CSV file
        level: Taxonomic level (phylum, family, genus)
    
    Returns:
        DataFrame with isolate counts, or None if file doesn't exist
    """
    if not isolate_file.exists():
        logger.warning(f"  Isolate file not found: {isolate_file.name} (skipping isolate tracking)")
        return None
    
    logger.info(f"Loading isolate data from: {isolate_file.name}")
    
    try:
        df = pd.read_csv(isolate_file)
        
        # Verify required columns
        if 'genome_source' not in df.columns or level not in df.columns:
            logger.warning(f"  Missing required columns in isolate file (skipping isolate tracking)")
            return None
        
        # Count isolates per taxon
        isolate_counts = df[df['genome_source'] == 'isolate'].groupby(level).size()
        total_counts = df.groupby(level).size()
        
        isolate_df = pd.DataFrame({
            'taxon': isolate_counts.index,
            'isolate_count': isolate_counts.values,
            'total_genomes': total_counts.loc[isolate_counts.index].values
        })
        
        isolate_df['isolate_percentage'] = (
            isolate_df['isolate_count'] / isolate_df['total_genomes'] * 100
        ).round(2)
        
        logger.info(f"  Loaded isolate data for {len(isolate_df):,} taxa")
        logger.info(f"  Total isolates: {isolate_df['isolate_count'].sum():,}")
        
        return isolate_df
        
    except Exception as e:
        logger.warning(f"  Error loading isolate data: {e} (skipping isolate tracking)")
        return None

