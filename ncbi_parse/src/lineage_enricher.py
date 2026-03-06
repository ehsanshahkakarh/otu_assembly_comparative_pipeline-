#!/usr/bin/env python3
"""
Lineage Enricher Module
========================

Adds taxonomic lineage information using taxonkit.
"""

import pandas as pd
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Dict, Tuple

logger = logging.getLogger('species_parser')


def add_lineage_information(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add lineage information for each species_taxid using taxonkit.

    Args:
        df: DataFrame with species_taxid column

    Returns:
        DataFrame with added columns:
        - species_name (extracted from lineage)
        - lineage
        - lineage_ranks
        - lineage_taxids
    """
    logger.info("Adding lineage information using taxonkit...")

    # Get unique species_taxid values
    unique_taxids = df['species_taxid'].dropna().unique().astype(int).tolist()
    logger.info(f"Retrieving lineages for {len(unique_taxids):,} unique species...")

    # Get lineages using taxonkit
    lineage_data = _get_lineages_from_taxids(unique_taxids)

    # Map lineage data to dataframe
    df['lineage'] = df['species_taxid'].map(
        lambda x: lineage_data.get(str(int(x)), ('', '', ''))[0] if pd.notna(x) else ''
    )
    df['lineage_ranks'] = df['species_taxid'].map(
        lambda x: lineage_data.get(str(int(x)), ('', '', ''))[1] if pd.notna(x) else ''
    )
    df['lineage_taxids'] = df['species_taxid'].map(
        lambda x: lineage_data.get(str(int(x)), ('', '', ''))[2] if pd.notna(x) else ''
    )

    # Extract species name from lineage
    df['species_name'] = df.apply(_extract_species_name, axis=1)

    # Log results
    lineages_found = (df['lineage'] != '').sum()
    species_names_found = (df['species_name'] != '').sum()
    logger.info(f"Lineages retrieved: {lineages_found:,}/{len(df):,} ({lineages_found/len(df)*100:.1f}%)")
    logger.info(f"Species names extracted: {species_names_found:,}/{len(df):,} ({species_names_found/len(df)*100:.1f}%)")

    return df


def _extract_species_name(row: pd.Series) -> str:
    """
    Extract species name from lineage based on ranks.

    Args:
        row: DataFrame row with lineage, lineage_ranks columns

    Returns:
        Species name or empty string if not found
    """
    lineage = row.get('lineage', '')
    ranks = row.get('lineage_ranks', '')

    if not lineage or not ranks:
        return ''

    lineage_parts = lineage.split(';')
    rank_parts = ranks.split(';')

    # Find the species rank
    for lineage_part, rank_part in zip(lineage_parts, rank_parts):
        if rank_part.strip().lower() == 'species':
            return lineage_part.strip()

    return ''


def _get_lineages_from_taxids(taxids: list) -> Dict[str, Tuple[str, str, str]]:
    """
    Get lineages for taxids using taxonkit lineage -R -t.
    
    Args:
        taxids: List of taxid integers
        
    Returns:
        Dictionary mapping taxid to (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}
    
    lineage_data = {}
    
    # Create temporary file with taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")
    
    try:
        # Run taxonkit lineage with -R (show rank) and -t (show taxid lineage)
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode == 0:
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        taxid = parts[0].strip()
                        lineage = parts[1].strip()
                        lineage_taxids = parts[2].strip()
                        lineage_ranks = parts[3].strip()
                        
                        if lineage and lineage != taxid:
                            lineage_data[taxid] = (lineage, lineage_ranks, lineage_taxids)
        else:
            logger.warning(f"Taxonkit warning: {result.stderr}")
            
    except Exception as e:
        logger.error(f"Error running taxonkit lineage: {e}")
    finally:
        Path(temp_filename).unlink(missing_ok=True)
    
    logger.info(f"Retrieved {len(lineage_data):,} lineages from taxonkit")
    return lineage_data

