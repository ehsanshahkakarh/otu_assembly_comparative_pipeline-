#!/usr/bin/env python3
"""
Optimized Species Searcher
===========================

Optimized species searching with:
- Pre-processed lineage data
- Vectorized pandas operations
- Compiled regex patterns
"""

import pandas as pd
import re
import logging
from typing import Set

logger = logging.getLogger(__name__)


def preprocess_species_data(species_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pre-process species data for faster searching.
    
    Optimizations:
    - Convert lineage_taxids to sets for O(1) lookup
    - Cache lineage strings
    
    Args:
        species_df: Species DataFrame
    
    Returns:
        DataFrame with pre-processed columns
    """
    logger.info("Pre-processing species data for faster searching...")
    
    # Create a copy to avoid modifying original
    df = species_df.copy()
    
    # Pre-process lineage_taxids into sets for faster matching
    logger.info("  Converting lineage_taxids to sets...")
    df['lineage_taxids_set'] = df['lineage_taxids'].apply(
        lambda x: set(str(x).split(';')) if pd.notna(x) and x != '' else set()
    )
    
    logger.info("  Pre-processing complete!")
    return df


def search_species_optimized(
    species_df: pd.DataFrame,
    census_name: str,
    census_taxid: str,
    possible_names: Set[str],
    level: str
) -> pd.DataFrame:
    """
    Optimized species search using pre-processed data and vectorized operations.
    
    Args:
        species_df: Pre-processed species DataFrame (with lineage_taxids_set)
        census_name: Census taxon name
        census_taxid: Census taxon ID (clean string)
        possible_names: Set of all possible names for this taxon
        level: Taxonomic level
    
    Returns:
        DataFrame of matched species
    """
    # Initialize match containers
    taxid_matched = pd.DataFrame()
    name_matched = pd.DataFrame()
    
    # PRIORITY 1: Direct taxid matching using pre-processed sets
    if census_taxid and census_taxid != '':
        # Use vectorized set operation
        taxid_mask = species_df['lineage_taxids_set'].apply(lambda x: census_taxid in x)
        taxid_matched = species_df[taxid_mask]
        
        if len(taxid_matched) > 0:
            logger.debug(f"  Taxid match: Found {len(taxid_matched)} species")
    
    # PRIORITY 2: Name matching in lineage (vectorized)
    if len(possible_names) > 0:
        # Build combined regex pattern for all names at once
        escaped_names = [re.escape(name) for name in possible_names]
        # Pattern: match name as whole word in lineage (bounded by ; or start/end)
        combined_pattern = '|'.join([f'(^{n};|;{n};|;{n}$|^{n}$)' for n in escaped_names])
        
        # Vectorized string matching
        name_mask = species_df['lineage'].str.contains(combined_pattern, regex=True, na=False)
        name_matched = species_df[name_mask]
        
        if len(name_matched) > 0:
            logger.debug(f"  Name match: Found {len(name_matched)} species")
    
    # Combine results (remove duplicates based on species_name)
    all_matched = pd.concat([taxid_matched, name_matched]).drop_duplicates(subset=['species_name'])

    return all_matched

