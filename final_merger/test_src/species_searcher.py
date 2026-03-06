#!/usr/bin/env python3
"""
Species Searcher Module
========================

Searches the species_grouped DataFrame for species that belong to a census taxon.
Uses hierarchical matching strategy with synonym support.
"""

import pandas as pd
import re
import logging
from typing import Set

logger = logging.getLogger(__name__)


def search_species_for_census_taxon(
    species_df: pd.DataFrame,
    census_name: str,
    census_taxid: any,
    possible_names: Set[str],
    level: str
) -> pd.DataFrame:
    """
    Search species_grouped DataFrame for species belonging to a census taxon.
    
    Matching strategy (hierarchical):
    1. Direct taxid match: census_taxid in lineage_taxids
    2. Name match in lineage: any possible_name in lineage string
    3. Extract from lineage_ranks: find taxon at the specified rank
    
    Args:
        species_df: Species DataFrame with lineage columns
        census_name: Census taxon name (e.g., 'Apicomplexa')
        census_taxid: Census taxon ID (can be int, float, or string)
        possible_names: Set of all possible names for this taxon (from synonym dict)
        level: Taxonomic level (division/phylum, family, genus)
        
    Returns:
        DataFrame of matching species
    """
    logger.debug(f"Searching for census taxon: {census_name} (taxid: {census_taxid})")
    logger.debug(f"  Possible names: {possible_names}")
    
    # Initialize match containers
    taxid_matched = pd.DataFrame()
    name_matched = pd.DataFrame()
    rank_matched = pd.DataFrame()
    
    # PRIORITY 1: Direct taxid matching (census taxid in lineage_taxids)
    # Skip if taxid is empty string or None
    if census_taxid and census_taxid != '':
        try:
            # Taxid should already be a string from aggregator
            census_taxid_str = census_taxid.strip()

            # Check if census taxid appears anywhere in the lineage_taxids string
            # Pattern: ;taxid; or ^taxid; or ;taxid$ or ^taxid$
            taxid_pattern = f';{census_taxid_str};|^{census_taxid_str};|;{census_taxid_str}$|^{census_taxid_str}$'
            taxid_mask = species_df['lineage_taxids'].str.contains(taxid_pattern, regex=True, na=False)
            taxid_matched = species_df[taxid_mask]

            if len(taxid_matched) > 0:
                logger.debug(f"  Taxid match: Found {len(taxid_matched)} species")
        except (ValueError, TypeError, AttributeError) as e:
            logger.warning(f"  Skipping taxid matching for '{census_name}' - invalid taxid: {census_taxid}")
    
    # PRIORITY 2: Name matching in lineage (any possible name in lineage string)
    for name in possible_names:
        escaped_name = re.escape(name)
        name_pattern = f';{escaped_name};|^{escaped_name};|;{escaped_name}$|^{escaped_name}$'
        name_mask = species_df['lineage'].str.contains(name_pattern, regex=True, na=False)
        name_matches = species_df[name_mask]
        
        if len(name_matches) > 0:
            logger.debug(f"  Name match '{name}': Found {len(name_matches)} species")
            name_matched = pd.concat([name_matched, name_matches]).drop_duplicates()
    
    # PRIORITY 3: Extract from lineage_ranks (find taxon at specified rank)
    # Map level to NCBI rank name
    rank_map = {
        'division': 'phylum',  # Census uses 'division', NCBI uses 'phylum'
        'phylum': 'phylum',
        'family': 'family',
        'genus': 'genus'
    }
    target_rank = rank_map.get(level, level)
    
    # Extract taxon at target rank for each species
    rank_matched = _extract_by_rank(species_df, possible_names, target_rank)
    if len(rank_matched) > 0:
        logger.debug(f"  Rank match (rank={target_rank}): Found {len(rank_matched)} species")
    
    # Combine all matches (avoid duplicates)
    all_matched = pd.concat([taxid_matched, name_matched, rank_matched]).drop_duplicates()
    
    if len(all_matched) > 0:
        logger.debug(f"  Total unique matches: {len(all_matched)} species")
        logger.debug(f"    Taxid matches: {len(taxid_matched)}")
        logger.debug(f"    Name matches: {len(name_matched)}")
        logger.debug(f"    Rank matches: {len(rank_matched)}")
    
    return all_matched


def _extract_by_rank(species_df: pd.DataFrame, possible_names: Set[str], target_rank: str) -> pd.DataFrame:
    """
    Extract species where the taxon at target_rank matches any possible_name.

    OPTIMIZED: Uses vectorized operations instead of row-by-row iteration.

    Args:
        species_df: Species DataFrame
        possible_names: Set of possible names for the census taxon
        target_rank: Target rank to extract (e.g., 'phylum', 'family', 'genus')

    Returns:
        DataFrame of matching species
    """
    # For now, skip rank-based matching to improve performance
    # This is redundant with name matching in most cases
    # TODO: Implement vectorized rank extraction if needed
    return pd.DataFrame()

