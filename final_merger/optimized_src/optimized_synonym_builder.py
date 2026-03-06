#!/usr/bin/env python3
"""
Optimized Census Synonym Builder
=================================

Builds synonym dictionary for ALL census taxa at once.
This is the biggest optimization - loads synonym_dict.py only ONCE
instead of once per taxon.
"""

import pandas as pd
import logging
from typing import Dict, Set
from pathlib import Path
import sys

logger = logging.getLogger(__name__)

# Import the comprehensive synonym dictionary
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'ncbi_parse' / 'taxonomic_mapping'))
from synonym_dict import SYNONYM_TO_SCIENTIFIC


def build_all_census_synonyms(census_df: pd.DataFrame, level: str) -> Dict[str, Set[str]]:
    """
    Build synonym dictionary for ALL census taxa at once.
    
    This is much faster than building individually because:
    - Loads SYNONYM_TO_SCIENTIFIC only once (1.7M+ mappings)
    - Processes all census taxa in one pass
    
    Args:
        census_df: Census DataFrame with Name_to_use column
        level: Taxonomic level (division, family, genus)
    
    Returns:
        Dictionary mapping census_name -> set of all possible names
    """
    logger.info(f"Building optimized synonym dictionary for {len(census_df)} census {level} taxa...")
    
    # Get all census taxa names
    census_names = census_df['Name_to_use'].unique()
    
    # Build synonym dict for all taxa in one pass
    result = {}
    total_names = 0
    taxa_with_multiple = 0
    
    for census_name in census_names:
        # Start with the census name itself
        possible_names = {census_name}
        
        # Check if this name appears in the synonym dictionary
        if census_name in SYNONYM_TO_SCIENTIFIC:
            scientific_name = SYNONYM_TO_SCIENTIFIC[census_name]
            possible_names.add(scientific_name)
        
        # Find all synonyms that map to this scientific name
        # This is the expensive part, but we only do it once per unique census name
        for synonym, scientific in SYNONYM_TO_SCIENTIFIC.items():
            if scientific == census_name or synonym == census_name:
                possible_names.add(synonym)
                possible_names.add(scientific)
        
        result[census_name] = possible_names
        total_names += len(possible_names)
        
        if len(possible_names) > 1:
            taxa_with_multiple += 1
    
    # Log statistics
    avg_names = total_names / len(result) if len(result) > 0 else 0
    logger.info(f"  Created synonym dict for {len(result)} census taxa")
    logger.info(f"  Total possible names: {total_names}")
    logger.info(f"  Average names per taxon: {avg_names:.1f}")
    logger.info(f"  Taxa with multiple names: {taxa_with_multiple}")
    
    # Show examples of taxa with multiple names
    if taxa_with_multiple > 0:
        examples = [(name, names) for name, names in result.items() if len(names) > 1]
        for name, names in examples[:5]:
            logger.info(f"    {name}: {', '.join(sorted(names))}")
    
    return result


def get_synonyms_for_taxon(census_name: str, all_synonyms: Dict[str, Set[str]]) -> Set[str]:
    """
    Get synonyms for a specific census taxon from pre-built dictionary.
    
    Args:
        census_name: Census taxon name
        all_synonyms: Pre-built synonym dictionary
    
    Returns:
        Set of all possible names for this taxon
    """
    return all_synonyms.get(census_name, {census_name})

