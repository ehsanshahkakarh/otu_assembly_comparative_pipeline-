#!/usr/bin/env python3
"""
Census Synonym Builder Module
==============================

Builds a comprehensive synonym dictionary for each census taxon using NCBI names.dmp.
This allows matching census taxa even when NCBI uses different names or synonyms.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Dict, Set
import sys

logger = logging.getLogger(__name__)

# Add taxonomic_mapping to path for synonym_dict import
TAXONOMIC_MAPPING_DIR = Path(__file__).resolve().parents[2] / 'ncbi_parse' / 'taxonomic_mapping'
sys.path.insert(0, str(TAXONOMIC_MAPPING_DIR))

try:
    from synonym_dict import SYNONYM_TO_SCIENTIFIC
    logger.info(f"✅ Loaded synonym dictionary with {len(SYNONYM_TO_SCIENTIFIC):,} mappings")
except ImportError as e:
    logger.warning(f"⚠️  Could not load synonym_dict: {e}")
    SYNONYM_TO_SCIENTIFIC = {}


def build_census_synonym_dict(census_df: pd.DataFrame, level: str) -> Dict[str, Set[str]]:
    """
    Build a dictionary mapping each census taxon to all its possible names.
    
    For each census taxon, we create a set of all possible names that could
    match it in the species_grouped file:
    - The census name itself
    - Scientific name (if different from census name)
    - All synonyms from NCBI names.dmp
    
    Args:
        census_df: Census DataFrame with Name_to_use column
        level: Taxonomic level (division, family, genus)
        
    Returns:
        Dictionary: {census_name: {set of all possible matching names}}
    
    Example:
        {
            'Apicomplexa': {'Apicomplexa', 'Sporozoa', ...},
            'Euryarchaeota': {'Euryarchaeota', 'Methanobacteriota', ...}
        }
    """
    logger.info(f"Building synonym dictionary for {len(census_df)} census {level} taxa...")
    
    census_synonym_dict = {}
    
    for _, row in census_df.iterrows():
        census_name = row['Name_to_use']
        
        # Start with the census name itself
        possible_names = {census_name}
        
        # Check if this name has a scientific name mapping
        if census_name in SYNONYM_TO_SCIENTIFIC:
            scientific_name = SYNONYM_TO_SCIENTIFIC[census_name]
            possible_names.add(scientific_name)
            logger.debug(f"  {census_name} → scientific: {scientific_name}")
        
        # Also check reverse: if census_name is a scientific name,
        # find all synonyms that map to it
        for synonym, scientific in SYNONYM_TO_SCIENTIFIC.items():
            if scientific == census_name:
                possible_names.add(synonym)
        
        census_synonym_dict[census_name] = possible_names
        
        if len(possible_names) > 1:
            logger.debug(f"  {census_name}: {len(possible_names)} possible names")
    
    # Log statistics
    total_names = sum(len(names) for names in census_synonym_dict.values())
    avg_names = total_names / len(census_synonym_dict) if census_synonym_dict else 0
    
    logger.info(f"  Created synonym dict for {len(census_synonym_dict)} census taxa")
    logger.info(f"  Total possible names: {total_names}")
    logger.info(f"  Average names per taxon: {avg_names:.1f}")
    
    # Show examples of taxa with multiple names
    multi_name_taxa = [(name, names) for name, names in census_synonym_dict.items() if len(names) > 1]
    if multi_name_taxa:
        logger.info(f"  Taxa with multiple names: {len(multi_name_taxa)}")
        for name, names in multi_name_taxa[:5]:
            logger.info(f"    {name}: {', '.join(sorted(names))}")
    
    return census_synonym_dict

