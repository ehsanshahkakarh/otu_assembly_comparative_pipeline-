#!/usr/bin/env python3
"""
Isolate Tracker Module
======================

Adds isolate genome information to matched results.
"""

import pandas as pd
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def add_isolate_information(
    matched_df: pd.DataFrame,
    isolate_df: Optional[pd.DataFrame],
    level: str,
    census_level: str
) -> pd.DataFrame:
    """
    Add isolate genome counts and percentages to matched results.
    
    Args:
        matched_df: DataFrame with matched census-NCBI results
        isolate_df: DataFrame with isolate counts per taxon (or None)
        level: NCBI taxonomic level (phylum, family, genus)
        census_level: Census level name (division, family, genus)
    
    Returns:
        DataFrame with added isolate_count and isolate_percentage columns
    """
    logger.info("Adding isolate information...")
    
    # Initialize isolate columns
    matched_df['isolate_count'] = 0
    matched_df['isolate_percentage'] = 0.0
    
    if isolate_df is None or isolate_df.empty:
        logger.info("  No isolate data available, all isolate counts set to 0")
        return matched_df
    
    # For each matched taxon, look up isolate counts
    isolate_added = 0
    
    for idx, row in matched_df.iterrows():
        taxon_name = row[census_level]
        total_genomes = row['ncbi_genome_count']
        
        if total_genomes > 0:
            # Look up isolate count for this taxon
            isolate_match = isolate_df[isolate_df['taxon'] == taxon_name]
            
            if not isolate_match.empty:
                isolate_count = isolate_match['isolate_count'].iloc[0]
                isolate_pct = (isolate_count / total_genomes * 100) if total_genomes > 0 else 0
                
                matched_df.at[idx, 'isolate_count'] = int(isolate_count)
                matched_df.at[idx, 'isolate_percentage'] = round(isolate_pct, 2)
                isolate_added += 1
    
    total_isolates = matched_df['isolate_count'].sum()
    logger.info(f"  Added isolate data for {isolate_added}/{len(matched_df)} taxa")
    logger.info(f"  Total isolates: {total_isolates:,}")
    
    return matched_df

