#!/usr/bin/env python3
"""
Metrics Calculator Module
==========================

Calculates novelty, overrepresentation, and coverage metrics.
"""

import pandas as pd
import logging

logger = logging.getLogger(__name__)


def calculate_metrics(matched_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate novelty factor and overrepresentation factor.

    Args:
        matched_df: DataFrame with matched census-NCBI results

    Returns:
        DataFrame with added metric columns
    """
    logger.info("Calculating metrics...")

    # Novelty factor: census_otu_count / ncbi_species_count
    # Higher values indicate more environmental diversity than genomic representation
    matched_df['novelty_factor'] = matched_df.apply(
        lambda row: round(row['census_otu_count'] / row['ncbi_species_count'], 3)
        if row['ncbi_species_count'] > 0 else float('inf'),
        axis=1
    )

    # Overrepresentation factor: ncbi_species_count / census_otu_count
    # Higher values indicate database bias toward cultured taxa
    matched_df['overrepresentation_factor'] = matched_df.apply(
        lambda row: round(row['ncbi_species_count'] / row['census_otu_count'], 3)
        if row['census_otu_count'] > 0 else float('inf'),
        axis=1
    )
    
    # Log summary statistics
    matched_only = matched_df[matched_df['match_status'] == 'matched']
    
    if len(matched_only) > 0:
        avg_novelty = matched_only[matched_only['novelty_factor'] != float('inf')]['novelty_factor'].mean()
        avg_overrep = matched_only[matched_only['overrepresentation_factor'] != float('inf')]['overrepresentation_factor'].mean()
        
        logger.info(f"  Average novelty factor: {avg_novelty:.2f}")
        logger.info(f"  Average overrepresentation factor: {avg_overrep:.2f}")
        
        # Count high novelty taxa (novelty > 10)
        high_novelty = (matched_only['novelty_factor'] > 10).sum()
        logger.info(f"  High novelty taxa (>10): {high_novelty}")
        
        # Count well-represented taxa (novelty < 2)
        well_represented = (matched_only['novelty_factor'] < 2).sum()
        logger.info(f"  Well-represented taxa (<2): {well_represented}")
    
    return matched_df

