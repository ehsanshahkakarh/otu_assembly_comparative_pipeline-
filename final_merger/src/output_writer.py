#!/usr/bin/env python3
"""
Output Writer Module
====================

Saves merged results to CSV files.
"""

import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def save_merged_output(
    merged_df: pd.DataFrame,
    output_file: Path,
    census_level: str
) -> Path:
    """
    Save merged results to CSV file with proper sorting.
    
    Args:
        merged_df: DataFrame with merged census-NCBI results
        output_file: Path to output CSV file
        census_level: Census level name (for logging)
    
    Returns:
        Path to saved output file
    """
    logger.info(f"Saving merged data to: {output_file.name}")
    
    # Sort: matched first (by genome count desc), then unmatched (by otu count desc)
    matched_df = merged_df[merged_df['match_status'] == 'matched'].sort_values(
        'ncbi_genome_count', ascending=False
    )
    unmatched_df = merged_df[merged_df['match_status'] != 'matched'].sort_values(
        'census_otu_count', ascending=False
    )
    sorted_df = pd.concat([matched_df, unmatched_df], ignore_index=True)
    
    # Save to CSV
    sorted_df.to_csv(output_file, index=False)
    
    # Log summary
    matched_count = len(matched_df)
    unmatched_count = len(unmatched_df)
    match_rate = (matched_count / len(sorted_df) * 100) if len(sorted_df) > 0 else 0
    
    logger.info(f"  Total entries: {len(sorted_df):,}")
    logger.info(f"  Matched: {matched_count:,} ({match_rate:.1f}%)")
    logger.info(f"  Census-only: {unmatched_count:,}")
    logger.info(f"  Output saved: {output_file}")
    
    return output_file

