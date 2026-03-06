#!/usr/bin/env python3
"""
Output Writer Module
====================

Saves species-grouped data to CSV file.
"""

import pandas as pd
import logging
from pathlib import Path
from datetime import datetime

logger = logging.getLogger('species_parser')


def save_species_output(df: pd.DataFrame, output_dir: Path) -> Path:
    """
    Save species-grouped data to CSV file.
    
    Args:
        df: DataFrame with species statistics and lineage information
        output_dir: Directory to save output file
        
    Returns:
        Path to saved output file
    """
    logger.info("Saving species-grouped output...")
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define column order
    column_order = [
        'species_taxid',
        'species_name',
        'total_genome_count',
        'isolate_genome_count',
        'uncultured_genome_count',
        'isolate_genome_percentage',
        'species_genome_percentage',
        'lineage',
        'lineage_ranks',
        'lineage_taxids'
    ]
    
    # Reorder columns (only include columns that exist)
    existing_columns = [col for col in column_order if col in df.columns]
    df_output = df[existing_columns].copy()
    
    # Sort by total_genome_count descending
    df_output = df_output.sort_values('total_genome_count', ascending=False)
    
    # Generate output filename with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = output_dir / f'species_grouped_{timestamp}.csv'
    
    # Save to CSV
    df_output.to_csv(output_file, index=False)
    
    logger.info(f"✅ Saved species output: {output_file.name}")
    logger.info(f"   Total species: {len(df_output):,}")
    logger.info(f"   Total genomes: {df_output['total_genome_count'].sum():,}")
    logger.info(f"   File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    return output_file

