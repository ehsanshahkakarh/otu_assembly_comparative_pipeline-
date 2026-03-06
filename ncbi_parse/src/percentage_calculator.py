#!/usr/bin/env python3
"""
Percentage Calculator Module
=============================

Calculates species-level genome percentages.
"""

import pandas as pd
import logging

logger = logging.getLogger('species_parser')


def add_species_percentages(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add species_genome_percentage column.
    
    This represents what percentage of the total dataset each species represents.
    
    Args:
        df: DataFrame with total_genome_count column
        
    Returns:
        DataFrame with added species_genome_percentage column
    """
    logger.info("Calculating species genome percentages...")
    
    # Calculate total genomes in dataset
    total_genomes = df['total_genome_count'].sum()
    
    # Calculate percentage for each species
    df['species_genome_percentage'] = (
        df['total_genome_count'] / total_genomes * 100
    ).round(2)
    
    logger.info(f"Total genomes in dataset: {total_genomes:,}")
    logger.info(f"Species genome percentages calculated")
    
    # Log top species by percentage
    logger.info(f"\nTop 5 species by genome percentage:")
    top_species = df.nlargest(5, 'species_genome_percentage')
    for _, row in top_species.iterrows():
        species_display = row.get('species_name', f"taxid:{row['species_taxid']}")[:50]
        logger.info(
            f"  {species_display}: "
            f"{row['species_genome_percentage']:.2f}% "
            f"({row['total_genome_count']:,} genomes)"
        )
    
    return df

