#!/usr/bin/env python3
"""
Species Grouper Module
=======================

Groups genomes by species_taxid and generates statistics.
"""

import pandas as pd
import logging

logger = logging.getLogger('species_parser')


def group_by_species(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group genomes by species_taxid and calculate statistics.

    Args:
        df: DataFrame with genome data including genome_source column

    Returns:
        DataFrame with one row per species_taxid containing:
        - species_taxid
        - total_genome_count
        - isolate_genome_count
        - uncultured_genome_count
        - isolate_genome_percentage
    """
    logger.info("Grouping genomes by species_taxid...")

    # Group by species_taxid and aggregate
    species_stats = df.groupby('species_taxid').agg({
        'assembly_accession': 'count',
        'genome_source': lambda x: (x == 'isolate').sum()
    }).rename(columns={
        'assembly_accession': 'total_genome_count',
        'genome_source': 'isolate_genome_count'
    })

    # Calculate uncultured count
    species_stats['uncultured_genome_count'] = (
        species_stats['total_genome_count'] - species_stats['isolate_genome_count']
    )

    # Calculate isolate percentage
    species_stats['isolate_genome_percentage'] = (
        species_stats['isolate_genome_count'] / species_stats['total_genome_count'] * 100
    ).round(2)

    # Reset index to make species_taxid a column
    species_stats.reset_index(inplace=True)

    # Log summary statistics
    logger.info(f"Species grouping results:")
    logger.info(f"  Total unique species: {len(species_stats):,}")
    logger.info(f"  Total genomes: {species_stats['total_genome_count'].sum():,}")
    logger.info(f"  Species with isolates: {(species_stats['isolate_genome_count'] > 0).sum():,}")
    logger.info(f"  Species with only uncultured: {(species_stats['isolate_genome_count'] == 0).sum():,}")
    logger.info(f"  Species with both: {((species_stats['isolate_genome_count'] > 0) & (species_stats['uncultured_genome_count'] > 0)).sum():,}")

    # Show top species by genome count
    logger.info(f"\nTop 10 species by genome count:")
    top_species = species_stats.nlargest(10, 'total_genome_count')
    for _, row in top_species.iterrows():
        logger.info(f"  Species taxid {row.name if hasattr(row, 'name') else row['species_taxid']}: {row['total_genome_count']:,} genomes ({row['isolate_genome_percentage']:.1f}% isolate)")

    return species_stats

