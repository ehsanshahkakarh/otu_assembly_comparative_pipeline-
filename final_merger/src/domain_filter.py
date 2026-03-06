#!/usr/bin/env python3
"""
Domain Filter Module
====================

Filters NCBI data by domain (Eukaryota, Bacteria, Archaea).
"""

import pandas as pd
import logging
from typing import List

logger = logging.getLogger(__name__)


def filter_by_domain(ncbi_df: pd.DataFrame, domains: List[str]) -> pd.DataFrame:
    """
    Filter NCBI species data to specific domains.

    Args:
        ncbi_df: NCBI species DataFrame with 'domain' column
        domains: List of domains to keep (e.g., ['Eukaryota'] or ['Bacteria', 'Archaea'])

    Returns:
        Filtered DataFrame containing only species from specified domains
    """
    logger.info(f"Filtering NCBI data to domains: {', '.join(domains)}")

    original_count = len(ncbi_df)
    original_genomes = ncbi_df['total_genome_count'].sum()

    # Filter NCBI data
    filtered_ncbi = ncbi_df[ncbi_df['domain'].isin(domains)].copy()

    filtered_count = len(filtered_ncbi)
    filtered_genomes = filtered_ncbi['total_genome_count'].sum()

    logger.info(f"  Filtered species: {filtered_count:,}/{original_count:,} ({filtered_count/original_count*100:.1f}%)")
    logger.info(f"  Filtered genomes: {filtered_genomes:,}/{original_genomes:,} ({filtered_genomes/original_genomes*100:.1f}%)")

    return filtered_ncbi

