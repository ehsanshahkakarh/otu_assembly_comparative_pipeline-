#!/usr/bin/env python3
"""
Genome Classifier Module
========================

Classifies genomes as isolate or uncultured based on organism name
and RefSeq exclusion notes.

Adapted from ncbi_parse/py_ncbi/ncbi_parser/src/genome_classifier.py
"""

import pandas as pd
import logging
from typing import List

logger = logging.getLogger('species_parser')

# Default uncultured patterns
DEFAULT_UNCULTURED_PATTERNS = [
    'uncultured', 'environmental', 'metagenome', 'unclassified',
    'unknown', 'unidentified', 'mixed culture', 'enrichment culture',
    'derived from metagenome', 'metagenome-assembled', 'mag',
    'single amplified genome', 'sag', 'environmental sample'
]


def classify_genome_source(
    df: pd.DataFrame,
    uncultured_patterns: List[str] = None
) -> pd.DataFrame:
    """
    Classify genomes as isolate or uncultured.

    Logic:
    1. Default: All genomes are 'isolate'
    2. Check organism_name for uncultured patterns → 'uncultured'
    3. Check excluded_from_refseq for uncultured patterns → 'uncultured'

    Special case: enrichment culture with strain name remains 'isolate'

    Args:
        df: DataFrame with organism_name and optionally excluded_from_refseq columns
        uncultured_patterns: List of patterns to match (default: standard patterns)

    Returns:
        DataFrame with added genome_source column
    """
    logger.info("Classifying genome sources...")

    if uncultured_patterns is None:
        uncultured_patterns = DEFAULT_UNCULTURED_PATTERNS

    # Initialize all as isolate
    classification = pd.Series(['isolate'] * len(df), index=df.index)

    # Prepare organism names (lowercase for matching)
    organism_names = df['organism_name'].fillna('').str.lower()

    # Check organism names for uncultured patterns
    for pattern in uncultured_patterns:
        mask = organism_names.str.contains(pattern, na=False)

        # Special case: enrichment culture with strain name should remain as isolate
        if pattern == 'enrichment culture':
            has_strain = (
                organism_names.str.contains(r'\bstrain\b', na=False) |
                organism_names.str.contains(r'\bisolate\b', na=False) |
                organism_names.str.contains(r'\bstr\.\b', na=False)
            )
            mask = mask & (~has_strain)

        classification.loc[mask] = 'uncultured'

    # Check excluded_from_refseq column if available
    if 'excluded_from_refseq' in df.columns:
        excluded_notes = df['excluded_from_refseq'].fillna('').str.lower()
        for pattern in uncultured_patterns:
            mask = excluded_notes.str.contains(pattern, na=False)

            # Apply same enrichment culture exception
            if pattern == 'enrichment culture':
                has_strain = (
                    organism_names.str.contains(r'\bstrain\b', na=False) |
                    organism_names.str.contains(r'\bisolate\b', na=False) |
                    organism_names.str.contains(r'\bstr\.\b', na=False)
                )
                mask = mask & (~has_strain)

            classification.loc[mask] = 'uncultured'

    # Add classification to dataframe
    df = df.copy()
    df['genome_source'] = classification

    # Log results
    isolate_count = (classification == 'isolate').sum()
    uncultured_count = (classification == 'uncultured').sum()
    logger.info(
        f"Classification: {isolate_count:,} isolate ({isolate_count/len(df)*100:.1f}%), "
        f"{uncultured_count:,} uncultured ({uncultured_count/len(df)*100:.1f}%)"
    )

    return df

