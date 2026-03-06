#!/usr/bin/env python3
"""
Pipeline module for Division Context Adder

Adds division context to unmapped entries in all three 18S census CSV outputs.

For entries that have NO parent lineage context (just their own name), this
adds the division from the raw census data as a minimal taxonomic context.

Example:
    Before: WIM80-lineage | WIM80-lineage | original_name | NA
    After:  WIM80-lineage | Evosea;WIM80-lineage | division;family | NA;NA
"""

import logging
from pathlib import Path

from .division_context_adder import load_division_mapping, process_csv_with_division_context


def run_division_context_adder():
    """Run the division context adder pipeline step."""
    
    logger = logging.getLogger(__name__)
    
    # Paths
    base_dir = Path(__file__).parent.parent.parent
    census_file = base_dir / 'metadata' / 'eukcensus_18S.clusters.97.tsv'
    csv_dir = base_dir / 'csv_outputs'
    
    # Load division mapping
    logger.info("=" * 80)
    logger.info("LOADING DIVISION MAPPING")
    logger.info("=" * 80)
    division_map = load_division_mapping(str(census_file))
    
    # CSV files to process
    csv_files = [
        'eukcensus_18S_by_family.csv',
        'eukcensus_18S_by_genus.csv',
        'eukcensus_18S_by_division.csv'
    ]
    
    logger.info("\n" + "=" * 80)
    logger.info("PROCESSING CSV FILES")
    logger.info("=" * 80)
    
    total_updated = 0
    
    for csv_file in csv_files:
        input_path = csv_dir / csv_file
        output_path = csv_dir / csv_file.replace('.csv', '_with_division_context.csv')
        
        if not input_path.exists():
            logger.warning(f"File not found: {input_path}")
            continue
        
        logger.info(f"\n--- Processing {csv_file} ---")
        
        updated = process_csv_with_division_context(
            str(input_path),
            str(output_path),
            division_map
        )
        
        total_updated += updated
        logger.info(f"✅ Saved to: {output_path.name}")
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total entries updated across all files: {total_updated}")
    logger.info(f"Output files saved with '_with_division_context' suffix")
    logger.info("\nDone! ✅")

