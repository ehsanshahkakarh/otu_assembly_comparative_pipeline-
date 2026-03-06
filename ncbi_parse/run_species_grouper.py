#!/usr/bin/env python3
"""
Species Parser - Main Entry Point
==================================

Groups NCBI assembly data by species_taxid with lineage information.

Usage:
    python run_species_grouper.py                    # Process full dataset
    python run_species_grouper.py --sample 10000     # Process sample for testing
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import Config, setup_logging
from src.data_loader import load_assembly_data
from src.genome_classifier import classify_genome_source
from src.species_grouper import group_by_species
from src.lineage_enricher import add_lineage_information
# from src.taxid_fallback_resolver import resolve_missing_lineages  # Disabled - not needed with updated taxdump
from src.percentage_calculator import add_species_percentages
from src.output_writer import save_species_output


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Group NCBI assembly data by species_taxid with lineage information"
    )
    parser.add_argument(
        '--sample',
        type=int,
        help='Process only N records for testing'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory (default: nev_parse_meth/output/)'
    )
    args = parser.parse_args()

    # Initialize configuration
    config = Config()

    # Override output directory if specified
    if args.output_dir:
        config.paths.output_dir = Path(args.output_dir)

    # Setup logging
    logger = setup_logging(config.paths.logs_dir)

    logger.info("=" * 70)
    logger.info("SPECIES PARSER - NCBI Assembly Data by Species")
    logger.info("=" * 70)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Assembly file: {config.paths.assembly_file}")
    logger.info(f"Output directory: {config.paths.output_dir}")
    if args.sample:
        logger.info(f"Sample size: {args.sample:,} records")

    try:
        # Step 1: Load assembly data
        logger.info("\n[STEP 1/7] Loading assembly data...")
        df = load_assembly_data(config.paths.assembly_file, sample_size=args.sample)

        # Step 2: Classify genome sources (isolate vs uncultured)
        logger.info("\n[STEP 2/7] Classifying genome sources...")
        df = classify_genome_source(df, config.filters.uncultured_patterns)

        # Step 3: Group by species_taxid
        logger.info("\n[STEP 3/7] Grouping by species_taxid...")
        species_stats = group_by_species(df)

        # Step 4: Add lineage information
        logger.info("\n[STEP 4/5] Adding lineage information...")
        species_stats = add_lineage_information(species_stats)

        # Step 5: Calculate percentages
        logger.info("\n[STEP 5/5] Calculating species percentages...")
        species_stats = add_species_percentages(species_stats)

        # Step 6: Save output
        logger.info("\n[STEP 6/6] Saving output...")
        output_file = save_species_output(species_stats, config.paths.output_dir)

        # Summary
        logger.info("\n" + "=" * 70)
        logger.info("SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Input records: {len(df):,}")
        logger.info(f"Unique species: {len(species_stats):,}")
        logger.info(f"Total genomes: {species_stats['total_genome_count'].sum():,}")
        logger.info(f"Isolate genomes: {species_stats['isolate_genome_count'].sum():,}")
        logger.info(f"Uncultured genomes: {species_stats['uncultured_genome_count'].sum():,}")
        logger.info(f"\nOutput file: {output_file}")
        logger.info(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 70)

        return output_file

    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    output_file = main()
    print(f"\n📊 Species parsing complete!")
    print(f"📁 Output saved to: {output_file}")

