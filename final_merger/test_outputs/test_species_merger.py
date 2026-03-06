#!/usr/bin/env python3
"""
Test Species-Based Merger
==========================

Test implementation of merger using species_grouped_*.csv as data source.

This script demonstrates the correct workflow:
1. Load census data (census taxa are the BACKBONE)
2. Build synonym dictionary for census taxa using names.dmp
3. Load species_grouped_*.csv file
4. For each census taxon, search species_grouped for matching species
5. Aggregate genome counts from matched species
6. Save merged results

Usage:
    python test_species_merger.py --level division
    python test_species_merger.py --level family
    python test_species_merger.py --level genus
"""

import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from config import Config
from src.data_loader import load_census_data
from test_src.species_data_loader import load_species_grouped_data
from test_src.census_synonym_builder import build_census_synonym_dict
from test_src.species_searcher import search_species_for_census_taxon
from test_src.species_aggregator import aggregate_all_census_taxa


def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Test species-based merger using species_grouped_*.csv"
    )
    parser.add_argument(
        '--level',
        choices=['division', 'family', 'genus'],
        default='division',
        help='Taxonomic level to process'
    )
    parser.add_argument(
        '--census-type',
        choices=['18S', '16S'],
        default='18S',
        help='Census type (18S or 16S)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info("=" * 80)
    logger.info("TEST SPECIES-BASED MERGER")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Census type: {args.census_type}")
    logger.info(f"Level: {args.level}")
    
    try:
        # Initialize configuration
        config = Config()
        
        # Step 1: Load census data
        logger.info("\n[STEP 1/5] Loading census data...")
        census_file = config.paths.get_census_file(args.census_type, 
                                                    'phylum' if args.level == 'division' else args.level)
        census_df = load_census_data(census_file)
        
        # Step 2: Build synonym dictionary for census taxa
        logger.info("\n[STEP 2/5] Building synonym dictionary for census taxa...")
        census_synonym_dict = build_census_synonym_dict(census_df, args.level)
        
        # Step 3: Load species_grouped data
        logger.info("\n[STEP 3/5] Loading species_grouped data...")
        species_file = config.paths.parse_dir / 'ncbi_parse' / 'metadata' / 'nev_parse_meth' / 'output' / 'species_grouped_20260223_014859.csv'
        species_df = load_species_grouped_data(species_file)
        
        # Step 4: Search and aggregate for each census taxon
        logger.info("\n[STEP 4/5] Searching species for each census taxon...")
        merged_df = aggregate_all_census_taxa(
            census_df,
            species_df,
            census_synonym_dict,
            args.level,
            search_species_for_census_taxon
        )
        
        # Step 5: Save output
        logger.info("\n[STEP 5/5] Saving output...")
        output_dir = Path(__file__).parent / 'test_outputs'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = output_dir / f'test_{args.census_type.lower()}_species_merged_{args.level}_{timestamp}.csv'
        merged_df.to_csv(output_file, index=False)
        
        # Summary
        logger.info("\n" + "=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Census taxa: {len(census_df):,}")
        logger.info(f"Matched taxa: {(merged_df['match_status'] == 'matched').sum():,}")
        logger.info(f"Census-only taxa: {(merged_df['match_status'] == 'census_only').sum():,}")
        logger.info(f"Total NCBI genomes: {merged_df['ncbi_genome_count'].sum():,}")
        logger.info(f"Total NCBI species: {merged_df['ncbi_species_count'].sum():,}")
        logger.info(f"\nOutput file: {output_file}")
        logger.info(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 80)
        
        return output_file
        
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    try:
        output_file = main()
        print(f"\n✅ Test merger completed successfully!")
        print(f"📁 Output saved to: {output_file}")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

