#!/usr/bin/env python3
"""
18S Species-Based Merger - Complete Pipeline
=============================================

Runs the complete species-based merger for 18S census data at all three levels:
- Division
- Family  
- Genus

Outputs results to 18S_output/ directory.

Usage:
    python run_18s_species_merger.py
"""

import sys
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


def setup_logging(output_dir: Path):
    """Setup logging configuration with file output."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = output_dir / 'logs' / f'18s_species_merger_{timestamp}.log'
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


def run_single_level(level: str, census_df, species_df, config, output_dir: Path, logger):
    """Run merger for a single taxonomic level."""
    logger.info("=" * 80)
    logger.info(f"PROCESSING 18S {level.upper()} LEVEL")
    logger.info("=" * 80)
    
    try:
        # Build synonym dictionary for census taxa
        logger.info(f"Building synonym dictionary for {level} taxa...")
        census_synonym_dict = build_census_synonym_dict(census_df, level)
        
        # Search and aggregate for each census taxon
        logger.info(f"Searching species for each census {level} taxon...")
        merged_df = aggregate_all_census_taxa(
            census_df,
            species_df,
            census_synonym_dict,
            level,
            search_species_for_census_taxon
        )
        
        # Save output
        logger.info("Saving output...")
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = output_dir / f'18s_species_merged_{level}_{timestamp}.csv'
        merged_df.to_csv(output_file, index=False)
        
        # Summary
        logger.info("\n" + "=" * 80)
        logger.info(f"SUMMARY - {level.upper()}")
        logger.info("=" * 80)
        logger.info(f"Census taxa: {len(census_df):,}")
        logger.info(f"Matched taxa: {(merged_df['match_status'] == 'matched').sum():,}")
        logger.info(f"Census-only taxa: {(merged_df['match_status'] == 'census_only').sum():,}")
        logger.info(f"Total NCBI genomes: {merged_df['ncbi_genome_count'].sum():,}")
        logger.info(f"Total NCBI species: {merged_df['ncbi_species_count'].sum():,}")
        logger.info(f"Output file: {output_file}")
        logger.info("=" * 80)
        
        return output_file
        
    except Exception as e:
        logger.error(f"Error processing {level} level: {e}", exc_info=True)
        raise


def main():
    """Main execution function."""
    # Setup output directory
    output_dir = Path(__file__).parent / '18S_output'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    logger = setup_logging(output_dir)
    
    logger.info("=" * 80)
    logger.info("18S SPECIES-BASED MERGER - COMPLETE PIPELINE")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Output directory: {output_dir}")
    
    try:
        # Initialize configuration
        config = Config()
        
        # Load species_grouped data once (used for all levels)
        logger.info("\nLoading species_grouped data...")
        species_file = config.paths.parse_dir / 'ncbi_parse' / 'metadata' / 'nev_parse_meth' / 'output' / 'species_grouped_20260223_014859.csv'
        species_df = load_species_grouped_data(species_file)
        
        # Process all three levels
        levels = ['division', 'family', 'genus']
        output_files = []
        
        for level in levels:
            logger.info(f"\n{'=' * 80}")
            logger.info(f"LEVEL: {level.upper()}")
            logger.info(f"{'=' * 80}\n")
            
            # Load census data for this level
            census_file = config.paths.get_census_file('18S', 'phylum' if level == 'division' else level)
            logger.info(f"Loading census data from: {census_file.name}")
            census_df = load_census_data(census_file)
            
            # Run merger for this level
            output_file = run_single_level(level, census_df, species_df, config, output_dir, logger)
            output_files.append((level, output_file))
        
        # Final summary
        logger.info("\n" + "=" * 80)
        logger.info("COMPLETE PIPELINE SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"\nOutput files:")
        for level, output_file in output_files:
            logger.info(f"  {level:10s}: {output_file.name}")
        logger.info("=" * 80)
        
        return output_files
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    try:
        output_files = main()
        print(f"\n✅ 18S species merger completed successfully!")
        print(f"📁 Output directory: 18S_output/")
        print(f"\nGenerated files:")
        for level, output_file in output_files:
            print(f"  - {output_file.name}")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

