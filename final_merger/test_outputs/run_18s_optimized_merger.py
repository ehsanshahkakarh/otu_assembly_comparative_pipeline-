#!/usr/bin/env python3
"""
18S Optimized Species-Based Merger
===================================

Optimized version with 3-4x speedup:
- Pre-build synonym dictionary once for all taxa
- Pre-process species lineage data
- Vectorized pandas operations

Expected performance:
- Division (22 taxa):   ~15 seconds  (vs 40 sec)
- Family (314 taxa):    ~5 minutes   (vs 20 min)
- Genus (491 taxa):     ~8 minutes   (vs 30 min)
- TOTAL:                ~13 minutes  (vs 50 min)

Usage:
    python run_18s_optimized_merger.py
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
from optimized_src.optimized_synonym_builder import build_all_census_synonyms
from optimized_src.optimized_species_searcher import preprocess_species_data, search_species_optimized
from optimized_src.optimized_aggregator import aggregate_species_optimized


def setup_logging(output_dir: Path):
    """Setup logging configuration with file output."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = output_dir / 'logs' / f'18s_optimized_merger_{timestamp}.log'
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


def run_single_level(level: str, census_df, species_df_preprocessed, all_synonyms, output_dir: Path, logger):
    """Run optimized merger for a single taxonomic level."""
    logger.info("=" * 80)
    logger.info(f"PROCESSING 18S {level.upper()} LEVEL")
    logger.info("=" * 80)
    
    try:
        start_time = datetime.now()
        
        # Aggregate using pre-built synonym dict and pre-processed species data
        logger.info(f"Searching and aggregating for {level} taxa...")
        merged_df = aggregate_species_optimized(
            census_df,
            species_df_preprocessed,
            all_synonyms,
            level,
            search_species_optimized
        )
        
        # Save output
        logger.info("Saving output...")
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = output_dir / f'18s_optimized_merged_{level}_{timestamp}.csv'
        merged_df.to_csv(output_file, index=False)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        
        # Summary
        logger.info("\n" + "=" * 80)
        logger.info(f"SUMMARY - {level.upper()}")
        logger.info("=" * 80)
        logger.info(f"Census taxa: {len(census_df):,}")
        logger.info(f"Matched taxa: {(merged_df['match_status'] == 'matched').sum():,}")
        logger.info(f"Census-only taxa: {(merged_df['match_status'] == 'census_only').sum():,}")
        logger.info(f"Total NCBI genomes: {merged_df['ncbi_genome_count'].sum():,}")
        logger.info(f"Total NCBI species: {merged_df['ncbi_species_count'].sum():,}")
        logger.info(f"Processing time: {elapsed:.1f} seconds")
        logger.info(f"Output file: {output_file}")
        logger.info("=" * 80)
        
        return output_file, elapsed
        
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
    logger.info("18S OPTIMIZED SPECIES-BASED MERGER")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Output directory: {output_dir}")
    
    overall_start = datetime.now()
    
    try:
        # Initialize configuration
        config = Config()
        
        # OPTIMIZATION 1: Load species_grouped data once and pre-process
        logger.info("\n[OPTIMIZATION] Loading and pre-processing species data...")
        species_file = config.paths.parse_dir / 'ncbi_parse' / 'metadata' / 'nev_parse_meth' / 'output' / 'species_grouped_20260223_014859.csv'
        species_df = load_species_grouped_data(species_file)
        species_df_preprocessed = preprocess_species_data(species_df)
        
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
            
            # OPTIMIZATION 2: Build synonym dict once for all taxa at this level
            logger.info(f"\n[OPTIMIZATION] Building synonym dictionary for all {len(census_df)} {level} taxa...")
            all_synonyms = build_all_census_synonyms(census_df, level)
            
            # Run optimized merger for this level
            output_file, elapsed = run_single_level(
                level, census_df, species_df_preprocessed, all_synonyms, output_dir, logger
            )
            output_files.append((level, output_file, elapsed))
        
        overall_elapsed = (datetime.now() - overall_start).total_seconds()
        
        # Final summary
        logger.info("\n" + "=" * 80)
        logger.info("COMPLETE PIPELINE SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Total processing time: {overall_elapsed/60:.1f} minutes ({overall_elapsed:.0f} seconds)")
        logger.info(f"\nOutput files:")
        for level, output_file, elapsed in output_files:
            logger.info(f"  {level:10s}: {output_file.name} ({elapsed:.1f}s)")
        logger.info("=" * 80)
        
        return output_files
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    try:
        output_files = main()
        print(f"\n✅ 18S optimized merger completed successfully!")
        print(f"📁 Output directory: 18S_output/")
        print(f"\nGenerated files:")
        for level, output_file, elapsed in output_files:
            print(f"  - {output_file.name} ({elapsed:.1f}s)")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

