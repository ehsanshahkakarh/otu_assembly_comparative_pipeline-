#!/usr/bin/env python3
"""
16S NCBI Merger - Main Entry Point
===================================

Merges 16S eukcensus (prokaryotic) data with NCBI genome data.

Usage:
    python run_16s_ncbi_merger.py
"""

import sys
from pathlib import Path
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import Config, setup_logging
from src.data_loader import load_census_data, load_ncbi_species_data
from src.domain_filter import filter_by_domain
from src.lineage_matcher import match_taxa_by_lineage
from src.metrics_calculator import calculate_metrics
from src.output_writer import save_merged_output
from tqdm import tqdm


def main():
    """Main execution function for 16S NCBI merger."""
    # Initialize configuration
    config = Config()
    census_type = '16S'
    
    # Setup logging
    logger = setup_logging(config.paths.logs_dir, census_type)
    
    logger.info("=" * 80)
    logger.info("16S NCBI MERGER - Prokaryotic Census + NCBI Genomes")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Census type: {census_type}")
    logger.info(f"Domains: {', '.join(config.domains.get_domains(census_type))}")
    
    try:
        # Load NCBI species data once (used for all levels)
        logger.info("\n" + "=" * 80)
        logger.info("LOADING NCBI SPECIES DATA")
        logger.info("=" * 80)
        species_file = config.paths.get_ncbi_species_file()
        ncbi_species_df = load_ncbi_species_data(species_file)

        # Filter by domain once
        logger.info("\nFiltering by domain...")
        domains = config.domains.get_domains(census_type)
        ncbi_species_df = filter_by_domain(ncbi_species_df, domains)

        # Process each taxonomic level with progress bar
        logger.info("\n" + "=" * 80)
        logger.info("PROCESSING TAXONOMIC LEVELS")
        logger.info("=" * 80)

        for level in tqdm(config.taxonomic.ncbi_levels, desc="Processing levels", unit="level"):
            # Map to census level name
            census_level = 'division' if level == 'phylum' else level

            logger.info("\n" + "=" * 80)
            logger.info(f"PROCESSING {level.upper()} LEVEL (census: {census_level})")
            logger.info("=" * 80)

            # Step 1: Load census data
            logger.info("\n[STEP 1/5] Loading census data...")
            census_file = config.paths.get_census_file(census_type, level)
            census_df = load_census_data(census_file)

            # Step 2: Match taxa by lineage (aggregation happens here)
            logger.info("\n[STEP 2/5] Matching taxa by lineage and aggregating...")
            matched_df = match_taxa_by_lineage(census_df, ncbi_species_df, level, census_level)

            # Step 3: Calculate metrics
            logger.info("\n[STEP 3/5] Calculating metrics...")
            matched_df = calculate_metrics(matched_df)

            # Step 4: Save output
            logger.info("\n[STEP 4/5] Saving output...")
            output_file = config.paths.get_output_file(census_type, level)
            save_merged_output(matched_df, output_file, census_level)

            logger.info(f"\n✅ {level.upper()} level complete!")
        
        # Final summary
        logger.info("\n" + "=" * 80)
        logger.info("16S NCBI MERGER COMPLETE")
        logger.info("=" * 80)
        logger.info("Output files:")
        for level in config.taxonomic.ncbi_levels:
            output_file = config.paths.get_output_file(census_type, level)
            logger.info(f"  - {output_file}")
        logger.info(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 80)
        
        return True
        
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    try:
        success = main()
        print("\n🎉 16S NCBI merger completed successfully!")
        print(f"📁 Output files saved to: new_merger/outputs/")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

