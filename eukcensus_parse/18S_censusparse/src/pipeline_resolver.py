#!/usr/bin/env python3
"""
Pipeline module for Systematic Resolver

Reads unmapped taxa (families and genera) from taxonkit parser and applies
systematic resolutions using the known_parents database.
"""

import sys
import logging
from pathlib import Path

# Handle imports for both module and script execution
try:
    from .config import setup_directory_paths, setup_taxonkit_environment
    from .resolution_builder import build_all_resolutions, save_resolutions
    from .resolution_applier import (
        load_resolutions,
        apply_resolutions_to_csv,
        create_final_unmapped_log
    )
    from .known_parents import get_statistics
except ImportError:
    # Running as script or imported from parent, use absolute imports
    from config import setup_directory_paths, setup_taxonkit_environment
    from resolution_builder import build_all_resolutions, save_resolutions
    from resolution_applier import (
        load_resolutions,
        apply_resolutions_to_csv,
        create_final_unmapped_log
    )
    from known_parents import get_statistics


def setup_resolver_logging(log_dir: Path):
    """Set up logging for the systematic resolver."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "systematic_resolver.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def run_systematic_resolver():
    """Run the systematic resolver pipeline step."""
    # Set up paths
    paths = setup_directory_paths()
    setup_resolver_logging(paths.log_dir)
    
    logging.info("=" * 80)
    logging.info("SYSTEMATIC RESOLVER (Families + Genera)")
    logging.info("=" * 80)
    
    # Show statistics about known_parents database
    stats = get_statistics()
    logging.info(f"\nKnown Parents Database:")
    logging.info(f"  Total families: {stats['total_families']}")
    logging.info(f"  Total genera: {stats['total_genera']}")
    logging.info(f"  Unique divisions: {stats['unique_divisions']}")
    
    # Define paths
    resolver_dir = Path(__file__).parent.parent / "systematic_resolver"
    outputs_dir = resolver_dir / "outputs"
    outputs_dir.mkdir(exist_ok=True)
    
    resolutions_file = outputs_dir / "systematic_resolutions.json"
    unmapped_log = paths.log_dir / "eukcensus_taxonkit_only_unmapped.log"
    
    # Step 1: Build resolutions from unmapped log
    logging.info("\n" + "=" * 80)
    logging.info("STEP 1: Building Resolutions")
    logging.info("=" * 80)
    
    if not unmapped_log.exists():
        logging.error(f"Unmapped log not found: {unmapped_log}")
        logging.error("Please run taxonkit parser first!")
        sys.exit(1)
    
    resolutions = build_all_resolutions(unmapped_log)
    save_resolutions(resolutions, resolutions_file)
    
    # Step 2: Apply resolutions to CSV files
    logging.info("\n" + "=" * 80)
    logging.info("STEP 2: Applying Resolutions to CSV Files")
    logging.info("=" * 80)
    
    resolutions = load_resolutions(resolutions_file)
    
    # Apply to each level
    for level in ['division', 'family', 'genus']:
        input_csv = paths.csv_output_dir / f"eukcensus_taxonkit_only_by_{level}.csv"
        output_csv = paths.csv_output_dir / f"eukcensus_18S_by_{level}.csv"
        
        if input_csv.exists():
            apply_resolutions_to_csv(input_csv, output_csv, resolutions, level)
        else:
            logging.warning(f"Input CSV not found: {input_csv}")
    
    # Step 3: Create final unmapped log
    logging.info("\n" + "=" * 80)
    logging.info("STEP 3: Creating Final Unmapped Log")
    logging.info("=" * 80)
    
    final_unmapped_log = paths.log_dir / "eukcensus_18S_unmapped_final.log"
    create_final_unmapped_log(
        paths.csv_output_dir / "eukcensus_18S_by_division.csv",
        paths.csv_output_dir / "eukcensus_18S_by_family.csv",
        paths.csv_output_dir / "eukcensus_18S_by_genus.csv",
        final_unmapped_log
    )
    
    # Step 4: Clean up intermediate taxonkit_only files
    logging.info("\n" + "=" * 80)
    logging.info("STEP 4: Cleaning Up Intermediate Files")
    logging.info("=" * 80)

    intermediate_files = [
        paths.csv_output_dir / "eukcensus_taxonkit_only_by_division.csv",
        paths.csv_output_dir / "eukcensus_taxonkit_only_by_family.csv",
        paths.csv_output_dir / "eukcensus_taxonkit_only_by_genus.csv",
    ]

    for intermediate_file in intermediate_files:
        if intermediate_file.exists():
            intermediate_file.unlink()
            logging.info(f"Deleted intermediate file: {intermediate_file.name}")

    # Also clean up intermediate unmapped log
    if unmapped_log.exists():
        unmapped_log.unlink()
        logging.info(f"Deleted intermediate log: {unmapped_log.name}")

    logging.info("\n" + "=" * 80)
    logging.info("SYSTEMATIC RESOLVER COMPLETE")
    logging.info("=" * 80)
    logging.info(f"\nFinal output files:")
    logging.info(f"  - {paths.csv_output_dir / 'eukcensus_18S_by_division.csv'}")
    logging.info(f"  - {paths.csv_output_dir / 'eukcensus_18S_by_family.csv'}")
    logging.info(f"  - {paths.csv_output_dir / 'eukcensus_18S_by_genus.csv'}")
    logging.info(f"  - {final_unmapped_log}")
    logging.info(f"\nIntermediate 'taxonkit_only' files have been removed.")


def run_resolver_pipeline():
    """Wrapper function for pipeline integration."""
    run_systematic_resolver()
