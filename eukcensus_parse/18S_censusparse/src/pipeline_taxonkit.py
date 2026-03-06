#!/usr/bin/env python3
"""
Pipeline module for Taxonkit Parser

Processes EukCensus 18S cluster data using ONLY taxonkit for NCBI taxonomy lookups.
NO systematic resolution - pure NCBI taxonomy mapping.
"""

import os
import sys
import time
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Handle imports for both module and script execution
try:
    from .config import setup_directory_paths, setup_logging, setup_taxonkit_environment
    from .level_processor import process_taxonomic_level, write_level_to_csv
    from .taxonkit_utils import get_taxids_for_names, get_lineages_for_taxids
    from .unmapped_logger import create_unmapped_log
except ImportError:
    # Running as script or imported from parent, use absolute imports
    from config import setup_directory_paths, setup_logging, setup_taxonkit_environment
    from level_processor import process_taxonomic_level, write_level_to_csv
    from taxonkit_utils import get_taxids_for_names, get_lineages_for_taxids
    from unmapped_logger import create_unmapped_log


def run_taxonkit_parser(input_file=None, output_prefix=None):
    """Run the taxonkit parser pipeline step."""
    # Set up directory paths
    paths = setup_directory_paths()

    # Use provided arguments or defaults
    if input_file is None:
        input_file = paths.metadata_dir / "eukcensus_18S.clusters.97.tsv"
    elif not os.path.isabs(input_file):
        input_file = paths.metadata_dir / input_file

    if output_prefix is None:
        output_prefix = "eukcensus_taxonkit_only"

    # Set up logging
    setup_logging(paths.log_dir, output_prefix)
    start_time = time.time()

    logging.info("=" * 80)
    logging.info("18S TAXONKIT PARSER - Clean NCBI-only version")
    logging.info("NO systematic resolution - pure taxonkit lookups")
    logging.info("=" * 80)

    # Output file paths
    division_output = paths.csv_output_dir / f"{output_prefix}_by_division.csv"
    family_output = paths.csv_output_dir / f"{output_prefix}_by_family.csv"
    genus_output = paths.csv_output_dir / f"{output_prefix}_by_genus.csv"
    unmapped_log = paths.log_dir / f"{output_prefix}_unmapped.log"

    # Set up taxonkit environment
    setup_taxonkit_environment()

    # Read input file
    logging.info(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    logging.info(f"Total clusters: {len(df):,}")

    # Initialize dictionaries to store grouped data
    division_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})
    family_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})
    genus_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})

    # Process each taxonomic level to aggregate data
    logging.info("\n" + "=" * 80)
    logging.info("AGGREGATING TAXONOMIC DATA")
    logging.info("=" * 80)

    process_taxonomic_level(df, 'division', division_data)
    process_taxonomic_level(df, 'family', family_data)
    process_taxonomic_level(df, 'genus', genus_data)

    # Collect all unique names for batch processing
    logging.info("\n" + "=" * 80)
    logging.info("GETTING TAXIDS AND LINEAGES")
    logging.info("=" * 80)

    all_names = set()
    all_names.update(division_data.keys())
    all_names.update(family_data.keys())
    all_names.update(genus_data.keys())

    logging.info(f"Total unique names across all levels: {len(all_names):,}")

    # Get taxids for all names at once
    all_taxid_results, all_failed_names = get_taxids_for_names(list(all_names), "all_ranks")

    # Set up environment for taxonkit
    env = setup_taxonkit_environment()

    # Get lineages for all taxids
    all_taxids = [taxid for taxid in all_taxid_results.values() if taxid != "NA"]
    taxid_to_lineage = get_lineages_for_taxids(all_taxids, env)

    # Write results to CSV files
    logging.info("\n" + "=" * 80)
    logging.info("WRITING OUTPUT FILES")
    logging.info("=" * 80)

    division_filtered = write_level_to_csv(
        division_output, 'division', division_data, all_taxid_results, taxid_to_lineage, env
    )
    family_filtered = write_level_to_csv(
        family_output, 'family', family_data, all_taxid_results, taxid_to_lineage, env
    )
    genus_filtered = write_level_to_csv(
        genus_output, 'genus', genus_data, all_taxid_results, taxid_to_lineage, env
    )

    # Create unmapped log
    # Note: We use all_taxid_results for all three levels since we did a combined lookup
    # The unmapped_logger will filter by the actual names in each data_dict
    create_unmapped_log(
        division_data,
        family_data,
        genus_data,
        all_taxid_results,  # division_to_taxid
        all_taxid_results,  # family_to_taxid
        all_taxid_results,  # genus_to_taxid
        taxid_to_lineage,
        paths.log_dir,
        output_prefix
    )

    # Summary
    elapsed_time = time.time() - start_time
    logging.info("\n" + "=" * 80)
    logging.info("TAXONKIT PARSER COMPLETE")
    logging.info("=" * 80)
    logging.info(f"Total time: {elapsed_time:.1f} seconds")
    logging.info(f"\nOutput files:")
    logging.info(f"  - {division_output}")
    logging.info(f"  - {family_output}")
    logging.info(f"  - {genus_output}")
    logging.info(f"  - {unmapped_log}")
    logging.info("\nNext step: Run systematic resolver to handle unmapped taxa")


def run_taxonkit_pipeline():
    """Wrapper function for pipeline integration."""
    run_taxonkit_parser()
