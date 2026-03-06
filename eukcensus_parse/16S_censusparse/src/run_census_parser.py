"""
Enhanced EukCensus 16S Parser - Main Entry Point
================================================

PERFORMANCE OPTIMIZED VERSION:
- Disabled individual organelle recovery during initial processing
- All organellar sequences now handled by vectorized approach
- This prevents thousands of individual subprocess calls
- Expected speedup: 10-100x for datasets with many organellar sequences

Processes EukCensus 16S cluster data and generates taxonomic summaries
at division (phylum), family, and genus levels with NCBI taxonomy integration.
"""

import sys
import os
import csv
import time
import logging
import pandas as pd
from tqdm import tqdm

from .config import PathConfig, setup_logging
from .level_processor import process_taxonomic_level
from .taxonkit_utils import get_taxids_using_taxonkit, get_lineages_using_taxonkit
from .lineage_processor import append_name_to_lineage
from .unmapped_logger import create_comprehensive_unmapped_log
from .cached_ai_resolver import CachedAIResolver
from .resolution_applier import apply_all_resolutions


def main():
    """
    Main entry point for the 16S census parser.
    
    PERFORMANCE OPTIMIZED VERSION:
    - Disabled individual organelle recovery during initial processing
    - All organellar sequences now handled by vectorized approach
    - This prevents thousands of individual subprocess calls
    - Expected speedup: 10-100x for datasets with many organellar sequences
    """
    # Set up directory paths using PathConfig
    config = PathConfig()
    metadata_dir = config.metadata_dir
    csv_output_dir = config.csv_output_dir
    log_dir = config.log_dir

    # Parse command line arguments for custom input/output
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        # If relative path provided, make it relative to metadata directory
        if not os.path.isabs(input_file):
            input_file = metadata_dir / input_file
    else:
        input_file = config.input_file

    if len(sys.argv) > 2:
        output_prefix = sys.argv[2]
    else:
        output_prefix = "eukcensus16S"

    # Set up logging with proper log directory
    setup_logging(log_dir, output_prefix)
    start_time = time.time()

    # Output file paths - using csv_16S directory with new naming
    phylum_output = csv_output_dir / f"{output_prefix}_by_division.csv"
    family_output = csv_output_dir / f"{output_prefix}_by_family.csv"
    genus_output = csv_output_dir / f"{output_prefix}_by_genus.csv"

    # Log the paths being used
    logging.info(f"Input file: {input_file}")
    logging.info(f"Output directory: {csv_output_dir}")
    logging.info(f"Log directory: {log_dir}")

    print(f"Processing file: {input_file}")

    # Read the TSV file
    try:
        logging.info("Loading input file...")

        # Use chunked reading with progress bar for large files
        chunk_size = 50000
        chunks = []

        logging.info(f"Loading file in chunks of {chunk_size:,}...")

        # Read chunks with progress bar
        chunk_iterator = pd.read_csv(input_file, sep='\t', chunksize=chunk_size)
        for chunk in tqdm(chunk_iterator, desc="Loading chunks", unit="chunk"):
            chunks.append(chunk)

        # Combine chunks if multiple chunks were read
        if len(chunks) > 1:
            logging.info(f"Combining {len(chunks)} chunks...")
            df = pd.concat(chunks, ignore_index=True)
        else:
            df = chunks[0]

        logging.info(f"Successfully loaded {len(df)} rows")
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        return

    # Check if required columns exist
    required_columns = ['centroid', 'members', 'size', 'phylum', 'familiy', 'genus']
    for col in required_columns:
        if col not in df.columns:
            logging.error(f"Required column '{col}' not found in the input file")
            return

    # Process each taxonomic level with enhanced logic
    logging.info("Processing taxonomic levels...")
    phylum_data = process_taxonomic_level(df, 'phylum', 'phylum')
    family_data = process_taxonomic_level(df, 'familiy', 'family')
    genus_data = process_taxonomic_level(df, 'genus', 'genus')

    # Get taxids for each taxonomic level using appropriate names for lookup
    logging.info("Getting taxids using enhanced taxonkit processing...")

    # OPTIMIZATION: Vectorized collection and processing of unique names
    logging.info("Vectorized optimization: Collecting unique names across all ranks...")

    all_unique_names = set()
    all_unique_names.update(phylum_data.keys())
    all_unique_names.update(family_data.keys())
    all_unique_names.update(genus_data.keys())

    logging.info(f"Total unique names across all ranks: {len(all_unique_names)}")

    # Single vectorized lookup for all unique names (with optimized organelle handling)
    all_names_to_taxid = get_taxids_using_taxonkit(list(all_unique_names), "all_ranks")

    # Map results to individual ranks
    phylum_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in phylum_data.keys()}
    family_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in family_data.keys()}
    genus_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in genus_data.keys()}

    # Collect all valid taxids with progress tracking
    logging.info("Collecting taxids for lineage retrieval...")
    all_taxids = set()

    taxid_sources = [
        ("phylum", phylum_to_taxid),
        ("family", family_to_taxid),
        ("genus", genus_to_taxid)
    ]

    for _, taxid_dict in tqdm(taxid_sources, desc="Collecting taxids", unit="source"):
        for taxid in taxid_dict.values():
            if taxid != "NA":
                all_taxids.add(taxid)

    logging.info(f"Collected {len(all_taxids)} unique taxids for lineage retrieval")

    # Get lineages for all taxids
    taxid_to_lineage = get_lineages_using_taxonkit(list(all_taxids))

    # Calculate totals for percentage calculations
    logging.info("Calculating database totals for percentage calculations...")

    # Calculate total OTU count and total size count across all taxonomic levels
    total_otu_count = sum(data['otu_count'] for data in phylum_data.values())
    total_size_count = sum(data['size_count'] for data in phylum_data.values())

    logging.info(f"📊 Total OTUs in database: {total_otu_count:,}")
    logging.info(f"📊 Total size count in database: {total_size_count:,}")

    # Write results to CSV files with enhanced information including percentages
    logging.info("Writing results to CSV files...")

    # Write phylum data (now called division to match 18S)
    logging.info(f"Writing division data to {phylum_output}")
    with open(phylum_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'otu_percentage', 'size_count', 'size_percentage', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_phyla = sorted(phylum_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)

        for phylum, data in tqdm(sorted_phyla, desc="Writing division data", unit="entry"):
            taxid = phylum_to_taxid.get(phylum, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            env = os.environ.copy()
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, phylum, taxid, env, taxid_to_lineage
            )

            # Calculate percentages
            otu_percentage = round((data['otu_count'] / total_otu_count * 100), 2) if total_otu_count > 0 else 0
            size_percentage = round((data['size_count'] / total_size_count * 100), 2) if total_size_count > 0 else 0

            writer.writerow([
                phylum,  # Original name preserved
                taxid,
                data['otu_count'],
                otu_percentage,
                data['size_count'],
                size_percentage,
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Write family data
    logging.info(f"Writing family data to {family_output}")
    with open(family_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'otu_percentage', 'size_count', 'size_percentage', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_families = sorted(family_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)

        for family, data in tqdm(sorted_families, desc="Writing family data", unit="entry"):
            taxid = family_to_taxid.get(family, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            env = os.environ.copy()
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, family, taxid, env, taxid_to_lineage
            )

            # Calculate percentages
            otu_percentage = round((data['otu_count'] / total_otu_count * 100), 2) if total_otu_count > 0 else 0
            size_percentage = round((data['size_count'] / total_size_count * 100), 2) if total_size_count > 0 else 0

            writer.writerow([
                family,  # Original name preserved
                taxid,
                data['otu_count'],
                otu_percentage,
                data['size_count'],
                size_percentage,
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Write genus data
    logging.info(f"Writing genus data to {genus_output}")
    with open(genus_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'otu_percentage', 'size_count', 'size_percentage', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_genera = sorted(genus_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)

        for genus, data in tqdm(sorted_genera, desc="Writing genus data", unit="entry"):
            taxid = genus_to_taxid.get(genus, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            env = os.environ.copy()
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, genus, taxid, env, taxid_to_lineage
            )

            # Calculate percentages
            otu_percentage = round((data['otu_count'] / total_otu_count * 100), 2) if total_otu_count > 0 else 0
            size_percentage = round((data['size_count'] / total_size_count * 100), 2) if total_size_count > 0 else 0

            writer.writerow([
                genus,  # Original name preserved
                taxid,
                data['otu_count'],
                otu_percentage,
                data['size_count'],
                size_percentage,
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Create comprehensive unmapped log
    create_comprehensive_unmapped_log(
        phylum_data, family_data, genus_data,
        phylum_to_taxid, family_to_taxid, genus_to_taxid,
        taxid_to_lineage, log_dir, output_prefix
    )

    # ========== AI RESOLVER STEP ==========
    # Resolve unmapped taxa using AI-assisted resolution
    logging.info("\n" + "="*80)
    logging.info("STEP 4: AI-ASSISTED RESOLUTION")
    logging.info("="*80)

    # Collect unmapped taxa for each rank
    unmapped_phyla = [name for name, taxid in phylum_to_taxid.items() if taxid == "NA"]
    unmapped_families = [name for name, taxid in family_to_taxid.items() if taxid == "NA"]
    unmapped_genera = [name for name, taxid in genus_to_taxid.items() if taxid == "NA"]

    logging.info(f"Unmapped taxa counts:")
    logging.info(f"  Phyla: {len(unmapped_phyla)}")
    logging.info(f"  Families: {len(unmapped_families)}")
    logging.info(f"  Genera: {len(unmapped_genera)}")

    # Initialize resolver
    resolver = CachedAIResolver()

    # Resolve each rank
    resolutions_by_rank = {}

    if unmapped_phyla:
        logging.info(f"\nResolving {len(unmapped_phyla)} unmapped phyla...")
        resolutions_by_rank['phylum'] = resolver.resolve_unmapped_taxa(unmapped_phyla, rank="phylum")

    if unmapped_families:
        logging.info(f"\nResolving {len(unmapped_families)} unmapped families...")
        resolutions_by_rank['family'] = resolver.resolve_unmapped_taxa(unmapped_families, rank="family")

    if unmapped_genera:
        logging.info(f"\nResolving {len(unmapped_genera)} unmapped genera...")
        resolutions_by_rank['genus'] = resolver.resolve_unmapped_taxa(unmapped_genera, rank="genus")

    # Apply resolutions to CSV files
    if any(resolutions_by_rank.values()):
        applied_counts = apply_all_resolutions(
            csv_output_dir, log_dir, resolutions_by_rank,
            output_prefix, total_otu_count, total_size_count
        )

        # Log resolution statistics
        total_resolved = sum(applied_counts.values())
        logging.info(f"\n✅ Successfully resolved and applied {total_resolved} taxa!")
    else:
        logging.info("\nNo resolutions to apply (all taxa already mapped or need AI research)")

    # Print cache statistics
    stats = resolver.get_cache_statistics()
    logging.info(f"\nAI Cache Statistics:")
    logging.info(f"  Total cached resolutions: {stats['total_resolutions']}")
    logging.info(f"  Validated: {stats['validated']}")
    logging.info(f"  Unvalidated: {stats['unvalidated']}")

    # ========== END AI RESOLVER STEP ==========

    # Calculate and log performance metrics
    end_time = time.time()
    processing_time = end_time - start_time
    total_entries = len(phylum_data) + len(family_data) + len(genus_data)

    # Calculate success rates
    phylum_success = len([t for t in phylum_to_taxid.values() if t != 'NA'])
    family_success = len([t for t in family_to_taxid.values() if t != 'NA'])
    genus_success = len([t for t in genus_to_taxid.values() if t != 'NA'])

    phylum_success_rate = (phylum_success / len(phylum_data) * 100) if len(phylum_data) > 0 else 0
    family_success_rate = (family_success / len(family_data) * 100) if len(family_data) > 0 else 0
    genus_success_rate = (genus_success / len(genus_data) * 100) if len(genus_data) > 0 else 0

    logging.info("Saving results to CSV files...")
    logging.info(f"Saved {len(phylum_data)} division entries to {phylum_output}")
    logging.info(f"Division summary: {phylum_success} with taxids, {len([t for t in phylum_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Saved {len(family_data)} family entries to {family_output}")
    logging.info(f"Family summary: {family_success} with taxids, {len([t for t in family_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Saved {len(genus_data)} genus entries to {genus_output}")
    logging.info(f"Genus summary: {genus_success} with taxids, {len([t for t in genus_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Processing complete in {processing_time:.2f} seconds")
    logging.info(f"Total entries processed: {total_entries}")
    logging.info(f"Performance: {total_entries/processing_time:.1f} entries/second")

    logging.info("Processing complete! Generated the following files:")
    logging.info(f"- {phylum_output}")
    logging.info(f"- {family_output}")
    logging.info(f"- {genus_output}")
    logging.info(f"- {log_dir / f'{output_prefix}_comprehensive_unmapped.log'}")
    logging.info(f"- {log_dir / 'eukcensus16S_processing.log'}")
    logging.info("Output files organized in csv_16S/ and logs/ directories")


if __name__ == "__main__":
    main()

