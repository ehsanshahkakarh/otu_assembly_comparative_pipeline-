#!/usr/bin/env python3
"""
GTDB-NCBI Accession-Based Merger
=================================

Merges GTDB taxonomy metadata with NCBI assembly metadata based on accession numbers.

This script uses configuration from config.yaml for all paths and settings.

Input files (configured in config.yaml):
- GTDB: Archaea taxonomy file
- GTDB: Bacteria taxonomy file
- NCBI: Assembly summary file

Output files (configured in config.yaml):
- archaea_mapped.csv - GTDB Archaea accessions that matched NCBI
- archaea_unmapped.csv - GTDB Archaea accessions NOT in NCBI
- bacteria_mapped.csv - GTDB Bacteria accessions that matched NCBI
- bacteria_unmapped.csv - GTDB Bacteria accessions NOT in NCBI
- ncbi_unmapped.csv - NCBI accessions NOT in GTDB (both domains)

Usage:
    python merge_gtdb_ncbi_accessions.py [--config CONFIG_PATH]
"""

import pandas as pd
from pathlib import Path
import logging
import sys
import argparse
from config_loader import load_config

# Global config and logger
config = None
logger = None


def setup_logging(config):
    """Setup logging based on configuration."""
    log_file = config.get_log_file('merger')
    log_level = config.get('logging', 'level', default='INFO')
    log_format = config.get('logging', 'format', default='%(asctime)s - %(levelname)s - %(message)s')

    logging.basicConfig(
        level=getattr(logging, log_level),
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode='w')
        ]
    )
    return logging.getLogger(__name__)


def load_gtdb_taxonomy(file_path: Path, domain: str) -> pd.DataFrame:
    """Load GTDB taxonomy file and extract accessions."""
    logger.info(f"Loading GTDB {domain} taxonomy from: {file_path}")

    df = pd.read_csv(file_path, sep='\t', header=None, names=['accession', 'taxonomy'])
    logger.info(f"  Loaded {len(df):,} {domain} records")

    # Clean accession: remove configured prefixes
    prefixes = config.get('processing', 'remove_prefixes', default=['RS_', 'GB_'])
    prefix_pattern = '^(' + '|'.join(prefixes) + ')'
    df['accession'] = df['accession'].str.replace(prefix_pattern, '', regex=True)
    logger.info(f"  Cleaned accessions (removed {', '.join(prefixes)} prefixes)")

    # Parse taxonomy string to extract all levels
    taxonomy_levels = config.get('processing', 'taxonomy_levels',
                                  default=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    # Map level names to GTDB prefixes
    level_prefixes = {
        'domain': 'd__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__',
        'family': 'f__', 'genus': 'g__', 'species': 's__'
    }

    for level in taxonomy_levels:
        if level in level_prefixes:
            prefix = level_prefixes[level]
            df[level] = df['taxonomy'].str.extract(f'{prefix}([^;]+)')

    # Add source domain
    df['gtdb_domain'] = domain

    return df


def load_ncbi_assembly(file_path: Path) -> pd.DataFrame:
    """Load NCBI assembly summary file."""
    logger.info(f"Loading NCBI assembly metadata from: {file_path}")

    df = pd.read_csv(file_path, sep='\t', skiprows=1, low_memory=False)

    # Clean column names (remove leading #)
    df.columns = df.columns.str.replace('^#', '', regex=True)

    logger.info(f"  Loaded {len(df):,} NCBI assembly records")

    # Get column names from config
    primary_col = config.get('processing', 'ncbi_columns', 'primary', default='assembly_accession')
    secondary_col = config.get('processing', 'ncbi_columns', 'secondary', default='gbrs_paired_asm')

    # Rename primary column for clarity
    if primary_col in df.columns:
        df = df.rename(columns={primary_col: 'accession_gca'})
    else:
        raise ValueError(f"Primary NCBI column '{primary_col}' not found in assembly file")

    # Fill NaN values in secondary column with empty string for easier processing
    if secondary_col in df.columns:
        df[secondary_col] = df[secondary_col].fillna('')
    else:
        logger.warning(f"Secondary NCBI column '{secondary_col}' not found, using empty values")
        df[secondary_col] = ''

    return df


def merge_and_split(gtdb_df: pd.DataFrame, ncbi_df: pd.DataFrame, domain: str) -> tuple:
    """
    Merge GTDB and NCBI on accession, return mapped and unmapped.
    Checks both GCA (assembly_accession) and GCF (gbrs_paired_asm) columns.

    Returns:
        (mapped_df, unmapped_df)
    """
    logger.info(f"\nMerging {domain} data...")

    # Get accession sets
    gtdb_accessions = set(gtdb_df['accession'])
    ncbi_gca_accessions = set(ncbi_df['accession_gca'])
    ncbi_gcf_accessions = set(ncbi_df['gbrs_paired_asm']) - {''}  # Remove empty strings

    # Combine both NCBI accession types
    all_ncbi_accessions = ncbi_gca_accessions | ncbi_gcf_accessions

    # Find common and unique
    common = gtdb_accessions & all_ncbi_accessions
    gtdb_only = gtdb_accessions - all_ncbi_accessions

    logger.info(f"  GTDB {domain} accessions: {len(gtdb_accessions):,}")
    logger.info(f"  NCBI GCA accessions: {len(ncbi_gca_accessions):,}")
    logger.info(f"  NCBI GCF accessions (paired): {len(ncbi_gcf_accessions):,}")
    logger.info(f"  Common with NCBI: {len(common):,} ({len(common)/len(gtdb_accessions)*100:.1f}%)")
    logger.info(f"  GTDB-only: {len(gtdb_only):,} ({len(gtdb_only)/len(gtdb_accessions)*100:.1f}%)")

    # Merge on accession - try both GCA and GCF columns
    # First, try matching on GCA (assembly_accession)
    mapped_gca = pd.merge(
        gtdb_df,
        ncbi_df,
        left_on='accession',
        right_on='accession_gca',
        how='inner'
    )
    mapped_gca['match_type'] = 'GCA'

    # Second, try matching on GCF (gbrs_paired_asm)
    # Only match non-empty GCF accessions
    ncbi_with_gcf = ncbi_df[ncbi_df['gbrs_paired_asm'] != ''].copy()
    mapped_gcf = pd.merge(
        gtdb_df,
        ncbi_with_gcf,
        left_on='accession',
        right_on='gbrs_paired_asm',
        how='inner'
    )
    mapped_gcf['match_type'] = 'GCF'

    # Combine both matches and remove duplicates (prefer GCA match if both exist)
    mapped = pd.concat([mapped_gca, mapped_gcf], ignore_index=True)
    mapped = mapped.drop_duplicates(subset=['accession'], keep='first')

    logger.info(f"  Matched via GCA: {len(mapped_gca):,}")
    logger.info(f"  Matched via GCF: {len(mapped_gcf):,}")
    logger.info(f"  Total unique matches: {len(mapped):,}")

    # Unmapped are those not in the common set
    unmapped = gtdb_df[gtdb_df['accession'].isin(gtdb_only)]

    return mapped, unmapped


def main():
    """Main execution."""
    global config, logger

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='GTDB-NCBI Accession-Based Merger')
    parser.add_argument('--config', type=str, default=None,
                       help='Path to config.yaml (default: ./config.yaml)')
    args = parser.parse_args()

    # Load configuration
    config_path = Path(args.config) if args.config else None
    config = load_config(config_path)

    # Setup logging
    logger = setup_logging(config)

    # Print header
    metadata = config.metadata
    logger.info("="*70)
    logger.info(f"{metadata.get('pipeline_name', 'GTDB-NCBI MAPPER')}")
    logger.info(f"Version: {metadata.get('version', 'unknown')}")
    logger.info("="*70)

    # Get input paths from config
    input_paths = config.get_input_paths()
    archaea_file = input_paths['archaea']
    bacteria_file = input_paths['bacteria']
    ncbi_file = input_paths['ncbi_assembly']

    # Verify input files exist
    for name, path in input_paths.items():
        if not path.exists():
            logger.error(f"Input file not found: {name} = {path}")
            logger.error("Please check config.yaml and ensure all input files exist")
            sys.exit(1)

    # Load data
    archaea_df = load_gtdb_taxonomy(archaea_file, 'Archaea')
    bacteria_df = load_gtdb_taxonomy(bacteria_file, 'Bacteria')
    ncbi_df = load_ncbi_assembly(ncbi_file)

    # Merge Archaea
    archaea_mapped, archaea_unmapped = merge_and_split(archaea_df, ncbi_df, 'Archaea')

    # Merge Bacteria
    bacteria_mapped, bacteria_unmapped = merge_and_split(bacteria_df, ncbi_df, 'Bacteria')

    # Find NCBI accessions not in GTDB
    # Check both GCA and GCF columns
    all_gtdb_accessions = set(archaea_df['accession']) | set(bacteria_df['accession'])

    secondary_col = config.get('processing', 'ncbi_columns', 'secondary', default='gbrs_paired_asm')

    # NCBI accessions that don't match GTDB in either GCA or GCF column
    ncbi_gca_only = set(ncbi_df['accession_gca']) - all_gtdb_accessions
    ncbi_gcf_only = set(ncbi_df[secondary_col]) - all_gtdb_accessions - {''}

    # An NCBI record is "unmapped" only if BOTH its GCA and GCF don't match GTDB
    ncbi_unmapped = ncbi_df[
        (ncbi_df['accession_gca'].isin(ncbi_gca_only)) &
        ((ncbi_df[secondary_col] == '') | (ncbi_df[secondary_col].isin(ncbi_gcf_only)))
    ]

    logger.info(f"\nNCBI records not in GTDB (neither GCA nor GCF matched): {len(ncbi_unmapped):,}")

    # Save outputs using config
    logger.info("\nSaving output files...")

    output_files = {
        'archaea_mapped': (archaea_mapped, 'raw', 'archaea_mapped'),
        'archaea_unmapped': (archaea_unmapped, 'raw', 'archaea_unmapped'),
        'bacteria_mapped': (bacteria_mapped, 'raw', 'bacteria_mapped'),
        'bacteria_unmapped': (bacteria_unmapped, 'raw', 'bacteria_unmapped'),
        'ncbi_unmapped': (ncbi_unmapped, 'raw', 'ncbi_unmapped')
    }

    for name, (df, category, config_key) in output_files.items():
        output_path = config.get_output_file(category, config_key)
        df.to_csv(output_path, index=False)
        logger.info(f"  ✅ {output_path.name} ({len(df):,} records)")

    logger.info("\n" + "="*70)
    logger.info("MERGE COMPLETE!")
    logger.info("="*70)
    logger.info(f"\nOutput directory: {config.get_path('output', 'raw')}")


if __name__ == '__main__':
    main()

