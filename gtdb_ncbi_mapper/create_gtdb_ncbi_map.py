#!/usr/bin/env python3
"""
Create Clean GTDB-NCBI Mapping Tables
======================================

Creates simplified mapping tables between GTDB and NCBI taxonomies.

This script uses configuration from config.yaml for all paths and settings.

Output files (configured in config.yaml):
1. gtdb_to_ncbi_full_map.csv - Complete mapping with all fields
2. gtdb_to_ncbi_simple_map.csv - Simplified accession + taxonomy mapping
3. gtdb_to_ncbi_taxonomy_map.csv - Taxonomy comparison only

Usage:
    python create_gtdb_ncbi_map.py [--config CONFIG_PATH]
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
    log_level = config.get('logging', 'level', default='INFO')
    log_format = config.get('logging', 'format', default='%(asctime)s - %(levelname)s - %(message)s')

    # Try to get mapping log file, fall back to stdout only
    try:
        log_file = config.get_log_file('mapping')
        handlers = [
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode='w')
        ]
    except (ValueError, KeyError):
        handlers = [logging.StreamHandler(sys.stdout)]

    logging.basicConfig(
        level=getattr(logging, log_level),
        format=log_format,
        handlers=handlers
    )
    return logging.getLogger(__name__)


def create_full_map(archaea_df: pd.DataFrame, bacteria_df: pd.DataFrame, output_dir: Path):
    """Create complete GTDB-NCBI mapping with all fields."""
    logger.info("Creating full GTDB-NCBI map...")
    
    # Combine both domains
    full_map = pd.concat([archaea_df, bacteria_df], ignore_index=True)
    
    # Save
    output_file = output_dir / 'gtdb_to_ncbi_full_map.csv'
    full_map.to_csv(output_file, index=False)
    logger.info(f"✅ Saved full map: {output_file}")
    logger.info(f"   Total records: {len(full_map):,}")
    
    return full_map


def create_simple_map(full_map: pd.DataFrame, output_dir: Path):
    """Create simplified accession + taxonomy mapping."""
    logger.info("\nCreating simple GTDB-NCBI map...")
    
    # Select key columns
    simple_map = full_map[[
        'accession',
        'match_type',
        'gtdb_domain',
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'organism_name',
        'taxid',
        'species_taxid',
        'assembly_level',
        'genome_size',
        'gc_percent'
    ]].copy()
    
    # Rename for clarity
    simple_map.rename(columns={
        'domain': 'ncbi_domain',
        'phylum': 'gtdb_phylum',
        'class': 'gtdb_class',
        'order': 'gtdb_order',
        'family': 'gtdb_family',
        'genus': 'gtdb_genus',
        'species': 'gtdb_species',
        'organism_name': 'ncbi_organism_name',
        'taxid': 'ncbi_taxid',
        'species_taxid': 'ncbi_species_taxid'
    }, inplace=True)
    
    # Save
    output_file = output_dir / 'gtdb_to_ncbi_simple_map.csv'
    simple_map.to_csv(output_file, index=False)
    logger.info(f"✅ Saved simple map: {output_file}")
    logger.info(f"   Total records: {len(simple_map):,}")
    
    return simple_map


def create_taxonomy_comparison(full_map: pd.DataFrame, output_dir: Path):
    """Create taxonomy-only comparison table."""
    logger.info("\nCreating taxonomy comparison map...")
    
    # Parse GTDB taxonomy string
    def parse_gtdb_taxonomy(tax_string):
        """Parse GTDB taxonomy string into components."""
        parts = tax_string.split(';')
        return {
            'gtdb_domain': parts[0].replace('d__', '') if len(parts) > 0 else '',
            'gtdb_phylum': parts[1].replace('p__', '') if len(parts) > 1 else '',
            'gtdb_class': parts[2].replace('c__', '') if len(parts) > 2 else '',
            'gtdb_order': parts[3].replace('o__', '') if len(parts) > 3 else '',
            'gtdb_family': parts[4].replace('f__', '') if len(parts) > 4 else '',
            'gtdb_genus': parts[5].replace('g__', '') if len(parts) > 5 else '',
            'gtdb_species': parts[6].replace('s__', '') if len(parts) > 6 else ''
        }
    
    # Parse taxonomy for all records
    taxonomy_data = full_map['taxonomy'].apply(parse_gtdb_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_data.tolist())
    
    # Combine with accession and NCBI organism name
    taxonomy_map = pd.concat([
        full_map[['accession', 'organism_name', 'taxid', 'match_type']],
        taxonomy_df
    ], axis=1)
    
    # Rename for clarity
    taxonomy_map.rename(columns={
        'organism_name': 'ncbi_organism_name',
        'taxid': 'ncbi_taxid'
    }, inplace=True)
    
    # Save
    output_file = output_dir / 'gtdb_to_ncbi_taxonomy_map.csv'
    taxonomy_map.to_csv(output_file, index=False)
    logger.info(f"✅ Saved taxonomy map: {output_file}")
    logger.info(f"   Total records: {len(taxonomy_map):,}")
    
    return taxonomy_map


def generate_summary_stats(full_map: pd.DataFrame, output_dir: Path):
    """Generate summary statistics."""
    logger.info("\n" + "="*70)
    logger.info("MAPPING SUMMARY STATISTICS")
    logger.info("="*70)
    
    # Overall stats
    logger.info(f"\nTotal mapped genomes: {len(full_map):,}")
    logger.info(f"  Archaea: {len(full_map[full_map['gtdb_domain'] == 'Archaea']):,}")
    logger.info(f"  Bacteria: {len(full_map[full_map['gtdb_domain'] == 'Bacteria']):,}")
    
    # Match type breakdown
    logger.info("\nMatch type breakdown:")
    match_counts = full_map['match_type'].value_counts()
    for match_type, count in match_counts.items():
        pct = (count / len(full_map)) * 100
        logger.info(f"  {match_type}: {count:,} ({pct:.1f}%)")
    
    # Assembly level breakdown
    logger.info("\nAssembly level breakdown:")
    assembly_counts = full_map['assembly_level'].value_counts()
    for level, count in assembly_counts.head(5).items():
        pct = (count / len(full_map)) * 100
        logger.info(f"  {level}: {count:,} ({pct:.1f}%)")
    
    # Unique phyla
    unique_phyla = full_map['phylum'].nunique()
    logger.info(f"\nUnique GTDB phyla represented: {unique_phyla}")
    
    # Top 10 phyla by genome count
    logger.info("\nTop 10 phyla by genome count:")
    top_phyla = full_map['phylum'].value_counts().head(10)
    for phylum, count in top_phyla.items():
        logger.info(f"  {phylum}: {count:,}")


def main():
    """Main execution."""
    global config, logger

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='GTDB-NCBI Mapping Table Generator')
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
    logger.info(f"{metadata.get('pipeline_name', 'GTDB-NCBI MAPPER')} - Mapping Table Generator")
    logger.info(f"Version: {metadata.get('version', 'unknown')}")
    logger.info("="*70)

    # Get paths from config
    input_dir = config.get_path('output', 'raw')
    output_dir = config.get_path('output', 'mapping_tables', create_dir=True)

    # Load mapped data
    logger.info("\nLoading mapped data...")
    archaea_file = config.get_output_file('raw', 'archaea_mapped')
    bacteria_file = config.get_output_file('raw', 'bacteria_mapped')

    # Check if input files exist
    if not archaea_file.exists() or not bacteria_file.exists():
        logger.error("Input files not found!")
        logger.error(f"  Archaea: {archaea_file}")
        logger.error(f"  Bacteria: {bacteria_file}")
        logger.error("\nPlease run merge_gtdb_ncbi_accessions.py first to generate these files.")
        sys.exit(1)

    archaea_df = pd.read_csv(archaea_file)
    bacteria_df = pd.read_csv(bacteria_file)

    logger.info(f"  Archaea: {len(archaea_df):,} records")
    logger.info(f"  Bacteria: {len(bacteria_df):,} records")

    # Create mapping tables
    full_map = create_full_map(archaea_df, bacteria_df, output_dir)
    simple_map = create_simple_map(full_map, output_dir)
    taxonomy_map = create_taxonomy_comparison(full_map, output_dir)

    # Generate summary stats
    generate_summary_stats(full_map, output_dir)

    logger.info("\n" + "="*70)
    logger.info("✅ MAPPING TABLES CREATED SUCCESSFULLY!")
    logger.info("="*70)
    logger.info(f"\nOutput directory: {output_dir}")
    logger.info("\nFiles created:")
    logger.info(f"  1. {config.get('output_files', 'mapping', 'full')} - Complete mapping (all fields)")
    logger.info(f"  2. {config.get('output_files', 'mapping', 'simple')} - Simplified mapping (key fields)")
    logger.info(f"  3. {config.get('output_files', 'mapping', 'taxonomy')} - Taxonomy comparison only")


if __name__ == '__main__':
    main()

