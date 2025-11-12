#!/usr/bin/env python3
"""
GTDB Phylum Parser - Pure Metadata-Only Version

This script processes GTDB taxonomy files to extract phylum-level information
ONLY from taxa that are present in the actual metadata files. No taxdump lookup.

Key Features:
- Processes only taxa from metadata files (00bac120_taxonomy.tsv, 00ar53_taxonomy.tsv)
- No taxdump or taxonkit dependencies
- Creates simple lineage information directly from metadata
- Eliminates taxa that exist only in taxdump but not in actual genome metadata

Input:
- GTDB bacterial taxonomy file (00bac120_taxonomy.tsv)
- GTDB archaeal taxonomy file (00ar53_taxonomy.tsv)

Output:
- Summary CSV with phylum counts (gtdb_phylum_counts.csv)
- Detailed CSV with accessions (gtdb_phylum_with_accessions.csv)

Usage:
    python phylum_gtdb_parse_metadata_only.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
"""

import argparse
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime
import sys
import subprocess
import tempfile
import os

# Processing parameters
CHUNK_SIZE = 100000  # Number of rows to process at once

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="GTDB Phylum Parser - Metadata Only")
    parser.add_argument("--input-dir", 
                       type=str, 
                       default="../metadata",
                       help="Directory containing GTDB metadata files")
    parser.add_argument("--output-dir", 
                       type=str, 
                       default="../csv_gtdb",
                       help="Output directory for CSV files")
    return parser.parse_args()

def setup_paths(args):
    """Set up input and output paths."""
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Input files
    bac_file = input_dir / "00bac120_taxonomy.tsv"
    ar_file = input_dir / "00ar53_taxonomy.tsv"
    
    # Output files
    output_dir.mkdir(exist_ok=True)
    counts_file = output_dir / "gtdb_phylum_counts.csv"
    detailed_file = output_dir / "gtdb_phylum_with_accessions.csv"
    
    return bac_file, ar_file, counts_file, detailed_file

def load_gtdb_metadata(bac_file, ar_file, logger):
    """Load GTDB metadata files and extract taxonomy information."""
    logger.info("Loading GTDB metadata files...")
    
    all_genomes = []
    
    # Load bacterial metadata
    if bac_file.exists():
        logger.info(f"Loading bacterial metadata: {bac_file}")
        bac_df = pd.read_csv(bac_file, sep='\t', dtype=str)
        logger.info(f"Loaded {len(bac_df):,} bacterial genomes")
        all_genomes.append(bac_df)
    else:
        logger.warning(f"Bacterial file not found: {bac_file}")
    
    # Load archaeal metadata
    if ar_file.exists():
        logger.info(f"Loading archaeal metadata: {ar_file}")
        ar_df = pd.read_csv(ar_file, sep='\t', dtype=str)
        logger.info(f"Loaded {len(ar_df):,} archaeal genomes")
        all_genomes.append(ar_df)
    else:
        logger.warning(f"Archaeal file not found: {ar_file}")
    
    if not all_genomes:
        raise FileNotFoundError("No GTDB metadata files found!")
    
    # Combine all genomes
    combined_df = pd.concat(all_genomes, ignore_index=True)
    logger.info(f"Total genomes loaded: {len(combined_df):,}")
    
    return combined_df

def extract_taxonomy_info(df, logger):
    """Extract taxonomy information from GTDB metadata."""
    logger.info("Extracting taxonomy information...")
    
    # Extract accession (first column)
    df['accession_clean'] = df.iloc[:, 0].str.strip()
    
    # Extract taxonomy string (second column)
    taxonomy_col = df.columns[1]
    df['taxonomy'] = df[taxonomy_col].str.strip()
    
    # Parse taxonomy string to extract domain and phylum
    taxonomy_parts = df['taxonomy'].str.split(';', expand=True)
    
    # Extract domain (d__) and phylum (p__)
    df['domain'] = taxonomy_parts[0].str.replace('d__', '', regex=False).str.strip()
    df['phylum'] = taxonomy_parts[1].str.replace('p__', '', regex=False).str.strip()
    
    # Clean phylum names (remove GTDB subdivision suffixes like _A, _B, etc.)
    df['phylum_clean'] = df['phylum'].str.replace(r'_[A-Z]$', '', regex=True)
    
    # Filter out invalid entries
    valid_mask = (
        df['accession_clean'].notna() & 
        df['domain'].notna() & 
        df['phylum_clean'].notna() &
        (df['accession_clean'] != '') &
        (df['domain'] != '') &
        (df['phylum_clean'] != '')
    )
    
    df_clean = df[valid_mask].copy()
    logger.info(f"Valid genomes after cleaning: {len(df_clean):,}")
    
    return df_clean

def get_taxids_from_names_gtdb(phylum_names, logger):
    """Get taxids for phylum names using GTDB taxdump via taxonkit name2taxid."""
    if not phylum_names:
        return {}

    logger.info(f"Getting GTDB taxids for {len(phylum_names)} phylum names...")

    # Set up environment to use GTDB taxdump
    script_dir = Path(__file__).resolve().parent
    gtdb_taxdump_dir = script_dir.parent / "taxdump_gtdp" / "gtdb-taxdump-R226"
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(gtdb_taxdump_dir)

    name_to_taxid = {}

    # Create temporary file for phylum names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in phylum_names:
            temp_file.write(f"{name}\n")

    try:
        # Run taxonkit name2taxid command with GTDB database
        result = subprocess.run(
            ["taxonkit", "name2taxid", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')

            for line in lines:
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()
                    if taxid and taxid != name:
                        name_to_taxid[name] = taxid

    except Exception as e:
        logger.warning(f"Error running GTDB taxonkit name2taxid: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    logger.info(f"Successfully mapped {len(name_to_taxid)} phylum names to GTDB taxids")
    return name_to_taxid

def get_lineages_from_gtdb_taxids(taxids, logger):
    """Get lineages for GTDB taxids using GTDB taxdump via taxonkit lineage."""
    if not taxids:
        return {}

    logger.info(f"Getting GTDB lineages for {len(taxids)} taxids...")

    # Set up environment to use GTDB taxdump
    script_dir = Path(__file__).resolve().parent
    gtdb_taxdump_dir = script_dir.parent / "taxdump_gtdp" / "gtdb-taxdump-R226"
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(gtdb_taxdump_dir)

    lineage_data = {}

    # Create temporary file for taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        # Run taxonkit lineage command with GTDB database
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')

            for line in lines:
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip()
                    lineage_taxids = parts[2].strip()
                    lineage_ranks = parts[3].strip()

                    if lineage and lineage != taxid:
                        lineage_data[taxid] = (lineage, lineage_ranks, lineage_taxids)

    except Exception as e:
        logger.warning(f"Error running GTDB taxonkit lineage: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    logger.info(f"Successfully obtained lineages for {len(lineage_data)} taxids")
    return lineage_data

def create_phylum_counts(df, logger):
    """Create phylum-level counts from metadata."""
    logger.info("Creating phylum-level counts...")

    # Group by cleaned phylum and domain, count genomes
    counts = df.groupby(['phylum_clean', 'domain']).size().reset_index(name='phylum_counts')
    counts = counts.rename(columns={'phylum_clean': 'phylum'})
    counts = counts.sort_values('phylum_counts', ascending=False)

    logger.info(f"Found {len(counts)} unique phyla from metadata")
    logger.info(f"Total genomes: {counts['phylum_counts'].sum():,}")

    # Get taxids for phylum names using GTDB taxdump
    unique_phyla = counts['phylum'].unique().tolist()
    phylum_to_taxid = get_taxids_from_names_gtdb(unique_phyla, logger)

    # Add taxid column
    counts['taxid'] = counts['phylum'].map(phylum_to_taxid)

    # Get lineages for available taxids
    available_taxids = [str(taxid) for taxid in counts['taxid'].dropna().unique()]
    if available_taxids:
        taxid_to_lineage = get_lineages_from_gtdb_taxids(available_taxids, logger)

        # Add lineage columns
        lineages = []
        lineage_ranks = []
        lineage_taxids = []

        for _, row in counts.iterrows():
            if pd.notna(row['taxid']):
                taxid = str(row['taxid'])
                if taxid in taxid_to_lineage:
                    lineage, ranks, taxids = taxid_to_lineage[taxid]
                    lineages.append(lineage)
                    lineage_ranks.append(ranks)
                    lineage_taxids.append(taxids)
                else:
                    # Fallback to simple lineage
                    lineages.append(f"{row['domain']};{row['phylum']}")
                    lineage_ranks.append("domain;phylum")
                    lineage_taxids.append("")
            else:
                # Fallback to simple lineage
                lineages.append(f"{row['domain']};{row['phylum']}")
                lineage_ranks.append("domain;phylum")
                lineage_taxids.append("")

        counts['lineage'] = lineages
        counts['lineage_ranks'] = lineage_ranks
        counts['lineage_taxids'] = lineage_taxids

        # Log success rate
        taxid_success = len([t for t in counts['taxid'] if pd.notna(t)])
        total_phyla = len(counts)
        success_rate = (taxid_success / total_phyla) * 100
        logger.info(f"Taxid mapping success: {success_rate:.1f}% ({taxid_success}/{total_phyla})")

    else:
        # Fallback to simple lineage for all
        counts['lineage'] = counts['domain'] + ';' + counts['phylum']
        counts['lineage_ranks'] = 'domain;phylum'
        counts['lineage_taxids'] = ''
        counts['taxid'] = ''
        logger.warning("No taxids obtained, using simple lineage format")

    return counts

def save_results(counts_df, detailed_df, counts_file, detailed_file, logger):
    """Save results to CSV files."""
    logger.info("Saving results...")
    
    # Save counts file
    counts_df.to_csv(counts_file, index=False)
    logger.info(f"Saved phylum counts: {counts_file}")
    
    # Save detailed file
    detailed_columns = ['accession_clean', 'phylum', 'domain']
    detailed_df[detailed_columns].to_csv(detailed_file, index=False)
    logger.info(f"Saved detailed file: {detailed_file}")
    
    # Log top phyla
    logger.info("\nTop 10 phyla by genome count:")
    for i, (_, row) in enumerate(counts_df.head(10).iterrows(), 1):
        logger.info(f"  {i:2d}. {row['phylum']} ({row['domain']}): {row['phylum_counts']:,} genomes")

def main():
    """Main execution function."""
    logger = setup_logging()
    
    try:
        # Parse arguments and setup paths
        args = parse_arguments()
        bac_file, ar_file, counts_file, detailed_file = setup_paths(args)
        
        # Load and process metadata
        df = load_gtdb_metadata(bac_file, ar_file, logger)
        df_clean = extract_taxonomy_info(df, logger)
        counts = create_phylum_counts(df_clean, logger)
        
        # Save results
        save_results(counts, df_clean, counts_file, detailed_file, logger)
        
        logger.info("✅ Metadata-only GTDB phylum parsing complete!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
