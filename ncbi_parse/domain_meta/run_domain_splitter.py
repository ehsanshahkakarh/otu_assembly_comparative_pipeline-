#!/usr/bin/env python3
"""
Domain Splitter Runner
=======================

Main script to split ncbi_parse output into domain-specific CSV files.
Updated to work with current ncbi_parse directory structure.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from domain_splitter import split_by_domain, get_domain_summary


def setup_logging(log_dir: Path) -> None:
    """Configure logging to file and console."""
    log_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f'domain_splitter_{timestamp}.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logger = logging.getLogger('domain_splitter')
    logger.info(f"Log file: {log_file}")
    return logger


def find_latest_species_file(output_dir: Path) -> Path:
    """Find the most recent species_grouped CSV file."""
    species_files = list(output_dir.glob('species_grouped_*.csv'))
    
    if not species_files:
        raise FileNotFoundError(
            f"No species_grouped_*.csv files found in {output_dir}"
        )
    
    # Sort by modification time, most recent first
    latest_file = max(species_files, key=lambda p: p.stat().st_mtime)
    return latest_file


def main():
    """Main execution function."""
    # Setup paths - UPDATED to point to current ncbi_parse structure
    base_dir = Path(__file__).parent
    ncbi_parse_dir = base_dir.parent  # Go up to ncbi_parse directory
    ncbi_parse_output_dir = ncbi_parse_dir / 'output'
    output_dir = base_dir / 'output'
    log_dir = base_dir / 'logs'
    
    # Setup logging
    logger = setup_logging(log_dir)
    
    logger.info("=" * 80)
    logger.info("Domain Splitter - Split species data by taxonomic domain")
    logger.info("=" * 80)
    logger.info(f"Working directory: {ncbi_parse_dir}")
    logger.info(f"Input directory: {ncbi_parse_output_dir}")
    
    try:
        # Find latest species_grouped file
        input_file = find_latest_species_file(ncbi_parse_output_dir)
        logger.info(f"\nInput file: {input_file}")
        logger.info(f"File date: {datetime.fromtimestamp(input_file.stat().st_mtime)}")
        
        # Split by domain
        logger.info("\n" + "=" * 80)
        logger.info("Splitting species data by domain...")
        logger.info("=" * 80)
        
        output_files = split_by_domain(input_file, output_dir)
        
        # Generate summary
        logger.info("\n" + "=" * 80)
        logger.info("Domain Summary Statistics")
        logger.info("=" * 80)
        
        summary_df = get_domain_summary(output_dir)
        
        # Save summary
        summary_file = output_dir / f'domain_summary_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"\n✅ Summary saved: {summary_file.name}\n")
        
        # Print summary table
        print("\n" + "=" * 80)
        print("DOMAIN SUMMARY")
        print("=" * 80)
        for _, row in summary_df.iterrows():
            print(f"\n{row['domain'].upper()}")
            print(f"  Species: {row['species_count']:,}")
            print(f"  Total genomes: {row['total_genomes']:,}")
            print(f"  Isolate genomes: {row['isolate_genomes']:,} ({row['isolate_percentage']:.1f}%)")
            print(f"  Uncultured genomes: {row['uncultured_genomes']:,}")
            print(f"  Avg genomes/species: {row['avg_genomes_per_species']:.1f}")
        
        print("\n" + "=" * 80)
        print("✅ Domain splitting complete!")
        print("=" * 80)
        print(f"\nOutput directory: {output_dir}")
        print(f"Files created:")
        for domain, filepath in output_files.items():
            print(f"  - {filepath.name}")
        print(f"  - {summary_file.name}")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()

