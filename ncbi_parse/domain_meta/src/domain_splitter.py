#!/usr/bin/env python3
"""
Domain Splitter Module
=======================

Splits species-grouped data from nev_parse_meth into separate CSV files by domain.
Extracts domain from lineage and creates domain-specific CSV files.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger('domain_splitter')


def extract_domain(lineage: str, lineage_ranks: str) -> str:
    """
    Extract domain from lineage string.
    
    Args:
        lineage: Semicolon-separated lineage string
        lineage_ranks: Semicolon-separated ranks string
        
    Returns:
        Domain name (Bacteria, Archaea, Eukaryota, Viruses) or 'Unknown'
    """
    if pd.isna(lineage) or not lineage:
        return 'Unknown'
    
    lineage_parts = lineage.split(';')
    
    # Handle viruses (acellular)
    if lineage_parts[0].strip() == 'Viruses':
        return 'Viruses'
    
    # Handle cellular organisms (domain is at index 1)
    if len(lineage_parts) > 1:
        domain = lineage_parts[1].strip()
        if domain in ['Bacteria', 'Archaea', 'Eukaryota']:
            return domain
    
    return 'Unknown'


def split_by_domain(input_file: Path, output_dir: Path) -> Dict[str, Path]:
    """
    Split species-grouped data by domain and save to separate CSV files.
    
    Args:
        input_file: Path to species_grouped CSV file from nev_parse_meth
        output_dir: Directory to save domain-specific CSV files
        
    Returns:
        Dictionary mapping domain names to output file paths
    """
    logger.info(f"Loading species-grouped data from: {input_file.name}")
    
    # Load the data
    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df):,} species records")
    
    # Extract domain for each species
    df['domain'] = df.apply(
        lambda row: extract_domain(row['lineage'], row['lineage_ranks']),
        axis=1
    )
    
    # Count species by domain
    domain_counts = df['domain'].value_counts()
    logger.info("\nSpecies distribution by domain:")
    for domain, count in domain_counts.items():
        logger.info(f"  {domain}: {count:,} species")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Split and save by domain
    output_files = {}
    for domain in domain_counts.index:
        # Filter data for this domain
        domain_df = df[df['domain'] == domain].copy()
        
        # Sort by total_genome_count descending
        domain_df = domain_df.sort_values('total_genome_count', ascending=False)
        
        # Calculate domain-specific statistics
        total_genomes = domain_df['total_genome_count'].sum()
        total_isolates = domain_df['isolate_genome_count'].sum()
        total_uncultured = domain_df['uncultured_genome_count'].sum()
        
        # Save to CSV
        output_file = output_dir / f'{domain.lower()}_species.csv'
        domain_df.to_csv(output_file, index=False)
        output_files[domain] = output_file
        
        logger.info(f"\n✅ Saved {domain} data: {output_file.name}")
        logger.info(f"   Species: {len(domain_df):,}")
        logger.info(f"   Total genomes: {total_genomes:,}")
        logger.info(f"   Isolate genomes: {total_isolates:,} ({total_isolates/total_genomes*100:.1f}%)")
        logger.info(f"   Uncultured genomes: {total_uncultured:,} ({total_uncultured/total_genomes*100:.1f}%)")
        logger.info(f"   File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    return output_files


def get_domain_summary(output_dir: Path) -> pd.DataFrame:
    """
    Generate summary statistics for all domains.
    
    Args:
        output_dir: Directory containing domain-specific CSV files
        
    Returns:
        DataFrame with summary statistics for each domain
    """
    summary_data = []
    
    for domain_file in output_dir.glob('*_species.csv'):
        domain = domain_file.stem.replace('_species', '').capitalize()
        df = pd.read_csv(domain_file)
        
        summary_data.append({
            'domain': domain,
            'species_count': len(df),
            'total_genomes': df['total_genome_count'].sum(),
            'isolate_genomes': df['isolate_genome_count'].sum(),
            'uncultured_genomes': df['uncultured_genome_count'].sum(),
            'isolate_percentage': (df['isolate_genome_count'].sum() / df['total_genome_count'].sum() * 100),
            'avg_genomes_per_species': df['total_genome_count'].mean()
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('total_genomes', ascending=False)
    
    return summary_df

