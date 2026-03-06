#!/usr/bin/env python3
"""
Configuration Module for NCBI Census Merger
============================================

Centralized configuration for merging census data (16S/18S) with NCBI genome data.
"""

import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class PathConfig:
    """File path configuration."""
    # base_dir should point to new_merger directory
    base_dir: Path = field(default_factory=lambda: Path(__file__).resolve().parent)
    
    # Parent directory (parse_repaa_table)
    parse_dir: Path = field(init=False)
    
    # Input directories
    ncbi_species_dir: Path = field(init=False)  # Species-level data from nev_parse_meth
    census_18s_dir: Path = field(init=False)
    census_16s_dir: Path = field(init=False)

    # Output directories
    output_dir: Path = field(init=False)
    logs_dir: Path = field(init=False)

    def __post_init__(self):
        """Set default paths based on base_dir."""
        self.parse_dir = self.base_dir.parent

        # Input directories
        self.ncbi_species_dir = self.parse_dir / 'ncbi_parse' / 'output'
        self.census_18s_dir = self.parse_dir / '18S_censusparse' / 'output'
        self.census_16s_dir = self.parse_dir / '16S_censusparse' / 'output'
        
        # Output directories
        self.output_dir = self.base_dir / 'outputs'
        self.logs_dir = self.base_dir / 'logs'
        
        # Create output directories if they don't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
    
    def get_ncbi_species_file(self) -> Path:
        """Get the most recent species_grouped_*.csv file from nev_parse_meth output."""
        import glob

        # Find all species_grouped files
        pattern = str(self.ncbi_species_dir / 'species_grouped_*.csv')
        species_files = glob.glob(pattern)

        if not species_files:
            raise FileNotFoundError(f"No species_grouped_*.csv files found in {self.ncbi_species_dir}")

        # Return the most recent one (sorted by filename which includes timestamp)
        most_recent = sorted(species_files)[-1]
        return Path(most_recent)
    
    def get_census_file(self, census_type: str, level: str) -> Path:
        """Get census file path for a taxonomic level."""
        if census_type == '18S':
            # 18S uses 'division' instead of 'phylum'
            census_level = 'division' if level == 'phylum' else level
            return self.census_18s_dir / f'eukcensus_18S_by_{census_level}.csv'
        elif census_type == '16S':
            # 16S uses 'division' instead of 'phylum'
            census_level = 'division' if level == 'phylum' else level
            return self.census_16s_dir / f'eukcensus16S_by_{census_level}.csv'
        else:
            raise ValueError(f"Invalid census_type: {census_type}. Must be '18S' or '16S'")
    
    def get_output_file(self, census_type: str, level: str) -> Path:
        """Get output file path for a taxonomic level."""
        # Convert phylum to division for output naming
        output_level = 'division' if level == 'phylum' else level
        return self.output_dir / f'{census_type.lower()}_ncbi_merged_{output_level}.csv'


@dataclass
class DomainConfig:
    """Domain filtering configuration."""
    
    # Domain filters for each census type
    domains_18s: List[str] = field(default_factory=lambda: ['Eukaryota'])
    domains_16s: List[str] = field(default_factory=lambda: ['Bacteria', 'Archaea'])
    
    def get_domains(self, census_type: str) -> List[str]:
        """Get domain list for a census type."""
        if census_type == '18S':
            return self.domains_18s
        elif census_type == '16S':
            return self.domains_16s
        else:
            raise ValueError(f"Invalid census_type: {census_type}. Must be '18S' or '16S'")


@dataclass
class TaxonomicConfig:
    """Taxonomic level configuration."""
    
    # NCBI taxonomic levels
    ncbi_levels: List[str] = field(default_factory=lambda: ['phylum', 'family', 'genus'])
    
    # Census level mapping (census name -> NCBI name)
    census_level_map: Dict[str, str] = field(default_factory=lambda: {
        'division': 'phylum',
        'family': 'family',
        'genus': 'genus'
    })


@dataclass
class Config:
    """Main configuration container."""
    paths: PathConfig = field(default_factory=PathConfig)
    domains: DomainConfig = field(default_factory=DomainConfig)
    taxonomic: TaxonomicConfig = field(default_factory=TaxonomicConfig)


def setup_logging(logs_dir: Path, census_type: str) -> logging.Logger:
    """
    Set up logging configuration.
    
    Args:
        logs_dir: Directory for log files
        census_type: '18S' or '16S'
    
    Returns:
        Configured logger instance
    """
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger(f'ncbi_merger_{census_type}')
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler
    from datetime import datetime
    log_file = logs_dir / f'{census_type.lower()}_ncbi_merger_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(console_format)
    logger.addHandler(file_handler)
    
    return logger

