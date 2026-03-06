#!/usr/bin/env python3
"""
Configuration Module for Species Parser
========================================

Centralized configuration for species-level NCBI assembly parsing.
"""

import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List


@dataclass
class PathConfig:
    """File path configuration."""
    # base_dir should point to nev_parse_meth directory
    base_dir: Path = field(default_factory=lambda: Path(__file__).resolve().parent)

    # Input paths
    assembly_file: Path = None

    # Output paths
    output_dir: Path = None
    logs_dir: Path = None

    def __post_init__(self):
        """Set default paths based on base_dir."""
        if self.assembly_file is None:
            # Look for assembly file in parent metadata directory
            self.assembly_file = self.base_dir.parent / '00assembly_summary_genbank.txt'
        if self.output_dir is None:
            self.output_dir = self.base_dir / 'output'
        if self.logs_dir is None:
            self.logs_dir = self.base_dir / 'logs'


@dataclass
class FilterConfig:
    """Filter and classification configuration."""

    # Uncultured patterns for genome source classification
    # Matches ncbi_parse/py_ncbi/ncbi_parser/src/genome_classifier.py
    uncultured_patterns: List[str] = field(default_factory=lambda: [
        'uncultured', 'environmental', 'metagenome', 'unclassified',
        'unknown', 'unidentified', 'mixed culture', 'enrichment culture',
        'derived from metagenome', 'metagenome-assembled', 'mag',
        'single amplified genome', 'sag', 'environmental sample'
    ])
    
    # TaxonKit batch size for processing
    taxonkit_batch_size: int = 1000
    
    # Taxonomic ranks to extract
    taxonomic_ranks: List[str] = field(default_factory=lambda: [
        'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
    ])


@dataclass
class Config:
    """Main configuration container."""
    paths: PathConfig = field(default_factory=PathConfig)
    filters: FilterConfig = field(default_factory=FilterConfig)


def setup_logging(logs_dir: Path) -> logging.Logger:
    """Set up logging configuration."""
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger('species_parser')
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
    log_file = logs_dir / f'species_parser_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(console_format)
    logger.addHandler(file_handler)
    
    return logger

