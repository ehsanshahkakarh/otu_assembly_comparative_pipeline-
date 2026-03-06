"""
16S Census Parser - Modular Package
====================================

A modular package for processing EukCensus 16S cluster data with enhanced
organelle handling and comprehensive taxa preservation.

Modules:
    config: Configuration and directory path setup
    organelle_handler: Organelle detection and cleaning
    taxon_cleaner: Taxonomic name cleaning and validation
    taxonkit_utils: NCBI taxonkit integration
    lineage_processor: Lineage handling and manipulation
    level_processor: Taxonomic level processing
    unmapped_logger: Unmapped taxa logging and analysis
    run_census_parser: Main entry point

Usage:
    python -m census_parser.run_census_parser
"""

__version__ = "1.0.0"
__author__ = "Enhanced EukCensus Analysis Team"

from . import config
from . import organelle_handler
from . import taxon_cleaner
from . import taxonkit_utils
from . import lineage_processor
from . import level_processor
from . import unmapped_logger
from . import run_census_parser

__all__ = [
    'config',
    'organelle_handler',
    'taxon_cleaner',
    'taxonkit_utils',
    'lineage_processor',
    'level_processor',
    'unmapped_logger',
    'run_census_parser',
]

