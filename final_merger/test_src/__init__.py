#!/usr/bin/env python3
"""
Test Source Modules for Species-Grouped Based Merger
======================================================

This test implementation uses species_grouped_*.csv as the data source
and searches it using census taxa as the backbone.
"""

from .species_data_loader import load_species_grouped_data
from .census_synonym_builder import build_census_synonym_dict
from .species_searcher import search_species_for_census_taxon
from .species_aggregator import aggregate_species_matches

__all__ = [
    'load_species_grouped_data',
    'build_census_synonym_dict',
    'search_species_for_census_taxon',
    'aggregate_species_matches',
]

