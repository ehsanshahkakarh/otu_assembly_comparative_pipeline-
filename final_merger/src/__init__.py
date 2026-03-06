"""
Source Modules for NCBI Census Merger
======================================

Modular components for merging census data with NCBI genome data.
"""

from .data_loader import load_census_data, load_ncbi_species_data
from .domain_filter import filter_by_domain
from .lineage_matcher import match_taxa_by_lineage
from .metrics_calculator import calculate_metrics
from .output_writer import save_merged_output

__all__ = [
    'load_census_data',
    'load_ncbi_species_data',
    'filter_by_domain',
    'match_taxa_by_lineage',
    'calculate_metrics',
    'save_merged_output'
]

