"""
Species Parser Source Modules
==============================

Modular components for species-level NCBI assembly parsing.
"""

from .data_loader import load_assembly_data
from .genome_classifier import classify_genome_source
from .species_grouper import group_by_species
from .lineage_enricher import add_lineage_information
from .taxid_fallback_resolver import resolve_missing_lineages
from .percentage_calculator import add_species_percentages
from .output_writer import save_species_output

__all__ = [
    'load_assembly_data',
    'classify_genome_source',
    'group_by_species',
    'add_lineage_information',
    'resolve_missing_lineages',
    'add_species_percentages',
    'save_species_output'
]

