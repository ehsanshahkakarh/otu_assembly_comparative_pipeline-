#!/usr/bin/env python3
"""
Taxon validation module for 18S Census Parser

Handles filtering and validation logic for eukaryotic taxonomic names.
"""


def should_filter_taxon(taxon_name: str) -> bool:
    """
    Check if a taxon name should be filtered out from final output.

    Modified to KEEP unclassified (.U.) entries as they represent significant biological diversity.
    Based on analysis showing .U. entries represent 20% of division, 59% of family, and 72% of genus data.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    # Previously filtered out .U. entries, but analysis shows this causes massive data loss:
    # - 20% loss at division level
    # - 59% loss at family level
    # - 72% loss at genus level
    # Therefore, we now KEEP all .U. entries to preserve biological diversity

    # Only filter out truly problematic entries (none currently defined)
    # Future filtering criteria can be added here if needed

    return False


def should_append_name_to_lineage(name_to_use: str) -> bool:
    """
    Check if the name_to_use should be appended to the lineage.

    Appends the original name for entries that contain:
    - Numbers (e.g., "Theileria1")
    - .U. patterns (e.g., "Eukaryota.U.family")
    - Underscores (e.g., "Embryophyceae_XX")

    Args:
        name_to_use: The original taxonomic name

    Returns:
        True if the name should be appended to lineage, False otherwise
    """
    if not name_to_use or name_to_use == "NA":
        return False

    # Check for numbers
    if any(char.isdigit() for char in name_to_use):
        return True

    # Check for .U. patterns
    if ".U." in name_to_use:
        return True

    # Check for underscores
    if "_" in name_to_use:
        return True

    return False

