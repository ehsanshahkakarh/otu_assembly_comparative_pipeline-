#!/usr/bin/env python3
"""
Taxon name cleaning module for 18S Census Parser

Handles name cleaning, pattern extraction, and taxonomic mapping for eukaryotic taxa.
"""

import re
from typing import Optional, Dict


def clean_taxon_name(taxon_name: str) -> str:
    """
    Clean a taxon name by replacing underscores with two spaces and removing trailing numbers.

    This handles EukCensus patterns like "_XX" by removing everything after the underscore
    and adding two spaces to help with taxonkit matching.

    Args:
        taxon_name: The taxon name to clean

    Returns:
        The cleaned taxon name
    """
    # Replace underscores with two spaces and remove trailing numbers
    cleaned = taxon_name.replace("_", "  ")
    return strip_trailing_numbers(cleaned)


def strip_trailing_numbers(taxon_name: str) -> str:
    """
    Remove trailing numbers from taxon names (e.g., "Theileria1" -> "Theileria").

    Args:
        taxon_name: The taxon name to process

    Returns:
        The taxon name with trailing numbers removed
    """
    # Handle taxa with numbers at the end (e.g., "Cryptosporidium15", "Eimeria1", "Plasmodium1")
    # This regex matches any word characters followed by one or more digits at the end
    match = re.match(r'^(.+?)(\d+)$', taxon_name.strip())
    if match:
        base_name = match.group(1).rstrip()  # Remove any trailing spaces
        return base_name

    # If no trailing numbers found, return as is
    return taxon_name


def extract_taxa_from_hyphenated(taxon_name: str) -> Optional[str]:
    """
    Extract taxa from hyphenated names based on different patterns.

    Args:
        taxon_name: The hyphenated taxon name

    Returns:
        The extracted taxon name or None if no valid pattern found
    """
    if "-" not in taxon_name:
        return None

    # Pattern 1: [taxa]-lineage → extract first part
    if taxon_name.endswith("-lineage") or "_X" in taxon_name and "-lineage" in taxon_name:
        base_name = taxon_name.split("-lineage")[0]
        # Remove any _X suffix
        if "_" in base_name:
            base_name = base_name.split("_")[0]
        return base_name

    # Pattern 2: [taxa]-Group → extract first part
    if "-Group" in taxon_name:
        return taxon_name.split("-Group")[0]

    # Pattern 3: X-[taxa]_XX → extract middle part (existing logic)
    parts = taxon_name.split("-")
    if len(parts) >= 2:
        second_part = parts[1]
        # Remove any _XXX, _XX, _X suffix
        if "_" in second_part:
            clean_part = second_part.split("_")[0]
        else:
            clean_part = second_part

        # Only return if it looks like a valid taxonomic name (not a research clade)
        if not any(keyword in clean_part.lower() for keyword in ["clade", "group", "relatives"]):
            return clean_part

    return None


def extract_genus(taxon_name: str) -> str:
    """
    Extract the genus part from a taxon name, handling trailing numbers.

    Handles special cases:
    - "Candidatus Genus" -> "Candidatus Genus" (keep both words)
    - "candidate division Name" -> "candidate division Name" (keep both words)
    - "Genus species" -> "Genus" (normal case)
    - "Genus_species" -> "Genus" (underscore-separated)

    Args:
        taxon_name: The taxon name to extract genus from

    Returns:
        The genus part of the taxon name
    """
    # First clean the name to remove organelle information
    if "." in taxon_name:
        # For names like "Genus_species.Mitochondria", get the part before the dot
        parts = taxon_name.split(".")
        name_part = parts[0]

        # If the name part contains an underscore, get the first part (genus)
        if "_" in name_part:
            genus = name_part.split("_")[0]
        else:
            genus = name_part

        # Strip trailing numbers from the genus
        return strip_trailing_numbers(genus)

    # For names like "Genus_species", get the first part
    if "_" in taxon_name:
        genus = taxon_name.split("_")[0]
        return strip_trailing_numbers(genus)

    # For names with spaces (after clean_taxon_name processing)
    if " " in taxon_name:
        parts = taxon_name.strip().split()

        # Special case: "Candidatus" is a prefix for uncultivated prokaryotes
        # Keep both "Candidatus" and the following genus name
        if len(parts) >= 2 and parts[0].lower() == 'candidatus':
            return f"{parts[0]} {parts[1]}"

        # Special case: "candidate division" is a prefix for environmental clades
        # Keep both "candidate division" and the following name
        if len(parts) >= 3 and parts[0].lower() == 'candidate' and parts[1].lower() == 'division':
            return f"{parts[0]} {parts[1]} {parts[2]}"

        # Normal case: return first word (genus)
        return strip_trailing_numbers(parts[0])

    # For names like "Genus" or "Genus1", return the genus part (without numbers)
    return strip_trailing_numbers(taxon_name)


def get_taxonomic_mapping() -> Dict[str, Optional[str]]:
    """
    Get mapping from outdated/informal taxonomic names to modern valid names.

    Returns:
        Dictionary mapping old names to new names (None means unmappable)
    """
    return {
        # Outdated taxonomic names to modern equivalents
        "Maxillopoda": "Copepoda",
        "Embryophyceae_XX": "Embryophyta",
        "Chytridiomycetaceae": "Chytridiomycetes",
        "Ophryoglenida": "Ophryoglenidae",

        # Informal groups to formal taxonomy
        "Blastocystis-Group": "Blastocystis",
        "Flamella-lineage": "Flamella",
        "Endostelium-lineage": "Endostelium",
        "Protaspa-lineage": "Protaspa",
        "Rhogostoma-lineage": "Rhogostoma",

        # Research clades to broader valid groups where possible
        "Neobodonidae": "Bodonidae",  # Neobodonidae is often considered part of Bodonidae
        "Vermamoebidae": "Amoebidae",  # Map to broader amoeba family
        "Tholoniidae": "Tholonia",  # Map to genus level
        "Nolandellidae": "Nolandella",  # Map to genus level
        "Paradinidae": "Paradinium",  # Map to genus level
        "Skeletonemaceae": "Skeletonema",  # Map to genus level

        # Remove suffixes for research annotations
        "Dino-Group-II_X": "Dinoflagellata",  # Map to broader group
        "Endomyxa-Ascetosporea_XX": "Endomyxa",
        "Eupetalomonads_X": "Eupetalomonas",
        "Filoretidae_X": "Filoreta",
        "Novel-clade-10_X": None,  # Keep as unmappable

        # Genus level mappings
        "Paradinida_XX": "Paradinium",
        "Cryptosporidium15": "Cryptosporidium",
        "Filosa-Thecofilosea_XXX": "Thecofilosea",
        "Craniata_XXX": "Craniata",
        "Mataza-lineage_X": "Mataza",
        "Novel-Clade-4_X": None,  # Keep as unmappable
        "Corallicolla": "Corallicola",  # Fix spelling
        "OLIGO5_XX": None,  # Environmental clade, unmappable

        # Environmental clades that should remain unmapped (return None)
        "MAST-12": None,
        "MAST-3E": None,
        "LKM74-lineage": None,
        "NC12A-lineage": None,
        "NC12B-lineage": None,
        "NPK2-lineage": None,
        "WIM80-lineage": None,
        "AND16-lineage": None,
        "LOS7N/I-lineage": None,
        "Mariager-Fjord-lineage": None,
        "Mb5C-lineage_X": None,
        "CCW10-lineage_X": None,
        "OLIGO2": None,
        "Nucleohelea": None,  # Environmental group
    }


def apply_taxonomic_mapping(taxon_name: str, mapping_dict: Dict[str, Optional[str]]) -> str:
    """
    Apply taxonomic mapping to convert outdated names to modern ones.

    Args:
        taxon_name: Original taxon name
        mapping_dict: Dictionary of mappings

    Returns:
        Mapped taxon name or original if no mapping exists
    """
    if taxon_name in mapping_dict:
        mapped_name = mapping_dict[taxon_name]
        if mapped_name is not None:
            return mapped_name
        else:
            return taxon_name
    return taxon_name

