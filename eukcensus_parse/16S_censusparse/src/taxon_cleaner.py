"""
Taxonomic Name Cleaning and Validation for 16S Census Parser
=============================================================

Handles cleaning, validation, and extraction of taxonomic names with
support for organelles, uncultured taxa, and rank-appropriate filtering.
"""

import pandas as pd
import re
from .organelle_handler import detect_organelle_type


def extract_meaningful_taxonomic_part(taxon_name):
    """
    Extract meaningful taxonomic information from complex names.

    Handles cases like:
    - uncultured_Alphaproteobacteria_bacterium -> Alphaproteobacteria
    - uncultured_Rickettsia_sp -> Rickettsia
    - marine_metagenome -> None (not taxonomically meaningful)

    Args:
        taxon_name: The taxon name to analyze

    Returns:
        Meaningful taxonomic part or None if no meaningful part found
    """
    if not taxon_name or pd.isna(taxon_name):
        return None

    # Skip purely environmental/technical terms
    environmental_terms = ['metagenome', 'environmental', 'sample', 'clone', 'specimen']
    if any(term in taxon_name.lower() for term in environmental_terms):
        return None

    # For uncultured names, try to extract the taxonomic part
    if 'uncultured' in taxon_name.lower():
        parts = taxon_name.replace('_', ' ').split()

        # Look for meaningful taxonomic terms (not generic terms)
        generic_terms = ['uncultured', 'bacterium', 'organism', 'eukaryote', 'sp', 'species']
        meaningful_parts = []

        for part in parts:
            if part.lower() not in generic_terms and len(part) > 2:
                # Check if this looks like a taxonomic name
                if part[0].isupper() and part[1:].islower():
                    meaningful_parts.append(part)
                elif part.endswith('bacteria') or part.endswith('proteobacteria'):
                    meaningful_parts.append(part)

        if meaningful_parts:
            # Return the most specific meaningful part
            return meaningful_parts[-1] if len(meaningful_parts) == 1 else ' '.join(meaningful_parts[:2])

    return None


def strip_trailing_numbers(taxon_name):
    """
    Remove trailing numbers from taxon names (e.g., "Theileria1" -> "Theileria").

    Args:
        taxon_name: The taxon name to process

    Returns:
        The taxon name with trailing numbers removed
    """
    if not taxon_name:
        return taxon_name
        
    # Handle taxa with numbers at the end (e.g., "Cryptosporidium15", "Eimeria1", "Plasmodium1")
    match = re.match(r'^(.+?)(\d+)$', taxon_name.strip())
    if match:
        base_name = match.group(1).rstrip()  # Remove any trailing spaces
        return base_name

    # If no trailing numbers found, return as is
    return taxon_name


def clean_organelle_taxon_name(taxon_name):
    """
    Enhanced cleaning for organelle-containing taxon names with improved Candidatus handling.

    Handles cases like:
    - Vitis_vinifera:plas.Chloroplast -> Vitis vinifera
    - uncultured_bacterium.Mitochondria -> uncultured bacterium
    - Genus_species.Plastid -> Genus species
    - Candidatus names -> preserve (no stripping needed with updated NCBI taxonomy)

    Args:
        taxon_name: The taxon name to clean

    Returns:
        Cleaned taxon name suitable for taxonomic lookup
    """
    if pd.isna(taxon_name) or taxon_name == "Unknown":
        return taxon_name

    # Detect and handle organelles
    is_organelle, organelle_type, base_name = detect_organelle_type(taxon_name)

    if is_organelle:
        # Quietly handle organelle detection
        taxon_name = base_name

    # Candidatus taxa are now preserved in NCBI taxonomy - no stripping needed
    # The NCBI taxonomic mapping scripts have been updated to handle Candidatus taxa properly

    # Handle .U. patterns (e.g., "Gammaproteobacteria.U.family" -> "Gammaproteobacteria")
    if '.U.' in taxon_name:
        # Extract the base taxonomic name before the .U. pattern
        base_name = taxon_name.split('.U.')[0]
        return strip_trailing_numbers(base_name)

    # Replace underscores with two spaces for species-level names and EukCensus patterns
    if '_' in taxon_name:
        # Check if this looks like a binomial species name
        parts = taxon_name.split('_')
        if len(parts) == 2 and not any(char.isdigit() for char in parts[1]):
            # Likely a species name like Genus_species - use single space for species
            cleaned = taxon_name.replace('_', ' ')
            return strip_trailing_numbers(cleaned)
        else:
            # Handle other underscore cases (like _XX patterns) - use two spaces
            cleaned = taxon_name.replace('_', '  ')
            return strip_trailing_numbers(cleaned)

    # Handle trailing numbers
    return strip_trailing_numbers(taxon_name)


def extract_genus_from_species(species_name):
    """
    Extract genus name from a species name.

    Handles special cases:
    - "Candidatus Genus" -> "Candidatus Genus" (keep both words)
    - "candidate division Name" -> "candidate division Name" (keep both words)
    - "Genus species" -> "Genus" (normal case)

    Args:
        species_name: Full species name (e.g., "Vitis vinifera" or "Candidatus Edwardsbacteria")

    Returns:
        Genus name (e.g., "Vitis" or "Candidatus Edwardsbacteria")
    """
    if not species_name or pd.isna(species_name):
        return species_name

    parts = species_name.strip().split()

    # Special case: "Candidatus" is a prefix for uncultivated prokaryotes
    # Keep both "Candidatus" and the following genus name
    if len(parts) >= 2 and parts[0].lower() == 'candidatus':
        return f"{parts[0]} {parts[1]}"

    # Special case: "candidate division" is a prefix for environmental clades
    # Keep both "candidate division" and the following name
    if len(parts) >= 3 and parts[0].lower() == 'candidate' and parts[1].lower() == 'division':
        return f"{parts[0]} {parts[1]} {parts[2]}"

    # Normal case: return first word (genus)
    if len(parts) >= 1:
        return parts[0]

    return species_name


def clean_taxon_name(taxon_name):
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


def should_filter_taxon(taxon_name):
    """
    Check if a taxon name should be filtered out.

    UPDATED: Removed .U. filtering to preserve unidentified taxa for downstream visualization.
    Now only filters out truly empty/null entries.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    if not taxon_name or pd.isna(taxon_name):
        return True

    # REMOVED: .U. pattern filtering - these entries are now preserved for visualization
    # REMOVED: unidentified pattern filtering - these entries are now preserved

    # Only filter out completely empty strings after stripping whitespace
    if isinstance(taxon_name, str) and taxon_name.strip() == "":
        return True

    return False


def validate_rank_appropriateness(taxon_name, target_rank, lineage_info):
    """
    Validate if a taxon entry is appropriate for the target taxonomic rank.

    Args:
        taxon_name: Original taxon name
        target_rank: Target rank ('genus', 'family', 'phylum')
        lineage_info: Tuple of (lineage, lineage_ranks, lineage_taxids)

    Returns:
        True if appropriate, False if should be filtered
    """
    import logging

    if not lineage_info or not lineage_info[1]:  # No rank information
        return True  # Let it through, will be handled by name-based logic

    lineage, ranks, taxids = lineage_info
    rank_list = ranks.split(';') if ranks else []

    # Check if the lineage contains our target rank
    if target_rank not in rank_list:
        logging.info(f"⚠️ Taxon '{taxon_name}' lacks {target_rank} rank in lineage, filtering out")
        return False

    # For genus-level parsing, filter out entries that are clearly family or higher
    if target_rank == 'genus':
        if 'species' in rank_list:
            # Species-level entry is OK for genus parsing (we extract genus)
            return True
        elif 'genus' in rank_list:
            # Genus-level entry is perfect
            return True
        else:
            # Higher-level entry (family, order, etc.) - might not be appropriate
            logging.info(f"⚠️ Taxon '{taxon_name}' appears to be higher than genus level, filtering out")
            return False

    # For family-level parsing
    elif target_rank == 'family':
        if 'species' in rank_list or 'genus' in rank_list:
            # Too specific for family-level parsing
            logging.info(f"⚠️ Taxon '{taxon_name}' is too specific for family-level parsing, filtering out")
            return False
        elif 'family' in rank_list:
            return True
        else:
            # Higher level might be OK
            return True

    # For phylum-level parsing
    elif target_rank == 'phylum':
        if any(rank in rank_list for rank in ['species', 'genus', 'family']):
            # Too specific for phylum-level parsing
            logging.info(f"⚠️ Taxon '{taxon_name}' is too specific for phylum-level parsing, filtering out")
            return False
        elif 'phylum' in rank_list:
            return True
        else:
            return True

    return True


def extract_appropriate_rank_name(taxon_name, target_rank, lineage_info=None):
    """
    Extract the appropriate taxonomic rank name from a taxon entry with DISABLED organelle recovery.

    For example, if parsing at genus level but entry is species-level,
    extract the genus portion.

    NOTE: Organelle recovery is now DISABLED here to prevent performance bottlenecks.
    Organellar sequences are handled by the vectorized approach in get_taxids_using_taxonkit().

    Args:
        taxon_name: Original taxon name
        target_rank: Target taxonomic rank ('genus', 'family', 'phylum')
        lineage_info: Optional lineage information tuple (lineage, ranks, taxids)

    Returns:
        Appropriate rank name or None if should be filtered
    """
    # PERFORMANCE FIX: Disable individual organelle recovery here
    # Organellar sequences will be handled by vectorized approach later
    # This prevents thousands of individual subprocess calls during initial processing

    cleaned_name = clean_organelle_taxon_name(taxon_name)

    # If we have lineage information, use it to validate rank appropriateness
    if lineage_info and lineage_info[1]:  # lineage_ranks exists
        lineage, ranks, taxids = lineage_info
        rank_list = ranks.split(';') if ranks else []

        # Check if the entry matches our target rank
        if target_rank in rank_list:
            # Find the position of our target rank
            try:
                target_idx = rank_list.index(target_rank)
                lineage_parts = lineage.split(';') if lineage else []
                if target_idx < len(lineage_parts):
                    return lineage_parts[target_idx]
            except (ValueError, IndexError):
                pass

    # Fallback to name-based extraction
    if target_rank == 'genus':
        # For genus level, extract genus from species names
        if ' ' in cleaned_name and len(cleaned_name.split()) >= 2:
            # Likely a species name, extract genus
            return extract_genus_from_species(cleaned_name)
        else:
            # Already genus-level or higher
            return cleaned_name

    elif target_rank == 'family':
        # For family level, we generally keep the name as-is unless it's clearly species
        if ' ' in cleaned_name and len(cleaned_name.split()) >= 2:
            # This is a species name, which is inappropriate for family-level parsing
            return None
        return cleaned_name

    elif target_rank == 'phylum':
        # For phylum level, filter out species and genus names
        if ' ' in cleaned_name:
            # This is likely a species name
            return None
        # Could still be genus or family, but we'll let taxonkit determine appropriateness
        return cleaned_name

    return cleaned_name

