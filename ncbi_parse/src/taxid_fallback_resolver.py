#!/usr/bin/env python3
"""
Taxid Fallback Resolver Module
================================

Resolves deleted/missing species_taxid values by using the taxid column
from the assembly file as a fallback, then extracting species-level taxid
from the lineage.
"""

import pandas as pd
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Dict, Tuple, Set

logger = logging.getLogger('species_parser')


def resolve_missing_lineages(
    species_df: pd.DataFrame,
    assembly_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Resolve missing lineages by using taxid column as fallback.
    
    Args:
        species_df: DataFrame with species_taxid and lineage columns
        assembly_df: Original assembly DataFrame with both taxid and species_taxid columns
        
    Returns:
        DataFrame with recovered lineages for previously missing entries
    """
    logger.info("Resolving missing lineages using taxid fallback...")
    
    # Identify species with missing lineages
    missing_mask = species_df['lineage'] == ''
    missing_count = missing_mask.sum()
    
    if missing_count == 0:
        logger.info("No missing lineages to resolve!")
        return species_df
    
    logger.info(f"Found {missing_count:,} species with missing lineages")
    missing_species_taxids = species_df[missing_mask]['species_taxid'].tolist()
    
    # Build mapping from species_taxid to alternative taxids from assembly file
    logger.info("Building species_taxid -> taxid mapping from assembly file...")
    taxid_mapping = _build_taxid_mapping(assembly_df, missing_species_taxids)
    
    # Get lineages for the alternative taxids
    logger.info("Retrieving lineages for alternative taxids...")
    alternative_lineages = _get_lineages_for_alternatives(taxid_mapping)
    
    # Extract species-level information from alternative lineages
    logger.info("Extracting species-level taxids from lineages...")
    species_level_data = _extract_species_level_data(alternative_lineages)
    
    # Update the species_df with recovered lineages
    recovered_count = 0
    for species_taxid, species_data in species_level_data.items():
        mask = species_df['species_taxid'] == species_taxid
        if mask.any():
            species_df.loc[mask, 'lineage'] = species_data['lineage']
            species_df.loc[mask, 'lineage_ranks'] = species_data['lineage_ranks']
            species_df.loc[mask, 'lineage_taxids'] = species_data['lineage_taxids']
            species_df.loc[mask, 'species_name'] = species_data['species_name']
            recovered_count += 1
    
    logger.info(f"✅ Recovered lineages for {recovered_count:,}/{missing_count:,} species ({recovered_count/missing_count*100:.1f}%)")
    
    # Log which ones couldn't be recovered
    still_missing = (species_df['lineage'] == '').sum()
    if still_missing > 0:
        logger.warning(f"⚠️  Still missing lineages for {still_missing:,} species")
    
    return species_df


def _build_taxid_mapping(
    assembly_df: pd.DataFrame,
    missing_species_taxids: list
) -> Dict[int, Set[int]]:
    """
    Build mapping from species_taxid to set of alternative taxids.
    
    Args:
        assembly_df: Assembly DataFrame with taxid and species_taxid columns
        missing_species_taxids: List of species_taxids with missing lineages
        
    Returns:
        Dictionary mapping species_taxid -> set of alternative taxids
    """
    mapping = {}
    
    for species_taxid in missing_species_taxids:
        # Find all genomes with this species_taxid
        mask = assembly_df['species_taxid'] == species_taxid
        if mask.any():
            # Get unique taxids from the taxid column
            alternative_taxids = assembly_df.loc[mask, 'taxid'].dropna().unique()
            # Convert to integers and store as set
            mapping[species_taxid] = set(int(x) for x in alternative_taxids if pd.notna(x))
    
    total_alternatives = sum(len(v) for v in mapping.values())
    logger.info(f"Found {total_alternatives:,} alternative taxids for {len(mapping):,} missing species")
    
    return mapping


def _get_lineages_for_alternatives(
    taxid_mapping: Dict[int, Set[int]]
) -> Dict[int, Dict[int, Tuple[str, str, str]]]:
    """
    Get lineages for all alternative taxids.
    
    Args:
        taxid_mapping: Dict mapping species_taxid -> set of alternative taxids
        
    Returns:
        Dict mapping species_taxid -> {alt_taxid: (lineage, ranks, taxids)}
    """
    # Collect all unique alternative taxids
    all_alt_taxids = set()
    for alt_taxids in taxid_mapping.values():
        all_alt_taxids.update(alt_taxids)
    
    logger.info(f"Querying taxonkit for {len(all_alt_taxids):,} alternative taxids...")
    
    # Get lineages using taxonkit
    lineage_results = _query_taxonkit(list(all_alt_taxids))
    
    # Organize by species_taxid
    result = {}
    for species_taxid, alt_taxids in taxid_mapping.items():
        result[species_taxid] = {}
        for alt_taxid in alt_taxids:
            if alt_taxid in lineage_results:
                result[species_taxid][alt_taxid] = lineage_results[alt_taxid]
    
    return result


def _query_taxonkit(taxids: list) -> Dict[int, Tuple[str, str, str]]:
    """
    Query taxonkit for lineages of given taxids.
    
    Args:
        taxids: List of taxids to query
        
    Returns:
        Dict mapping taxid -> (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}
    
    lineage_data = {}
    
    # Create temporary file with taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")
    
    try:
        # Run taxonkit lineage with -R (show rank) and -t (show taxid lineage)
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if result.returncode == 0:
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        taxid = int(parts[0].strip())
                        lineage = parts[1].strip()
                        lineage_taxids = parts[2].strip()
                        lineage_ranks = parts[3].strip()

                        # Only store if lineage is not empty
                        if lineage and lineage != str(taxid):
                            lineage_data[taxid] = (lineage, lineage_ranks, lineage_taxids)
        else:
            logger.warning(f"Taxonkit warning: {result.stderr}")

    except Exception as e:
        logger.error(f"Error running taxonkit lineage: {e}")
    finally:
        Path(temp_filename).unlink(missing_ok=True)

    logger.info(f"Retrieved {len(lineage_data):,} valid lineages from alternative taxids")
    return lineage_data


def _extract_species_level_data(
    alternative_lineages: Dict[int, Dict[int, Tuple[str, str, str]]]
) -> Dict[int, Dict[str, str]]:
    """
    Extract species-level taxid and lineage from alternative taxid lineages.

    For each species_taxid, find the species-level taxid from the lineage
    and create a unified species-level lineage.

    Args:
        alternative_lineages: Dict mapping species_taxid -> {alt_taxid: (lineage, ranks, taxids)}

    Returns:
        Dict mapping species_taxid -> {lineage, lineage_ranks, lineage_taxids, species_name}
    """
    species_data = {}

    for species_taxid, alt_lineages in alternative_lineages.items():
        if not alt_lineages:
            continue

        # Try to find a species-level lineage from any alternative
        species_level_found = False

        for alt_taxid, (lineage, lineage_ranks, lineage_taxids) in alt_lineages.items():
            # Split the lineage components
            lineage_parts = lineage.split(';')
            rank_parts = lineage_ranks.split(';')
            taxid_parts = lineage_taxids.split(';')

            # Find the species rank
            species_idx = None
            for i, rank in enumerate(rank_parts):
                if rank.strip().lower() == 'species':
                    species_idx = i
                    break

            if species_idx is not None:
                # Extract up to species level
                species_lineage = ';'.join(lineage_parts[:species_idx + 1])
                species_ranks = ';'.join(rank_parts[:species_idx + 1])
                species_taxids = ';'.join(taxid_parts[:species_idx + 1])
                species_name = lineage_parts[species_idx].strip()

                species_data[species_taxid] = {
                    'lineage': species_lineage,
                    'lineage_ranks': species_ranks,
                    'lineage_taxids': species_taxids,
                    'species_name': species_name
                }
                species_level_found = True
                break

        # If no species rank found, use the full lineage of the first alternative
        if not species_level_found and alt_lineages:
            first_alt = next(iter(alt_lineages.values()))
            lineage, lineage_ranks, lineage_taxids = first_alt

            # Try to extract last meaningful name
            lineage_parts = lineage.split(';')
            species_name = lineage_parts[-1].strip() if lineage_parts else ''

            species_data[species_taxid] = {
                'lineage': lineage,
                'lineage_ranks': lineage_ranks,
                'lineage_taxids': lineage_taxids,
                'species_name': species_name
            }

    logger.info(f"Extracted species-level data for {len(species_data):,} species")
    return species_data


