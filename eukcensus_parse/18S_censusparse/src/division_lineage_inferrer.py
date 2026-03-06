#!/usr/bin/env python3
"""
Division Lineage Inferrer Module
==================================

For unmapped family/genus names that have NA taxid, this module attempts to
infer a partial lineage from the 'division' field in the raw census data.

The division field often contains higher-level taxonomic information (e.g., 
"Evosea", "Dinophyceae", "Stramenopiles") that can be used to place unmapped
entries in a broader taxonomic context.

Strategy:
1. For entries with taxid=NA but non-empty division field
2. Look up the division name in NCBI taxonomy using taxonkit
3. If found, use that lineage as a partial context
4. Append the unmapped name to the division lineage
"""

import logging
import subprocess
import os
from typing import Dict, Tuple, Optional

logger = logging.getLogger(__name__)


def infer_lineage_from_division(
    division_name: str,
    unmapped_name: str,
    rank_level: str = "family",
    data_dir: str = "~/.taxonkit"
) -> Optional[Tuple[str, str, str]]:
    """
    Infer a partial lineage for an unmapped name using its division.

    Args:
        division_name: The division field from raw census data (e.g., "Evosea")
        unmapped_name: The unmapped family/genus name
        rank_level: The rank of the unmapped name ("family" or "genus")
        data_dir: Path to taxonkit data directory

    Returns:
        Tuple of (lineage, lineage_ranks, lineage_taxids) or None if division not found
    """
    if not division_name or division_name == "NA":
        return None

    # Try to get lineage for the division name
    division_lineage = _get_division_lineage(division_name, data_dir)

    if not division_lineage:
        return None

    lineage, lineage_ranks, lineage_taxids = division_lineage

    # Append the unmapped name to the division lineage
    full_lineage = f"{lineage};{unmapped_name}"
    full_ranks = f"{lineage_ranks};{rank_level}"
    full_taxids = f"{lineage_taxids};NA"

    logger.info(f"Inferred lineage for '{unmapped_name}' from division '{division_name}'")

    return (full_lineage, full_ranks, full_taxids)


def _get_division_lineage(division_name: str, data_dir: str = "~/.taxonkit") -> Optional[Tuple[str, str, str]]:
    """
    Get NCBI lineage for a division name using taxonkit.

    Args:
        division_name: Name to look up (e.g., "Evosea", "Dinophyceae")
        data_dir: Path to taxonkit data directory

    Returns:
        Tuple of (lineage, lineage_ranks, lineage_taxids) or None if not found
    """
    try:
        # Use taxonkit name2taxid to get taxid from name
        result = subprocess.run(
            ["taxonkit", "name2taxid", "--data-dir", data_dir],
            input=division_name + "\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10
        )

        if result.returncode != 0 or not result.stdout.strip():
            logger.debug(f"No taxid found for division '{division_name}'")
            return None

        # Parse the output: "name\ttaxid"
        parts = result.stdout.strip().split('\t')
        if len(parts) < 2:
            return None

        taxid = parts[1].strip()

        if not taxid or taxid == "":
            return None

        # Now get the full lineage for this taxid
        lineage_result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", "--data-dir", data_dir],
            input=taxid + "\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10
        )

        if lineage_result.returncode != 0 or not lineage_result.stdout.strip():
            return None

        # Parse lineage output: "taxid\tlineage\tlineage_taxids\tlineage_ranks"
        lineage_parts = lineage_result.stdout.strip().split('\t')
        if len(lineage_parts) < 4:
            return None

        lineage = lineage_parts[1].strip()
        lineage_taxids = lineage_parts[2].strip()
        lineage_ranks = lineage_parts[3].strip()

        if lineage and lineage != taxid:
            logger.debug(f"Found lineage for division '{division_name}' (taxid {taxid})")
            return (lineage, lineage_ranks, lineage_taxids)

        return None

    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout while looking up division '{division_name}'")
        return None
    except Exception as e:
        logger.error(f"Error looking up division '{division_name}': {e}")
        return None


def batch_infer_lineages(
    unmapped_entries: Dict[str, Dict],
    raw_census_data: Dict[str, str]
) -> Dict[str, Tuple[str, str, str]]:
    """
    Batch process unmapped entries to infer lineages from divisions.
    
    Args:
        unmapped_entries: Dict mapping name -> {metadata}
        raw_census_data: Dict mapping name -> division from raw census file
        
    Returns:
        Dict mapping name -> (lineage, lineage_ranks, lineage_taxids)
    """
    inferred_lineages = {}
    
    for name, metadata in unmapped_entries.items():
        division = raw_census_data.get(name)
        
        if not division or division == "NA":
            continue
        
        rank_level = metadata.get('rank', 'family')
        
        lineage_data = infer_lineage_from_division(division, name, rank_level)
        
        if lineage_data:
            inferred_lineages[name] = lineage_data
    
    logger.info(f"Inferred lineages for {len(inferred_lineages)}/{len(unmapped_entries)} unmapped entries")
    
    return inferred_lineages

