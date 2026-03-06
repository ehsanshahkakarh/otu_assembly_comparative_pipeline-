#!/usr/bin/env python3
"""
Resolution Builder for Systematic Resolver

Reads unmapped taxa (families and genera) from taxonkit parser log and builds
systematic resolutions using the known_parents database.
"""

import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, Tuple, List

from .known_parents import KNOWN_PARENTS, get_parent_info


def get_lineage_from_taxonkit(taxid: str, env: dict) -> Tuple[str, str, str]:
    """
    Get lineage information for a taxid using taxonkit.
    
    Args:
        taxid: NCBI taxonomy ID
        env: Environment dictionary for taxonkit
        
    Returns:
        Tuple of (lineage, lineage_ranks, lineage_taxids)
    """
    try:
        # Use taxonkit lineage to get full lineage
        result = subprocess.run(
            ['taxonkit', 'lineage', '-t', '-r'],
            input=taxid,
            capture_output=True,
            text=True,
            env=env,
            check=True
        )
        
        if result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 3:
                lineage = parts[1]
                lineage_taxids = parts[2]
                
                # Get ranks
                rank_result = subprocess.run(
                    ['taxonkit', 'lineage', '-t', '-r', '-R'],
                    input=taxid,
                    capture_output=True,
                    text=True,
                    env=env,
                    check=True
                )
                
                if rank_result.stdout.strip():
                    rank_parts = rank_result.stdout.strip().split('\t')
                    if len(rank_parts) >= 4:
                        lineage_ranks = rank_parts[3]
                        return lineage, lineage_ranks, lineage_taxids
        
        return "", "", ""
        
    except Exception as e:
        logging.warning(f"Failed to get lineage for taxid {taxid}: {e}")
        return "", "", ""


def build_resolution(taxon_name: str, parent_taxid: str, parent_name: str, rank: str, env: dict) -> Dict:
    """
    Build a systematic resolution for a taxon by appending it to parent lineage.

    Args:
        taxon_name: Name of the taxon (family or genus) to resolve
        parent_taxid: NCBI taxid of the parent taxon
        parent_name: Name of the parent taxon
        rank: Taxonomic rank ('family' or 'genus')
        env: Environment dictionary for taxonkit

    Returns:
        Dictionary with resolution information
    """
    # Get parent lineage from taxonkit
    parent_lineage, parent_ranks, parent_taxids = get_lineage_from_taxonkit(parent_taxid, env)

    if not parent_lineage:
        logging.error(f"Failed to get lineage for parent {parent_name} (taxid: {parent_taxid})")
        return None

    # Append taxon name to parent lineage
    resolved_lineage = f"{parent_lineage};{taxon_name}"
    resolved_ranks = f"{parent_ranks};{rank}"
    resolved_taxids = f"{parent_taxids};NA"

    resolution = {
        'taxon_name': taxon_name,
        'rank': rank,
        'parent_taxid': parent_taxid,
        'parent_name': parent_name,
        'lineage': resolved_lineage,
        'lineage_ranks': resolved_ranks,
        'lineage_taxids': resolved_taxids,
        'resolution_method': 'parent_lookup_append'
    }

    logging.info(f"Built resolution for {taxon_name} → {parent_name}")
    return resolution


def build_all_resolutions(env: dict, unmapped_log_path: Path = None) -> Dict[str, Dict]:
    """
    Build resolutions for all taxa (families and genera) in the known_parents database.

    Args:
        env: Environment dictionary for taxonkit
        unmapped_log_path: Optional path to unmapped log to filter taxa

    Returns:
        Dictionary mapping taxon names to resolution information
    """
    logging.info("Building systematic resolutions...")

    # If unmapped log provided, read it to get list of unmapped taxa
    unmapped_taxa = set()
    if unmapped_log_path and unmapped_log_path.exists():
        logging.info(f"Reading unmapped taxa from {unmapped_log_path}")
        with open(unmapped_log_path, 'r') as f:
            for line in f:
                if line.startswith('FAMILY |') or line.startswith('GENUS |'):
                    parts = line.split('|')
                    if len(parts) >= 2:
                        taxon_name = parts[1].strip()
                        unmapped_taxa.add(taxon_name)
        logging.info(f"Found {len(unmapped_taxa)} unmapped taxa in log")

    resolutions = {}

    for taxon_name, (parent_taxid, parent_name, notes, rank) in KNOWN_PARENTS.items():
        # If unmapped log provided, only resolve taxa that are actually unmapped
        if unmapped_taxa and taxon_name not in unmapped_taxa:
            logging.debug(f"Skipping {taxon_name} - not in unmapped list")
            continue

        resolution = build_resolution(taxon_name, parent_taxid, parent_name, rank, env)
        if resolution:
            resolutions[taxon_name] = resolution

    logging.info(f"Built {len(resolutions)} systematic resolutions")
    return resolutions


def save_resolutions(resolutions: Dict[str, Dict], output_file: Path):
    """
    Save resolutions to JSON file.
    
    Args:
        resolutions: Dictionary of resolutions
        output_file: Path to output JSON file
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(resolutions, f, indent=2)
    
    logging.info(f"Saved {len(resolutions)} resolutions to {output_file}")

