#!/usr/bin/env python3
"""
Lineage processing module for 18S Census Parser

Handles lineage manipulation and CSV field cleaning.
"""

import subprocess
import logging
from typing import Tuple, Optional, Dict

from .taxon_validator import should_append_name_to_lineage


def clean_csv_field(field) -> str:
    """
    Clean a field for CSV output by removing newlines and other problematic characters.

    Args:
        field: The field to clean

    Returns:
        Cleaned field safe for CSV
    """
    if field is None:
        return ""

    # Convert to string and remove newlines, carriage returns, and tabs
    cleaned = str(field).replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')

    # Remove extra whitespace
    cleaned = ' '.join(cleaned.split())

    return cleaned


def append_name_to_lineage(
    lineage: str,
    lineage_ranks: str,
    lineage_taxids: str,
    name_to_use: str,
    taxid: str,
    env: Optional[dict] = None,
    taxid_to_lineage_cache: Optional[Dict[str, Tuple[str, str, str]]] = None
) -> Tuple[str, str, str]:
    """
    Append the name_to_use to the end of lineage components if it meets criteria.
    If we have a taxid but no lineage, try to retrieve the lineage first.

    Args:
        lineage: Original lineage string
        lineage_ranks: Original lineage ranks string
        lineage_taxids: Original lineage taxids string
        name_to_use: The original taxonomic name
        taxid: The taxid for this entry
        env: Environment for taxonkit subprocess calls
        taxid_to_lineage_cache: Dictionary to update with newly retrieved lineages

    Returns:
        Tuple of (updated_lineage, updated_ranks, updated_taxids)
    """
    if not should_append_name_to_lineage(name_to_use):
        return lineage, lineage_ranks, lineage_taxids

    # If we have a taxid but no lineage, try to get the lineage
    if taxid != "NA" and not lineage and env:
        try:
            result = subprocess.run(
                ["taxonkit", "lineage", "-R", "-t"],
                input=taxid,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )

            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 4:
                    lineage = parts[1]
                    lineage_taxids = parts[2]
                    lineage_ranks = parts[3]

                    # Update the cache so the log reflects the newly retrieved lineage
                    if taxid_to_lineage_cache is not None:
                        taxid_to_lineage_cache[taxid] = (lineage, lineage_ranks, lineage_taxids)
                        logging.info(f"Retrieved missing lineage for taxid {taxid}: {lineage}")

        except Exception as e:
            logging.warning(f"Failed to retrieve lineage for taxid {taxid}: {e}")

    # Append the name_to_use to lineage components
    if lineage:
        updated_lineage = f"{lineage};{name_to_use}"
    else:
        updated_lineage = name_to_use

    if lineage_ranks:
        updated_ranks = f"{lineage_ranks};original_name"
    else:
        updated_ranks = "original_name"

    if lineage_taxids:
        updated_taxids = f"{lineage_taxids};{taxid}" if taxid != "NA" else f"{lineage_taxids};NA"
    else:
        updated_taxids = taxid if taxid != "NA" else "NA"

    return updated_lineage, updated_ranks, updated_taxids

