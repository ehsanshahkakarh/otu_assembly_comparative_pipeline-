"""
Lineage Processing and Manipulation for 16S Census Parser
=========================================================

Handles lineage string manipulation and appending original names
to lineage information for special cases.
"""

import subprocess
import logging


def should_append_name_to_lineage(name_to_use):
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


def append_name_to_lineage(lineage, lineage_ranks, lineage_taxids, name_to_use, taxid, env=None, taxid_to_lineage_cache=None):
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

