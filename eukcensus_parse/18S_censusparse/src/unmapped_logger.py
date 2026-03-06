#!/usr/bin/env python3
"""
Unmapped logger module for 18S Taxonkit Parser (Clean NCBI-only version)

Creates a simple log of families that failed taxonkit lookup.
NO resolution filtering - just reports what taxonkit couldn't find.
"""

import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Tuple

from .taxon_cleaner import clean_taxon_name


def create_unmapped_log(
    division_data: Dict[str, Dict[str, int]],
    family_data: Dict[str, Dict[str, int]],
    genus_data: Dict[str, Dict[str, int]],
    division_to_taxid: Dict[str, str],
    family_to_taxid: Dict[str, str],
    genus_to_taxid: Dict[str, str],
    taxid_to_lineage: Dict[str, Tuple[str, str, str]],
    log_dir: Path,
    output_prefix: str
) -> Path:
    """
    Create a simple log of all taxonomic names that failed taxonkit lookup.
    
    This is the PRE-RESOLUTION log - it shows what NCBI taxonomy doesn't have.

    Args:
        division_data, family_data, genus_data: Data dictionaries for each rank
        division_to_taxid, family_to_taxid, genus_to_taxid: Taxid mapping dictionaries
        taxid_to_lineage: Lineage information dictionary
        log_dir: Directory to store log files
        output_prefix: Prefix for output files

    Returns:
        Path to the created log file
    """
    log_file = log_dir / f"{output_prefix}_unmapped_from_taxonkit.log"
    logging.info(f"Creating unmapped log (taxonkit failures only): {log_file}")

    with open(log_file, 'w') as f:
        f.write("# Unmapped Names Log - 18S Taxonkit Parser\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# This log contains all taxonomic names that failed taxonkit lookup\n")
        f.write("# These names are NOT found in NCBI taxonomy database\n")
        f.write("# Format: Rank | Original_Name | Cleaned_Name | OTU_Count | Size_Count | Taxid | Reason\n\n")

        # Summary statistics
        f.write("=== SUMMARY STATISTICS ===\n\n")

        # Calculate statistics for each rank
        # NOTE: all_taxid_results is shared across all ranks, so we need to filter by actual names in each rank
        rank_stats = {}
        for rank_name, data_dict, taxid_dict in [
            ('division', division_data, division_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            total = len(data_dict)
            # Count mapped taxa: those that have a taxid AND that taxid has a lineage
            mapped = 0
            for name in data_dict.keys():
                taxid = taxid_dict.get(name, "NA")
                if taxid != "NA" and taxid in taxid_to_lineage:
                    mapped += 1

            unmapped = total - mapped

            rank_stats[rank_name] = {
                'total': total,
                'mapped': mapped,
                'unmapped': unmapped
            }

            f.write(f"{rank_name.capitalize()}:\n")
            f.write(f"  Total unique {rank_name}s: {total}\n")
            f.write(f"  Found in NCBI: {mapped}\n")
            f.write(f"  NOT in NCBI: {unmapped}\n\n")

        # Overall statistics
        total_all = sum(stats['total'] for stats in rank_stats.values())
        mapped_all = sum(stats['mapped'] for stats in rank_stats.values())
        unmapped_all = total_all - mapped_all

        f.write(f"Overall:\n")
        f.write(f"  Total unique taxa: {total_all}\n")
        f.write(f"  Found in NCBI: {mapped_all}\n")
        f.write(f"  NOT in NCBI: {unmapped_all}\n\n")

        # Detailed unmapped entries by rank
        for rank_name, data_dict, taxid_dict in [
            ('division', division_data, division_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            f.write(f"=== {rank_name.upper()} LEVEL UNMAPPED NAMES ===\n")

            unmapped_entries = []
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    # Determine reason for failure
                    if taxid == "NA":
                        reason = "NO_TAXID_FOUND"
                    else:
                        reason = "TAXID_NO_LINEAGE"

                    unmapped_entries.append({
                        'original_name': orig_name,
                        'cleaned_name': clean_taxon_name(orig_name),
                        'otu_count': data['otu_count'],
                        'size_count': data.get('size_count', 0),
                        'taxid': taxid,
                        'reason': reason
                    })

            f.write(f"Total unmapped {rank_name} entries: {len(unmapped_entries)}\n\n")

            # Sort by OTU count (descending)
            unmapped_entries.sort(key=lambda x: x['otu_count'], reverse=True)

            for entry in unmapped_entries:
                f.write(f"{rank_name.upper()} | {entry['original_name']} | {entry['cleaned_name']} | {entry['otu_count']} | {entry['size_count']} | {entry['taxid']} | {entry['reason']}\n")

            f.write(f"\n")

        f.write("\n=== NEXT STEPS ===\n")
        f.write("These unmapped families can be resolved using the systematic_resolver:\n")
        f.write("  python systematic_resolver/run_systematic_resolver.py\n")
        f.write("\nThe systematic resolver will:\n")
        f.write("1. Read this unmapped log\n")
        f.write("2. Apply known parent taxid mappings\n")
        f.write("3. Generate lineages by appending family names to parent lineages\n")
        f.write("4. Create systematic_resolutions.json for merging back into CSV files\n")

    logging.info(f"Unmapped log written to {log_file}")
    logging.info(f"Found {rank_stats['family']['unmapped']} unmapped families for systematic resolution")
    
    return log_file

