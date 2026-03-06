#!/usr/bin/env python3
"""
Level processing module for 18S Taxonkit Parser (Clean NCBI-only version)

Handles taxonomic level aggregation and CSV writing using ONLY taxonkit lookups.
NO systematic resolution integration - pure NCBI taxonomy.
"""

import csv
import logging
from pathlib import Path
from typing import Dict, Tuple
from tqdm import tqdm

from .taxon_validator import should_filter_taxon
from .lineage_processor import append_name_to_lineage, clean_csv_field


def process_taxonomic_level(
    df,
    level: str,
    data_dict: Dict[str, Dict[str, int]]
) -> None:
    """
    Process a single taxonomic level from the dataframe.

    Args:
        df: Pandas DataFrame with cluster data
        level: Taxonomic level name ('division', 'family', or 'genus')
        data_dict: Dictionary to store aggregated data
    """
    logging.info(f"Processing {level} level...")

    # Filter only truly empty/null entries, keep "Unknown" and .U. entries
    level_df = df[~df[level].isna()]
    filtered_count = len(df) - len(level_df)
    logging.info(f"Filtered out {filtered_count} null/empty entries from {level} (keeping Unknown and .U. entries)")

    # Group by taxon and count occurrences (OTU clusters) and sum sizes (sequence counts)
    grouped = level_df.groupby(level).agg({
        'centroid': 'count',  # Count OTU occurrences
        'size': 'sum'  # Sum of sequence counts (cluster sizes)
    }).reset_index()

    # Store in data dictionary with progress bar
    for _, row in tqdm(grouped.iterrows(), total=len(grouped), desc=f"Processing {level} taxa", leave=False):
        taxon = row[level]
        if not should_filter_taxon(taxon):
            data_dict[taxon]['otu_count'] = row['centroid']  # Number of OTU occurrences
            data_dict[taxon]['size_count'] = row['size']  # Sum of sequence counts

    logging.info(f"Processed {len(data_dict)} unique {level} entries")


def write_level_to_csv(
    output_file: Path,
    level_name: str,
    data_dict: Dict[str, Dict[str, int]],
    taxid_dict: Dict[str, str],
    taxid_to_lineage: Dict[str, Tuple[str, str, str]],
    env: dict
) -> int:
    """
    Write taxonomic level data to CSV file using ONLY taxonkit data.
    
    NO systematic resolution - if taxonkit doesn't find it, it stays as NA.

    Args:
        output_file: Path to output CSV file
        level_name: Name of the taxonomic level
        data_dict: Dictionary with aggregated data
        taxid_dict: Dictionary mapping names to taxids
        taxid_to_lineage: Dictionary mapping taxids to lineage tuples
        env: Environment for taxonkit subprocess calls

    Returns:
        Number of filtered entries
    """
    logging.info(f"Writing {level_name} data to {output_file}")

    filtered_count = 0

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Name_to_use', 'taxid', 'otu_count', 'size_count', 'lineage', 'lineage_ranks', 'lineage_taxids']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Sort by OTU count (descending) for better readability
        sorted_items = sorted(data_dict.items(), key=lambda x: x[1]['otu_count'], reverse=True)

        for name, data in tqdm(sorted_items, desc=f"Writing {level_name} data", unit="entry"):
            # Skip if should be filtered
            if should_filter_taxon(name):
                filtered_count += 1
                continue

            # Get taxid and lineage from taxonkit results ONLY
            taxid = taxid_dict.get(name, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, name, taxid, env, taxid_to_lineage
            )

            # Clean ALL fields for CSV output (including taxid!)
            clean_taxid = clean_csv_field(taxid)
            lineage = clean_csv_field(lineage)
            lineage_ranks = clean_csv_field(lineage_ranks)
            lineage_taxids = clean_csv_field(lineage_taxids)

            writer.writerow({
                'Name_to_use': name,
                'taxid': clean_taxid,
                'otu_count': data['otu_count'],
                'size_count': data.get('size_count', 0),
                'lineage': lineage,
                'lineage_ranks': lineage_ranks,
                'lineage_taxids': lineage_taxids
            })

    logging.info(f"Wrote {len(data_dict) - filtered_count} entries to {output_file}")
    if filtered_count > 0:
        logging.info(f"Filtered out {filtered_count} entries")

    return filtered_count

