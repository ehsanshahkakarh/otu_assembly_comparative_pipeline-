#!/usr/bin/env python3
"""
Division Context Adder Module
==============================

For unmapped entries that have NO parent lineage context (just their own name),
this module adds the division information from the raw census data as a 
minimal taxonomic context.

Example:
    Before: WIM80-lineage | WIM80-lineage | original_name | NA
    After:  WIM80-lineage | Evosea;WIM80-lineage | phylum;family | 2605435;NA

This is NOT a full lineage lookup - it just adds the division as a parent
to provide some taxonomic context for otherwise completely isolated entries.
"""

import logging
from typing import Dict, Tuple, Optional

logger = logging.getLogger(__name__)


def load_division_mapping(census_file_path: str) -> Dict[str, str]:
    """
    Load mapping of family/genus names to their division from raw census file.
    
    Args:
        census_file_path: Path to eukcensus_18S.clusters.97.tsv
        
    Returns:
        Dict mapping name -> division
    """
    logger.info(f"Loading division mapping from {census_file_path}")
    
    division_map = {}
    
    with open(census_file_path, 'r') as f:
        # Skip header
        header = f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                division = parts[3]
                family = parts[4]
                genus = parts[5]
                
                # Map family and genus to their division
                if family and family != "":
                    division_map[family] = division
                if genus and genus != "":
                    division_map[genus] = division
    
    logger.info(f"Loaded division info for {len(division_map)} names")
    return division_map


def add_division_context(
    name: str,
    current_lineage: str,
    current_ranks: str,
    current_taxids: str,
    division: str,
    rank_level: str = "family"
) -> Tuple[str, str, str]:
    """
    Add division as parent context if the entry has no parent lineage.
    
    Args:
        name: The unmapped name (e.g., "WIM80-lineage")
        current_lineage: Current lineage string
        current_ranks: Current ranks string
        current_taxids: Current taxids string
        division: Division from raw census (e.g., "Evosea")
        rank_level: Rank of the unmapped name ("family" or "genus")
        
    Returns:
        Tuple of (updated_lineage, updated_ranks, updated_taxids)
    """
    # Check if this entry has NO parent context (lineage is just the name itself)
    if current_lineage == name or current_lineage == "":
        # Add division as parent
        new_lineage = f"{division};{name}"
        new_ranks = f"division;{rank_level}"
        new_taxids = f"NA;NA"  # We don't have taxids for either
        
        logger.info(f"Added division context '{division}' to '{name}'")
        return (new_lineage, new_ranks, new_taxids)
    
    # Entry already has parent context, return unchanged
    return (current_lineage, current_ranks, current_taxids)


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


def process_csv_with_division_context(
    input_csv: str,
    output_csv: str,
    division_map: Dict[str, str]
) -> int:
    """
    Process CSV file and add division context to entries that need it.

    Args:
        input_csv: Path to input CSV file
        output_csv: Path to output CSV file
        division_map: Dict mapping name -> division

    Returns:
        Number of entries updated
    """
    import csv

    logger.info(f"Processing {input_csv}")

    updated_count = 0
    cleaned_count = 0

    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()

        for row in reader:
            name = row['Name_to_use']
            taxid = row['taxid']
            lineage = row['lineage']
            ranks = row['lineage_ranks']
            taxids = row['lineage_taxids']

            # Clean all fields to remove any newlines or problematic characters
            # This fixes any corrupted data from the input file
            taxid_cleaned = clean_csv_field(taxid)
            if taxid_cleaned != taxid:
                logger.warning(f"Cleaned corrupted taxid for '{name}': '{taxid}' -> '{taxid_cleaned}'")
                taxid = taxid_cleaned
                row['taxid'] = taxid
                cleaned_count += 1

            # Clean other fields as well
            row['lineage'] = clean_csv_field(lineage)
            row['lineage_ranks'] = clean_csv_field(ranks)
            row['lineage_taxids'] = clean_csv_field(taxids)

            # Only process unmapped entries (taxid == NA)
            if taxid == 'NA':
                division = division_map.get(name)

                if division and division != "":
                    new_lineage, new_ranks, new_taxids = add_division_context(
                        name, row['lineage'], row['lineage_ranks'], row['lineage_taxids'], division
                    )

                    if new_lineage != row['lineage']:
                        row['lineage'] = new_lineage
                        row['lineage_ranks'] = new_ranks
                        row['lineage_taxids'] = new_taxids
                        updated_count += 1

            writer.writerow(row)

    if cleaned_count > 0:
        logger.warning(f"Cleaned {cleaned_count} corrupted fields from input file")
    logger.info(f"Updated {updated_count} entries with division context")
    logger.info(f"Output saved to {output_csv}")

    return updated_count

