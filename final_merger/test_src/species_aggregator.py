#!/usr/bin/env python3
"""
Species Aggregator Module
==========================

Aggregates matched species data to calculate genome counts, species counts,
and isolate statistics for a census taxon.
"""

import pandas as pd
import logging

logger = logging.getLogger(__name__)


def aggregate_species_matches(matched_species: pd.DataFrame) -> dict:
    """
    Aggregate matched species to calculate summary statistics.
    
    Args:
        matched_species: DataFrame of species that match a census taxon
        
    Returns:
        Dictionary with aggregated statistics:
        - total_genomes: Total genome count
        - isolate_genomes: Isolate genome count
        - uncultured_genomes: Uncultured genome count
        - isolate_percentage: Percentage of genomes that are isolates
        - species_count: Number of unique species
        - matched_species_count: Number of species matched (same as species_count)
    """
    if len(matched_species) == 0:
        return {
            'total_genomes': 0,
            'isolate_genomes': 0,
            'uncultured_genomes': 0,
            'isolate_percentage': 0.0,
            'species_count': 0,
            'matched_species_count': 0
        }
    
    # Sum genome counts
    total_genomes = matched_species['total_genome_count'].sum()
    isolate_genomes = matched_species['isolate_genome_count'].sum()
    uncultured_genomes = matched_species['uncultured_genome_count'].sum()
    
    # Calculate isolate percentage
    isolate_percentage = (isolate_genomes / total_genomes * 100) if total_genomes > 0 else 0.0
    
    # Count unique species
    species_count = len(matched_species)
    
    logger.debug(f"  Aggregated {species_count} species:")
    logger.debug(f"    Total genomes: {total_genomes:,}")
    logger.debug(f"    Isolate genomes: {isolate_genomes:,} ({isolate_percentage:.1f}%)")
    logger.debug(f"    Uncultured genomes: {uncultured_genomes:,}")
    
    return {
        'total_genomes': int(total_genomes),
        'isolate_genomes': int(isolate_genomes),
        'uncultured_genomes': int(uncultured_genomes),
        'isolate_percentage': round(isolate_percentage, 2),
        'species_count': species_count,
        'matched_species_count': species_count
    }


def aggregate_all_census_taxa(
    census_df: pd.DataFrame,
    species_df: pd.DataFrame,
    census_synonym_dict: dict,
    level: str,
    search_function
) -> pd.DataFrame:
    """
    Aggregate species data for all census taxa.

    Args:
        census_df: Census DataFrame
        species_df: Species DataFrame
        census_synonym_dict: Dictionary mapping census names to possible names
        level: Taxonomic level
        search_function: Function to search species for a census taxon

    Returns:
        DataFrame with aggregated results for all census taxa
    """
    logger.info(f"Aggregating species data for {len(census_df)} census {level} taxa...")

    results = []
    matched_count = 0
    total_taxa = len(census_df)

    for idx, census_row in census_df.iterrows():
        # Progress tracking
        if (idx + 1) % 10 == 0 or (idx + 1) == total_taxa:
            logger.info(f"  Processing {idx + 1}/{total_taxa} taxa...")
        census_name = census_row['Name_to_use']

        # Convert taxid to string, handle NaN and clean whitespace
        census_taxid_raw = census_row.get('taxid', None)
        if pd.notna(census_taxid_raw):
            # Clean any whitespace/newlines first
            taxid_str = str(census_taxid_raw).strip().split('\n')[0].split()[0]
            # Convert to int first to remove .0, then to string
            census_taxid = str(int(float(taxid_str)))
        else:
            census_taxid = ''  # Empty string for missing taxids

        census_otus = census_row['otu_count']
        census_size = census_row['size_count']
        
        # Get possible names for this census taxon
        possible_names = census_synonym_dict.get(census_name, {census_name})
        
        # Search for matching species
        matched_species = search_function(
            species_df,
            census_name,
            census_taxid,
            possible_names,
            level
        )
        
        # Aggregate matched species
        aggregated = aggregate_species_matches(matched_species)
        
        # Track match status
        if aggregated['species_count'] > 0:
            matched_count += 1
            match_status = 'matched'
        else:
            match_status = 'census_only'
        
        # Calculate novelty and overrepresentation factors
        # Novelty factor: census_otu_count / ncbi_species_count
        # Higher values = more environmental diversity than genomic representation
        if aggregated['species_count'] > 0:
            novelty_factor = round(census_otus / aggregated['species_count'], 3)
        else:
            novelty_factor = float('inf')

        # Overrepresentation factor: ncbi_species_count / census_otu_count
        # Higher values = database bias toward cultured taxa
        if census_otus > 0:
            overrepresentation_factor = round(aggregated['species_count'] / census_otus, 3)
        else:
            overrepresentation_factor = float('inf')

        # Combine census and NCBI data
        result = {
            level: census_name,
            'census_taxid': census_taxid,
            'census_otu_count': census_otus,
            'census_size_count': census_size,
            'ncbi_genome_count': aggregated['total_genomes'],
            'ncbi_isolate_count': aggregated['isolate_genomes'],
            'ncbi_uncultured_count': aggregated['uncultured_genomes'],
            'ncbi_isolate_percentage': aggregated['isolate_percentage'],
            'ncbi_species_count': aggregated['species_count'],
            'novelty_factor': novelty_factor,
            'overrepresentation_factor': overrepresentation_factor,
            'match_status': match_status
        }
        results.append(result)
    
    # Create output DataFrame
    output_df = pd.DataFrame(results)

    match_rate = (matched_count / len(census_df) * 100) if len(census_df) > 0 else 0
    logger.info(f"  Matched: {matched_count}/{len(census_df)} ({match_rate:.1f}%)")
    logger.info(f"  Census-only: {len(census_df) - matched_count}")

    # Log metrics summary
    matched_only = output_df[output_df['match_status'] == 'matched']
    if len(matched_only) > 0:
        # Filter out inf values for averaging
        finite_novelty = matched_only[matched_only['novelty_factor'] != float('inf')]['novelty_factor']
        finite_overrep = matched_only[matched_only['overrepresentation_factor'] != float('inf')]['overrepresentation_factor']

        if len(finite_novelty) > 0:
            avg_novelty = finite_novelty.mean()
            logger.info(f"  Average novelty factor: {avg_novelty:.2f}")

            # Count high novelty taxa (novelty > 10)
            high_novelty = (finite_novelty > 10).sum()
            logger.info(f"  High novelty taxa (>10): {high_novelty}")

            # Count well-represented taxa (novelty < 2)
            well_represented = (finite_novelty < 2).sum()
            logger.info(f"  Well-represented taxa (<2): {well_represented}")

        if len(finite_overrep) > 0:
            avg_overrep = finite_overrep.mean()
            logger.info(f"  Average overrepresentation factor: {avg_overrep:.2f}")

    return output_df

