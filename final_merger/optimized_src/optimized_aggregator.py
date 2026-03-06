#!/usr/bin/env python3
"""
Optimized Species Aggregator
=============================

Aggregates matched species with optimized synonym lookup.
"""

import pandas as pd
import logging
from typing import Dict, Set, Callable

logger = logging.getLogger(__name__)


def aggregate_species_optimized(
    census_df: pd.DataFrame,
    species_df: pd.DataFrame,
    all_synonyms: Dict[str, Set[str]],
    level: str,
    search_function: Callable
) -> pd.DataFrame:
    """
    Aggregate species data for all census taxa using pre-built synonym dictionary.
    
    Args:
        census_df: Census DataFrame
        species_df: Pre-processed species DataFrame
        all_synonyms: Pre-built synonym dictionary for ALL census taxa
        level: Taxonomic level
        search_function: Function to search species for a census taxon
    
    Returns:
        DataFrame with merged census and NCBI data
    """
    logger.info(f"Aggregating species data for {len(census_df)} census {level} taxa...")
    
    results = []
    matched_count = 0
    
    for idx, census_row in census_df.iterrows():
        # Progress logging
        if (idx + 1) % 10 == 0:
            logger.info(f"  Processing {idx + 1}/{len(census_df)} taxa...")
        
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
        
        # Get synonyms from pre-built dictionary (FAST - no reload!)
        possible_names = all_synonyms.get(census_name, {census_name})
        
        # Search for matching species
        matched_species = search_function(
            species_df,
            census_name,
            census_taxid,
            possible_names,
            level
        )
        
        # Aggregate genome counts
        if len(matched_species) > 0:
            matched_count += 1
            match_status = 'matched'

            aggregated = {
                'total_genomes': matched_species['total_genome_count'].sum(),
                'isolate_genomes': matched_species['isolate_genome_count'].sum(),
                'uncultured_genomes': matched_species['uncultured_genome_count'].sum(),
                'species_count': len(matched_species)
            }

            # Calculate isolate percentage
            if aggregated['total_genomes'] > 0:
                aggregated['isolate_percentage'] = round(
                    (aggregated['isolate_genomes'] / aggregated['total_genomes']) * 100, 2
                )
            else:
                aggregated['isolate_percentage'] = 0.0
        else:
            match_status = 'census_only'
            aggregated = {
                'total_genomes': 0,
                'isolate_genomes': 0,
                'uncultured_genomes': 0,
                'isolate_percentage': 0.0,
                'species_count': 0
            }
        
        # Calculate novelty and overrepresentation factors
        if aggregated['species_count'] > 0:
            novelty_factor = round(census_otus / aggregated['species_count'], 3)
        else:
            novelty_factor = float('inf')
        
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
    
    logger.info(f"  Processing {len(census_df)}/{len(census_df)} taxa...")
    
    # Create output DataFrame
    output_df = pd.DataFrame(results)

    # Sort by novelty factor (descending - highest novelty first)
    # Separate finite and infinite values
    df_finite = output_df[output_df['novelty_factor'] != float('inf')].copy()
    df_inf = output_df[output_df['novelty_factor'] == float('inf')].copy()

    # Sort finite values by novelty (descending)
    df_finite_sorted = df_finite.sort_values('novelty_factor', ascending=False)

    # Concatenate: finite (sorted) + infinite
    output_df = pd.concat([df_finite_sorted, df_inf])

    match_rate = (matched_count / len(census_df) * 100) if len(census_df) > 0 else 0
    logger.info(f"  Matched: {matched_count}/{len(census_df)} ({match_rate:.1f}%)")
    logger.info(f"  Census-only: {len(census_df) - matched_count}")
    
    # Log metrics summary
    matched_only = output_df[output_df['match_status'] == 'matched']
    if len(matched_only) > 0:
        finite_novelty = matched_only[matched_only['novelty_factor'] != float('inf')]['novelty_factor']
        finite_overrep = matched_only[matched_only['overrepresentation_factor'] != float('inf')]['overrepresentation_factor']
        
        if len(finite_novelty) > 0:
            avg_novelty = finite_novelty.mean()
            logger.info(f"  Average novelty factor: {avg_novelty:.2f}")
            high_novelty = (finite_novelty > 10).sum()
            logger.info(f"  High novelty taxa (>10): {high_novelty}")
            well_represented = (finite_novelty < 2).sum()
            logger.info(f"  Well-represented taxa (<2): {well_represented}")
        
        if len(finite_overrep) > 0:
            avg_overrep = finite_overrep.mean()
            logger.info(f"  Average overrepresentation factor: {avg_overrep:.2f}")
    
    return output_df

