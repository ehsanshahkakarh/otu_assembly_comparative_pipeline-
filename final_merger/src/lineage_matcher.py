#!/usr/bin/env python3
"""
Lineage Matcher Module
======================

Performs hierarchical matching between census and NCBI species-level data using:
1. Direct taxid matching (exact taxid match)
2. Hierarchical taxid matching (census taxid found in NCBI lineage_taxids)

Then aggregates genome counts and species counts after matching.
"""

import pandas as pd
import re
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def match_taxa_by_lineage(
    census_df: pd.DataFrame,
    ncbi_species_df: pd.DataFrame,
    level: str,
    census_level: str
) -> pd.DataFrame:
    """
    Match census taxa to NCBI species using hierarchical taxid matching,
    then aggregate genome counts and species counts.

    Matching strategy:
    1. Direct taxid match: census taxid == NCBI species_taxid
    2. Hierarchical taxid match: census taxid found in NCBI lineage_taxids string
       Pattern: ;taxid; or ^taxid; or ;taxid$ or ^taxid$

    Args:
        census_df: Census DataFrame with Name_to_use, taxid, otu_count, size_count, lineage, lineage_taxids
        ncbi_species_df: NCBI species DataFrame with species_taxid, total_genome_count,
                        isolate_genome_count, uncultured_genome_count, lineage, lineage_taxids
        level: NCBI taxonomic level (phylum, family, genus) - used for logging only
        census_level: Census level name (division, family, genus)

    Returns:
        DataFrame with matched results containing census and aggregated NCBI data
    """
    logger.info(f"Matching {len(census_df):,} census {census_level} entries to NCBI species...")
    logger.info(f"  NCBI species pool: {len(ncbi_species_df):,} species")
    logger.info(f"  Using hierarchical taxid matching: direct taxid → hierarchical taxid")

    # Prepare NCBI data (fill NaN with empty string)
    ncbi_species_df = ncbi_species_df.copy()
    ncbi_species_df['lineage'] = ncbi_species_df['lineage'].fillna('')
    ncbi_species_df['lineage_taxids'] = ncbi_species_df['lineage_taxids'].fillna('')

    # Calculate totals for percentages
    total_census_otus = census_df['otu_count'].sum()
    total_census_size = census_df['size_count'].sum()
    total_ncbi_genomes = ncbi_species_df['total_genome_count'].sum()
    total_ncbi_species = len(ncbi_species_df)  # Each row is a species

    logger.info(f"  Census totals: {total_census_otus:,} OTUs, {total_census_size:,} size")
    logger.info(f"  NCBI totals: {total_ncbi_genomes:,} genomes, {total_ncbi_species:,} species")

    # Process each census entry
    results = []
    matched_count = 0
    match_type_counts = {'direct_taxid': 0, 'hierarchical_taxid': 0}

    for _, census_row in tqdm(census_df.iterrows(), total=len(census_df),
                               desc=f"  Matching {census_level}", unit="taxa", leave=False):
        taxon_name = census_row['Name_to_use']
        census_taxid = census_row.get('taxid', None)
        census_otus = census_row['otu_count']
        census_size = census_row['size_count']

        # Calculate census percentages
        otu_percentage = (census_otus / total_census_otus * 100) if total_census_otus > 0 else 0
        size_percentage = (census_size / total_census_size * 100) if total_census_size > 0 else 0

        # Initialize match containers
        direct_taxid_matched = pd.DataFrame()
        hierarchical_taxid_matched = pd.DataFrame()

        # PRIORITY 1: Direct taxid matching (census taxid == NCBI species_taxid)
        if census_taxid and pd.notna(census_taxid):
            try:
                # Clean and convert taxid
                if isinstance(census_taxid, str):
                    census_taxid_clean = census_taxid.split('\n')[0].strip()
                    census_taxid_str = census_taxid_clean
                else:
                    census_taxid_str = str(int(census_taxid))

                # Match against species_taxid
                direct_taxid_mask = ncbi_species_df['species_taxid'].astype(str) == census_taxid_str
                direct_taxid_matched = ncbi_species_df[direct_taxid_mask]
            except (ValueError, TypeError) as e:
                logger.warning(f"  Skipping taxid matching for '{taxon_name}' - invalid taxid: {census_taxid}")
                direct_taxid_matched = pd.DataFrame()

        # PRIORITY 2: Hierarchical taxid matching (census taxid in NCBI lineage_taxids)
        if census_taxid and pd.notna(census_taxid):
            try:
                # Clean and convert taxid
                if isinstance(census_taxid, str):
                    census_taxid_clean = census_taxid.split('\n')[0].strip()
                    census_taxid_str = census_taxid_clean
                else:
                    census_taxid_str = str(int(census_taxid))

                # Check if census taxid appears anywhere in the lineage_taxids string
                # Pattern: ;taxid; or ^taxid; or ;taxid$ or ^taxid$
                taxid_pattern = f';{census_taxid_str};|^{census_taxid_str};|;{census_taxid_str}$|^{census_taxid_str}$'
                hierarchical_taxid_mask = ncbi_species_df['lineage_taxids'].str.contains(taxid_pattern, regex=True, na=False)
                hierarchical_taxid_matched = ncbi_species_df[hierarchical_taxid_mask]
            except (ValueError, TypeError) as e:
                logger.warning(f"  Skipping hierarchical taxid matching for '{taxon_name}' - invalid taxid: {census_taxid}")
                hierarchical_taxid_matched = pd.DataFrame()

        # Combine all matches (avoid duplicates, maintain priority order)
        all_matched_indices = set()
        matched_species_list = []

        # Add matches in priority order
        for match_df in [direct_taxid_matched, hierarchical_taxid_matched]:
            if not match_df.empty:
                # Only add rows we haven't seen yet
                new_indices = set(match_df.index) - all_matched_indices
                if new_indices:
                    matched_species_list.append(match_df.loc[list(new_indices)])
                    all_matched_indices.update(new_indices)

        # Combine all unique matches
        if matched_species_list:
            matched_species = pd.concat(matched_species_list, ignore_index=True)
        else:
            matched_species = pd.DataFrame()

        if len(matched_species) > 0:
            # Aggregate counts from all matched species
            total_species = len(matched_species)  # Count unique species
            total_genomes = matched_species['total_genome_count'].sum()
            total_isolates = matched_species['isolate_genome_count'].sum()
            total_uncultured = matched_species['uncultured_genome_count'].sum()

            # Calculate NCBI percentages
            genome_pct_db = (total_genomes / total_ncbi_genomes * 100) if total_ncbi_genomes > 0 else 0
            species_pct = (total_species / total_ncbi_species * 100) if total_ncbi_species > 0 else 0

            # Count match types (avoid double-counting)
            direct_taxid_count = len(direct_taxid_matched)
            hierarchical_taxid_count = len(hierarchical_taxid_matched) - len(hierarchical_taxid_matched[hierarchical_taxid_matched.index.isin(direct_taxid_matched.index)])

            # Track match type statistics
            if direct_taxid_count > 0:
                match_type_counts['direct_taxid'] += 1
            elif hierarchical_taxid_count > 0:
                match_type_counts['hierarchical_taxid'] += 1

            # Get domain (most common)
            domain_counts = matched_species['domain'].value_counts()
            domain = domain_counts.index[0] if len(domain_counts) > 0 else 'Unknown'

            match_status = 'matched'
            matched_count += 1
        else:
            total_species = total_genomes = total_isolates = total_uncultured = 0
            direct_taxid_count = hierarchical_taxid_count = 0
            genome_pct_db = species_pct = 0
            domain = 'Unknown'
            match_status = 'census_only'

        # Format census_taxid as integer (remove .0 from float)
        if census_taxid and pd.notna(census_taxid):
            try:
                if isinstance(census_taxid, str):
                    census_taxid_clean = census_taxid.split('\n')[0].strip()
                    census_taxid_formatted = int(census_taxid_clean)
                else:
                    census_taxid_formatted = int(census_taxid)
            except (ValueError, TypeError):
                census_taxid_formatted = str(census_taxid)  # Keep as-is if can't convert
        else:
            census_taxid_formatted = 'N/A'

        result = {
            census_level: taxon_name,
            'census_taxid': census_taxid_formatted,
            'census_otu_count': census_otus,
            'census_size_count': census_size,
            'otu_percentage': round(otu_percentage, 2),
            'size_percentage': round(size_percentage, 2),
            'ncbi_genome_count': total_genomes,
            'ncbi_species_count': total_species,
            'genome_pct_db': round(genome_pct_db, 2),
            'species_pct': round(species_pct, 2),
            'direct_name_matches': 0,  # Not used with species-level data
            'direct_taxid_matches': direct_taxid_count,
            'hierarchical_taxid_matches': hierarchical_taxid_count,
            'lineage_name_matches': 0,  # Not used with species-level data
            'total_matches': len(matched_species),
            'match_status': match_status,
            'domain': domain,
            'isolate_count': total_isolates,
            'isolate_percentage': round((total_isolates / total_genomes * 100) if total_genomes > 0 else 0, 2)
        }
        results.append(result)

    # Create output dataframe
    output_df = pd.DataFrame(results)

    match_rate = (matched_count / len(census_df) * 100) if len(census_df) > 0 else 0
    logger.info(f"  Matched: {matched_count}/{len(census_df)} ({match_rate:.1f}%)")
    logger.info(f"  Census-only: {len(census_df) - matched_count}")
    logger.info(f"  Match type breakdown:")
    logger.info(f"    Direct taxid matches: {match_type_counts['direct_taxid']}")
    logger.info(f"    Hierarchical taxid matches: {match_type_counts['hierarchical_taxid']}")

    # Log aggregated totals
    total_matched_genomes = output_df[output_df['match_status'] == 'matched']['ncbi_genome_count'].sum()
    total_matched_species = output_df[output_df['match_status'] == 'matched']['ncbi_species_count'].sum()
    logger.info(f"  Aggregated totals: {total_matched_genomes:,} genomes, {total_matched_species:,} species")

    return output_df

