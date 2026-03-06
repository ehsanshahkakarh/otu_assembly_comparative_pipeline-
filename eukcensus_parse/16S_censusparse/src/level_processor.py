"""
Taxonomic Level Processing for 16S Census Parser
=================================================

Processes specific taxonomic levels (division/phylum, family, genus)
with enhanced filtering and rank validation.
"""

import pandas as pd
from collections import defaultdict
from tqdm import tqdm

from .taxon_cleaner import (
    should_filter_taxon,
    extract_appropriate_rank_name,
    clean_organelle_taxon_name
)


def process_taxonomic_level(df, col_name, target_rank):
    """
    Process a specific taxonomic level with enhanced filtering and rank validation.

    Args:
        df: Input dataframe
        col_name: Column name for this taxonomic level
        target_rank: Target taxonomic rank

    Returns:
        Dictionary of processed data for this level, keyed by ORIGINAL names
    """
    print(f"📊 Processing {target_rank} level...")

    # Filter out only NaN values - preserve all identified taxa including .U. entries
    level_df = df[~df[col_name].isna()]
    filtered_count = len(df) - len(level_df)

    # Initialize data dictionary - KEY CHANGE: Use original names as keys
    # Added size_count to track total sequence count (sum of cluster sizes)
    data_dict = defaultdict(lambda: {'otu_count': 0, 'size_count': 0, 'cleaned_name': '', 'appropriate_name': ''})

    # Process each entry with progress bar
    for _, row in tqdm(level_df.iterrows(), total=len(level_df), desc=f"Processing {target_rank} entries", leave=False, ncols=80):
        original_taxon = row[col_name]

        if should_filter_taxon(original_taxon):
            continue

        # Extract appropriate rank name for taxonkit lookup
        appropriate_name = extract_appropriate_rank_name(original_taxon, target_rank)

        if appropriate_name is None:
            continue  # Filtered out

        # Get the size value from the row (number of sequences in this cluster)
        cluster_size = row.get('size', 1)  # Default to 1 if size column missing
        if pd.isna(cluster_size):
            cluster_size = 1

        # Store data using ORIGINAL name as key
        data_dict[original_taxon]['otu_count'] += 1  # Count of clusters/OTUs
        data_dict[original_taxon]['size_count'] += int(cluster_size)  # Sum of sequence counts
        data_dict[original_taxon]['cleaned_name'] = clean_organelle_taxon_name(original_taxon)
        data_dict[original_taxon]['appropriate_name'] = appropriate_name

    print(f"✅ Processed {len(data_dict)} unique {target_rank} entries (filtered {filtered_count} null entries, preserved .U. and unidentified taxa)")

    return dict(data_dict)

