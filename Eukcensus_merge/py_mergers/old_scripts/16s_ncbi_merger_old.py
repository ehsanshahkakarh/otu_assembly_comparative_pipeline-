#!/usr/bin/env python3
"""
16S Census + NCBI Merger - Vectorized Version
=============================================

Fast, vectorized merger using pandas operations for lineage-based matching.
Matches the perfect 18S merger structure and performance.

INPUT FILES:
-----------
NCBI Data:
  - ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_family_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_genus_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_*_with_accessions.csv (for isolate analysis)

16S Census Data:
  - 16S_censusparse/csv_16S/eukcensus16S_by_division.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_family.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_genus.csv

OUTPUT FILES:
------------
  - merged_output/16s_merged/results/16s_ncbi_merged_clean_*.csv
  - merged_output/16s_merged/analysis_summary/16s_ncbi_merger_clean_summary.csv

LOGIC:
------
1. Vectorized Lineage Matching: Use pandas string operations to find census taxa in NCBI lineages
2. Vectorized Aggregation: Use pandas groupby and sum operations for count aggregation
3. Isolate Analysis: Vectorized isolate percentage calculations
4. Clean Output: Same format as 18S version
5. Prokaryotic Focus: Filters for Bacteria and Archaea domains, includes organellar sequences

Author: Vectorized Merger Team
Date: 2025

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys

def setup_logging():
    """Setup minimal logging."""
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent.parent.parent
    logs_dir = base_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(logs_dir / '16s_ncbi_merger.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

class SixteenSNCBIMerger:
    """Vectorized 16S Census + NCBI Merger matching 18S structure."""

    def __init__(self):
        """
        Initialize the merger with data directories.

        Args:
            ncbi_data_dir: Directory containing NCBI CSV files
            census_data_dir: Directory containing 16S census CSV files
            output_dir: Directory for output files
        """
        self.ncbi_data_dir = Path(ncbi_data_dir)
        self.census_data_dir = Path(census_data_dir)
        self.output_dir = Path(output_dir) / "16s_merged"  # Organized subdirectory
        self.output_dir.mkdir(exist_ok=True)

        # Taxonomic levels to process
        self.taxonomic_levels = ['phylum', 'family', 'genus']

        # Total OTUs in the 16S census dataset (287,469 rows - 1 header = 287,468 OTUs)
        self.total_otus = 287468

        # Total sequences will be calculated from the census data size_count column
        self.total_sequences = None

        # Data containers
        self.ncbi_data = {}
        self.census_data = {}
        self.merged_data = {}
        self.classified_data = {}  # For isolate/uncultured classification data
        self.total_genome_counts = {}  # For total available genomes per taxa
        
    def load_ncbi_data(self) -> None:
        """Load NCBI taxonomic count data from new unified parser files."""
        logger.info("Loading NCBI data from new unified parser files...")

        for level in self.taxonomic_levels:
            # Load new unified parser file with comprehensive data
            ncbi_file = self.ncbi_data_dir / f"ncbi_{level}_counts.csv"

            if ncbi_file.exists():
                df = pd.read_csv(ncbi_file)

                # Check for new parser columns (all should be present)
                expected_cols = [
                    level, 'domain',
                    f'{level}_genome_count', f'{level}_genome_percentage',
                    f'{level}_species_count', f'{level}_species_percentage',
                    'taxid', 'lineage', 'lineage_ranks', 'lineage_taxids'
                ]

                missing_cols = [col for col in expected_cols if col not in df.columns]
                if missing_cols:
                    logger.warning(f"Missing columns in {ncbi_file}: {missing_cols}")
                    logger.warning(f"Available columns: {list(df.columns)}")
                    # Still try to use what's available, fallback to legacy if needed
                    if len(missing_cols) > 4:  # Too many missing columns
                        logger.warning(f"Too many missing columns, falling back to legacy loading")
                        self._load_ncbi_data_legacy(level)
                        continue

                # Rename columns for consistency with merger expectations
                df = df.rename(columns={
                    level: 'taxon_name',
                    f'{level}_genome_count': 'ncbi_genome_count',
                    f'{level}_species_count': 'ncbi_species_count',
                    f'{level}_genome_percentage': 'ncbi_genome_percentage',
                    f'{level}_species_percentage': 'ncbi_species_percentage'
                })

                # Store the data
                self.ncbi_data[level] = df
                logger.info(f"Loaded new NCBI {level} data: {len(df)} entries")
                logger.info(f"  - Genome count range: {df['ncbi_genome_count'].min():,} to {df['ncbi_genome_count'].max():,}")
                logger.info(f"  - Species count range: {df['ncbi_species_count'].min():,} to {df['ncbi_species_count'].max():,}")
                logger.info(f"  - Includes pre-calculated percentages and full lineage information")
            else:
                logger.warning(f"New NCBI file not found: {ncbi_file}")
                # Fallback to legacy loading method
                logger.info(f"Attempting legacy loading for {level}")
                self._load_ncbi_data_legacy(level)

    def _load_ncbi_data_legacy(self, level: str) -> None:
        """Fallback method for loading NCBI data from separate files."""
        logger.info(f"Using legacy loading method for NCBI {level} data...")

        # Load species count data
        species_file = self.ncbi_data_dir / f"ncbi_{level}_species_counts.csv"
        # Load genome count data
        genome_file = self.ncbi_data_dir / f"ncbi_{level}_counts.csv"

        if species_file.exists() and genome_file.exists():
            # Load species counts
            species_df = pd.read_csv(species_file)
            species_count_col = f"{level}_species_count"
            if species_count_col not in species_df.columns:
                logger.error(f"No {species_count_col} column found in {species_file}")
                return

            # Load genome counts
            genome_df = pd.read_csv(genome_file)
            genome_count_col = f"{level}_count"
            if genome_count_col not in genome_df.columns:
                logger.error(f"No {genome_count_col} column found in {genome_file}")
                return

            # Merge species and genome data on taxon name
            species_df = species_df.rename(columns={level: 'taxon_name'})
            genome_df = genome_df.rename(columns={level: 'taxon_name'})

            # Merge the dataframes
            df = pd.merge(species_df, genome_df[['taxon_name', genome_count_col]],
                         on='taxon_name', how='left')

            # Standardize column names
            df = df.rename(columns={
                species_count_col: 'ncbi_species_count',
                genome_count_col: 'ncbi_genome_count'
            })

            # Fill missing genome counts with 0
            df['ncbi_genome_count'] = df['ncbi_genome_count'].fillna(0).astype(int)

            self.ncbi_data[level] = df
            logger.info(f"Loaded legacy NCBI {level} data: {len(df)} entries")
        else:
            logger.warning(f"Legacy NCBI files not found: {species_file} or {genome_file}")

    def load_classified_accession_data(self) -> None:
        """Load accession data for isolate/uncultured analysis from new unified parser files."""
        logger.info("Loading accession data for isolate analysis from new parser files...")

        for level in self.taxonomic_levels:
            # Load new unified accession files (already includes genome_source classification)
            accession_file = self.ncbi_data_dir / f"ncbi_{level}_with_accessions.csv"

            if accession_file.exists():
                # Load accession data (already classified by new parser)
                accession_df = pd.read_csv(accession_file)

                # Check if genome_source column exists (should be in new format)
                if 'genome_source' not in accession_df.columns:
                    logger.warning(f"genome_source column not found in {accession_file}")
                    logger.warning(f"Available columns: {list(accession_df.columns)}")
                    logger.warning(f"Skipping isolate analysis for {level}")
                    continue

                # Calculate total genome counts from accession file (actual available genomes)
                total_counts = accession_df.groupby(level).size().reset_index(name='total_genome_count')

                # Calculate isolate counts from accession data
                isolate_df = accession_df[accession_df['genome_source'] == 'isolate']
                isolate_counts = isolate_df.groupby(level).size().reset_index(name='isolate_count')

                # Merge total and isolate counts
                metrics_df = pd.merge(total_counts, isolate_counts, on=level, how='left')
                metrics_df['isolate_count'] = metrics_df['isolate_count'].fillna(0).astype(int)

                # Rename taxon column for consistency with merger expectations
                metrics_df = metrics_df.rename(columns={level: 'taxon_name'})

                self.classified_data[level] = metrics_df
                self.total_genome_counts[level] = dict(zip(metrics_df['taxon_name'], metrics_df['total_genome_count']))

                # Store the total genome count across all taxa for genome percentage calculations
                total_genomes_all_taxa = accession_df.shape[0]  # Total rows in accession file
                setattr(self, f'total_genomes_{level}', total_genomes_all_taxa)

                logger.info(f"Loaded {level} accession data: {len(metrics_df)} taxa")
                logger.info(f"Total available genomes for {level}: {metrics_df['total_genome_count'].sum():,}")
                logger.info(f"Total genomes across all {level} taxa: {total_genomes_all_taxa:,}")
                logger.info(f"Isolate genomes for {level}: {metrics_df['isolate_count'].sum():,}")
            else:
                logger.warning(f"Accession file not found: {accession_file}")
                logger.warning(f"Isolate analysis will not be available for {level} level")

    def calculate_isolate_species_counts(self, level: str) -> dict:
        """
        Calculate isolate species counts and total species counts using species subset accession files.

        Args:
            level: Taxonomic level ('phylum', 'family', 'genus')

        Returns:
            dict with isolate_species_counts and total_species_per_taxon from species subset files
        """
        # Load the species subset accession file that contains isolate/uncultured classifications
        species_subset_accession_file = self.ncbi_data_dir / f"ncbi_{level}_species_subset_with_accessions_classified.csv"
        if not species_subset_accession_file.exists():
            logger.warning(f"Species subset accession file not found: {species_subset_accession_file}")
            # Fallback to regular species counts file for total species counts only
            species_subset_file = self.ncbi_data_dir / f"ncbi_{level}_species_counts.csv"
            if species_subset_file.exists():
                species_subset_df = pd.read_csv(species_subset_file)
                species_count_col = f"{level}_species_count"
                total_species_per_taxon = dict(zip(species_subset_df[level], species_subset_df[species_count_col]))
                logger.info(f"Loaded species subset totals for {level}: {len(total_species_per_taxon)} taxa (no isolate data)")
                return {
                    'isolate_species_counts': {},
                    'total_species_per_taxon': total_species_per_taxon
                }
            else:
                logger.warning(f"No species subset files found for {level}")
                return {}

        # Load the species subset accession file
        species_subset_df = pd.read_csv(species_subset_accession_file)
        logger.info(f"Loaded species subset accession file for {level}: {len(species_subset_df)} entries")

        # Total number of unique species in the database (number of rows in species subset file)
        total_species_in_database = len(species_subset_df)

        # Count total species per taxon (this represents the total species count from species subset)
        total_species_per_taxon = species_subset_df.groupby(level).size().to_dict()

        # Count isolate species per taxon (species that are classified as isolates)
        isolate_species_df = species_subset_df[species_subset_df['genome_source'] == 'isolate']
        isolate_species_counts = isolate_species_df.groupby(level).size().to_dict()

        logger.info(f"Calculated species counts for {level}: {len(total_species_per_taxon)} taxa total, {len(isolate_species_counts)} taxa with isolates")
        logger.info(f"Total species in database for {level}: {total_species_in_database}")

        return {
            'isolate_species_counts': isolate_species_counts,
            'total_species_per_taxon': total_species_per_taxon,
            'total_species_in_database': total_species_in_database
        }

    def load_census_data(self) -> None:
        """Load 16S census data from enhanced parser files."""
        logger.info("Loading 16S census data from enhanced parser files...")

        for level in self.taxonomic_levels:
            # Map phylum to division for 16S census files
            census_level = "division" if level == "phylum" else level
            census_file = self.census_data_dir / f"eukcensus16S_by_{census_level}.csv"
            if census_file.exists():
                df = pd.read_csv(census_file)

                # Check for enhanced parser columns
                expected_cols = ['Name_to_use', 'otu_count', 'otu_percentage', 'size_count', 'size_percentage']
                missing_cols = [col for col in expected_cols if col not in df.columns]

                if missing_cols:
                    logger.warning(f"Missing enhanced columns in {census_file}: {missing_cols}")
                    # Use available columns
                    available_cols = [col for col in expected_cols if col in df.columns]
                    logger.info(f"Using available columns: {available_cols}")

                # Standardize column names (preserve enhanced percentages if available)
                rename_dict = {
                    'Name_to_use': 'taxon_name',
                    'otu_count': 'census_otu_count',
                    'size_count': 'census_size_count'
                }

                # Preserve enhanced percentage columns if they exist
                if 'otu_percentage' in df.columns:
                    rename_dict['otu_percentage'] = 'census_otu_percentage'
                if 'size_percentage' in df.columns:
                    rename_dict['size_percentage'] = 'census_size_percentage'

                df = df.rename(columns=rename_dict)
                self.census_data[level] = df

                logger.info(f"Loaded enhanced census {level} data: {len(df)} entries")
                if 'census_otu_percentage' in df.columns:
                    logger.info(f"  - OTU percentage range: {df['census_otu_percentage'].min():.2f}% to {df['census_otu_percentage'].max():.2f}%")
                if 'census_size_percentage' in df.columns:
                    logger.info(f"  - Size percentage range: {df['census_size_percentage'].min():.2f}% to {df['census_size_percentage'].max():.2f}%")
            else:
                logger.warning(f"Census {level} file not found: {census_file}")

    def calculate_total_sequences(self) -> None:
        """
        Calculate total sequences from the original 16S census data.
        This reads the raw census file to get the sum of all size_count values.
        """
        if self.total_sequences is not None:
            return  # Already calculated

        logger.info("Calculating total sequences from census data...")

        # Read the original census file to get total sequences
        original_census_file = self.census_data_dir / "eukcensus_16S.clusters.97.tsv"
        if original_census_file.exists():
            try:
                # Read the original census file
                df = pd.read_csv(original_census_file, sep='\t')
                if 'size' in df.columns:
                    self.total_sequences = df['size'].sum()
                    logger.info(f"Total sequences calculated: {self.total_sequences:,}")
                else:
                    logger.warning("Size column not found in original census file, using default")
                    self.total_sequences = 1000000  # Default fallback
            except Exception as e:
                logger.error(f"Error reading original census file: {e}")
                self.total_sequences = 1000000  # Default fallback
        else:
            logger.warning(f"Original census file not found: {original_census_file}")
            self.total_sequences = 1000000  # Default fallback

    def clean_taxon_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean and standardize taxonomic names for better matching.
        
        Args:
            df: DataFrame with taxon_name column
            
        Returns:
            DataFrame with cleaned names
        """
        df = df.copy()
        
        # Remove common prefixes and suffixes
        df['taxon_name_clean'] = df['taxon_name'].str.replace(r'^Candidatus\s+', '', regex=True)
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'\.Mitochondria$', '', regex=True)
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'\.Chloroplast$', '', regex=True)
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'\.Apicoplast$', '', regex=True)
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r':plas\.Chloroplast$', '', regex=True)
        
        # Handle special cases
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'_[A-Z]$', '', regex=True)  # Remove GTDB-style suffixes
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'_\d+$', '', regex=True)   # Remove numeric suffixes
        
        # Strip whitespace
        df['taxon_name_clean'] = df['taxon_name_clean'].str.strip()
        
        return df
    
    def filter_prokaryotic_data(self, df: pd.DataFrame, data_source: str) -> pd.DataFrame:
        """
        Filter data to include only prokaryotic organisms (Bacteria and Archaea).
        
        Args:
            df: Input DataFrame
            data_source: 'ncbi' or 'census'
            
        Returns:
            Filtered DataFrame
        """
        if data_source == 'ncbi' and 'domain' in df.columns:
            # Filter NCBI data by domain
            prokaryotic_domains = ['Bacteria', 'Archaea']
            df_filtered = df[df['domain'].isin(prokaryotic_domains)].copy()
            logger.info(f"Filtered NCBI data: {len(df_filtered)} prokaryotic entries from {len(df)} total")
            
        elif data_source == 'census':
            # Filter census data by lineage information
            if 'lineage' in df.columns:
                # Keep entries that have Bacteria or Archaea in lineage, or organellar sequences
                mask = (
                    df['lineage'].str.contains('Bacteria', na=False) |
                    df['lineage'].str.contains('Archaea', na=False) |
                    df['taxon_name'].str.contains('Mitochondria|Chloroplast|Apicoplast', na=False)
                )
                df_filtered = df[mask].copy()
                logger.info(f"Filtered census data: {len(df_filtered)} prokaryotic/organellar entries from {len(df)} total")
            else:
                # If no lineage info, assume all are prokaryotic (16S data)
                df_filtered = df.copy()
                logger.info(f"No lineage filtering applied to census data: {len(df_filtered)} entries")
        else:
            df_filtered = df.copy()
            
        return df_filtered

    def add_isolate_metrics(self, df: pd.DataFrame, level: str) -> pd.DataFrame:
        """
        Add isolate count and percentage metrics to the merged DataFrame.

        Args:
            df: Merged DataFrame
            level: Taxonomic level ('phylum', 'family', 'genus')

        Returns:
            DataFrame with added isolate metrics
        """
        if level not in self.classified_data:
            logger.warning(f"No classified data available for {level} - skipping isolate metrics")
            return df

        classified_df = self.classified_data[level]

        # Merge with classified data to get isolate counts
        # Use left merge to preserve all taxa in the current dataset
        df = pd.merge(df, classified_df[['taxon_name', 'total_genome_count', 'isolate_count']],
                     on='taxon_name', how='left')

        # Fill missing values with 0
        df['total_genome_count'] = df['total_genome_count'].fillna(0).astype(int)
        df['isolate_count'] = df['isolate_count'].fillna(0).astype(int)

        # CONSISTENT DENOMINATOR APPROACH:
        # - Genome-related metrics: use total_genomes_all_taxa (total genomes in database)
        # - Species-related metrics: use total_species_in_database (total species in database)

        total_genomes_all_taxa = getattr(self, f'total_genomes_{level}', 1)

        # GENOME METRICS (denominator: total genomes in database)
        df['isolate_genome_percentage'] = np.where(
            total_genomes_all_taxa > 0,
            (df['isolate_count'] / total_genomes_all_taxa * 100).round(2),
            0
        )

        df['dataset_genome_percentage'] = np.where(
            total_genomes_all_taxa > 0,
            (df['ncbi_genome_count'] / total_genomes_all_taxa * 100).round(2),
            0
        )

        # SPECIES METRICS (denominator: total species in database)
        isolate_species_data = self.calculate_isolate_species_counts(level)
        if isolate_species_data:
            isolate_species_counts = isolate_species_data.get('isolate_species_counts', {})
            total_species_in_database = isolate_species_data.get('total_species_in_database', 1)

            # Add isolate species count column
            df['isolate_species_count'] = df['taxon_name'].map(isolate_species_counts).fillna(0).astype(int)

            # Calculate species percentages using consistent denominator (total species in database)
            df['isolate_species_percentage'] = np.where(
                total_species_in_database > 0,
                (df['isolate_species_count'] / total_species_in_database * 100).round(2),
                0
            )

            df['dataset_species_percentage'] = np.where(
                total_species_in_database > 0,
                (df['ncbi_species_count'] / total_species_in_database * 100).round(2),
                0
            )

        logger.info(f"Added isolate metrics for {level} level")
        return df

    def merge_taxonomic_level(self, level: str) -> pd.DataFrame:
        """
        Merge NCBI and census data for a specific taxonomic level.
        
        Args:
            level: Taxonomic level ('phylum', 'family', 'genus')
            
        Returns:
            Merged DataFrame
        """
        logger.info(f"Merging {level} level data...")
        
        # Get data for this level
        ncbi_df = self.ncbi_data.get(level)
        census_df = self.census_data.get(level)
        
        if ncbi_df is None and census_df is None:
            logger.warning(f"No data available for {level} level")
            return pd.DataFrame()
        
        # Initialize result DataFrame
        merged_df = pd.DataFrame()
        
        if ncbi_df is not None:
            # Filter and clean NCBI data
            ncbi_df = self.filter_prokaryotic_data(ncbi_df, 'ncbi')
            ncbi_df = self.clean_taxon_names(ncbi_df)
            
        if census_df is not None:
            # Filter and clean census data
            census_df = self.filter_prokaryotic_data(census_df, 'census')
            census_df = self.clean_taxon_names(census_df)
        
        # Perform the merge
        if ncbi_df is not None and census_df is not None:
            # Full outer join on cleaned names
            merged_df = pd.merge(
                ncbi_df, census_df,
                left_on='taxon_name_clean',
                right_on='taxon_name_clean',
                how='outer',
                suffixes=('_ncbi', '_census')
            )

            # Create unified taxon name - handle the suffix columns correctly
            if 'taxon_name_ncbi' in merged_df.columns and 'taxon_name_census' in merged_df.columns:
                merged_df['taxon_name'] = merged_df['taxon_name_ncbi'].fillna(merged_df['taxon_name_census'])
            elif 'taxon_name_ncbi' in merged_df.columns:
                merged_df['taxon_name'] = merged_df['taxon_name_ncbi']
            elif 'taxon_name_census' in merged_df.columns:
                merged_df['taxon_name'] = merged_df['taxon_name_census']
            else:
                # Fallback - use the original taxon_name if it exists
                if 'taxon_name' not in merged_df.columns:
                    merged_df['taxon_name'] = merged_df['taxon_name_clean']

        elif ncbi_df is not None:
            # Only NCBI data available
            merged_df = ncbi_df.copy()

        elif census_df is not None:
            # Only census data available
            merged_df = census_df.copy()
        
        # Calculate coverage metrics
        merged_df = self.calculate_coverage_metrics(merged_df)

        # Add isolate metrics
        merged_df = self.add_isolate_metrics(merged_df, level)

        # Group by taxonomic name and aggregate counts (in case of duplicates)
        if not merged_df.empty and 'taxon_name' in merged_df.columns:
            # Check if there are any duplicates to group
            duplicate_count = merged_df.duplicated(subset=['taxon_name']).sum()
            if duplicate_count > 0:
                logger.info(f"Found {duplicate_count} duplicate entries to aggregate")

                agg_dict = {
                    'census_otu_count': 'sum',
                    'census_size_count': 'sum',
                    'ncbi_species_count': 'sum',
                    'ncbi_genome_count': 'sum',
                    'taxon_name_clean': 'first',
                    'in_ncbi': 'max',
                    'in_census': 'max',
                    'in_both': 'max',
                    'total_genome_count': 'first',
                    'isolate_count': 'first',
                    'isolate_genome_percentage': 'first',
                    'dataset_genome_percentage': 'mean',
                    'isolate_species_count': 'first',
                    'isolate_species_percentage': 'first',
                    'dataset_species_percentage': 'first'
                }

                # Only include columns that exist
                agg_dict = {k: v for k, v in agg_dict.items() if k in merged_df.columns}

                merged_df = merged_df.groupby('taxon_name').agg(agg_dict).reset_index()

        # Calculate coverage percentage first for sorting
        merged_df['Coverage_Percentage'] = np.where(
            (merged_df['census_otu_count'] > 0) & (merged_df['ncbi_species_count'] > 0),
            (merged_df['ncbi_species_count'] / merged_df['census_otu_count'] * 100).round(2),
            0
        )

        # Sort by Coverage Percentage in descending order
        merged_df = merged_df.sort_values('Coverage_Percentage', ascending=False, na_position='last')

        # Reorder columns for better readability
        merged_df = self.reorder_columns(merged_df, level)

        logger.info(f"Merged {level} data: {len(merged_df)} total entries")

        return merged_df

    def calculate_coverage_metrics(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate coverage and representation metrics.

        Args:
            df: Merged DataFrame

        Returns:
            DataFrame with additional metrics
        """
        df = df.copy()

        # Fill NaN values with 0 for calculations
        df['ncbi_species_count'] = df['ncbi_species_count'].fillna(0)
        df['census_otu_count'] = df['census_otu_count'].fillna(0)

        # Presence indicators
        df['in_ncbi'] = (df['ncbi_species_count'] > 0).astype(int)
        df['in_census'] = (df['census_otu_count'] > 0).astype(int)
        df['in_both'] = ((df['ncbi_species_count'] > 0) & (df['census_otu_count'] > 0)).astype(int)

        # Coverage ratios (avoid division by zero)
        df['ncbi_to_census_ratio'] = np.where(
            df['census_otu_count'] > 0,
            df['ncbi_species_count'] / df['census_otu_count'],
            np.inf
        )

        df['census_to_ncbi_ratio'] = np.where(
            df['ncbi_species_count'] > 0,
            df['census_otu_count'] / df['ncbi_species_count'],
            np.inf
        )

        # Representation bias (log2 ratio)
        df['representation_bias'] = np.where(
            (df['ncbi_species_count'] > 0) & (df['census_otu_count'] > 0),
            np.log2(df['ncbi_species_count'] / df['census_otu_count']),
            np.nan
        )

        return df

    def reorder_columns(self, df: pd.DataFrame, level: str) -> pd.DataFrame:
        """
        Reorder columns for clean output format.

        Args:
            df: DataFrame to reorder
            level: Taxonomic level

        Returns:
            DataFrame with essential columns only
        """
        # Calculate coverage percentage (for matched taxa only)
        df['Coverage_Percentage'] = np.where(
            (df['census_otu_count'] > 0) & (df['ncbi_species_count'] > 0),
            (df['ncbi_species_count'] / df['census_otu_count'] * 100).round(2),
            0
        )

        # Create clean DataFrame with essential columns only

        # Create clean DataFrame with renamed columns
        clean_df = pd.DataFrame()

        # Map existing columns to clean names (matching original parser column names)
        if 'taxon_name' in df.columns:
            clean_df[level.capitalize()] = df['taxon_name']
        if 'census_otu_count' in df.columns:
            clean_df['otu_count'] = df['census_otu_count'].fillna(0).astype(int)
            # Use enhanced OTU percentage if available, otherwise calculate
            if 'census_otu_percentage' in df.columns:
                clean_df['otu_percentage'] = df['census_otu_percentage'].fillna(0).round(2)
                logger.info(f"Using enhanced OTU percentages from parser")
            else:
                clean_df['otu_percentage'] = (df['census_otu_count'].fillna(0) / self.total_otus * 100).round(2)
                logger.info(f"Calculating OTU percentages in merger")
        if 'census_size_count' in df.columns:
            clean_df['size_count'] = df['census_size_count'].fillna(0).astype(int)
            # Use enhanced Size percentage if available, otherwise calculate
            if 'census_size_percentage' in df.columns:
                clean_df['size_percentage'] = df['census_size_percentage'].fillna(0).round(2)
                logger.info(f"Using enhanced Size percentages from parser")
            else:
                if self.total_sequences and self.total_sequences > 0:
                    clean_df['size_percentage'] = (df['census_size_count'].fillna(0) / self.total_sequences * 100).round(2)
                    logger.info(f"Calculating Size percentages in merger")
                else:
                    clean_df['size_percentage'] = 0.0
        # Add NCBI data with original parser column names (level-specific)
        if 'ncbi_species_count' in df.columns:
            clean_df[f'{level}_species_count'] = df['ncbi_species_count'].fillna(0).astype(int)
        if 'ncbi_genome_count' in df.columns:
            clean_df[f'{level}_genome_count'] = df['ncbi_genome_count'].fillna(0).astype(int)

        # Add enhanced NCBI percentages with original parser column names (level-specific)
        if 'ncbi_genome_percentage' in df.columns:
            clean_df[f'{level}_genome_percentage'] = df['ncbi_genome_percentage'].fillna(0).round(2)
            logger.info(f"Using enhanced NCBI genome percentages from parser")
        if 'ncbi_species_percentage' in df.columns:
            clean_df[f'{level}_species_percentage'] = df['ncbi_species_percentage'].fillna(0).round(2)
            logger.info(f"Using enhanced NCBI species percentages from parser")
        if 'Coverage_Percentage' in df.columns:
            clean_df['Coverage_Percentage'] = df['Coverage_Percentage']

        # Add isolate metrics with cleaner names
        if 'isolate_count' in df.columns:
            clean_df['Isolate_Count'] = df['isolate_count'].fillna(0).astype(int)
        if 'isolate_genome_percentage' in df.columns:
            clean_df['Isolate_Genome_Percentage'] = df['isolate_genome_percentage']
        if 'dataset_genome_percentage' in df.columns:
            clean_df['Dataset_Genome_Percentage'] = df['dataset_genome_percentage']

        # Add species metrics with cleaner names
        if 'isolate_species_count' in df.columns:
            clean_df['Isolate_Species_Count'] = df['isolate_species_count'].fillna(0).astype(int)
        if 'isolate_species_percentage' in df.columns:
            clean_df['Isolate_Species_Percentage'] = df['isolate_species_percentage']
        if 'dataset_species_percentage' in df.columns:
            clean_df['Dataset_Species_Percentage'] = df['dataset_species_percentage']

        # Removed Total_Genome_Count column as requested

        # Add visualization compatibility columns
        if 'in_both' in df.columns:
            clean_df['in_both'] = df['in_both']
        # Note: Removed redundant taxon_name_clean and taxon_name columns
        # The main taxon column (e.g., 'Phylum') already serves this purpose

        return clean_df

    def generate_summary_statistics(self, merged_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Generate summary statistics across all taxonomic levels."""
        logger.info("Generating summary statistics...")

        summary_stats = []

        for level, df in merged_data.items():
            if df.empty:
                continue

            total_taxa = len(df)
            # Matched taxa = taxa that exist in BOTH census and NCBI (in_both = 1)
            matched_taxa = len(df[df['in_both'] == 1])

            # Coverage percentage = matched taxa / total taxa * 100
            coverage_pct = (matched_taxa / total_taxa * 100) if total_taxa > 0 else 0

            summary_stats.append({
                'Taxonomic_Level': level.capitalize(),
                'Total_Taxa': total_taxa,
                'Matched_Taxa': matched_taxa,
                'Coverage_Percentage': f"{coverage_pct:.1f}%"
            })

        return pd.DataFrame(summary_stats)

    def save_results(self, merged_data: Dict[str, pd.DataFrame], summary_stats: pd.DataFrame) -> None:
        """
        Save merged results and summary statistics.

        Args:
            merged_data: Dictionary of merged DataFrames
            summary_stats: Summary statistics DataFrame
        """
        logger.info("Saving results...")

        # Save individual taxonomic level files
        for level, df in merged_data.items():
            if not df.empty:
                # Apply taxonomic filters if enabled
                original_count = len(df)
                level_col = level.capitalize()

                if level_col in df.columns:
                    # Apply Candidatus filter if enabled
                    if 'EXCLUDE_CANDIDAT_TAXA' in globals() and EXCLUDE_CANDIDAT_TAXA:
                        df_before = len(df)
                        df = df[~df[level_col].str.lower().str.startswith('candidat', na=False)]
                        candidat_removed = df_before - len(df)
                        if candidat_removed > 0:
                            logger.info(f"Candidatus filter: {level} - removed {candidat_removed} entries")

                    # Apply .U. filter if enabled
                    if 'EXCLUDE_U_TAXA' in globals() and EXCLUDE_U_TAXA:
                        df_before = len(df)
                        df = df[~df[level_col].str.contains(r'\.U\.', na=False, regex=True)]
                        u_removed = df_before - len(df)
                        if u_removed > 0:
                            logger.info(f".U. filter: {level} - removed {u_removed} entries")

                    total_removed = original_count - len(df)
                    if total_removed > 0:
                        logger.info(f"Total filtered: {level} - removed {total_removed} entries, kept {len(df)}")
                else:
                    logger.warning(f"Cannot apply taxonomic filters: {level_col} column not found")

                output_file = self.output_dir / f"16s_ncbi_merged_{level}.csv"
                df.to_csv(output_file, index=False)
                logger.info(f"Saved {level} merged data: {output_file}")

        # Save summary statistics
        summary_file = self.output_dir / "16s_ncbi_merger_summary.csv"
        summary_stats.to_csv(summary_file, index=False)
        logger.info(f"Saved summary statistics: {summary_file}")

        # Save detailed log
        self.save_detailed_log(merged_data, summary_stats)

    def save_detailed_log(self, merged_data: Dict[str, pd.DataFrame], summary_stats: pd.DataFrame) -> None:
        """
        Save detailed analysis log.

        Args:
            merged_data: Dictionary of merged DataFrames
            summary_stats: Summary statistics DataFrame
        """
        log_file = self.output_dir / "16s_ncbi_merger_analysis.txt"

        with open(log_file, 'w') as f:
            f.write("16S Census + NCBI Merger Analysis Report\n")
            f.write("=" * 50 + "\n\n")

            # Overall summary
            f.write("SUMMARY STATISTICS\n")
            f.write("-" * 20 + "\n")
            for _, row in summary_stats.iterrows():
                f.write(f"\n{row['Taxonomic_Level'].upper()} LEVEL:\n")
                f.write(f"  Total taxa: {row['Total_Taxa']}\n")
                f.write(f"  Matched taxa: {row['Matched_Taxa']}\n")
                f.write(f"  Coverage: {row['Coverage_Percentage']}\n")

            # Detailed breakdowns
            f.write("\n\nDETAILED ANALYSIS\n")
            f.write("-" * 20 + "\n")

            for level, df in merged_data.items():
                if df.empty:
                    continue

                f.write(f"\n{level.upper()} LEVEL DETAILS:\n")

                # Isolate metrics summary
                if 'isolate_genome_percentage' in df.columns and 'dataset_genome_percentage' in df.columns:
                    avg_isolate_genome_pct = df['isolate_genome_percentage'].mean()
                    avg_dataset_genome_pct = df['dataset_genome_percentage'].mean()
                    total_isolates = df['isolate_count'].sum() if 'isolate_count' in df.columns else 0
                    total_genomes = df['total_genome_count'].sum() if 'total_genome_count' in df.columns else 0

                    f.write(f"\nISOLATE METRICS:\n")
                    f.write(f"  Average isolate genome percentage: {avg_isolate_genome_pct:.2f}%\n")
                    f.write(f"  Average dataset genome percentage: {avg_dataset_genome_pct:.2f}%\n")
                    f.write(f"  Total isolate genomes: {total_isolates:,}\n")
                    f.write(f"  Total available genomes: {total_genomes:,}\n")

                # Top represented taxa
                if 'representation_bias' in df.columns:
                    top_represented = df.nlargest(10, 'representation_bias')
                    f.write(f"\nTop 10 over-represented in NCBI:\n")
                    for _, row in top_represented.iterrows():
                        if not pd.isna(row['representation_bias']):
                            f.write(f"  {row['taxon_name']}: {row['representation_bias']:.2f}\n")

                    bottom_represented = df.nsmallest(10, 'representation_bias')
                    f.write(f"\nTop 10 under-represented in NCBI:\n")
                    for _, row in bottom_represented.iterrows():
                        if not pd.isna(row['representation_bias']):
                            f.write(f"  {row['taxon_name']}: {row['representation_bias']:.2f}\n")

        logger.info(f"Saved detailed analysis: {log_file}")

    def run_merger(self) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
        """
        Run the complete merger pipeline.

        Returns:
            Tuple of (merged_data, summary_stats)
        """
        logger.info("Starting 16S census + NCBI merger...")

        # Load data
        self.load_ncbi_data()
        self.load_census_data()
        self.calculate_total_sequences()  # Calculate total sequences for Size_Percentage
        self.load_classified_accession_data()

        # Merge data for each taxonomic level
        for level in self.taxonomic_levels:
            merged_df = self.merge_taxonomic_level(level)
            if not merged_df.empty:
                self.merged_data[level] = merged_df

        # Generate summary statistics
        summary_stats = self.generate_summary_statistics(self.merged_data)

        # Save results
        self.save_results(self.merged_data, summary_stats)

        logger.info("16S census + NCBI merger completed successfully!")

        return self.merged_data, summary_stats


def main():
    """
    Main function to run the 16S + NCBI merger.
    """
    # Get the script directory and workspace root
    script_dir = Path(__file__).parent
    workspace_root = script_dir.parent.parent.parent.parent.parent / "parse_repaa_table"  # Go up to metadata_proj/parse_repaa_table

    # Default paths relative to workspace root
    ncbi_data_dir = workspace_root / "ncbi_parse" / "csv_ncbi"
    census_data_dir = workspace_root / "16S_censusparse" / "csv_16S"
    output_dir = script_dir / "merged_output"

    # Allow command line arguments
    if len(sys.argv) >= 2:
        ncbi_data_dir = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        census_data_dir = Path(sys.argv[2])
    if len(sys.argv) >= 4:
        output_dir = Path(sys.argv[3])

    # Create and run merger
    merger = SixteenSNCBIMerger(str(ncbi_data_dir), str(census_data_dir), str(output_dir))
    _, summary_stats = merger.run_merger()

    # Print summary to console
    print("\n16S + NCBI MERGER SUMMARY")
    print("=" * 30)
    print(summary_stats.to_string(index=False))


if __name__ == "__main__":
    main()
