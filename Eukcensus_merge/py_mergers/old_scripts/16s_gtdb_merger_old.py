#!/usr/bin/env python3
"""
16S Census + GTDB Merger
========================

Merges GTDB genomic data with 16S census data to create comprehensive comparison tables.
Handles taxonomic matching across different data sources and provides coverage analysis.

INPUT FILES:
-----------
GTDB Data:
  - gtdb_parse/csv_gtdb/gtdb_phylum_species_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_phylum_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_family_species_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_family_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_genus_species_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_genus_counts.csv

16S Census Data:
  - 16S_censusparse/csv_16S/eukcensus16S_by_division.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_family.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_genus.csv
  - 16S_censusparse/csv_16S/eukcensus_16S.clusters.97.tsv (totals calculation)

OUTPUT FILES:
------------
Merged Results:
  - merged_output/16s_merged/16s_gtdb_merged_phylum.csv
  - merged_output/16s_merged/16s_gtdb_merged_family.csv
  - merged_output/16s_merged/16s_gtdb_merged_genus.csv

Summary & Analysis:
  - merged_output/16s_merged/16s_gtdb_merger_summary.csv
  - merged_output/16s_merged/16s_gtdb_merger_analysis.txt
  - logs/16s_gtdb_merger.log

MERGING LOGIC:
-------------
1. Direct Taxonomic Matching: Merges GTDB and census data on cleaned taxonomic names
   using full outer join to preserve all entries
2. Prokaryotic Focus: Filters for Bacteria and Archaea domains, includes organellar sequences
3. Name Standardization: Removes Candidatus prefixes, organellar suffixes, and GTDB-style
   suffixes for better matching
4. Coverage Calculation: Coverage = (GTDB_Species_Count / Census_OTU_Count) × 100
5. Representation Analysis: Calculates representation bias using log2 ratios
6. Candidatus Filtering: Optional filtering of Candidatus taxa (configurable)
7. Comprehensive Output: Includes presence indicators and percentage calculations

Author: Enhanced EukCensus Parser Team
Date: 2024
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import logging
from typing import Dict, Tuple

# ============================================================================
# CANDIDAT FILTER CONFIGURATION - Simple output filter
# ============================================================================
# To EXCLUDE candidat taxa from output: uncomment the next line
EXCLUDE_CANDIDAT_TAXA = True
# To INCLUDE candidat taxa in output: comment out the line above
# ============================================================================

# Setup logging
script_dir = Path(__file__).parent
log_dir = script_dir / 'logs'
log_dir.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_dir / '16s_gtdb_merger.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class SixteenSGTDBMerger:
    """
    Merges GTDB genomic data with 16S census data for prokaryotic organisms.
    """
    
    def __init__(self, gtdb_data_dir: str, census_data_dir: str, output_dir: str = "merged_output"):
        """
        Initialize the merger with data directories.
        
        Args:
            gtdb_data_dir: Directory containing GTDB CSV files
            census_data_dir: Directory containing 16S census CSV files  
            output_dir: Directory for output files
        """
        self.gtdb_data_dir = Path(gtdb_data_dir)
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
        self.gtdb_data = {}
        self.census_data = {}
        self.merged_data = {}
        
    def load_gtdb_data(self) -> None:
        """Load GTDB taxonomic count data."""
        logger.info("Loading GTDB data...")

        for level in self.taxonomic_levels:
            # Load species count data
            species_file = self.gtdb_data_dir / f"gtdb_{level}_species_counts.csv"
            # Load genome count data
            genome_file = self.gtdb_data_dir / f"gtdb_{level}_counts.csv"

            if species_file.exists() and genome_file.exists():
                # Load species counts
                species_df = pd.read_csv(species_file)
                species_count_col = f"{level}_species_count"
                if species_count_col not in species_df.columns:
                    logger.error(f"No {species_count_col} column found in {species_file}")
                    logger.info(f"Available columns in {species_file}: {list(species_df.columns)}")
                    continue

                # Load genome counts
                genome_df = pd.read_csv(genome_file)
                
                # Try different possible column names for genome counts
                possible_genome_cols = [f"{level}_count", f"{level}_counts", f"genome_count", f"count"]
                genome_count_col = None
                
                for col_name in possible_genome_cols:
                    if col_name in genome_df.columns:
                        genome_count_col = col_name
                        break
                
                if genome_count_col is None:
                    logger.error(f"No genome count column found in {genome_file}")
                    logger.info(f"Available columns in {genome_file}: {list(genome_df.columns)}")
                    continue

                # Merge species and genome data on taxon name
                species_df = species_df.rename(columns={level: 'taxon_name'})
                genome_df = genome_df.rename(columns={level: 'taxon_name'})

                # Merge the dataframes
                df = pd.merge(species_df, genome_df[['taxon_name', genome_count_col]],
                             on='taxon_name', how='left')

                # Standardize column names
                df = df.rename(columns={
                    species_count_col: 'gtdb_species_count',
                    genome_count_col: 'gtdb_genome_count'
                })

                # Fill missing genome counts with 0
                df['gtdb_genome_count'] = df['gtdb_genome_count'].fillna(0).astype(int)

                self.gtdb_data[level] = df
                logger.info(f"Using {species_count_col} and {genome_count_col} columns for {level}")
                logger.info(f"Loaded GTDB {level} data: {len(df)} entries")
            else:
                logger.warning(f"GTDB files not found: {species_file} or {genome_file}")
                
    def load_census_data(self) -> None:
        """Load 16S census data."""
        logger.info("Loading 16S census data...")
        
        for level in self.taxonomic_levels:
            # Map phylum to division for 16S census files
            census_level = "division" if level == "phylum" else level
            census_file = self.census_data_dir / f"eukcensus16S_by_{census_level}.csv"
            if census_file.exists():
                df = pd.read_csv(census_file)
                # Standardize column names
                df = df.rename(columns={
                    'Name_to_use': 'taxon_name',
                    'otu_count': 'census_otu_count',
                    'size_count': 'census_size_count'  # Add size_count column processing
                })
                self.census_data[level] = df
                logger.info(f"Loaded census {level} data: {len(df)} entries")
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
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'_[A-Z]$', '', regex=True)  # Remove GTDB suffixes
        df['taxon_name_clean'] = df['taxon_name_clean'].str.replace(r'_\d+$', '', regex=True)   # Remove numeric suffixes
        
        # Strip whitespace
        df['taxon_name_clean'] = df['taxon_name_clean'].str.strip()
        
        return df
    
    def filter_prokaryotic_data(self, df: pd.DataFrame, data_source: str) -> pd.DataFrame:
        """
        Filter data to include only prokaryotic organisms (Bacteria and Archaea).
        
        Args:
            df: Input DataFrame
            data_source: 'gtdb' or 'census'
            
        Returns:
            Filtered DataFrame
        """
        if data_source == 'gtdb' and 'domain' in df.columns:
            # Filter GTDB data by domain
            prokaryotic_domains = ['Bacteria', 'Archaea']
            df_filtered = df[df['domain'].isin(prokaryotic_domains)].copy()
            logger.info(f"Filtered GTDB data: {len(df_filtered)} prokaryotic entries from {len(df)} total")
            
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
    
    def merge_taxonomic_level(self, level: str) -> pd.DataFrame:
        """
        Merge GTDB and census data for a specific taxonomic level.
        
        Args:
            level: Taxonomic level ('phylum', 'family', 'genus')
            
        Returns:
            Merged DataFrame
        """
        logger.info(f"Merging {level} level data...")
        
        # Get data for this level
        gtdb_df = self.gtdb_data.get(level)
        census_df = self.census_data.get(level)
        
        if gtdb_df is None and census_df is None:
            logger.warning(f"No data available for {level} level")
            return pd.DataFrame()
        
        # Initialize result DataFrame
        merged_df = pd.DataFrame()
        
        if gtdb_df is not None:
            # Filter and clean GTDB data
            gtdb_df = self.filter_prokaryotic_data(gtdb_df, 'gtdb')
            gtdb_df = self.clean_taxon_names(gtdb_df)
            
        if census_df is not None:
            # Filter and clean census data
            census_df = self.filter_prokaryotic_data(census_df, 'census')
            census_df = self.clean_taxon_names(census_df)
        
        # Perform the merge
        if gtdb_df is not None and census_df is not None:
            # Full outer join on cleaned names
            merged_df = pd.merge(
                gtdb_df, census_df,
                left_on='taxon_name_clean',
                right_on='taxon_name_clean',
                how='outer',
                suffixes=('_gtdb', '_census')
            )
            
            # Create unified taxon name
            merged_df['taxon_name'] = merged_df['taxon_name_gtdb'].fillna(merged_df['taxon_name_census'])
            
        elif gtdb_df is not None:
            # Only GTDB data available
            merged_df = gtdb_df.copy()

        elif census_df is not None:
            # Only census data available
            merged_df = census_df.copy()
        
        # Calculate coverage metrics
        merged_df = self.calculate_coverage_metrics(merged_df)

        # Group by taxonomic name and aggregate counts (in case of duplicates)
        if not merged_df.empty and 'taxon_name' in merged_df.columns:
            # Check if there are any duplicates to group
            duplicate_count = merged_df.duplicated(subset=['taxon_name']).sum()
            if duplicate_count > 0:
                logger.info(f"Found {duplicate_count} duplicate entries to aggregate")

                agg_dict = {
                    'census_otu_count': 'sum',
                    'census_size_count': 'sum',  # Add size_count aggregation
                    'gtdb_species_count': 'sum',
                    'gtdb_genome_count': 'sum',
                    'taxon_name_clean': 'first',
                    'in_gtdb': 'max',
                    'in_census': 'max',
                    'in_both': 'max'
                }

                # Only include columns that exist
                agg_dict = {k: v for k, v in agg_dict.items() if k in merged_df.columns}

                merged_df = merged_df.groupby('taxon_name').agg(agg_dict).reset_index()

        # Calculate coverage percentage first for sorting
        merged_df['Coverage_Percentage'] = np.where(
            (merged_df['census_otu_count'] > 0) & (merged_df['gtdb_species_count'] > 0),
            (merged_df['gtdb_species_count'] / merged_df['census_otu_count'] * 100).round(2),
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
        df['gtdb_species_count'] = df['gtdb_species_count'].fillna(0)
        df['census_otu_count'] = df['census_otu_count'].fillna(0)

        # Presence indicators
        df['in_gtdb'] = (df['gtdb_species_count'] > 0).astype(int)
        df['in_census'] = (df['census_otu_count'] > 0).astype(int)
        df['in_both'] = ((df['gtdb_species_count'] > 0) & (df['census_otu_count'] > 0)).astype(int)

        # Coverage ratios (avoid division by zero)
        df['gtdb_to_census_ratio'] = np.where(
            df['census_otu_count'] > 0,
            df['gtdb_species_count'] / df['census_otu_count'],
            np.inf
        )

        df['census_to_gtdb_ratio'] = np.where(
            df['gtdb_species_count'] > 0,
            df['census_otu_count'] / df['gtdb_species_count'],
            np.inf
        )

        # Representation bias (log2 ratio)
        df['representation_bias'] = np.where(
            (df['gtdb_species_count'] > 0) & (df['census_otu_count'] > 0),
            np.log2(df['gtdb_species_count'] / df['census_otu_count']),
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
            (df['census_otu_count'] > 0) & (df['gtdb_species_count'] > 0),
            (df['gtdb_species_count'] / df['census_otu_count'] * 100).round(2),
            0
        )

        # Create clean DataFrame with essential columns only
        clean_df = pd.DataFrame()

        # Map existing columns to clean names
        if 'taxon_name' in df.columns:
            clean_df[level.capitalize()] = df['taxon_name']
        if 'census_otu_count' in df.columns:
            clean_df['Census_OTU_Count'] = df['census_otu_count'].fillna(0).astype(int)
            # Calculate OTU percentage
            clean_df['OTU_Percentage'] = (df['census_otu_count'].fillna(0) / self.total_otus * 100).round(2)
        if 'census_size_count' in df.columns:
            clean_df['Census_Size_Count'] = df['census_size_count'].fillna(0).astype(int)
            # Calculate Size percentage
            if self.total_sequences and self.total_sequences > 0:
                clean_df['Size_Percentage'] = (df['census_size_count'].fillna(0) / self.total_sequences * 100).round(2)
            else:
                clean_df['Size_Percentage'] = 0.0
        if 'gtdb_species_count' in df.columns:
            clean_df['GTDB_Species_Count'] = df['gtdb_species_count'].fillna(0).astype(int)
        if 'gtdb_genome_count' in df.columns:
            clean_df['GTDB_Genome_Count'] = df['gtdb_genome_count'].fillna(0).astype(int)
        if 'Coverage_Percentage' in df.columns:
            clean_df['Coverage_Percentage'] = df['Coverage_Percentage']

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
            # Matched taxa = taxa that exist in BOTH census and GTDB (in_both = 1)
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
                # Apply candidat filter if enabled
                if 'EXCLUDE_CANDIDAT_TAXA' in globals() and EXCLUDE_CANDIDAT_TAXA:
                    original_count = len(df)
                    df_filtered = df[~df['taxon_name'].str.lower().str.startswith('candidat', na=False)]
                    filtered_count = len(df_filtered)
                    if original_count != filtered_count:
                        logger.info(f"Candidatus filter: {level} - removed {original_count - filtered_count} entries, kept {filtered_count}")
                    df = df_filtered

                output_file = self.output_dir / f"16s_gtdb_merged_{level}.csv"
                df.to_csv(output_file, index=False)
                logger.info(f"Saved {level} merged data: {output_file}")

        # Save summary statistics
        summary_file = self.output_dir / "16s_gtdb_merger_summary.csv"
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
        log_file = self.output_dir / "16s_gtdb_merger_analysis.txt"

        with open(log_file, 'w') as f:
            f.write("16S Census + GTDB Merger Analysis Report\n")
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

                # Top represented taxa
                if 'representation_bias' in df.columns:
                    top_represented = df.nlargest(10, 'representation_bias')
                    f.write(f"\nTop 10 over-represented in GTDB:\n")
                    for _, row in top_represented.iterrows():
                        if not pd.isna(row['representation_bias']):
                            f.write(f"  {row['taxon_name']}: {row['representation_bias']:.2f}\n")

                    bottom_represented = df.nsmallest(10, 'representation_bias')
                    f.write(f"\nTop 10 under-represented in GTDB:\n")
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
        logger.info("Starting 16S census + GTDB merger...")

        # Load data
        self.load_gtdb_data()
        self.load_census_data()
        self.calculate_total_sequences()  # Calculate total sequences for Size_Percentage

        # Merge data for each taxonomic level
        for level in self.taxonomic_levels:
            merged_df = self.merge_taxonomic_level(level)
            if not merged_df.empty:
                self.merged_data[level] = merged_df

        # Generate summary statistics
        summary_stats = self.generate_summary_statistics(self.merged_data)

        # Save results
        self.save_results(self.merged_data, summary_stats)

        logger.info("16S census + GTDB merger completed successfully!")

        return self.merged_data, summary_stats


def main():
    """
    Main function to run the 16S + GTDB merger.
    """
    # Get the script directory and workspace root
    script_dir = Path(__file__).parent
    workspace_root = script_dir.parent.parent.parent.parent.parent / "parse_repaa_table"  # Go up to metadata_proj/parse_repaa_table

    # Default paths relative to workspace root
    gtdb_data_dir = workspace_root / "gtdb_parse" / "csv_gtdb"
    census_data_dir = workspace_root / "16S_censusparse" / "csv_16S"
    output_dir = script_dir / "merged_output"

    # Allow command line arguments
    if len(sys.argv) >= 2:
        gtdb_data_dir = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        census_data_dir = Path(sys.argv[2])
    if len(sys.argv) >= 4:
        output_dir = Path(sys.argv[3])

    # Create and run merger
    merger = SixteenSGTDBMerger(str(gtdb_data_dir), str(census_data_dir), str(output_dir))
    _, summary_stats = merger.run_merger()

    # Print summary to console
    print("\n16S + GTDB MERGER SUMMARY")
    print("=" * 30)
    print(summary_stats.to_string(index=False))


if __name__ == "__main__":
    main()
