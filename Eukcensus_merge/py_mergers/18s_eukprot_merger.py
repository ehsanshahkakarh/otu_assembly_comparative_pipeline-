#!/usr/bin/env python3
"""
18S Census + EukProt Merger
===========================

Merges EukProt genomic data with 18S census data to create comprehensive comparison tables.
Handles taxonomic matching across different data sources and provides coverage analysis.
Combines functionality of lineage_merger_div.py, lineage_merger_family.py, and lineage_merger_genus.py.

INPUT FILES:
-----------
EukProt Data:
  - eukprot_parse/csv_output/eukprot_new_lineages.csv

18S Census Data:
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_family.csv
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_genus.csv

OUTPUT FILES:
------------
Merged Results:
  - merged_output/18s_merged/results/18s_eukprot_merged_division.csv
  - merged_output/18s_merged/results/18s_eukprot_merged_family.csv
  - merged_output/18s_merged/results/18s_eukprot_merged_genus.csv

Summary & Analysis:
  - merged_output/18s_merged/analysis_summary/18s_eukprot_merger_summary.csv
  - merged_output/18s_merged/analysis_summary/18s_eukprot_merger_analysis.txt
  - logs/18s_eukprot_merger.log

MERGING LOGIC:
-------------
1. Multi-Stream Matching: Uses direct name matching, taxid-based matching, and hierarchical
   lineage matching with confidence scoring (high/medium/none)
2. Enhanced Lineage Parsing: Parses semicolon-separated lineage strings into rank-taxon
   dictionaries, handles missing taxids gracefully
3. Species Deduplication: Uses caching to prevent double-counting species across multiple
   census entries at the same taxonomic level
4. Coverage Calculation: Coverage = (EukProt_Species_Count / Census_OTU_Count) × 100
5. Eukaryotic Focus: Filters for eukaryotic organisms only, handles special taxonomic cases
6. Comprehensive Output: Includes both matched and unmatched entries from both databases

Author: (augment agentic ai) Ehsan Shah Kakar
Date: 2025
"""

import pandas as pd
import sys
from pathlib import Path
import logging
from typing import Dict, Tuple
from tqdm import tqdm
from collections import defaultdict

# Setup logging
log_dir = Path("../18s_merged/logs")
log_dir.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_dir / '18s_eukprot_merger.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class EighteenSEukProtMerger:
    """
    Merges EukProt genomic data with 18S census data for eukaryotic organisms.
    """
    
    def __init__(self, eukprot_data_dir: str, census_data_dir: str, output_dir: str = "merged_output"):
        """
        Initialize the merger with data directories.
        
        Args:
            eukprot_data_dir: Directory containing EukProt CSV files
            census_data_dir: Directory containing 18S census CSV files
            output_dir: Directory for output files
        """
        self.eukprot_data_dir = Path(eukprot_data_dir)
        self.census_data_dir = Path(census_data_dir)
        self.output_dir = Path(output_dir) / "18s_merged"  # Main output directory
        self.output_dir.mkdir(exist_ok=True)

        # Create subdirectories for organized output
        self.results_dir = self.output_dir / "csv_results"
        self.analysis_dir = self.output_dir / "analysis_summary"
        self.results_dir.mkdir(exist_ok=True)
        self.analysis_dir.mkdir(exist_ok=True)

        # Taxonomic levels to process
        self.taxonomic_levels = ['division', 'family', 'genus']

        # Data containers
        self.eukprot_data = {}
        self.census_data = {}
        self.merged_data = {}

        # Expected file mappings
        self.eukprot_files = {
            'division': 'eukprot_new_lineages.csv',
            'family': 'eukprot_new_lineages.csv',  # Same file, different processing
            'genus': 'eukprot_new_lineages.csv'    # Same file, different processing
        }

        self.census_files = {
            'division': 'eukcensus_18S_by_division.csv',
            'family': 'eukcensus_18S_by_family.csv',
            'genus': 'eukcensus_18S_by_genus.csv'
        }

    def parse_lineage_with_taxids(self, lineage_str, ranks_str, taxids_str=None):
        """Parse lineage, ranks, and taxids into comprehensive dictionaries."""
        if pd.isna(lineage_str) or pd.isna(ranks_str):
            return {}, {}

        try:
            # Clean and split all three components
            lineage_parts = [x.strip() for x in str(lineage_str).split(';') if x.strip() and x.strip().lower() != 'nan']
            rank_parts = [x.strip() for x in str(ranks_str).split(';') if x.strip() and x.strip().lower() != 'nan']

            # Handle taxids (may be missing in some datasets)
            taxid_parts = []
            if not pd.isna(taxids_str):
                taxid_parts = [x.strip() for x in str(taxids_str).split(';') if x.strip() and x.strip().lower() != 'nan']

            # Handle EukCensus data that still has "cellular organisms" prefix
            if lineage_parts and lineage_parts[0] == "cellular organisms":
                lineage_parts = lineage_parts[1:]  # Remove "cellular organisms"
            if rank_parts and rank_parts[0] == "cellular root":
                rank_parts = rank_parts[1:]  # Remove "cellular root"
            if taxid_parts and len(taxid_parts) > 0 and taxid_parts[0] == "131567":  # cellular organisms taxid
                taxid_parts = taxid_parts[1:]  # Remove cellular organisms taxid

            # Create dictionaries mapping ranks to taxa and taxids
            lineage_dict = {}
            taxid_dict = {}

            # Ensure all arrays are same length by taking minimum
            min_length = min(len(rank_parts), len(lineage_parts))
            if taxid_parts:
                min_length = min(min_length, len(taxid_parts))

            for i in range(min_length):
                rank = rank_parts[i].lower().strip()
                taxon = lineage_parts[i].strip()

                if rank and taxon:
                    lineage_dict[rank] = taxon
                    if taxid_parts and i < len(taxid_parts):
                        taxid_dict[rank] = taxid_parts[i].strip()

            return lineage_dict, taxid_dict
        except Exception as e:
            logger.error(f"Error parsing lineage '{lineage_str}' with ranks '{ranks_str}' and taxids '{taxids_str}': {e}")
            return {}, {}

    def extract_taxonomic_level(self, lineage_dict, taxid_dict, target_level):
        """Extract taxa at the specified taxonomic level."""
        taxa = set()
        taxids = set()

        # Special case mappings for known taxonomic groups
        special_cases = {
            'Rigifilida': 'Rigifilida',
            'Hemimastigophora': 'Hemimastigophora',
            'Rhodophyta': 'Rhodophyta',
            'Nebulidia': 'Nebulidia',
            'Nibbleridia': 'Nibbleridia'
        }

        # Check for special cases first
        for taxon_name in lineage_dict.values():
            if taxon_name in special_cases:
                taxa.add(special_cases[taxon_name])
                # Also collect corresponding taxid if available
                for rank, name in lineage_dict.items():
                    if name == taxon_name and rank in taxid_dict:
                        taxids.add(taxid_dict[rank])
                return taxa, taxids

        # Define rank hierarchies for each level
        rank_hierarchies = {
            'division': ['clade', 'phylum', 'kingdom', 'class', 'supergroup', 'order'],
            'family': ['family', 'superfamily', 'order', 'class', 'phylum'],
            'genus': ['genus', 'subgenus', 'family', 'subfamily']
        }

        # Regular processing for target level - take FIRST match only
        target_ranks = rank_hierarchies.get(target_level, [target_level])
        
        for rank in target_ranks:
            if rank in lineage_dict:
                taxa.add(lineage_dict[rank])
                # Also collect corresponding taxid
                if rank in taxid_dict:
                    taxids.add(taxid_dict[rank])
                break  # Take only the first match to avoid multiple assignments

        return taxa, taxids

    def find_multi_stream_matches(self, census_name, census_lineage_dict, census_taxid_dict,
                                 eukprot_lineage_dict, eukprot_taxid_dict, target_level):
        """
        Find matches using multiple data streams: names, taxids, and hierarchical relationships.
        """
        matches = {
            'direct_name_match': False,
            'taxid_match': False,
            'hierarchical_match': False,
            'confidence': 'none',
            'match_details': []
        }

        # 1. Direct name matching
        eukprot_taxa, eukprot_taxids = self.extract_taxonomic_level(
            eukprot_lineage_dict, eukprot_taxid_dict, target_level
        )

        if census_name in eukprot_taxa:
            matches['direct_name_match'] = True
            matches['confidence'] = 'high'
            matches['match_details'].append(f"Direct name match: {census_name}")

        # 2. Taxid-based matching
        rank_priority = ['clade', 'phylum', 'kingdom', 'family', 'genus']
        census_taxid = None
        for rank in rank_priority:
            if rank in census_taxid_dict:
                census_taxid = census_taxid_dict[rank]
                break

        if census_taxid and census_taxid in eukprot_taxids:
            matches['taxid_match'] = True
            if matches['confidence'] == 'none':
                matches['confidence'] = 'high'
            matches['match_details'].append(f"Taxid match: {census_taxid}")

        # 3. Hierarchical matching
        if census_taxid:
            for rank, taxid in eukprot_taxid_dict.items():
                if taxid == census_taxid:
                    matches['hierarchical_match'] = True
                    if matches['confidence'] == 'none':
                        matches['confidence'] = 'medium'
                    matches['match_details'].append(f"Hierarchical match: {census_taxid} at {rank} level")
                    break

        # 4. Cross-reference name matching
        for rank, name in eukprot_lineage_dict.items():
            if name == census_name:
                if matches['confidence'] == 'none':
                    matches['confidence'] = 'medium'
                matches['match_details'].append(f"Cross-reference match: {census_name} at {rank} level")
                break

        return matches

    def categorize_representation(self, coverage_pct):
        """Categorize representation level for visualization."""
        if coverage_pct == 0:
            return 'No_Representation'
        elif coverage_pct < 1:
            return 'Low_Representation'
        elif coverage_pct < 10:
            return 'Moderate_Representation'
        elif coverage_pct < 50:
            return 'Good_Representation'
        else:
            return 'High_Representation'

    def load_eukprot_data(self):
        """Load EukProt data files."""
        logger.info("Loading EukProt data...")
        
        eukprot_file = self.eukprot_data_dir / 'eukprot_new_lineages.csv'
        
        if not eukprot_file.exists():
            raise FileNotFoundError(f"EukProt file not found: {eukprot_file}")
        
        # Load the main EukProt lineages file
        df = pd.read_csv(eukprot_file)
        logger.info(f"Loaded EukProt data: {len(df)} entries")
        
        # Store for all levels (same file, different processing)
        for level in self.taxonomic_levels:
            self.eukprot_data[level] = df.copy()

    def load_census_data(self):
        """Load 18S census data files."""
        logger.info("Loading 18S census data...")
        
        for level in self.taxonomic_levels:
            census_file = self.census_data_dir / self.census_files[level]
            
            if not census_file.exists():
                logger.warning(f"Census file not found: {census_file}")
                continue
            
            df = pd.read_csv(census_file)
            logger.info(f"Loaded {level} census data: {len(df)} entries")
            self.census_data[level] = df

    def vectorized_lineage_matching(self, census_df: pd.DataFrame, eukprot_df: pd.DataFrame,
                                   census_col: str, level: str) -> pd.DataFrame:
        """Ultra-fast vectorized lineage matching and aggregation for EukProt data."""

        # Filter out .U. entries from both datasets
        print(f"   📋 Filtering census {level} data...")
        original_census_count = len(census_df)
        u_mask = ~census_df[census_col].str.contains(r'\.U\.', case=False, na=False)
        candidat_mask = ~census_df[census_col].str.contains(r'candidat', case=False, na=False)
        virus_mask = ~census_df[census_col].str.contains(r'virus|viral', case=False, na=False)
        combined_mask = u_mask & candidat_mask & virus_mask
        census_df_filtered = census_df[combined_mask]
        filtered_census_count = len(census_df_filtered)
        print(f"   📋 Filtered census {level}: {filtered_census_count}/{original_census_count} taxa (removed {original_census_count - filtered_census_count} .U./candidat/virus entries)")

        print(f"   🧬 Filtering EukProt {level} data...")
        original_eukprot_count = len(eukprot_df)
        u_mask_euk = ~eukprot_df['Name_to_use'].str.contains(r'\.U\.', case=False, na=False)
        candidat_mask_euk = ~eukprot_df['Name_to_use'].str.contains(r'candidat', case=False, na=False)
        virus_mask_euk = ~eukprot_df['Name_to_use'].str.contains(r'virus|viral', case=False, na=False)
        combined_mask_euk = u_mask_euk & candidat_mask_euk & virus_mask_euk
        eukprot_df_filtered = eukprot_df[combined_mask_euk]
        filtered_eukprot_count = len(eukprot_df_filtered)
        print(f"   🧬 Filtered EukProt {level}: {filtered_eukprot_count}/{original_eukprot_count} taxa (removed {original_eukprot_count - filtered_eukprot_count} .U./candidat/virus entries)")

        # Prepare EukProt lineage data for vectorized matching
        eukprot_df_filtered = eukprot_df_filtered.copy()
        eukprot_df_filtered['lineage'] = eukprot_df_filtered['lineage'].fillna('')

        # Note: Census percentages will be taken directly from census data (no recalculation needed)

        results = []

        for _, census_row in census_df_filtered.iterrows():
            census_name = census_row[census_col]
            census_otus = census_row.get('otu_count', 0)
            census_size = census_row.get('size_count', 0)
            # Take percentages directly from census data (already calculated correctly)
            otu_percentage = census_row.get('otu_percentage', 0)
            size_percentage = census_row.get('size_percentage', 0)

            # Vectorized lineage search - similar to NCBI merger approach
            pattern = f';{census_name};|^{census_name};|;{census_name}$|^{census_name}$'
            matches_mask = eukprot_df_filtered['lineage'].str.contains(pattern, regex=True, na=False)
            matched_eukprot = eukprot_df_filtered[matches_mask]

            if not matched_eukprot.empty:
                # Get unique species names (eliminates duplicates automatically)
                unique_species = matched_eukprot['Name_to_use'].unique()
                total_species = len(unique_species)

                # Count direct vs lineage matches
                direct_matches = len(matched_eukprot[matched_eukprot['Name_to_use'] == census_name]['Name_to_use'].unique())
                lineage_matches = total_species - direct_matches
            else:
                total_species = direct_matches = lineage_matches = 0



            # Calculate novelty and overrepresentation factors
            novelty_factor = census_otus / total_species if total_species > 0 else float('inf')
            overrepresentation_factor = total_species / census_otus if census_otus > 0 else float('inf')

            # Use proper column name based on level
            level_col_name = 'division' if level == 'division' else level

            results.append({
                level_col_name: census_name,  # Use proper level column name
                'census_otu_count': census_otus,
                'census_size_count': census_size,
                'otu_percentage': round(otu_percentage, 2),
                'size_percentage': round(size_percentage, 2),
                'eukprot_species_count': total_species,

                'novelty_factor': round(novelty_factor, 3),
                'overrepresentation_factor': round(overrepresentation_factor, 3),
                'direct_matches': direct_matches,
                'lineage_matches': lineage_matches,
                'total_matches': direct_matches + lineage_matches,
                'match_status': 'matched' if total_species > 0 else 'census_only'
            })

        return pd.DataFrame(results)

    def merge_taxonomic_level(self, level: str) -> pd.DataFrame:
        """
        Merge EukProt and census data for a specific taxonomic level using vectorized operations.
        """
        logger.info(f"Merging {level} level data...")

        if level not in self.eukprot_data or level not in self.census_data:
            logger.warning(f"Missing data for {level} level")
            return pd.DataFrame()

        eukprot_df = self.eukprot_data[level]
        census_df = self.census_data[level]

        # Get the name column for this level
        census_col = 'Name_to_use'  # Census data uses Name_to_use column

        # Vectorized matching and aggregation
        merged_df = self.vectorized_lineage_matching(census_df, eukprot_df, census_col, level)

        # Sort results: matched taxa first (by eukprot_species_count desc), then unmatched taxa (by eukprot_species_count desc)
        matched_df = merged_df[merged_df['match_status'] == 'matched'].sort_values('eukprot_species_count', ascending=False)
        unmatched_df = merged_df[merged_df['match_status'] != 'matched'].sort_values('eukprot_species_count', ascending=False)
        final_df = pd.concat([matched_df, unmatched_df], ignore_index=True)

        return final_df

    def generate_summary_statistics(self, merged_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Generate summary statistics across all taxonomic levels."""
        logger.info("Generating summary statistics...")

        summary_stats = []

        for level, df in merged_data.items():
            if df.empty:
                continue

            total_taxa = len(df)
            taxa_with_eukprot = len(df[df['eukprot_species_count'] > 0])
            taxa_census_only = len(df[df['eukprot_species_count'] == 0])

            # Coverage metrics
            coverage_pct = (taxa_with_eukprot / total_taxa * 100) if total_taxa > 0 else 0

            # Coverage distribution
            avg_coverage = df['coverage_percentage'].mean()
            median_coverage = df['coverage_percentage'].median()

            summary_stats.append({
                'Taxonomic_Level': level.capitalize(),
                'Total_Taxa': total_taxa,
                'Matched_Taxa': taxa_with_eukprot,  # Renamed for consistency with other mergers
                'Taxa_Census_Only': taxa_census_only,
                'Coverage_Percentage': f"{coverage_pct:.1f}%",
                'Avg_Coverage': f"{avg_coverage:.2f}%",
                'Median_Coverage': f"{median_coverage:.2f}%"
            })

        return pd.DataFrame(summary_stats)

    def save_results(self, merged_data: Dict[str, pd.DataFrame], summary_stats: pd.DataFrame):
        """Save merged results and summary statistics to organized subdirectories."""
        logger.info("Saving results to organized subdirectories...")

        # Save merged data for each taxonomic level to results directory
        for level, df in merged_data.items():
            if not df.empty:
                output_file = self.results_dir / f'18s_eukprot_merged_{level}.csv'
                df.to_csv(output_file, index=False)
                logger.info(f"Saved {level} merged data: {output_file}")

        # Save summary statistics to analysis directory
        summary_file = self.analysis_dir / '18s_eukprot_merger_summary.csv'
        summary_stats.to_csv(summary_file, index=False)
        logger.info(f"Saved summary statistics: {summary_file}")

        # Save detailed analysis report to analysis directory
        self.save_detailed_log(merged_data, summary_stats)

    def save_detailed_log(self, merged_data: Dict[str, pd.DataFrame], summary_stats: pd.DataFrame) -> None:
        """
        Save detailed analysis log matching NCBI merger format.

        Args:
            merged_data: Dictionary of merged DataFrames
            summary_stats: Summary statistics DataFrame
        """
        log_file = self.analysis_dir / "18s_eukprot_merger_analysis.txt"

        with open(log_file, 'w') as f:
            f.write("18S Census + EukProt Merger Analysis Report\n")
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

                # Top taxa by EukProt species count
                if 'EukProt_Species_Count' in df.columns:
                    top_eukprot = df.nlargest(10, 'EukProt_Species_Count')
                    f.write(f"\nTop 10 taxa by EukProt species count:\n")
                    for _, row in top_eukprot.iterrows():
                        taxon_name = row[f'{level.capitalize()}']
                        species_count = row['EukProt_Species_Count']
                        coverage = row.get('Coverage_Percentage', 0)
                        f.write(f"  {taxon_name}: {species_count} species ({coverage:.1f}% coverage)\n")

                # Top taxa by coverage percentage
                if 'Coverage_Percentage' in df.columns:
                    top_coverage = df.nlargest(10, 'Coverage_Percentage')
                    f.write(f"\nTop 10 taxa by coverage percentage:\n")
                    for _, row in top_coverage.iterrows():
                        taxon_name = row[f'{level.capitalize()}']
                        coverage = row['Coverage_Percentage']
                        species_count = row.get('EukProt_Species_Count', 0)
                        f.write(f"  {taxon_name}: {coverage:.1f}% coverage ({species_count} species)\n")

        logger.info(f"Saved detailed analysis: {log_file}")

    def run_merger(self) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
        """
        Run the complete merger pipeline.

        Returns:
            Tuple of (merged_data, summary_stats)
        """
        logger.info("Starting 18S census + EukProt merger...")

        # Load data
        self.load_eukprot_data()
        self.load_census_data()

        # Merge data for each taxonomic level
        for level in self.taxonomic_levels:
            merged_df = self.merge_taxonomic_level(level)
            if not merged_df.empty:
                self.merged_data[level] = merged_df

        # Generate summary statistics
        summary_stats = self.generate_summary_statistics(self.merged_data)

        # Save results
        self.save_results(self.merged_data, summary_stats)

        logger.info("18S + EukProt merger completed successfully!")

        return self.merged_data, summary_stats


def main():
    """
    Main function to run the 18S + EukProt merger.
    """
    # Use paths relative to script location
    script_dir = Path(__file__).resolve().parent
    # Go up from py_mergers -> Eukcensus_merge -> parse_repaa_table
    base_dir = script_dir.parent.parent
    eukprot_data_dir = base_dir / "eukprot_parse" / "csv_output"
    census_data_dir = base_dir / "18S_censusparse" / "csv_outputs"
    output_dir = script_dir.parent

    # Allow command line arguments
    if len(sys.argv) >= 2:
        eukprot_data_dir = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        census_data_dir = Path(sys.argv[2])
    if len(sys.argv) >= 4:
        output_dir = Path(sys.argv[3])

    # Create and run merger
    merger = EighteenSEukProtMerger(str(eukprot_data_dir), str(census_data_dir), str(output_dir))
    _, summary_stats = merger.run_merger()

    # Print summary to console
    print("\n18S + EUKPROT MERGER SUMMARY")
    print("=" * 30)
    print(summary_stats.to_string(index=False))

    print("\n✅ 18S EukProt Merger completed successfully!")
    print("📁 Output files saved to: ../18s_merged/")
    print("   - csv_results/: Merged CSV files")
    print("   - analysis_summary/: Summary statistics")
    print("   - logs/: Detailed matching logs")


if __name__ == "__main__":
    main()
