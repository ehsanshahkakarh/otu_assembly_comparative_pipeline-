#!/usr/bin/env python3
"""
16S Census + GTDB Merger - Vectorized Version
=============================================

Fast, vectorized merger using pandas operations for lineage-based matching.
Matches the perfect 18S merger structure and performance.

INPUT FILES:
-----------
GTDB Data:
  - gtdb_parse/csv_gtdb/gtdb_phylum_counts.csv
  - gtdb_parse/csv_gtdb/gtdb_family_counts.csv  
  - gtdb_parse/csv_gtdb/gtdb_genus_counts.csv

16S Census Data:
  - 16S_censusparse/csv_16S/eukcensus16S_by_division.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_family.csv
  - 16S_censusparse/csv_16S/eukcensus16S_by_genus.csv

OUTPUT FILES:
------------
  - merged_output/16s_merged/results/16s_gtdb_merged_clean_*.csv
  - merged_output/16s_merged/analysis_summary/16s_gtdb_merger_clean_summary.csv

LOGIC:
------
1. Vectorized Lineage Matching: Use pandas string operations to find census taxa in GTDB lineages
2. Vectorized Aggregation: Use pandas groupby and sum operations for count aggregation
3. Clean Output: Same format as 18S version
4. Prokaryotic Focus: Filters for Bacteria and Archaea domains

Author: Vectorized Merger Team
Date: 2025
"""

import pandas as pd
from pathlib import Path
import logging
import sys

def setup_logging():
    """Setup minimal logging."""
    # Use relative path to 16s_merged logs directory
    logs_dir = Path("../16s_merged/logs")
    logs_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(logs_dir / '16s_gtdb_merger.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

class SixteenSGTDBMerger:
    """Vectorized 16S Census + GTDB Merger matching 18S structure."""
    
    def __init__(self):
        """Initialize vectorized merger."""
        self.logger = setup_logging()
        self.gtdb_data = {}
        self.census_data = {}
        self.matches_log = []
        self.levels = ['phylum', 'family', 'genus']

        # Set up paths relative to script location
        script_dir = Path(__file__).resolve().parent
        # Go up from py_mergers -> Eukcensus_merge -> parse_repaa_table
        self.base_dir = script_dir.parent.parent
        self.output_dir = script_dir.parent  # Eukcensus_merge directory

        # Total counts for percentage calculations
        self.total_gtdb_genomes = 0
        self.total_gtdb_species = 0
    
    def load_data(self):
        """Load all required data files."""
        self.logger.info("Loading GTDB and Census data...")
        
        # Load GTDB data
        for level in self.levels:
            gtdb_file = self.base_dir / "gtdb_parse/csv_gtdb" / f"gtdb_{level}_counts.csv"
            if gtdb_file.exists():
                df = pd.read_csv(gtdb_file)
                
                # Filter for prokaryotes (Bacteria + Archaea) - already prokaryotic in GTDB
                if 'domain' in df.columns:
                    df = df[df['domain'].isin(['Bacteria', 'Archaea'])]
                
                # Remove candidat entries for clean output
                if level in df.columns:
                    df = df[~df[level].str.contains('candidat', case=False, na=False)]
                
                self.gtdb_data[level] = df

                # Calculate level-specific totals for percentage calculations (from prokaryote-only data)
                prokaryote_only_df = df[df['domain'].isin(['Bacteria', 'Archaea'])] if 'domain' in df.columns else df
                total_counts = prokaryote_only_df[f'{level}_counts'].sum()
                setattr(self, f'total_gtdb_counts_{level}', total_counts)
                print(f"   📊 GTDB {level} prokaryote database totals: {total_counts:,} genomes")

                self.logger.info(f"Loaded GTDB {level} data: {len(df)} entries")
            else:
                self.logger.error(f"GTDB file not found: {gtdb_file}")
                
        # Load Census data (same as NCBI merger)
        census_files = {
            'phylum': self.base_dir / "16S_censusparse/csv_16S/eukcensus16S_by_division.csv",
            'family': self.base_dir / "16S_censusparse/csv_16S/eukcensus16S_by_family.csv",
            'genus': self.base_dir / "16S_censusparse/csv_16S/eukcensus16S_by_genus.csv"
        }
        
        for level, file_path in census_files.items():
            if file_path.exists():
                df = pd.read_csv(file_path)
                
                # Rename columns to match expected format
                df = df.rename(columns={
                    'otu_count': 'census_otu_count',
                    'size_count': 'census_size_count'
                })
                
                # Filter for prokaryotes - remove eukaryotic entries
                if len(df.columns) > 0:
                    first_col = df.columns[0]
                    # Remove obvious eukaryotic taxa
                    euk_patterns = ['Eukaryota', 'Opisthokonta', 'Streptophyta', 'Chlorophyta', 
                                   'Stramenopiles', 'Alveolata', 'Rhizaria', 'Discoba']
                    for pattern in euk_patterns:
                        df = df[~df[first_col].str.contains(pattern, case=False, na=False)]
                
                self.census_data[level] = df
                self.logger.info(f"Loaded Census {level} data: {len(df)} entries")
            else:
                self.logger.error(f"Census file not found: {file_path}")

        # GTDB totals are calculated during data loading
    


    def vectorized_lineage_matching(self, census_df: pd.DataFrame, gtdb_df: pd.DataFrame,
                                   census_col: str, gtdb_col: str, level: str) -> pd.DataFrame:
        """Ultra-fast vectorized lineage matching and aggregation."""

        # Prepare lineage data
        gtdb_df = gtdb_df.copy()

        results = []
        for _, census_row in census_df.iterrows():
            census_name = census_row[census_col]

            # Vectorized lineage matching - escape special regex characters
            import re
            escaped_name = re.escape(census_name)
            pattern = f';{escaped_name};|^{escaped_name};|;{escaped_name}$|^{escaped_name}$'
            matches_mask = gtdb_df['lineage'].str.contains(pattern, regex=True, na=False)
            matched_gtdb = gtdb_df[matches_mask]

            if not matched_gtdb.empty:
                # Vectorized aggregation - GTDB only has counts, not separate genome/species
                total_counts = matched_gtdb[f'{level}_counts'].sum()

                # Calculate percentages using correct total denominators from GTDB database
                total_gtdb_counts = getattr(self, f'total_gtdb_counts_{level}', 1)

                genome_pct_db = (total_counts / total_gtdb_counts * 100) if total_gtdb_counts > 0 else 0
                species_pct = genome_pct_db  # GTDB doesn't separate genome/species, so use same value

                # Count match types
                direct_matches = int(census_name in matched_gtdb[gtdb_col].values)
                lineage_matches = len(matched_gtdb) - direct_matches

                # Log the matches
                for _, match_row in matched_gtdb.iterrows():
                    match_type = 'direct_match' if match_row[gtdb_col] == census_name else 'lineage_match'
                    self.matches_log.append({
                        'level': level,
                        'census_taxon': census_name,
                        'gtdb_taxon': match_row[gtdb_col],
                        'match_type': match_type,
                        'gtdb_counts': match_row[f'{level}_counts'],
                        'lineage': match_row['lineage']
                    })
            else:
                total_counts = direct_matches = lineage_matches = 0
                genome_pct_db = species_pct = 0



            # Calculate novelty and overrepresentation factors
            census_otus = census_row['census_otu_count']
            novelty_factor = census_otus / total_counts if total_counts > 0 else float('inf')
            overrepresentation_factor = total_counts / census_otus if census_otus > 0 else float('inf')

            results.append({
                level: census_name,
                'census_otu_count': census_row['census_otu_count'],
                'census_size_count': census_row['census_size_count'],
                'otu_percentage': census_row['otu_percentage'],
                'size_percentage': census_row['size_percentage'],
                'gtdb_genome_count': total_counts,  # GTDB counts represent genomes
                'gtdb_species_count': total_counts,  # GTDB doesn't separate, so use same value
                'isolate_count': 0,  # GTDB doesn't have isolate classification
                'genome_pct_db': round(genome_pct_db, 2),
                'species_pct': round(species_pct, 2),
                'isolate_percentage': 0.0,  # GTDB doesn't have isolate classification

                'novelty_factor': round(novelty_factor, 3),
                'overrepresentation_factor': round(overrepresentation_factor, 3),
                'direct_matches': direct_matches,
                'lineage_matches': lineage_matches,
                'total_matches': direct_matches + lineage_matches,
                'match_status': 'matched' if total_counts > 0 else 'census_only'
            })

        return pd.DataFrame(results)

    def add_gtdb_only_taxa(self, merged_df: pd.DataFrame, gtdb_df: pd.DataFrame,
                          census_col: str, gtdb_col: str, level: str) -> pd.DataFrame:
        """Add GTDB-only taxa that weren't matched to any census taxon."""

        # Find all GTDB taxa that were matched
        all_matched = set()
        for _, row in merged_df.iterrows():
            if row['total_matches'] > 0:
                census_name = row[level]  # Use level name since we changed the column structure
                pattern = f';{census_name};|^{census_name};|;{census_name}$|^{census_name}$'
                matches_mask = gtdb_df['lineage'].str.contains(pattern, regex=True, na=False)
                matched_taxa = gtdb_df.loc[matches_mask, gtdb_col].tolist()
                all_matched.update(matched_taxa)

        # Find unmatched GTDB taxa
        unmatched_gtdb = gtdb_df[~gtdb_df[gtdb_col].isin(all_matched)]

        gtdb_only_rows = []

        for _, gtdb_row in unmatched_gtdb.iterrows():
            gtdb_name = gtdb_row[gtdb_col]
            counts = gtdb_row[f'{level}_counts']

            # Calculate percentages using correct total denominators from GTDB database
            total_gtdb_counts = getattr(self, f'total_gtdb_counts_{level}', 1)

            genome_pct_db = (counts / total_gtdb_counts * 100) if total_gtdb_counts > 0 else 0
            species_pct = genome_pct_db  # GTDB doesn't separate genome/species

            if counts > 0:  # Only add taxa with actual data
                gtdb_only_rows.append({
                    level: gtdb_name,  # Use level name instead of census_col
                    'census_otu_count': 0,
                    'census_size_count': 0,
                    'otu_percentage': 0,
                    'size_percentage': 0,
                    'gtdb_genome_count': counts,
                    'gtdb_species_count': counts,  # GTDB doesn't separate, so use same value
                    'isolate_count': 0,
                    'genome_pct_db': round(genome_pct_db, 2),
                    'species_pct': round(species_pct, 2),
                    'isolate_percentage': 0.0,
                    'coverage_percentage': 0,
                    'direct_matches': 0,
                    'lineage_matches': 0,
                    'total_matches': 0,
                    'match_status': 'gtdb_only'
                })

        if gtdb_only_rows:
            gtdb_only_df = pd.DataFrame(gtdb_only_rows)
            return pd.concat([merged_df, gtdb_only_df], ignore_index=True)

        return merged_df

    def process_level(self, level: str) -> pd.DataFrame:
        """Process one taxonomic level with vectorized operations."""
        gtdb_df = self.gtdb_data[level]
        census_df = self.census_data[level]

        # Column names
        gtdb_col = level
        census_col = census_df.columns[0]

        # Vectorized matching and aggregation
        merged_df = self.vectorized_lineage_matching(census_df, gtdb_df, census_col, gtdb_col, level)

        # Add GTDB-only taxa
        final_df = self.add_gtdb_only_taxa(merged_df, gtdb_df, census_col, gtdb_col, level)

        # Sort results: matched taxa first (by gtdb_genome_count desc), then unmatched taxa (by gtdb_genome_count desc)
        matched_df = final_df[final_df['match_status'] == 'matched'].sort_values('gtdb_genome_count', ascending=False)
        unmatched_df = final_df[final_df['match_status'] != 'matched'].sort_values('gtdb_genome_count', ascending=False)
        final_df = pd.concat([matched_df, unmatched_df], ignore_index=True)

        return final_df

    def save_matches_log(self):
        """Save detailed matches log to CSV file."""
        if self.matches_log:
            matches_df = pd.DataFrame(self.matches_log)
            output_dir = self.output_dir / "16s_merged/logs"
            output_dir.mkdir(parents=True, exist_ok=True)
            matches_df.to_csv(output_dir / "16s_gtdb_detailed_matches.csv", index=False)
            self.logger.info(f"Saved detailed matches log: {len(matches_df)} entries")

    def run_merger(self):
        """Run the complete merger process."""
        self.logger.info("Starting 16S GTDB vectorized merger...")

        # Load data
        self.load_data()

        # Create output directories
        results_dir = self.output_dir / "16s_merged/csv_results"
        summary_dir = self.output_dir / "16s_merged/analysis_summary"
        results_dir.mkdir(parents=True, exist_ok=True)
        summary_dir.mkdir(parents=True, exist_ok=True)

        # Process each level
        summary_data = []
        for level in self.levels:
            if level in self.gtdb_data and level in self.census_data:
                self.logger.info(f"Processing {level} level...")

                merged_df = self.process_level(level)

                # Save results
                output_file = results_dir / f"16s_gtdb_merged_clean_{level}.csv"
                merged_df.to_csv(output_file, index=False)
                self.logger.info(f"Saved {level} results: {len(merged_df)} entries")

                # Collect summary stats
                matched_count = len(merged_df[merged_df['match_status'] == 'matched'])
                census_only_count = len(merged_df[merged_df['match_status'] == 'census_only'])
                gtdb_only_count = len(merged_df[merged_df['match_status'] == 'gtdb_only'])

                summary_data.append({
                    'level': level,
                    'total_entries': len(merged_df),
                    'matched_taxa': matched_count,
                    'census_only_taxa': census_only_count,
                    'gtdb_only_taxa': gtdb_only_count,
                    'match_rate': f"{matched_count/len(merged_df)*100:.1f}%" if len(merged_df) > 0 else "0%"
                })

        # Save summary
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(summary_dir / "16s_gtdb_merger_clean_summary.csv", index=False)
            self.logger.info("Saved merger summary")

        # Save matches log
        self.save_matches_log()

        self.logger.info("16S GTDB merger completed successfully!")

def main():
    """Main execution function."""
    logger = setup_logging()

    try:
        merger = SixteenSGTDBMerger()
        merger.run_merger()

        print("\n✅ 16S GTDB Merger completed successfully!")
        print("📁 Output files saved to: ../16s_merged/")
        print("   - csv_results/: Merged CSV files")
        print("   - analysis_summary/: Summary statistics")
        print("   - logs/: Detailed matching logs")

    except Exception as e:
        logger.error(f"Merger failed: {e}")
        print(f"\n❌ Merger failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
