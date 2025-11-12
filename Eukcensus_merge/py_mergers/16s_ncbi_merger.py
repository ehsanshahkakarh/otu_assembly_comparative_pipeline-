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
1. Enhanced Matching Strategy:
   - Primary: Taxid-based matching (most reliable for handling name changes/synonyms)
   - Secondary: Vectorized lineage matching using pandas string operations
   - Combines both methods to maximize match accuracy
2. Vectorized Aggregation: Use pandas groupby and sum operations for count aggregation
3. Isolate Analysis: Vectorized isolate percentage calculations
4. Clean Output: Same format as 18S version with additional taxid match tracking
5. Prokaryotic Focus: Filters for Bacteria and Archaea domains, includes organellar sequences
6. Intra-database Aggregation: Handles cases where same taxid appears with different names

Author: Vectorized Merger Team
Date: 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys
import re

def setup_logging():
    """Setup minimal logging."""
    # Use relative path to 16s_merged logs directory
    logs_dir = Path("../16s_merged/logs")
    logs_dir.mkdir(parents=True, exist_ok=True)

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
        """Initialize vectorized merger."""
        self.logger = setup_logging()
        self.ncbi_data = {}
        self.census_data = {}
        self.isolate_data = {}
        self.matches_log = []
        self.levels = ['phylum', 'family', 'genus']

        # Set up paths relative to script location
        script_dir = Path(__file__).resolve().parent
        # Go up from py_mergers -> Eukcensus_merge -> parse_repaa_table
        self.base_dir = script_dir.parent.parent
        self.output_dir = script_dir.parent  # Eukcensus_merge directory

        # Total counts for percentage calculations
        self.total_ncbi_genomes = 0
        self.total_ncbi_species = 0
    
    def load_data(self):
        """Load all required data files."""
        self.logger.info("Loading NCBI and Census data...")

        # Load NCBI data
        for level in self.levels:
            ncbi_file = self.base_dir / "ncbi_parse/csv_ncbi" / f"ncbi_{level}_counts.csv"
            if ncbi_file.exists():
                df = pd.read_csv(ncbi_file)

                # Filter for prokaryotes (Bacteria + Archaea) - 16S specific
                if 'domain' in df.columns:
                    df = df[df['domain'].isin(['Bacteria', 'Archaea'])]

                # Remove candidat and .U. entries for clean output
                if level in df.columns:
                    df = df[~df[level].str.contains('candidat', case=False, na=False)]
                    df = df[~df[level].str.contains(r'\.U\.', case=False, na=False)]

                self.ncbi_data[level] = df

                # Calculate level-specific totals for percentage calculations (from prokaryote-only data)
                prokaryote_only_df = df[df['domain'].isin(['Bacteria', 'Archaea'])] if 'domain' in df.columns else df
                total_genomes = prokaryote_only_df[f'{level}_genome_count'].sum()
                total_species = prokaryote_only_df[f'{level}_species_count'].sum()
                setattr(self, f'total_ncbi_genomes_{level}', total_genomes)
                setattr(self, f'total_ncbi_species_{level}', total_species)
                print(f"   📊 NCBI {level} prokaryote database totals: {total_genomes:,} genomes, {total_species:,} species")

                self.logger.info(f"Loaded NCBI {level} data: {len(df)} entries")
            else:
                self.logger.error(f"NCBI file not found: {ncbi_file}")

        # Load Census data
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

        # Load isolate data
        self.load_isolate_data()
    
    def load_isolate_data(self):
        """Load isolate classification data."""
        for level in self.levels:
            isolate_file = self.base_dir / "ncbi_parse/csv_ncbi" / f"ncbi_{level}_with_accessions.csv"
            if isolate_file.exists():
                df = pd.read_csv(isolate_file)
                
                # Filter for prokaryotes
                if 'domain' in df.columns:
                    df = df[df['domain'].isin(['Bacteria', 'Archaea'])]
                
                # Calculate isolate counts
                if 'genome_source' in df.columns:
                    isolate_counts = df[df['genome_source'] == 'isolate'].groupby(level).size()
                    total_counts = df.groupby(level).size()
                    
                    isolate_df = pd.DataFrame({
                        level: isolate_counts.index,
                        'isolate_count': isolate_counts.values,
                        'total_genomes': total_counts.loc[isolate_counts.index].values
                    })
                    isolate_df['isolate_percentage'] = (isolate_df['isolate_count'] / isolate_df['total_genomes'] * 100).round(2)
                    self.isolate_data[level] = isolate_df
    


    def vectorized_lineage_matching(self, census_df: pd.DataFrame, ncbi_df: pd.DataFrame,
                                   census_col: str, ncbi_col: str, level: str) -> pd.DataFrame:
        """Ultra-fast vectorized lineage and taxid matching with aggregation."""

        # Prepare lineage data
        ncbi_df = ncbi_df.copy()
        isolate_df = self.isolate_data.get(level, pd.DataFrame())

        results = []
        for _, census_row in census_df.iterrows():
            census_name = census_row[census_col]
            census_taxid = census_row.get('taxid', None)

            # Method 1: Direct taxid matching (most reliable)
            taxid_matched_ncbi = pd.DataFrame()
            if census_taxid and 'taxid' in ncbi_df.columns:
                taxid_matches_mask = ncbi_df['taxid'] == census_taxid
                taxid_matched_ncbi = ncbi_df[taxid_matches_mask]

            # Method 2: Vectorized lineage matching
            escaped_name = re.escape(census_name)
            pattern = f';{escaped_name};|^{escaped_name};|;{escaped_name}$|^{escaped_name}$'
            lineage_matches_mask = ncbi_df['lineage'].str.contains(pattern, regex=True, na=False)
            lineage_matched_ncbi = ncbi_df[lineage_matches_mask]

            # Combine matches (taxid matches take priority, avoid duplicates)
            if not taxid_matched_ncbi.empty and not lineage_matched_ncbi.empty:
                # Remove lineage matches that are already covered by taxid matches
                lineage_only = lineage_matched_ncbi[~lineage_matched_ncbi.index.isin(taxid_matched_ncbi.index)]
                matched_ncbi = pd.concat([taxid_matched_ncbi, lineage_only], ignore_index=True)
            elif not taxid_matched_ncbi.empty:
                matched_ncbi = taxid_matched_ncbi
            elif not lineage_matched_ncbi.empty:
                matched_ncbi = lineage_matched_ncbi
            else:
                matched_ncbi = pd.DataFrame()

            if not matched_ncbi.empty:
                # Vectorized aggregation
                total_species = matched_ncbi[f'{level}_species_count'].sum()
                total_genomes = matched_ncbi[f'{level}_genome_count'].sum()

                # Calculate percentages using correct total denominators from NCBI database
                total_ncbi_genomes = getattr(self, f'total_ncbi_genomes_{level}', 1)
                total_ncbi_species = getattr(self, f'total_ncbi_species_{level}', 1)

                genome_pct_db = (total_genomes / total_ncbi_genomes * 100) if total_ncbi_genomes > 0 else 0
                species_pct = (total_species / total_ncbi_species * 100) if total_ncbi_species > 0 else 0

                # Count match types
                direct_matches = int(census_name in matched_ncbi[ncbi_col].values)
                taxid_matches = len(taxid_matched_ncbi) if not taxid_matched_ncbi.empty else 0
                lineage_matches = len(matched_ncbi) - direct_matches - taxid_matches

                # Log the matches
                for _, match_row in matched_ncbi.iterrows():
                    if match_row[ncbi_col] == census_name:
                        match_type = 'direct_match'
                    elif not taxid_matched_ncbi.empty and match_row.name in taxid_matched_ncbi.index:
                        match_type = 'taxid_match'
                    else:
                        match_type = 'lineage_match'

                    self.matches_log.append({
                        'level': level,
                        'census_taxon': census_name,
                        'census_taxid': census_taxid,
                        'ncbi_taxon': match_row[ncbi_col],
                        'ncbi_taxid': match_row.get('taxid', 'N/A'),
                        'match_type': match_type,
                        'ncbi_species_count': match_row[f'{level}_species_count'],
                        'ncbi_genome_count': match_row[f'{level}_genome_count'],
                        'lineage': match_row['lineage']
                    })
            else:
                total_species = total_genomes = direct_matches = taxid_matches = lineage_matches = 0
                genome_pct_db = species_pct = 0

            # Vectorized isolate aggregation - sum isolate counts from all matched NCBI taxa
            isolate_count = isolate_pct = 0
            if not isolate_df.empty and not matched_ncbi.empty:
                matched_isolate_taxa = isolate_df[isolate_df[level].isin(matched_ncbi[ncbi_col])]
                if not matched_isolate_taxa.empty:
                    isolate_count = matched_isolate_taxa['isolate_count'].sum()
                    total_matched_genomes = matched_isolate_taxa['total_genomes'].sum()
                    isolate_pct = (isolate_count / total_matched_genomes * 100).round(2) if total_matched_genomes > 0 else 0



            # Determine domain from matched NCBI taxa (for prokaryotes: Bacteria or Archaea)
            domain = 'Unknown'
            if not matched_ncbi.empty and 'domain' in matched_ncbi.columns:
                # Get the most common domain from matched taxa
                domain_counts = matched_ncbi['domain'].value_counts()
                if len(domain_counts) > 0:
                    domain = domain_counts.index[0]  # Most frequent domain

            # Calculate novelty and overrepresentation factors
            census_otus = census_row['census_otu_count']
            novelty_factor = census_otus / total_species if total_species > 0 else float('inf')
            overrepresentation_factor = total_species / census_otus if census_otus > 0 else float('inf')

            results.append({
                level: census_name,
                'census_otu_count': census_row['census_otu_count'],
                'census_size_count': census_row['census_size_count'],
                'otu_percentage': census_row['otu_percentage'],
                'size_percentage': census_row['size_percentage'],
                'ncbi_genome_count': total_genomes,
                'ncbi_species_count': total_species,
                'isolate_count': isolate_count,
                'genome_pct_db': round(genome_pct_db, 2),
                'species_pct': round(species_pct, 2),
                'isolate_percentage': isolate_pct,

                'novelty_factor': round(novelty_factor, 3),
                'overrepresentation_factor': round(overrepresentation_factor, 3),
                'direct_matches': direct_matches,
                'taxid_matches': taxid_matches,
                'lineage_matches': lineage_matches,
                'total_matches': direct_matches + taxid_matches + lineage_matches,
                'match_status': 'matched' if total_species > 0 else 'census_only',
                'domain': domain
            })

        return pd.DataFrame(results)

    def add_ncbi_only_taxa(self, merged_df: pd.DataFrame, ncbi_df: pd.DataFrame,
                          census_col: str, ncbi_col: str, level: str) -> pd.DataFrame:
        """Add NCBI-only taxa that weren't matched to any census taxon."""

        # Find all NCBI taxa that were matched (both lineage and taxid matches)
        all_matched = set()
        for _, row in merged_df.iterrows():
            if row['total_matches'] > 0:
                census_name = row[level]

                # Method 1: Find lineage-based matches
                escaped_name = re.escape(census_name)
                pattern = f';{escaped_name};|^{escaped_name};|;{escaped_name}$|^{escaped_name}$'
                lineage_matches_mask = ncbi_df['lineage'].str.contains(pattern, regex=True, na=False)
                lineage_matched_taxa = ncbi_df.loc[lineage_matches_mask, ncbi_col].tolist()
                all_matched.update(lineage_matched_taxa)

                # Method 2: Find taxid-based matches (if census data has taxid info)
                # Get the census taxid from the original census data
                census_df = self.census_data[level]
                census_row_match = census_df[census_df.iloc[:, 0] == census_name]
                if not census_row_match.empty and 'taxid' in census_row_match.columns:
                    census_taxid = census_row_match['taxid'].iloc[0]
                    if census_taxid and 'taxid' in ncbi_df.columns:
                        taxid_matches_mask = ncbi_df['taxid'] == census_taxid
                        taxid_matched_taxa = ncbi_df.loc[taxid_matches_mask, ncbi_col].tolist()
                        all_matched.update(taxid_matched_taxa)

        # Find unmatched NCBI taxa
        unmatched_ncbi = ncbi_df[~ncbi_df[ncbi_col].isin(all_matched)]

        ncbi_only_rows = []
        isolate_df = self.isolate_data.get(level, pd.DataFrame())

        for _, ncbi_row in unmatched_ncbi.iterrows():
            ncbi_name = ncbi_row[ncbi_col]
            genome_count = ncbi_row[f'{level}_genome_count']
            species_count = ncbi_row[f'{level}_species_count']

            # Calculate percentages using correct total denominators from NCBI database
            total_ncbi_genomes = getattr(self, f'total_ncbi_genomes_{level}', 1)
            total_ncbi_species = getattr(self, f'total_ncbi_species_{level}', 1)

            genome_pct_db = (genome_count / total_ncbi_genomes * 100) if total_ncbi_genomes > 0 else 0
            species_pct = (species_count / total_ncbi_species * 100) if total_ncbi_species > 0 else 0

            # Get isolate data
            isolate_count = isolate_pct = 0
            if not isolate_df.empty:
                isolate_match = isolate_df[isolate_df[level] == ncbi_name]
                if not isolate_match.empty:
                    isolate_count = isolate_match['isolate_count'].iloc[0]
                    isolate_pct = isolate_match['isolate_percentage'].iloc[0]

            if genome_count > 0 or species_count > 0:  # Only add taxa with actual data
                # Get domain from NCBI data
                domain = ncbi_row.get('domain', 'Unknown')

                ncbi_only_rows.append({
                    level: ncbi_name,  # Use level name instead of census_col
                    'census_otu_count': 0,
                    'census_size_count': 0,
                    'otu_percentage': 0,
                    'size_percentage': 0,
                    'ncbi_genome_count': genome_count,
                    'ncbi_species_count': species_count,
                    'isolate_count': isolate_count,
                    'genome_pct_db': round(genome_pct_db, 2),
                    'species_pct': round(species_pct, 2),
                    'isolate_percentage': isolate_pct,
                    'coverage_percentage': 0,
                    'direct_matches': 0,
                    'taxid_matches': 0,
                    'lineage_matches': 0,
                    'total_matches': 0,
                    'match_status': 'ncbi_only',
                    'domain': domain
                })

        if ncbi_only_rows:
            ncbi_only_df = pd.DataFrame(ncbi_only_rows)
            return pd.concat([merged_df, ncbi_only_df], ignore_index=True)

        return merged_df

    def process_level(self, level: str) -> pd.DataFrame:
        """Process one taxonomic level with vectorized operations."""
        ncbi_df = self.ncbi_data[level]
        census_df = self.census_data[level]

        # Column names
        ncbi_col = level
        census_col = census_df.columns[0]

        # Vectorized matching and aggregation
        merged_df = self.vectorized_lineage_matching(census_df, ncbi_df, census_col, ncbi_col, level)

        # Add NCBI-only taxa
        final_df = self.add_ncbi_only_taxa(merged_df, ncbi_df, census_col, ncbi_col, level)

        # Sort results: matched taxa first (by ncbi_genome_count desc), then unmatched taxa (by ncbi_genome_count desc)
        matched_df = final_df[final_df['match_status'] == 'matched'].sort_values('ncbi_genome_count', ascending=False)
        unmatched_df = final_df[final_df['match_status'] != 'matched'].sort_values('ncbi_genome_count', ascending=False)
        final_df = pd.concat([matched_df, unmatched_df], ignore_index=True)

        return final_df

    def save_matches_log(self):
        """Save detailed matches log to CSV file."""
        if self.matches_log:
            matches_df = pd.DataFrame(self.matches_log)
            output_dir = self.output_dir / "16s_merged/logs"
            output_dir.mkdir(parents=True, exist_ok=True)
            matches_df.to_csv(output_dir / "16s_ncbi_detailed_matches.csv", index=False)
            self.logger.info(f"Saved detailed matches log: {len(matches_df)} entries")

    def run_merger(self):
        """Run the complete merger process."""
        self.logger.info("Starting 16S NCBI vectorized merger...")

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
            if level in self.ncbi_data and level in self.census_data:
                self.logger.info(f"Processing {level} level...")

                merged_df = self.process_level(level)

                # Save results
                output_file = results_dir / f"16s_ncbi_merged_clean_{level}.csv"
                merged_df.to_csv(output_file, index=False)
                self.logger.info(f"Saved {level} results: {len(merged_df)} entries")

                # Collect summary stats
                matched_count = len(merged_df[merged_df['match_status'] == 'matched'])
                census_only_count = len(merged_df[merged_df['match_status'] == 'census_only'])
                ncbi_only_count = len(merged_df[merged_df['match_status'] == 'ncbi_only'])

                summary_data.append({
                    'level': level,
                    'total_entries': len(merged_df),
                    'matched_taxa': matched_count,
                    'census_only_taxa': census_only_count,
                    'ncbi_only_taxa': ncbi_only_count,
                    'match_rate': f"{matched_count/len(merged_df)*100:.1f}%" if len(merged_df) > 0 else "0%"
                })

        # Save summary
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(summary_dir / "16s_ncbi_merger_clean_summary.csv", index=False)
            self.logger.info("Saved merger summary")

        # Save matches log
        self.save_matches_log()

        self.logger.info("16S NCBI merger completed successfully!")

def main():
    """Main execution function."""
    logger = setup_logging()

    try:
        merger = SixteenSNCBIMerger()
        merger.run_merger()

        print("\n✅ 16S NCBI Merger completed successfully!")
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
