#!/usr/bin/env python3
"""
18S Census + NCBI Merger - Vectorized Version (Updated: 2025-01-13)
===================================================================

Fast, vectorized merger using pandas operations for lineage-based matching.
Both genome_pct_db and species_pct calculated as percentage of Eukaryota subset only.

INPUT FILES:
-----------
NCBI Data:
  - ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_family_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_genus_counts.csv
  - ncbi_parse/csv_ncbi/ncbi_*_with_accessions.csv (for isolate analysis)

18S Census Data:
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_family.csv
  - 18S_censusparse/csv_outputs/eukcensus_18S_by_genus.csv

OUTPUT FILES:
------------
  - merged_output/18s_merged/results/18s_ncbi_merged_clean_*.csv
  - merged_output/18s_merged/analysis_summary/18s_ncbi_merger_clean_summary.csv

LOGIC:
------
1. Filter NCBI data to Eukaryota domain only
2. Vectorized Lineage Matching: Use pandas string operations to find census taxa in NCBI lineages
3. Vectorized Aggregation: Use pandas groupby and sum operations for count aggregation
4. Isolate Analysis: Vectorized isolate percentage calculations
5. Clean Output: Same format as clean version

Date: 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys

def setup_logging():
    """Setup minimal logging."""
    # Use relative path to 18s_merged logs directory
    logs_dir = Path("../18s_merged/logs")
    logs_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.WARNING,  # Only warnings and errors
        format='%(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(logs_dir / '18s_ncbi_merger_vectorized.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )

logger = logging.getLogger(__name__)

class Vectorized18SNCBIMerger:
    """Ultra-fast vectorized 18S Census + NCBI merger."""
    
    def __init__(self, ncbi_data_dir: str, census_data_dir: str, output_dir: str):
        self.ncbi_data_dir = Path(ncbi_data_dir)
        self.census_data_dir = Path(census_data_dir)
        self.output_dir = Path(output_dir) / "18s_merged"
        
        # Create output directories
        self.results_dir = self.output_dir / "csv_results"
        self.analysis_dir = self.output_dir / "analysis_summary"
        self.logs_dir = self.output_dir / "logs"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        # Taxonomic levels
        self.levels = ['phylum', 'family', 'genus']
        self.census_levels = ['division', 'family', 'genus']

        # Initialize matches log
        self.matches_log = []
        
    def load_data(self):
        """Load all data files efficiently."""
        print("📂 Loading data files...")

        # Load NCBI data with filtering
        self.ncbi_data = {}
        for level in self.levels:
            file_path = self.ncbi_data_dir / f"ncbi_{level}_counts.csv"
            if file_path.exists():
                df = pd.read_csv(file_path)
                # Filter to keep only Eukaryota domain
                eukaryota_mask = df['domain'].str.contains('Eukaryota', case=False, na=False)
                df_filtered = df[eukaryota_mask]

                self.ncbi_data[level] = df_filtered
                print(f"   🧬 Filtered NCBI {level}: {len(df_filtered)}/{len(df)} taxa (kept only Eukaryota domain)")

                # Calculate total denominators for percentage calculations (from Eukaryota only)
                # Both genome_pct_db and species_pct use only Eukaryota subset
                total_genomes_euk = df_filtered[f'{level}_genome_count'].sum()
                total_species_euk = df_filtered[f'{level}_species_count'].sum()
                setattr(self, f'total_ncbi_genomes_{level}', total_genomes_euk)
                setattr(self, f'total_ncbi_species_{level}', total_species_euk)
                print(f"   📊 NCBI {level} Eukaryota totals: {total_genomes_euk:,} genomes, {total_species_euk:,} species")
        
        # Load census data with filtering
        self.census_data = {}
        for i, level in enumerate(self.levels):
            census_level = self.census_levels[i]
            file_path = self.census_data_dir / f"eukcensus_18S_by_{census_level}.csv"
            if file_path.exists():
                df = pd.read_csv(file_path)
                # Keep all census data including .U. entries (no filtering)
                census_col = df.columns[0]  # First column is taxon name

                self.census_data[level] = df
                print(f"   📋 Loaded census {census_level}: {len(df)} taxa")

                # Calculate total denominators for census percentages
                total_census_otus = df['otu_count'].sum()
                total_census_size = df['size_count'].sum()
                setattr(self, f'total_census_otus_{level}', total_census_otus)
                setattr(self, f'total_census_size_{level}', total_census_size)
                print(f"   📊 Census {census_level} database totals: {total_census_otus:,} OTUs, {total_census_size:,} size count")
        
        # Load isolate data
        self.isolate_data = {}
        for level in self.levels:
            file_path = self.ncbi_data_dir / f"ncbi_{level}_with_accessions.csv"
            if file_path.exists():
                df = pd.read_csv(file_path)
                if 'genome_source' in df.columns:
                    # Vectorized isolate calculations
                    isolate_counts = df[df['genome_source'] == 'isolate'].groupby(level).size()
                    total_counts = df.groupby(level).size()
                    
                    isolate_df = pd.DataFrame({
                        'taxon': isolate_counts.index,
                        'isolate_count': isolate_counts.values,
                        'total_genomes': total_counts.loc[isolate_counts.index].values
                    })
                    isolate_df['isolate_percentage'] = (isolate_df['isolate_count'] / isolate_df['total_genomes'] * 100).round(2)
                    self.isolate_data[level] = isolate_df
    
    def vectorized_lineage_matching(self, census_df: pd.DataFrame, ncbi_df: pd.DataFrame,
                                   census_col: str, ncbi_col: str, level: str) -> pd.DataFrame:
        """Ultra-fast vectorized lineage matching and aggregation."""

        # Prepare lineage data
        ncbi_df = ncbi_df.copy()
        ncbi_df['lineage'] = ncbi_df['lineage'].fillna('')

        results = []
        isolate_df = self.isolate_data.get(level, pd.DataFrame())

        # Get totals for percentage calculations
        total_census_otus = census_df['otu_count'].sum()
        total_census_size = census_df['size_count'].sum()

        for _, census_row in census_df.iterrows():
            census_name = census_row[census_col]
            census_otus = census_row.get('otu_count', 0)
            census_size = census_row.get('size_count', 0)

            # Calculate percentages using the totals from the current dataset
            otu_percentage = (census_otus / total_census_otus * 100) if total_census_otus > 0 else 0
            size_percentage = (census_size / total_census_size * 100) if total_census_size > 0 else 0

            # Vectorized lineage search
            pattern = f';{census_name};|^{census_name};|;{census_name}$|^{census_name}$'
            matches_mask = ncbi_df['lineage'].str.contains(pattern, regex=True, na=False)
            matched_ncbi = ncbi_df[matches_mask]
            
            if len(matched_ncbi) > 0:
                # Vectorized aggregation
                total_species = matched_ncbi[f'{level}_species_count'].sum()
                total_genomes = matched_ncbi[f'{level}_genome_count'].sum()

                # Calculate percentages using Eukaryota-only totals from NCBI database
                total_ncbi_genomes = getattr(self, f'total_ncbi_genomes_{level}', 1)
                total_ncbi_species = getattr(self, f'total_ncbi_species_{level}', 1)

                genome_pct_db = (total_genomes / total_ncbi_genomes * 100) if total_ncbi_genomes > 0 else 0
                species_pct = (total_species / total_ncbi_species * 100) if total_ncbi_species > 0 else 0

                # Count match types
                direct_matches = int(census_name in matched_ncbi[ncbi_col].values)
                lineage_matches = len(matched_ncbi) - direct_matches

                # Log the matches
                for _, match_row in matched_ncbi.iterrows():
                    match_type = 'direct_match' if match_row[ncbi_col] == census_name else 'lineage_match'
                    self.matches_log.append({
                        'level': level,
                        'census_taxon': census_name,
                        'ncbi_taxon': match_row[ncbi_col],
                        'match_type': match_type,
                        'ncbi_species_count': match_row[f'{level}_species_count'],
                        'ncbi_genome_count': match_row[f'{level}_genome_count'],
                        'lineage': match_row['lineage']
                    })
            else:
                total_species = total_genomes = direct_matches = lineage_matches = 0
                genome_pct_db = species_pct = 0
            
            # Vectorized isolate aggregation - sum isolate counts from all matched NCBI taxa
            isolate_count = isolate_pct = 0
            if not isolate_df.empty and not matched_ncbi.empty:
                # Get all matched NCBI taxa names
                matched_taxa_names = matched_ncbi[ncbi_col].unique()
                # Sum isolate counts from all matched taxa
                isolate_matches = isolate_df[isolate_df['taxon'].isin(matched_taxa_names)]
                if not isolate_matches.empty:
                    isolate_count = isolate_matches['isolate_count'].sum()
                    # Use the actual NCBI genome count as denominator for consistency
                    isolate_pct = (isolate_count / total_genomes * 100).round(2) if total_genomes > 0 else 0
            


            # Determine domain from matched NCBI taxa (for eukaryotes: mostly Eukaryota)
            domain = 'Unknown'
            if len(matched_ncbi) > 0 and 'domain' in matched_ncbi.columns:
                # Get the most common domain from matched taxa
                domain_counts = matched_ncbi['domain'].value_counts()
                if len(domain_counts) > 0:
                    domain = domain_counts.index[0]  # Most frequent domain

            # Calculate novelty and overrepresentation factors
            novelty_factor = census_otus / total_species if total_species > 0 else float('inf')
            overrepresentation_factor = total_species / census_otus if census_otus > 0 else float('inf')



            results.append({
                level: census_name,  # Use level name instead of census_col
                'census_otu_count': census_otus,
                'census_size_count': census_size,
                'otu_percentage': round(otu_percentage, 2),
                'size_percentage': round(size_percentage, 2),
                'ncbi_genome_count': total_genomes,
                'ncbi_species_count': total_species,
                'isolate_count': isolate_count,
                'genome_pct_db': round(genome_pct_db, 2),
                'species_pct': round(species_pct, 2),
                'isolate_percentage': isolate_pct,
                'novelty_factor': round(novelty_factor, 3),
                'overrepresentation_factor': round(overrepresentation_factor, 3),
                'direct_matches': direct_matches,
                'lineage_matches': lineage_matches,
                'total_matches': direct_matches + lineage_matches,
                'match_status': 'matched' if total_species > 0 else 'census_only',
                'domain': domain,
                'coverage_percentage': 0.0  # Add this field to match ncbi_only_rows structure
            })
        
        return pd.DataFrame(results)
    
    def add_ncbi_only_taxa(self, merged_df: pd.DataFrame, ncbi_df: pd.DataFrame, 
                          census_col: str, ncbi_col: str, level: str) -> pd.DataFrame:
        """Add NCBI-only taxa that weren't matched to any census taxon."""
        
        # Find all NCBI taxa that were matched
        all_matched = set()
        for _, row in merged_df.iterrows():
            if row['total_matches'] > 0:
                census_name = row[level]  # Use level name since we changed the column structure
                pattern = f';{census_name};|^{census_name};|;{census_name}$|^{census_name}$'
                matches_mask = ncbi_df['lineage'].str.contains(pattern, regex=True, na=False)
                matched_taxa = ncbi_df.loc[matches_mask, ncbi_col].tolist()
                all_matched.update(matched_taxa)
        
        # Add unmatched NCBI taxa
        isolate_df = self.isolate_data.get(level, pd.DataFrame())
        ncbi_only_rows = []
        
        for _, ncbi_row in ncbi_df.iterrows():
            ncbi_name = ncbi_row[ncbi_col]
            if ncbi_name not in all_matched:
                species_count = ncbi_row[f'{level}_species_count']
                genome_count = ncbi_row[f'{level}_genome_count']
                
                # Get isolate info for this specific NCBI taxon
                isolate_count = isolate_pct = 0
                if not isolate_df.empty:
                    isolate_match = isolate_df[isolate_df['taxon'] == ncbi_name]
                    if not isolate_match.empty:
                        isolate_count = isolate_match.iloc[0]['isolate_count']
                        isolate_pct = isolate_match.iloc[0]['isolate_percentage']
                
                # Calculate percentages using Eukaryota-only totals from NCBI database
                total_ncbi_genomes = getattr(self, f'total_ncbi_genomes_{level}', 1)
                total_ncbi_species = getattr(self, f'total_ncbi_species_{level}', 1)

                genome_pct_db = (genome_count / total_ncbi_genomes * 100) if total_ncbi_genomes > 0 else 0
                species_pct = (species_count / total_ncbi_species * 100) if total_ncbi_species > 0 else 0

                # Get domain from NCBI data
                domain = ncbi_row.get('domain', 'Unknown')

                ncbi_only_rows.append({
                    level: ncbi_name,  # Use level name instead of census_col
                    'census_otu_count': 0,
                    'census_size_count': 0,
                    'otu_percentage': 0.0,
                    'size_percentage': 0.0,
                    'ncbi_genome_count': genome_count,
                    'ncbi_species_count': species_count,
                    'isolate_count': isolate_count,
                    'genome_pct_db': round(genome_pct_db, 2),
                    'species_pct': round(species_pct, 2),
                    'isolate_percentage': isolate_pct,
                    'coverage_percentage': 0,
                    'direct_matches': 0,
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
            log_file = self.logs_dir / "18s_ncbi_detailed_matches.csv"
            matches_df.to_csv(log_file, index=False)
            print(f"📋 Matches log saved: {log_file}")

    def run_merger(self):
        """Run the complete vectorized merger."""
        print("🚀 Starting vectorized 18S + NCBI merger...")
        
        # Load all data
        self.load_data()
        
        # Process each level
        summary_stats = []
        
        for level in self.levels:
            if level in self.ncbi_data and level in self.census_data:
                print(f"⚡ Processing {level} level...")

                # Vectorized processing
                result_df = self.process_level(level)
                
                # Save results with same filename as clean version
                output_file = self.results_dir / f"18s_ncbi_merged_clean_{level}.csv"
                result_df.to_csv(output_file, index=False)
                
                # Calculate stats (matching 16S format)
                total_entries = len(result_df)
                matched_taxa = len(result_df[result_df['match_status'] == 'matched'])
                census_only_taxa = len(result_df[result_df['match_status'] == 'census_only'])
                ncbi_only_taxa = len(result_df[result_df['match_status'] == 'ncbi_only'])
                match_rate = f"{matched_taxa/total_entries*100:.1f}%" if total_entries > 0 else "0%"

                summary_stats.append({
                    'level': level,
                    'total_entries': total_entries,
                    'matched_taxa': matched_taxa,
                    'census_only_taxa': census_only_taxa,
                    'ncbi_only_taxa': ncbi_only_taxa,
                    'match_rate': match_rate
                })

                print(f"   ✅ {matched_taxa}/{total_entries} taxa matched ({match_rate} match rate)")
                print(f"      Census-only: {census_only_taxa}, NCBI-only: {ncbi_only_taxa}")
        
        # Save summary with same filename
        summary_df = pd.DataFrame(summary_stats)
        summary_file = self.analysis_dir / "18s_ncbi_merger_clean_summary.csv"
        summary_df.to_csv(summary_file, index=False)

        # Save detailed matches log
        self.save_matches_log()

        print("🎉 Vectorized merger completed!")
        print(f"📁 Results: {self.results_dir}")
        print(f"📋 Logs: {self.logs_dir}")
        
        return summary_stats

if __name__ == "__main__":
    # Setup logging
    setup_logging()

    # Use paths relative to script location
    script_dir = Path(__file__).resolve().parent
    # Go up from py_mergers -> Eukcensus_merge -> parse_repaa_table
    base_dir = script_dir.parent.parent
    merger = Vectorized18SNCBIMerger(
        ncbi_data_dir=str(base_dir / "ncbi_parse" / "csv_ncbi"),
        census_data_dir=str(base_dir / "18S_censusparse" / "csv_outputs"),
        output_dir=str(script_dir.parent)  # Output to Eukcensus_merge directory
    )
    
    results = merger.run_merger()

    print("\n✅ 18S NCBI Merger completed successfully!")
    print("📁 Output files saved to: ./18s_merged/")
    print("   - csv_results/: Merged CSV files")
    print("   - analysis_summary/: Summary statistics")
    print("   - logs/: Detailed matching logs")
