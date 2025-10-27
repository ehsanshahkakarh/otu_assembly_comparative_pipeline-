#!/usr/bin/env python3
"""
Family Genera Extractor
=======================

Takes family-level source data (novel and overrepresented) and extracts all genera
that belong to those families, showing what makes each family novel or overrepresented.

For Novel families: Focus on census data to show environmental diversity
For Overrepresented families: Focus on NCBI data to show database diversity

Usage:
    python extract_family_genera.py

Output:
    - 16s_family_novel_genera_breakdown.csv
    - 16s_family_overrepresented_genera_breakdown.csv  
    - 18s_family_novel_genera_breakdown.csv
    - 18s_family_overrepresented_genera_breakdown.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class FamilyGeneraExtractor:
    """Extract genera for novel and overrepresented families."""
    
    def __init__(self):
        self.base_dir = Path(".")
        self.source_data_dir = self.base_dir / "source_data"
        self.output_dir = self.source_data_dir / "family_genera_breakdown"
        self.output_dir.mkdir(exist_ok=True)
        
        # Original metadata directories
        self.census_16s_metadata = self.base_dir / ".." / ".." / "16S_censusparse" / "metadata" / "eukcensus_16S.clusters.97.tsv"
        self.census_18s_metadata = self.base_dir / ".." / ".." / "18S_censusparse" / "metadata" / "eukcensus_18S.clusters.97.tsv"
        self.ncbi_dir = self.base_dir / ".." / ".." / "ncbi_parse" / "csv_ncbi"
        self.merged_16s_dir = self.base_dir / ".." / "Eukcensus_merge" / "16s_merged" / "csv_results"
        self.merged_18s_dir = self.base_dir / ".." / "Eukcensus_merge" / "18s_merged" / "csv_results"
        
    def load_family_source_data(self, gene: str, factor_type: str) -> pd.DataFrame:
        """Load family source data for a specific gene and factor type."""
        filename = f"{gene}_family_{factor_type}_source_data.csv"
        filepath = self.source_data_dir / filename
        
        if not filepath.exists():
            logger.warning(f"Source data file not found: {filepath}")
            return pd.DataFrame()
            
        df = pd.read_csv(filepath)
        logger.info(f"Loaded {len(df)} {factor_type} families for {gene}")
        return df
    
    def load_genus_data(self, gene: str, source: str) -> pd.DataFrame:
        """Load genus-level data from census or NCBI sources."""
        if gene == "16s":
            if source == "census":
                # Load 16S census genus data
                filepath = self.census_16s_dir / "eukcensus16S_by_genus.csv"
            elif source == "ncbi":
                # Load 16S NCBI merged genus data
                filepath = self.merged_16s_dir / "16s_ncbi_merged_clean_genus.csv"
            elif source == "gtdb":
                # Load 16S GTDB merged genus data  
                filepath = self.merged_16s_dir / "16s_gtdb_merged_clean_genus.csv"
        elif gene == "18s":
            if source == "census":
                # Load 18S census genus data
                filepath = self.census_18s_dir / "eukcensus_18S_by_genus.csv"
            elif source == "ncbi":
                # Load 18S NCBI merged genus data
                filepath = self.merged_18s_dir / "18s_ncbi_merged_clean_genus.csv"
            elif source == "eukprot":
                # Load 18S EukProt merged genus data
                filepath = self.merged_18s_dir / "18s_eukprot_merged_genus.csv"
        
        if not filepath.exists():
            logger.warning(f"Genus data file not found: {filepath}")
            return pd.DataFrame()
            
        df = pd.read_csv(filepath)
        logger.info(f"Loaded {len(df)} genera from {gene} {source} data")
        return df

    def load_original_metadata(self, gene):
        """Load original metadata files with complete taxonomic information."""
        if gene == "16s":
            filepath = self.census_16s_metadata
            if not filepath.exists():
                logger.error(f"16S metadata file not found: {filepath}")
                return pd.DataFrame()

            logger.info(f"Loading 16S original metadata from: {filepath}")
            df = pd.read_csv(filepath, sep='\t', low_memory=False)

            # Rename columns for consistency (16S has 'familiy' typo)
            if 'familiy' in df.columns:
                df = df.rename(columns={'familiy': 'family'})

            logger.info(f"Loaded {len(df):,} records from 16S metadata")
            return df

        elif gene == "18s":
            filepath = self.census_18s_metadata
            if not filepath.exists():
                logger.error(f"18S metadata file not found: {filepath}")
                return pd.DataFrame()

            logger.info(f"Loading 18S original metadata from: {filepath}")
            df = pd.read_csv(filepath, sep='\t', low_memory=False)

            # 18S uses 'division' instead of 'phylum'
            if 'division' in df.columns:
                df = df.rename(columns={'division': 'phylum'})

            logger.info(f"Loaded {len(df):,} records from 18S metadata")
            return df

        else:
            logger.error(f"Unknown gene: {gene}")
            return pd.DataFrame()

    def extract_genera_from_metadata(self, metadata, family_name):
        """Extract genera belonging to a specific family from original metadata."""
        # Direct family match
        direct_match = metadata[metadata['family'].str.contains(family_name, case=False, na=False)]

        if not direct_match.empty:
            # Group by genus and sum sizes, include domain information
            genus_summary = direct_match.groupby('genus').agg({
                'size': 'sum',
                'phylum': 'first',
                'family': 'first'
            }).reset_index()

            # Add domain information based on phylum
            genus_summary['domain'] = genus_summary['phylum'].apply(self.get_domain_from_phylum)

            logger.info(f"Found {len(genus_summary)} genera for '{family_name}' via direct family match")
            return genus_summary

        # If no direct match, try hierarchical search in other taxonomic levels
        # This handles cases where family might be mentioned in higher-level classifications
        hierarchical_match = metadata[
            metadata['phylum'].str.contains(family_name, case=False, na=False) |
            metadata['genus'].str.contains(family_name, case=False, na=False)
        ]

        if not hierarchical_match.empty:
            genus_summary = hierarchical_match.groupby('genus').agg({
                'size': 'sum',
                'phylum': 'first',
                'family': 'first'
            }).reset_index()

            # Add domain information based on phylum
            genus_summary['domain'] = genus_summary['phylum'].apply(self.get_domain_from_phylum)

            logger.info(f"Found {len(genus_summary)} genera for '{family_name}' via hierarchical search")
            return genus_summary

        logger.info(f"No genera found for '{family_name}' in metadata")
        return pd.DataFrame()

    def get_domain_from_phylum(self, phylum):
        """Classify domain based on phylum name."""
        if pd.isna(phylum):
            return 'Unknown'

        phylum_str = str(phylum).lower()

        # Archaea phyla
        archaea_phyla = [
            'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'korarchaeota',
            'nanoarchaeota', 'aigarchaeota', 'lokiarchaeota', 'thorarchaeota',
            'odinarchaeota', 'heimdallarchaeota', 'candidatus asgardarchaeota'
        ]

        for archaea_phylum in archaea_phyla:
            if archaea_phylum in phylum_str:
                return 'Archaea'

        # Everything else is likely Bacteria for 16S data
        return 'Bacteria'
    
    def extract_family_lineage_info(self, genus_data: pd.DataFrame, taxon_name: str) -> pd.DataFrame:
        """Extract genera that belong to a specific taxon using lineage information."""
        if 'lineage' not in genus_data.columns:
            logger.warning("No lineage column found in genus data")
            return pd.DataFrame()

        # Filter genera that have the taxon name in their lineage
        # Use case-insensitive matching and handle NaN values
        family_genera = genus_data[
            genus_data['lineage'].str.contains(taxon_name, case=False, na=False)
        ].copy()

        logger.info(f"Found {len(family_genera)} genera for taxon '{taxon_name}' using lineage matching")
        return family_genera

    def parse_lineage_hierarchy(self, genus_data: pd.DataFrame, target_taxon: str) -> pd.DataFrame:
        """Parse lineage hierarchy to find genera belonging to target taxon at any level."""
        if 'lineage' not in genus_data.columns:
            return pd.DataFrame()

        # Split lineages and check each level
        matching_genera = []

        for idx, row in genus_data.iterrows():
            lineage = row.get('lineage', '')
            if pd.isna(lineage) or lineage == '':
                continue

            # Split lineage by semicolon and check each taxonomic level
            lineage_parts = [part.strip() for part in str(lineage).split(';')]

            # Check if target taxon appears in any part of the lineage
            for part in lineage_parts:
                if target_taxon.lower() in part.lower():
                    matching_genera.append(row)
                    break

        result_df = pd.DataFrame(matching_genera)
        logger.info(f"Found {len(result_df)} genera for '{target_taxon}' using hierarchical lineage parsing")
        return result_df
    
    def extract_family_genera_from_merged(self, merged_data: pd.DataFrame, family_name: str) -> pd.DataFrame:
        """Extract genera from merged data that belong to a specific family."""
        # For merged data, we need to look at the taxonomic hierarchy
        # This is more complex as we need to map family -> genus relationships
        
        # For now, use a simple approach - look for family name in any available columns
        family_genera = pd.DataFrame()
        
        # Check if there's a family column or lineage information
        search_columns = []
        if 'family' in merged_data.columns:
            search_columns.append('family')
        if 'lineage' in merged_data.columns:
            search_columns.append('lineage')
        
        for col in search_columns:
            matches = merged_data[
                merged_data[col].str.contains(family_name, case=False, na=False)
            ].copy()
            if len(matches) > 0:
                family_genera = pd.concat([family_genera, matches], ignore_index=True)
        
        # Remove duplicates if any
        if len(family_genera) > 0:
            family_genera = family_genera.drop_duplicates()
        
        return family_genera
    
    def process_novel_families(self, gene: str, database: str = None) -> pd.DataFrame:
        """Process novel families - focus on census data and group by genus."""
        logger.info(f"Processing novel families for {gene}")

        # Load novel family source data
        novel_families = self.load_family_source_data(gene, "novelty")
        if novel_families.empty:
            return pd.DataFrame()

        # Load original metadata instead of parsed genus data
        metadata = self.load_original_metadata(gene)
        if metadata.empty:
            return pd.DataFrame()

        # Collect all genera from all novel families
        all_genera = {}  # genus_name -> {data, families}

        for _, family_row in novel_families.iterrows():
            family_name = family_row['Taxon']
            logger.info(f"Processing novel family: {family_name}")

            # Extract genera for this family from original metadata
            family_genera = self.extract_genera_from_metadata(metadata, family_name)

            if not family_genera.empty:
                for _, genus_row in family_genera.iterrows():
                    genus_name = genus_row.get('genus', 'Unknown')

                    if genus_name not in all_genera:
                        all_genera[genus_name] = {
                            'genus_data': genus_row,
                            'families': [],
                            'total_census_otu': genus_row.get('size', 0),  # 'size' is sequence count in metadata
                            'total_census_size': genus_row.get('size', 0),
                            'otu_percentage': 0  # Will calculate if needed
                        }

                    # Add family info with domain from genus data
                    genus_domain = genus_row.get('domain', 'Unknown')
                    genus_phylum = genus_row.get('phylum', 'Unknown')

                    all_genera[genus_name]['families'].append({
                        'family': family_name,
                        'domain': genus_domain,
                        'phylum': genus_phylum,
                        'novelty_factor': family_row.get('Novelty_Factor', 0)
                    })

        # Convert to DataFrame with grouped results
        results = []
        for genus_name, genus_info in all_genera.items():
            families_list = '; '.join([f['family'] for f in genus_info['families']])
            avg_novelty = sum([f['novelty_factor'] for f in genus_info['families']]) / len(genus_info['families'])

            # Get domain info from the first family (they should all be the same domain)
            domain = genus_info['families'][0].get('domain', 'Unknown')
            phylum = genus_info['families'][0].get('phylum', 'Unknown')

            result_row = {
                'Genus': genus_name,
                'Domain': domain,
                'Phylum': phylum,
                'Genus_Census_OTU_Count': genus_info['total_census_otu'],
                'Genus_Census_Size_Count': genus_info['total_census_size'],
                'Genus_OTU_Percentage': genus_info['otu_percentage'],
                'Associated_Novel_Families': families_list,
                'Number_of_Families': len(genus_info['families']),
                'Average_Family_Novelty_Factor': round(avg_novelty, 3),
                'Data_Source': 'Census'
            }
            results.append(result_row)

        # Sort by associated novel families first, then by OTU count descending within each family
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values(['Associated_Novel_Families', 'Genus_Census_OTU_Count'],
                                              ascending=[True, False])

        return results_df

    def process_overrepresented_families(self, gene: str, database: str = None) -> pd.DataFrame:
        """Process overrepresented families - focus on NCBI/database data and group by genus."""
        logger.info(f"Processing overrepresented families for {gene}")

        # Load overrepresented family source data
        overrep_families = self.load_family_source_data(gene, "overrepresented")
        if overrep_families.empty:
            return pd.DataFrame()

        # Load NCBI/database genus data (merged data)
        if gene == "16s":
            if database == "gtdb":
                db_genera = self.load_genus_data(gene, "gtdb")
            else:
                db_genera = self.load_genus_data(gene, "ncbi")
        else:  # 18s
            if database == "eukprot":
                db_genera = self.load_genus_data(gene, "eukprot")
            else:
                db_genera = self.load_genus_data(gene, "ncbi")

        if db_genera.empty:
            return pd.DataFrame()

        # Also try to load census data for lineage information
        census_genera = self.load_genus_data(gene, "census")

        # Collect all genera from all overrepresented families
        all_genera = {}  # genus_name -> {data, families}

        for _, family_row in overrep_families.iterrows():
            family_name = family_row['Taxon']
            logger.info(f"Processing overrepresented family: {family_name}")

            # Try to find genera using census data with lineage info first
            family_genera = pd.DataFrame()
            if not census_genera.empty and 'lineage' in census_genera.columns:
                family_genera = self.parse_lineage_hierarchy(census_genera, family_name)

                # Match these genera with NCBI data
                if not family_genera.empty:
                    census_genus_names = family_genera['Name_to_use'].tolist()
                    # Find matching genera in NCBI data
                    ncbi_matches = db_genera[db_genera['genus'].isin(census_genus_names)]

                    # Merge census and NCBI data
                    for _, census_genus in family_genera.iterrows():
                        genus_name = census_genus.get('Name_to_use', 'Unknown')
                        ncbi_match = ncbi_matches[ncbi_matches['genus'] == genus_name]

                        if genus_name not in all_genera:
                            # Get database data if available - handle different column names
                            if database == "eukprot":
                                db_species = ncbi_match['eukprot_species_count'].iloc[0] if len(ncbi_match) > 0 else 0
                                db_genomes = 0  # EukProt doesn't have genome count
                            elif database == "gtdb":
                                db_species = ncbi_match['gtdb_species_count'].iloc[0] if len(ncbi_match) > 0 else 0
                                db_genomes = ncbi_match['gtdb_genome_count'].iloc[0] if len(ncbi_match) > 0 else 0
                            else:  # ncbi
                                db_species = ncbi_match['ncbi_species_count'].iloc[0] if len(ncbi_match) > 0 else 0
                                db_genomes = ncbi_match['ncbi_genome_count'].iloc[0] if len(ncbi_match) > 0 else 0

                            census_otu = ncbi_match['census_otu_count'].iloc[0] if len(ncbi_match) > 0 else census_genus.get('otu_count', 0)

                            all_genera[genus_name] = {
                                'families': [],
                                'ncbi_species_count': db_species,
                                'ncbi_genome_count': db_genomes,
                                'census_otu_count': census_otu,
                                'otu_percentage': census_genus.get('otu_percentage', 0)
                            }

                        # Add family info
                        all_genera[genus_name]['families'].append({
                            'family': family_name,
                            'domain': family_row.get('Domain', 'Unknown'),
                            'phylum': family_row.get('Phylum', family_row.get('Division', 'Unknown')),
                            'overrep_factor': family_row.get('Overrepresentation_Factor', 0)
                        })

        # Convert to DataFrame with grouped results
        results = []
        for genus_name, genus_info in all_genera.items():
            families_list = '; '.join([f['family'] for f in genus_info['families']])
            avg_overrep = sum([f['overrep_factor'] for f in genus_info['families']]) / len(genus_info['families'])

            result_row = {
                'Genus': genus_name,
                'Genus_DB_Species_Count': genus_info['ncbi_species_count'],
                'Genus_DB_Genome_Count': genus_info['ncbi_genome_count'],
                'Genus_Census_OTU_Count': genus_info['census_otu_count'],
                'Genus_OTU_Percentage': genus_info['otu_percentage'],
                'Associated_Overrepresented_Families': families_list,
                'Number_of_Families': len(genus_info['families']),
                'Average_Family_Overrepresentation_Factor': round(avg_overrep, 3),
                'Data_Source': 'NCBI/Database'
            }
            results.append(result_row)

        # Sort by NCBI species count descending
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values('Genus_DB_Species_Count', ascending=False)

        return results_df

    def run_analysis(self):
        """Run simplified analysis - only novel families from census data."""
        logger.info("Starting Simplified Family Genera Extraction Analysis")
        logger.info("Focus: Novel families from census data only")

        # Process 16S NCBI Novel families
        logger.info("=" * 50)
        logger.info("Processing 16S NCBI Novel Families")
        logger.info("=" * 50)

        novel_16s_ncbi = self.process_novel_families("16s", "ncbi")
        if not novel_16s_ncbi.empty:
            # Save combined file
            output_file = self.output_dir / "16s_ncbi_novel_genera_grouped.csv"
            novel_16s_ncbi.to_csv(output_file, index=False)
            logger.info(f"✅ Saved 16S NCBI novel genera (combined): {output_file} ({len(novel_16s_ncbi)} genera)")

            # Separate by domain
            bacteria_df = novel_16s_ncbi[novel_16s_ncbi['Domain'] == 'Bacteria']
            archaea_df = novel_16s_ncbi[novel_16s_ncbi['Domain'] == 'Archaea']

            if not bacteria_df.empty:
                bacteria_file = self.output_dir / "16s_ncbi_novel_genera_bacteria.csv"
                bacteria_df.to_csv(bacteria_file, index=False)
                logger.info(f"✅ Saved 16S NCBI novel genera (Bacteria): {bacteria_file} ({len(bacteria_df)} genera)")

            if not archaea_df.empty:
                archaea_file = self.output_dir / "16s_ncbi_novel_genera_archaea.csv"
                archaea_df.to_csv(archaea_file, index=False)
                logger.info(f"✅ Saved 16S NCBI novel genera (Archaea): {archaea_file} ({len(archaea_df)} genera)")
        else:
            logger.warning("No 16S NCBI novel genera found")

        # Process 18S NCBI Novel families
        logger.info("=" * 50)
        logger.info("Processing 18S NCBI Novel Families")
        logger.info("=" * 50)

        novel_18s_ncbi = self.process_novel_families("18s", "ncbi")
        if not novel_18s_ncbi.empty:
            output_file = self.output_dir / "18s_ncbi_novel_genera_grouped.csv"
            novel_18s_ncbi.to_csv(output_file, index=False)
            logger.info(f"✅ Saved 18S NCBI novel genera: {output_file} ({len(novel_18s_ncbi)} genera)")
        else:
            logger.warning("No 18S NCBI novel genera found")

        logger.info("=" * 50)
        logger.info("✅ Simplified Analysis Complete!")
        logger.info(f"📁 Output files saved to: {self.output_dir}")
        logger.info("📊 Focus: Novel families with genera from census data")
        logger.info("=" * 50)

def main():
    """Main function to run the family genera extraction."""
    extractor = FamilyGeneraExtractor()
    extractor.run_analysis()

if __name__ == "__main__":
    main()
