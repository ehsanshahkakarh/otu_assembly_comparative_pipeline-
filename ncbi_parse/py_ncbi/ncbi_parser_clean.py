#!/usr/bin/env python3
"""
Clean NCBI Taxonomic Parser

A modern, unified parser that replaces all the redundant phylum/family/genus specific scripts.
Eliminates code duplication and provides a clean, efficient interface for NCBI data processing.

Features:
- Single script handles all taxonomic levels (phylum, family, genus)
- Generates both genome counts and species subset counts
- Built-in genome source classification (isolate vs uncultured)
- Memory efficient processing with progress tracking
- Clean output formats compatible with existing workflows

Input:
- NCBI assembly summary file (00assembly_summary_genbank.txt)
- Taxonomic mapping file (taxid_to_{level}.csv)

Output:
- Counts CSV: ncbi_{level}_counts.csv (genome counts + species counts)
- Accessions CSV: ncbi_{level}_with_accessions.csv (with classification)

Usage:
    python ncbi_parser_clean.py --level phylum
    python ncbi_parser_clean.py --level family --output-dir custom_output
    python ncbi_parser_clean.py --level genus --input-dir custom_input
"""

import pandas as pd
import argparse
import sys
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, Optional, List, Tuple
from tqdm import tqdm

class NCBIParser:
    """Clean, unified NCBI taxonomic parser."""
    
    def __init__(self, level: str, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None):
        self.level = level
        self.script_dir = Path(__file__).resolve().parent
        
        # Set up directories with smart defaults
        self.input_dir = input_dir or self._find_input_dir()
        self.output_dir = output_dir or (self.script_dir.parent / "csv_ncbi")
        self.output_dir.mkdir(exist_ok=True)
        
        # Set up file paths
        self.assembly_file = self._find_assembly_file()
        self.mapping_file = self._find_mapping_file()
        
        print(f"NCBI {level.upper()} PARSER")
        print(f"Input: {self.assembly_file.name}")
        print(f"Output: {self.output_dir.name}")
    
    def _find_input_dir(self) -> Path:
        """Find the best input directory."""
        metadata_dir = self.script_dir.parent / "metadata"
        return metadata_dir if metadata_dir.exists() else self.script_dir
    
    def _find_assembly_file(self) -> Path:
        """Find the assembly file with fallback logic."""
        candidates = [
            self.input_dir / "00assembly_summary_genbank.txt",
            self.script_dir.parent / "metadata" / "00assembly_summary_genbank.txt",
            self.script_dir / "00assembly_summary_genbank.txt"
        ]
        
        for candidate in candidates:
            if candidate.exists():
                return candidate
        
        raise FileNotFoundError("Assembly file not found in any expected location")
    
    def _find_mapping_file(self) -> Path:
        """Find the taxonomic mapping file."""
        mapping_file = self.script_dir.parent / "taxonomic_mapping" / f"taxid_to_{self.level}.csv"
        if not mapping_file.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
        return mapping_file
    
    def load_data(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Load and merge assembly and mapping data."""
        print("📖 Loading data...")

        # Read assembly file with proper header handling
        with open(self.assembly_file, 'r') as f:
            lines = f.readlines()

        # Find the header line (starts with #assembly_accession)
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#assembly_accession'):
                header_line = line.strip().lstrip('#')
                data_start = i + 1
                break

        if header_line is None:
            raise ValueError("Could not find header line in assembly file")

        # Read the data with the correct header
        assembly_df = pd.read_csv(self.assembly_file, sep='\t', low_memory=False,
                                 skiprows=data_start, names=header_line.split('\t'))

        mapping_df = pd.read_csv(self.mapping_file)
        print(f"Loaded {len(assembly_df):,} assemblies, {len(mapping_df):,} mappings")

        # Fix data type mismatch - ensure both taxid columns are integers
        assembly_df['taxid'] = pd.to_numeric(assembly_df['taxid'], errors='coerce')
        mapping_df['taxid'] = pd.to_numeric(mapping_df['taxid'], errors='coerce')

        # Remove any rows with invalid taxids
        assembly_df = assembly_df.dropna(subset=['taxid'])
        mapping_df = mapping_df.dropna(subset=['taxid'])

        merged_df = pd.merge(assembly_df, mapping_df, on='taxid', how='inner')
        print(f"Merged: {len(merged_df):,} entries")

        return merged_df, assembly_df
    
    def classify_genome_source(self, df: pd.DataFrame) -> pd.Series:
        """
        Classify genomes as isolate or uncultured.

        Logic:
        1. Default: All genomes are 'isolate'
        2. Check organism_name for uncultured patterns → 'uncultured'
        3. Check excluded_from_refseq for uncultured patterns → 'uncultured'

        Uncultured patterns: uncultured, environmental, metagenome, unclassified,
                           unknown, unidentified, mixed culture, enrichment culture
        """
        print("Classifying genome sources...")

        # Initialize all as isolate
        classification = pd.Series(['isolate'] * len(df), index=df.index)

        # Uncultured patterns (case-insensitive)
        uncultured_patterns = [
            'uncultured', 'environmental', 'metagenome', 'unclassified',
            'unknown', 'unidentified', 'mixed culture', 'enrichment culture',
            'derived from metagenome', 'metagenome-assembled', 'mag',
            'single amplified genome', 'sag', 'environmental sample'
        ]

        # Check organism names
        organism_names = df['organism_name'].fillna('').str.lower()
        for pattern in uncultured_patterns:
            mask = organism_names.str.contains(pattern, na=False)

            # Special case: enrichment culture with strain name should remain as isolate
            if pattern == 'enrichment culture':
                # Check if organism name contains strain indicators
                has_strain = (
                    organism_names.str.contains(r'\bstrain\b', na=False) |
                    organism_names.str.contains(r'\bisolate\b', na=False) |
                    organism_names.str.contains(r'\bstr\.\b', na=False)
                )
                # Only classify as uncultured if no strain name present
                mask = mask & (~has_strain)

            classification.loc[mask] = 'uncultured'

        # Check excluded_from_refseq column if available
        if 'excluded_from_refseq' in df.columns:
            excluded_notes = df['excluded_from_refseq'].fillna('').str.lower()
            for pattern in uncultured_patterns:
                mask = excluded_notes.str.contains(pattern, na=False)

                # Apply same enrichment culture exception for excluded_from_refseq
                if pattern == 'enrichment culture':
                    has_strain = (
                        organism_names.str.contains(r'\bstrain\b', na=False) |
                        organism_names.str.contains(r'\bisolate\b', na=False) |
                        organism_names.str.contains(r'\bstr\.\b', na=False)
                    )
                    mask = mask & (~has_strain)

                classification.loc[mask] = 'uncultured'

        isolate_count = (classification == 'isolate').sum()
        uncultured_count = (classification == 'uncultured').sum()
        print(f"Classification: {isolate_count:,} isolate ({isolate_count/len(df)*100:.1f}%), {uncultured_count:,} uncultured ({uncultured_count/len(df)*100:.1f}%)")

        return classification
    
    def extract_species_taxid(self, df: pd.DataFrame) -> pd.Series:
        """Extract species taxid using taxonkit if not available."""
        if 'species_taxid' in df.columns:
            return df['species_taxid']

        print("Extracting species_taxid...")
        unique_taxids = df['taxid'].unique()
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
            temp_filename = temp_file.name
            for taxid in unique_taxids:
                temp_file.write(f"{taxid}\n")
        
        try:
            # Run taxonkit to get species taxids
            result = subprocess.run([
                "taxonkit", "reformat", "--taxid-field", "1", 
                "--format", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}", 
                temp_filename
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                print(f"Taxonkit warning: {result.stderr}")
                return df['taxid']  # Fallback to taxid
            
            # Parse results
            taxid_to_species = {}
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 8:
                        taxid = int(parts[0])
                        species_taxid = parts[7] if parts[7] else parts[0]
                        taxid_to_species[taxid] = int(species_taxid)
            
            # Map back to dataframe
            species_taxids = df['taxid'].map(taxid_to_species).fillna(df['taxid'])
            print(f"Extracted {len(taxid_to_species):,} species taxids")
            return species_taxids
            
        except Exception as e:
            print(f"Error with taxonkit: {e}")
            return df['taxid']  # Fallback to taxid
        finally:
            Path(temp_filename).unlink(missing_ok=True)
    
    def create_species_subset(self, df: pd.DataFrame, species_col: str) -> pd.DataFrame:
        """Create species subset with isolate preference."""
        print("Creating species subset...")

        # Add priority columns for sorting
        df_copy = df.copy()

        # Genome source priority (isolate = 0, uncultured = 1)
        df_copy['source_priority'] = (df_copy['genome_source'] == 'uncultured').astype(int)

        # Assembly level priority
        level_priority = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}
        df_copy['level_priority'] = df_copy['assembly_level'].map(level_priority).fillna(4)

        # Sort by species, then by priorities, then take first of each species
        df_sorted = df_copy.sort_values([species_col, 'source_priority', 'level_priority'])
        species_subset = df_sorted.groupby(species_col).first().reset_index()

        # Remove the priority columns
        species_subset = species_subset.drop(['source_priority', 'level_priority'], axis=1)

        print(f"Species subset: {len(species_subset):,} representatives")
        return species_subset

    def get_taxids_from_names(self, names: List[str]) -> Dict[str, str]:
        """Get taxids for taxonomic names using taxonkit name2taxid."""
        if not names:
            return {}

        name_to_taxid = {}

        # Create temporary file for names
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
            temp_filename = temp_file.name
            for name in names:
                temp_file.write(f"{name}\n")

        try:
            result = subprocess.run(
                ["taxonkit", "name2taxid", temp_filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            name = parts[0].strip()
                            taxid = parts[1].strip()
                            if taxid and taxid != "":
                                name_to_taxid[name] = taxid
        except Exception as e:
            print(f"Error running taxonkit name2taxid: {e}")
        finally:
            try:
                os.unlink(temp_filename)
            except:
                pass

        print(f"Name->taxid mappings: {len(name_to_taxid)}")
        return name_to_taxid

    def get_lineages_from_taxids(self, taxids: List[str]) -> Dict[str, Tuple[str, str, str]]:
        """Get lineages for taxids using taxonkit lineage -R -t."""
        if not taxids:
            return {}

        lineage_data = {}

        # Create temporary file for taxids
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
            temp_filename = temp_file.name
            for taxid in taxids:
                temp_file.write(f"{taxid}\n")

        try:
            result = subprocess.run(
                ["taxonkit", "lineage", "-R", "-t", temp_filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = line.split('\t')
                        if len(parts) >= 4:
                            taxid = parts[0].strip()
                            lineage = parts[1].strip()
                            lineage_taxids = parts[2].strip()
                            lineage_ranks = parts[3].strip()

                            if lineage and lineage != taxid:
                                lineage_data[taxid] = (lineage, lineage_ranks, lineage_taxids)
        except Exception as e:
            print(f"Error running taxonkit lineage: {e}")
        finally:
            try:
                os.unlink(temp_filename)
            except:
                pass

        print(f"Lineages: {len(lineage_data)}")
        return lineage_data
    
    def generate_counts(self, df: pd.DataFrame, species_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Generate count summaries."""
        print(f"📊 Generating {self.level} counts...")

        # Debug: Check if columns exist
        required_cols = [self.level, 'domain', 'assembly_accession', 'taxid']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"❌ Missing columns in df: {missing_cols}")
            print(f"Available columns: {list(df.columns)}")
            raise KeyError(f"Missing required columns: {missing_cols}")

        # Total genome counts
        genome_counts = df.groupby([self.level, 'domain']).agg({
            'assembly_accession': 'count'
        }).rename(columns={
            'assembly_accession': f'{self.level}_genome_count'
        }).reset_index()

        # Species counts - COUNT UNIQUE SPECIES, NOT ASSEMBLIES
        # FIXED: Previously counted assemblies in species_df, now counts unique species_taxid in full dataset
        species_counts = df.groupby([self.level, 'domain']).agg({
            'species_taxid': 'nunique'
        }).rename(columns={
            'species_taxid': f'{self.level}_species_count'
        }).reset_index()

        # Debug: Verify species counts are reasonable
        total_unique_species = df['species_taxid'].nunique()
        total_species_in_counts = species_counts[f'{self.level}_species_count'].sum()
        print(f"🔍 Debug: Total unique species in dataset: {total_unique_species:,}")
        print(f"🔍 Debug: Sum of species counts across {self.level}s: {total_species_in_counts:,}")
        if total_species_in_counts > total_unique_species:
            print(f"⚠️  Warning: Species counts sum ({total_species_in_counts:,}) > unique species ({total_unique_species:,})")
            print("   This suggests species are being counted multiple times across taxonomic groups")

        # Merge counts
        counts_df = pd.merge(genome_counts, species_counts, on=[self.level, 'domain'], how='outer').fillna(0)

        # Convert to int
        for col in [f'{self.level}_genome_count', f'{self.level}_species_count']:
            counts_df[col] = counts_df[col].astype(int)

        # Sort by genome count
        counts_df = counts_df.sort_values(f'{self.level}_genome_count', ascending=False)

        # Add lineage information using taxonkit
        print(f"Adding lineage information...")

        # Get unique taxonomic names for taxonkit processing
        unique_names = counts_df[self.level].unique().tolist()

        # Get taxids and lineages
        name_to_taxid = self.get_taxids_from_names(unique_names)

        # Add taxid column
        counts_df['taxid'] = counts_df[self.level].map(name_to_taxid)

        # Get lineage information for valid taxids
        valid_taxids = [taxid for taxid in name_to_taxid.values() if taxid]
        lineage_data = self.get_lineages_from_taxids(valid_taxids)

        # Add lineage columns
        lineages, lineage_ranks, lineage_taxids = [], [], []

        for _, row in counts_df.iterrows():
            taxid = row['taxid']
            if taxid and taxid in lineage_data:
                lineage, ranks, taxids = lineage_data[taxid]
                lineages.append(lineage)
                lineage_ranks.append(ranks)
                lineage_taxids.append(taxids)
            else:
                lineages.append("")
                lineage_ranks.append("")
                lineage_taxids.append("")

        counts_df['lineage'] = lineages
        counts_df['lineage_ranks'] = lineage_ranks
        counts_df['lineage_taxids'] = lineage_taxids

        # Calculate percentage columns
        print(f"📊 Calculating percentage columns...")

        # Calculate total genome count and total species count
        total_genome_count = counts_df[f'{self.level}_genome_count'].sum()
        total_species_count = counts_df[f'{self.level}_species_count'].sum()

        # Add genome percentage column
        counts_df[f'{self.level}_genome_percentage'] = (
            counts_df[f'{self.level}_genome_count'] / total_genome_count * 100
        ).round(2)

        # Add species percentage column
        counts_df[f'{self.level}_species_percentage'] = (
            counts_df[f'{self.level}_species_count'] / total_species_count * 100
        ).round(2)

        # Reorder columns for better readability: count, percentage, count, percentage
        base_columns = [self.level, 'domain']
        count_columns = [
            f'{self.level}_genome_count',
            f'{self.level}_genome_percentage',
            f'{self.level}_species_count',
            f'{self.level}_species_percentage'
        ]
        lineage_columns = ['taxid', 'lineage', 'lineage_ranks', 'lineage_taxids']

        # Reorder the dataframe columns
        counts_df = counts_df[base_columns + count_columns + lineage_columns]

        print(f"📊 Total genomes in database: {total_genome_count:,}")
        print(f"📊 Total unique species in database: {total_species_count:,}")
        print(f"Generated {len(counts_df)} {self.level} groups with lineage and percentages")
        return counts_df, species_counts
    
    def save_outputs(self, merged_df: pd.DataFrame, counts_df: pd.DataFrame):
        """Save simplified output files."""
        # 1. Unified counts file (includes species counts)
        counts_file = self.output_dir / f"ncbi_{self.level}_counts.csv"
        counts_df.to_csv(counts_file, index=False)

        # 2. Accessions file with classification
        accessions_file = self.output_dir / f"ncbi_{self.level}_with_accessions.csv"
        accessions_df = merged_df[[
            'assembly_accession', self.level, 'domain', 'organism_name',
            'assembly_level', 'genome_source', 'taxid'
        ]].copy()
        accessions_df = accessions_df.rename(columns={'assembly_accession': 'accession_clean'})
        accessions_df.to_csv(accessions_file, index=False)

        print(f"Saved: {counts_file.name}, {accessions_file.name}")
        print(f"Generated {len(counts_df)} {self.level} groups")
        print(f"Top 3 by genome count:")
        for _, row in counts_df.head(3).iterrows():
            print(f"  {row[self.level]}: {row[f'{self.level}_genome_count']:,} genomes, {row[f'{self.level}_species_count']:,} species")
    
    def run(self):
        """Run the complete parsing workflow."""
        try:
            # Load data
            merged_df, _ = self.load_data()

            # Add genome source classification
            merged_df['genome_source'] = self.classify_genome_source(merged_df)

            # Extract species taxid
            merged_df['species_taxid'] = self.extract_species_taxid(merged_df)

            # Create species subset
            species_df = self.create_species_subset(merged_df, 'species_taxid')

            # Generate counts
            counts_df, _ = self.generate_counts(merged_df, species_df)

            # Save outputs
            self.save_outputs(merged_df, counts_df)
            
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Clean NCBI taxonomic parser")
    parser.add_argument("--level", required=True, choices=['phylum', 'family', 'genus'],
                       help="Taxonomic level to process")
    parser.add_argument("--input-dir", type=Path, help="Input directory (default: auto-detect)")
    parser.add_argument("--output-dir", type=Path, help="Output directory (default: ../csv_ncbi)")
    
    args = parser.parse_args()
    
    # Create and run parser
    ncbi_parser = NCBIParser(args.level, args.input_dir, args.output_dir)
    ncbi_parser.run()

if __name__ == "__main__":
    main()
