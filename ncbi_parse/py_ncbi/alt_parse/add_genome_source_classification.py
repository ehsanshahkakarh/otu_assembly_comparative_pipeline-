#!/usr/bin/env python3
"""
Add genome source classification (isolate vs uncultured) to NCBI accession files.

This script adds a 'genome_source' column to the existing *_with_accessions.csv files
to classify genomes as either 'isolate' or 'uncultured' based on:
1. Accession patterns (MAG/SAG identifiers)
2. Assembly metadata (if available)
3. Taxonomic patterns (uncultured/environmental/metagenome indicators)

Author: AI Assistant
Date: 2025-06-24
"""

import pandas as pd
import re
import logging
from pathlib import Path
from typing import List, Optional
import json
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class GenomeSourceClassifier:
    """Classify genomes as isolate vs uncultured based on multiple criteria."""

    def __init__(self, metadata_dir: Optional[str] = None):
        """
        Initialize the classifier.

        Args:
            metadata_dir: Path to directory containing assembly metadata files
        """
        self.metadata_dir = Path(metadata_dir).resolve() if metadata_dir else None
        self.assembly_metadata = {}

        # Patterns indicating uncultured genomes (metagenome-based, non-isolate sequences)
        self.uncultured_patterns = [
            # Metagenome-Assembled Genomes
            r'\bMAG\b',  # Metagenome-Assembled Genome
            r'\bSAG\b',  # Single-cell Assembled Genome
            r'metagenome.assembled',
            r'single.cell.assembled',

            # Direct uncultured indicators
            r'\buncultured\b',
            r'\bunculture\b',
            r'environmental.sample',
            r'environmental.clone',
            r'environmental.sequence',

            # Metagenome indicators (key for excluded_from_refseq column)
            r'\bmetagenome\b',
            r'metagenomic',
            r'meta.genome',
            r'derived.from.metagenome',  # Common in excluded_from_refseq
            r'single.cell',
            r'single.amplified.genome',
            r'derived.from.single.cell',  # Common in excluded_from_refseq

            # Assembly/binning identifiers
            r'bin\s*\d+',  # bin123, bin 123
            r'contig\s*\d+',
            r'scaffold\s*\d+',
            r'assembly\s*\d+',

            # Environmental prefixes/suffixes
            r'_ENV_',
            r'_UNC_',
            r'_ENVIR_',
            r'_META_',
            r'_METAG_',

            # Candidatus (uncultured bacteria)
            r'\bcandidatus\b',
            r'\bCa\.\s',  # Candidatus abbreviation

            # Marine/environmental specific
            r'marine\s+group',
            r'deep.sea',
            r'hydrothermal',
            r'hot.spring',
            r'sediment',

            # Uncultured taxonomic groups
            r'uncultured\s+\w+',
            r'unidentified\s+\w+',
            r'unknown\s+\w+',

            # Assembly quality indicators (often in excluded_from_refseq)
            r'contaminated',
            r'fragmented.assembly',
            r'genus.undefined',

            # Genome project indicators
            r'genome.project',
            r'whole.genome.shotgun',
            r'WGS',
        ]

        # Create regex pattern string for pandas operations
        self.uncultured_pattern_str = '|'.join(self.uncultured_patterns)

        # Also keep compiled version for non-pandas operations
        self.uncultured_regex = re.compile(
            self.uncultured_pattern_str,
            re.IGNORECASE
        )

        # Load assembly metadata if available
        if self.metadata_dir and self.metadata_dir.exists():
            self._load_assembly_metadata()
    
    def _load_assembly_metadata(self):
        """Load assembly metadata from JSON files or summary tables."""
        logger.info(f"Loading assembly metadata from {self.metadata_dir}")
        
        # Look for assembly summary files
        summary_files = list(self.metadata_dir.glob("*assembly_summary*.txt"))
        json_files = list(self.metadata_dir.glob("*.json"))
        
        if summary_files:
            self._load_assembly_summary(summary_files[0])
        
        if json_files:
            self._load_json_metadata(json_files)
    
    def _load_assembly_summary(self, summary_file: Path):
        """Load NCBI assembly summary file."""
        try:
            # NCBI assembly summary files are tab-separated
            df = pd.read_csv(summary_file, sep='\t', comment='#', low_memory=False)
            
            # Key columns for classification
            key_columns = [
                'assembly_accession', 'asm_name', 'isolate', 
                'assembly_level', 'genome_rep', 'seq_rel_date'
            ]
            
            # Only keep columns that exist
            available_columns = [col for col in key_columns if col in df.columns]
            
            if 'assembly_accession' in df.columns:
                self.assembly_metadata = df[available_columns].set_index('assembly_accession').to_dict('index')
                logger.info(f"Loaded metadata for {len(self.assembly_metadata)} assemblies")
            
        except Exception as e:
            logger.warning(f"Could not load assembly summary: {e}")
    
    def _load_json_metadata(self, json_files: List[Path]):
        """Load individual JSON metadata files."""
        for json_file in json_files[:100]:  # Limit to avoid memory issues
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                # Extract accession and relevant metadata
                if 'reports' in data:
                    for report in data['reports']:
                        accession = report.get('accession')
                        if accession:
                            self.assembly_metadata[accession] = {
                                'assembly_type': report.get('assembly_info', {}).get('assembly_type'),
                                'assembly_level': report.get('assembly_info', {}).get('assembly_level'),
                                'biosample_attributes': report.get('biosample', {}).get('attributes', [])
                            }
            except Exception as e:
                logger.debug(f"Could not load {json_file}: {e}")
    
    def classify_genomes_vectorized(self, df: pd.DataFrame) -> pd.Series:
        """
        Vectorized classification of genomes as 'isolate' or 'uncultured'.

        Args:
            df: DataFrame with 'accession_clean' and optional taxonomic/metadata columns

        Returns:
            Series with genome source classifications
        """
        # Initialize all as 'isolate'
        classifications = pd.Series(['isolate'] * len(df), index=df.index)

        # Get accession and taxonomic columns
        accessions = df['accession_clean'].fillna('')

        # Find taxonomic column
        taxon_col = None
        for col in ['family', 'genus', 'phylum', 'species', 'organism_name']:
            if col in df.columns:
                taxon_col = col
                break

        taxon_names = df[taxon_col].fillna('') if taxon_col else pd.Series([''] * len(df))

        # Check for excluded_from_refseq column (contains metagenome information)
        excluded_col = None
        for col in ['excluded_from_refseq', 'assembly_notes', 'notes']:
            if col in df.columns:
                excluded_col = col
                break

        excluded_notes = df[excluded_col].fillna('') if excluded_col else pd.Series([''] * len(df))

        # Log which columns are being used for classification
        logger.info(f"Classification columns: accession_clean=✓, taxonomic={taxon_col or 'None'}, excluded_notes={excluded_col or 'None'}")

        # Vectorized pattern matching for accessions
        accession_uncultured = accessions.str.contains(
            self.uncultured_pattern_str, case=False, na=False, regex=True
        )

        # Vectorized pattern matching for taxonomic names
        taxon_uncultured = taxon_names.str.contains(
            self.uncultured_pattern_str, case=False, na=False, regex=True
        )

        # Note: excluded_from_refseq is handled by dedicated function _classify_from_excluded_refseq

        # Apply excluded_from_refseq classification first (highest priority for metagenome detection)
        excluded_classifications = self._classify_from_excluded_refseq(excluded_notes)
        excluded_mask = excluded_classifications.notna()
        classifications.loc[excluded_mask] = excluded_classifications.loc[excluded_mask]

        # Apply metadata-based classification if available (second priority)
        if self.assembly_metadata:
            metadata_classifications = self._classify_from_metadata_vectorized(accessions)
            metadata_mask = metadata_classifications.notna() & ~excluded_mask  # Don't override excluded_from_refseq
            classifications.loc[metadata_mask] = metadata_classifications.loc[metadata_mask]

            # For entries without metadata or excluded_from_refseq, use pattern-based classification
            no_special_classification_mask = ~excluded_mask & ~metadata_mask
            pattern_uncultured = (accession_uncultured | taxon_uncultured) & no_special_classification_mask
            classifications.loc[pattern_uncultured] = 'uncultured'
        else:
            # Use pattern-based classification for entries not classified by excluded_from_refseq
            no_excluded_mask = ~excluded_mask
            pattern_uncultured = (accession_uncultured | taxon_uncultured) & no_excluded_mask
            classifications.loc[pattern_uncultured] = 'uncultured'

        return classifications

    def _classify_from_excluded_refseq(self, excluded_notes: pd.Series) -> pd.Series:
        """
        Classify genomes based on excluded_from_refseq column content.

        This function specifically handles the NCBI excluded_from_refseq column
        which contains key information about metagenome-derived assemblies.

        Args:
            excluded_notes: Series containing excluded_from_refseq values

        Returns:
            Series with classifications ('uncultured', 'isolate', or None)
        """
        classifications = pd.Series([None] * len(excluded_notes), index=excluded_notes.index)

        # Check for metagenome indicators in excluded_from_refseq
        metagenome_mask = excluded_notes.str.contains(
            r'derived.from.metagenome|metagenome|single.cell|contaminated.*metagenome',
            case=False, na=False, regex=True
        )

        # Mark as uncultured if metagenome indicators found
        classifications.loc[metagenome_mask] = 'uncultured'

        return classifications

    def _classify_from_metadata_vectorized(self, accessions: pd.Series) -> pd.Series:
        """Vectorized classification based on assembly metadata."""
        classifications = pd.Series([None] * len(accessions), index=accessions.index)

        for idx, accession in accessions.items():
            if accession in self.assembly_metadata:
                metadata = self.assembly_metadata[accession]

                # Check assembly type
                assembly_type = metadata.get('assembly_type', '').lower()
                if 'metagenome' in assembly_type:
                    classifications.iloc[idx] = 'uncultured'
                    continue

                # Check biosample attributes if available
                biosample_attrs = metadata.get('biosample_attributes', [])
                for attr in biosample_attrs:
                    name = attr.get('name', '').lower()
                    value = attr.get('value', '').lower()

                    if name == 'sample_type' and 'metagenome' in value:
                        classifications.iloc[idx] = 'uncultured'
                        break
                    if name == 'isolation_source' and any(term in value for term in ['metagenome', 'environmental']):
                        classifications.iloc[idx] = 'uncultured'
                        break
                    if name == 'sample_type' and 'isolate' in value:
                        classifications.iloc[idx] = 'isolate'
                        break

        return classifications
    
    def process_accession_file(self, input_file: Path, output_file: Optional[Path] = None):
        """
        Process an accession file and add genome source classification.

        Args:
            input_file: Path to input CSV file
            output_file: Path to output CSV file (default: adds _classified suffix)
        """
        # Resolve paths properly
        input_file = Path(input_file).resolve()
        if not output_file:
            output_file = input_file.parent / f"{input_file.stem}_classified{input_file.suffix}"
        else:
            output_file = Path(output_file).resolve()

        logger.info(f"Processing {input_file}")

        # Check if input file exists
        if not input_file.exists():
            logger.error(f"Input file not found: {input_file}")
            return None

        # Read the input file
        try:
            df = pd.read_csv(input_file)
        except Exception as e:
            logger.error(f"Error reading {input_file}: {e}")
            return None

        # Check required columns
        if 'accession_clean' not in df.columns:
            logger.error(f"Required column 'accession_clean' not found in {input_file}")
            logger.info(f"Available columns: {list(df.columns)}")
            return None

        # Vectorized classification
        logger.info(f"Classifying {len(df)} genomes using vectorized operations...")
        df['genome_source'] = self.classify_genomes_vectorized(df)

        # Generate summary statistics
        source_counts = df['genome_source'].value_counts()
        logger.info(f"Classification summary for {input_file.name}:")
        for source, count in source_counts.items():
            percentage = (count / len(df)) * 100
            logger.info(f"  {source}: {count:,} ({percentage:.1f}%)")

        # Additional statistics if excluded_from_refseq column was used
        if 'excluded_from_refseq' in df.columns:
            excluded_with_metagenome = df['excluded_from_refseq'].str.contains(
                r'metagenome', case=False, na=False
            ).sum()
            logger.info(f"  Entries with 'metagenome' in excluded_from_refseq: {excluded_with_metagenome:,}")

            uncultured_from_excluded = df[
                (df['genome_source'] == 'uncultured') &
                df['excluded_from_refseq'].str.contains(r'metagenome', case=False, na=False)
            ]
            logger.info(f"  Uncultured classified from excluded_from_refseq: {len(uncultured_from_excluded):,}")

        # Save the result
        try:
            df.to_csv(output_file, index=False)
            logger.info(f"Saved classified data to {output_file}")
        except Exception as e:
            logger.error(f"Error saving to {output_file}: {e}")
            return None

        return df

    def test_assembly_file_classification(self, assembly_file: Path, output_file: Optional[Path] = None):
        """
        Test classification directly on NCBI assembly file to verify excluded_from_refseq function.

        Args:
            assembly_file: Path to NCBI assembly summary file
            output_file: Path to output CSV file
        """
        logger.info(f"Testing classification on assembly file: {assembly_file}")

        # Read assembly file
        try:
            df = pd.read_csv(assembly_file, sep='\t', skiprows=1, low_memory=False)
            df = df.rename(columns={'#assembly_accession': 'assembly_accession'})
        except Exception as e:
            logger.error(f"Error reading assembly file: {e}")
            return None

        logger.info(f"Loaded {len(df)} assembly entries")

        # Prepare columns for classification
        df['accession_clean'] = df['assembly_accession']

        # Check available columns
        logger.info(f"Available columns: {list(df.columns)}")

        # Run classification
        df['genome_source'] = self.classify_genomes_vectorized(df)

        # Generate summary
        source_counts = df['genome_source'].value_counts()
        logger.info(f"Classification results:")
        for source, count in source_counts.items():
            percentage = (count / len(df)) * 100
            logger.info(f"  {source}: {count:,} ({percentage:.1f}%)")

        # Show some examples of uncultured entries
        uncultured_df = df[df['genome_source'] == 'uncultured']
        if len(uncultured_df) > 0:
            logger.info(f"Sample uncultured entries:")
            for _, row in uncultured_df.head(5).iterrows():
                logger.info(f"  {row['accession_clean']}: {row.get('excluded_from_refseq', 'N/A')}")

        # Save if output file specified
        if output_file:
            test_df = df[['accession_clean', 'organism_name', 'excluded_from_refseq', 'genome_source']].copy()
            test_df.to_csv(output_file, index=False)
            logger.info(f"Saved test results to {output_file}")

        return df

def main():
    """Main function to process all accession files."""
    parser = argparse.ArgumentParser(
        description="Add genome source classification to NCBI accession files"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=".",
        help="Directory containing *_with_accessions.csv files (default: current directory)"
    )
    parser.add_argument(
        "--metadata-dir",
        type=str,
        help="Directory containing assembly metadata files"
    )
    parser.add_argument(
        "--output-suffix",
        type=str,
        default="_classified",
        help="Suffix to add to output files (default: _classified)"
    )
    parser.add_argument(
        "--files",
        nargs="+",
        help="Specific files to process (default: all *_with_accessions.csv files)"
    )
    parser.add_argument(
        "--test-assembly-file",
        type=str,
        help="Test mode: classify directly from NCBI assembly file to verify excluded_from_refseq function"
    )

    args = parser.parse_args()

    # Handle test mode
    if args.test_assembly_file:
        assembly_file = Path(args.test_assembly_file).resolve()
        if not assembly_file.exists():
            logger.error(f"Assembly file not found: {assembly_file}")
            return

        classifier = GenomeSourceClassifier()
        output_file = assembly_file.parent / "test_classification_results.csv"
        classifier.test_assembly_file_classification(assembly_file, output_file)
        return

    # Set up paths with proper resolution
    input_dir = Path(args.input_dir).resolve()
    metadata_dir = Path(args.metadata_dir).resolve() if args.metadata_dir else None

    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Metadata directory: {metadata_dir}")

    # Check if input directory exists
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        return

    # Auto-detect csv_ncbi subdirectory if running from parent directory
    if input_dir.name != "csv_ncbi" and (input_dir / "csv_ncbi").exists():
        csv_ncbi_dir = input_dir / "csv_ncbi"
        logger.info(f"Auto-detected csv_ncbi directory: {csv_ncbi_dir}")
        input_dir = csv_ncbi_dir

    # Initialize classifier
    classifier = GenomeSourceClassifier(
        metadata_dir if metadata_dir and metadata_dir.exists() else None
    )

    # Find accession files
    if args.files:
        # Process specific files
        accession_files = [input_dir / f for f in args.files]
        # Verify files exist
        accession_files = [f for f in accession_files if f.exists()]
    else:
        # Find all accession files
        accession_files = list(input_dir.glob("*_with_accessions.csv"))

    if not accession_files:
        logger.error(f"No accession files found in {input_dir}")
        if args.files:
            logger.error(f"Specified files: {args.files}")
        else:
            logger.error("Looking for files matching pattern: *_with_accessions.csv")
            # Show what files are actually in the directory
            all_csv_files = list(input_dir.glob("*.csv"))
            if all_csv_files:
                logger.info(f"Available CSV files in {input_dir}:")
                for f in all_csv_files[:10]:  # Show first 10 files
                    logger.info(f"  {f.name}")
                if len(all_csv_files) > 10:
                    logger.info(f"  ... and {len(all_csv_files) - 10} more files")
            else:
                logger.error(f"No CSV files found in {input_dir}")
        return

    logger.info(f"Found {len(accession_files)} accession files to process")

    # Process each file
    success_count = 0
    for file_path in accession_files:
        try:
            result = classifier.process_accession_file(file_path)
            if result is not None:
                success_count += 1
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")

    logger.info(f"Classification complete! Successfully processed {success_count}/{len(accession_files)} files")

if __name__ == "__main__":
    main()
