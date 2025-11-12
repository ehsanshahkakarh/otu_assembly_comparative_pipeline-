#!/usr/bin/env python3
"""
Unified NCBI Taxonomic Parser

This script processes NCBI assembly data to extract taxonomic information at any level.
Replaces the redundant phylum/family/genus specific parsers with a single unified approach.
Includes built-in genome source classification (isolate vs uncultured) and generates both
total genome counts and species subset counts in a comprehensive output.

Input:
- NCBI assembly summary file (00assembly_summary_genbank.txt)
- Taxonomic mapping file (taxid_to_{level}.csv)

Output:
- Comprehensive counts CSV with both genome and species counts (ncbi_{level}_comprehensive_counts.csv)
- Classified accessions CSV with all genomes (ncbi_{level}_with_accessions_classified.csv)
- Species subset classified accessions CSV (ncbi_{level}_species_subset_with_accessions_classified.csv)

Usage:
    python unified_ncbi_parser.py --level phylum [--species-subset] [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
    python unified_ncbi_parser.py --level family [--species-subset]
    python unified_ncbi_parser.py --level genus [--species-subset]

Note: --species-subset flag now only affects which accession file is prioritized;
      comprehensive counts always include both genome and species counts.
"""

import pandas as pd
import os
import argparse
import sys
import subprocess
import tempfile
import psutil
import logging
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from tqdm import tqdm

def get_memory_usage():
    """Get current memory usage in GB."""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return memory_info.rss / (1024 ** 3)  # Convert to GB

def setup_logging(script_name: str) -> Path:
    """Set up logging to error_log directory."""
    script_dir = Path(__file__).resolve().parent
    log_dir = script_dir / "error_log"
    log_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f"{script_name}_{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    return log_file

def get_taxids_from_names(names: List[str]) -> Dict[str, str]:
    """Get taxids for taxonomic names using taxonkit name2taxid."""
    if not names:
        return {}

    print(f"Getting taxids for {len(names)} taxonomic names...")
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

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines):
                if not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()
                    if taxid and taxid != "":
                        name_to_taxid[name] = taxid
                        if i < 3:  # Show first few examples
                            print(f"  Example: {name} -> {taxid}")
    except Exception as e:
        print(f"Error running taxonkit name2taxid: {e}")
    finally:
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(name_to_taxid)} name->taxid mappings")
    return name_to_taxid

def get_lineages_from_taxids(taxids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """Get lineages for taxids using taxonkit lineage -R -t."""
    if not taxids:
        return {}

    print(f"Getting lineages for {len(taxids)} taxids...")
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

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines):
                if not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip()
                    lineage_taxids = parts[2].strip()
                    lineage_ranks = parts[3].strip()
                    
                    if lineage and lineage != taxid:
                        lineage_data[taxid] = (lineage, lineage_ranks, lineage_taxids)
                        if i < 3:  # Show first few examples
                            print(f"  Example: Taxid {taxid}: {lineage}")
    except Exception as e:
        print(f"Error running taxonkit lineage: {e}")
    finally:
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(lineage_data)} lineages")
    return lineage_data

def create_isolate_preferential_species_subset(df: pd.DataFrame, species_column: str) -> pd.DataFrame:
    """Create species subset using isolate-preferential strategy."""
    print("Creating isolate-preferential species subset...")
    
    if 'genome_source' not in df.columns:
        print("⚠️  No genome_source column found, using first occurrence per species")
        return df.drop_duplicates(subset=[species_column], keep='first')
    
    # Separate isolate and uncultured genomes
    isolate_df = df[df['genome_source'] == 'isolate']
    uncultured_df = df[df['genome_source'] == 'uncultured']
    
    # Get one isolate per species (if available)
    isolate_subset = isolate_df.drop_duplicates(subset=[species_column], keep='first')
    
    # Get species that don't have isolates
    species_with_isolates = set(isolate_subset[species_column])
    uncultured_only = uncultured_df[~uncultured_df[species_column].isin(species_with_isolates)]
    uncultured_subset = uncultured_only.drop_duplicates(subset=[species_column], keep='first')
    
    # Combine isolate-preferred subset
    species_subset = pd.concat([isolate_subset, uncultured_subset], ignore_index=True)
    
    print(f"✅ Species subset: {len(species_subset)} genomes from {species_subset[species_column].nunique()} species")
    print(f"   - Isolate genomes: {len(isolate_subset)}")
    print(f"   - Uncultured genomes: {len(uncultured_subset)}")
    
    return species_subset

class GenomeSourceClassifier:
    """Classify genomes as isolate vs uncultured based on multiple criteria."""

    def __init__(self):
        """Initialize the classifier with comprehensive uncultured patterns."""
        # Patterns indicating uncultured genomes (from add_genome_source_classification.py)
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

    def classify_genomes_vectorized(self, df: pd.DataFrame) -> pd.Series:
        """
        Vectorized classification of genomes as 'isolate' or 'uncultured'.

        Priority order:
        1. excluded_from_refseq column (highest priority for metagenome detection)
        2. Pattern-based classification (accession + taxonomic names)
        """
        # Initialize all as 'isolate'
        classifications = pd.Series(['isolate'] * len(df), index=df.index)

        # Get accession and organism name columns
        accessions = df['assembly_accession'].fillna('') if 'assembly_accession' in df.columns else pd.Series([''] * len(df))
        organism_names = df['organism_name'].fillna('') if 'organism_name' in df.columns else pd.Series([''] * len(df))

        # Check for excluded_from_refseq column (contains metagenome information)
        excluded_notes = df['excluded_from_refseq'].fillna('') if 'excluded_from_refseq' in df.columns else pd.Series([''] * len(df))

        # Priority 1: excluded_from_refseq classification (highest priority)
        excluded_classifications = self._classify_from_excluded_refseq(excluded_notes)
        excluded_mask = excluded_classifications.notna()
        classifications.loc[excluded_mask] = excluded_classifications.loc[excluded_mask]

        # Priority 2: Pattern-based classification for entries not classified by excluded_from_refseq
        no_excluded_mask = ~excluded_mask

        # Vectorized pattern matching for accessions
        accession_uncultured = accessions.str.contains(
            self.uncultured_pattern_str, case=False, na=False, regex=True
        )

        # Vectorized pattern matching for organism names
        organism_uncultured = organism_names.str.contains(
            self.uncultured_pattern_str, case=False, na=False, regex=True
        )

        # Apply pattern-based classification
        pattern_uncultured = (accession_uncultured | organism_uncultured) & no_excluded_mask
        classifications.loc[pattern_uncultured] = 'uncultured'

        return classifications

    def _classify_from_excluded_refseq(self, excluded_notes: pd.Series) -> pd.Series:
        """
        Classify genomes based on excluded_from_refseq column content.
        This handles the NCBI excluded_from_refseq column which contains
        key information about metagenome-derived assemblies.
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

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Unified NCBI taxonomic parser for any taxonomic level")
    parser.add_argument("--level", required=True, choices=['phylum', 'family', 'genus'],
                       help="Taxonomic level to process")
    parser.add_argument("--species-subset", action='store_true',
                       help="Create species subset (one genome per species)")
    parser.add_argument("--input-dir", help="Directory containing input files (default: ../metadata if exists, else script directory)")
    parser.add_argument("--output-dir", help="Directory for output files (default: ../ncbi_parse - safe for comparison)")
    parser.add_argument("--mapping-file", help="Path to taxid mapping file")
    return parser.parse_args()

def setup_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path]:
    """Set up input and output file paths."""
    script_dir = Path(__file__).resolve().parent
    level = args.level

    # Default paths
    if args.input_dir:
        input_dir = Path(args.input_dir)
    else:
        metadata_dir = script_dir.parent / "metadata"
        input_dir = metadata_dir if metadata_dir.exists() else script_dir

    # Default output to ncbi_parse directory (parent of py_ncbi) to avoid overriding existing csv_ncbi files
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = script_dir.parent  # This is the ncbi_parse directory
    output_dir.mkdir(exist_ok=True)
    
    # Input files
    assembly_file = input_dir / "00assembly_summary_genbank.txt"
    if not assembly_file.exists():
        metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
        if metadata_assembly.exists():
            assembly_file = metadata_assembly
    
    mapping_file = Path(args.mapping_file) if args.mapping_file else \
                   script_dir.parent / "taxonomic_mapping" / f"taxid_to_{level}.csv"
    
    # Output files - always generate comprehensive counts with both genome and species counts
    output_file = output_dir / f"ncbi_{level}_comprehensive_counts.csv"
    
    return assembly_file, mapping_file, output_file, output_dir

def load_assembly_file(file_path: Path) -> pd.DataFrame:
    """Load NCBI assembly file with memory optimization."""
    try:
        dtype_dict = {
            'taxid': 'str',
            'assembly_accession': 'str',
            'organism_name': 'str',
            'excluded_from_refseq': 'str'
        }
        df = pd.read_csv(file_path, sep="\t", skiprows=1, low_memory=False, dtype=dtype_dict)
        return df.rename(columns={"#assembly_accession": "assembly_accession"})
    except Exception as e:
        print(f"❌ Error loading {file_path}: {str(e)}")
        raise

def load_taxid_mapping(file_path: Path, level: str) -> pd.DataFrame:
    """Load taxid mapping file with memory optimization."""
    try:
        dtype_dict = {
            'taxid': 'str',
            level: 'str',
            'domain': 'str'
        }
        return pd.read_csv(file_path, dtype=dtype_dict)
    except Exception as e:
        print(f"❌ Error loading mapping file {file_path}: {str(e)}")
        raise

def main() -> int:
    """Main function to run the unified parser."""
    args = parse_arguments()
    level = args.level

    try:
        # Set up logging
        script_name = f"{level}_ncbi_parser{'_species_subset' if args.species_subset else ''}"
        log_file = setup_logging(script_name)
        logging.info(f"Starting unified NCBI parser for {level} level")

        # Set up paths
        assembly_file, mapping_file, output_file, output_dir = setup_paths(args)
        print(f"📂 Input assembly file: {assembly_file}")
        print(f"📂 Input mapping file: {mapping_file}")
        print(f"📂 Output file: {output_file}")
        print(f"📂 Log file: {log_file}")

        # Load assembly file
        print("Loading assembly file...")
        print(f"💾 Memory usage before loading: {get_memory_usage():.2f} GB")
        with tqdm(desc="Loading assembly data", unit="rows") as pbar:
            df = load_assembly_file(assembly_file)
            pbar.update(len(df))
        print(f"✅ Loaded {len(df)} assembly entries")
        print(f"💾 Memory usage after loading assembly: {get_memory_usage():.2f} GB")

        # Load taxid mapping file
        print("Loading taxid mapping file...")
        with tqdm(desc="Loading taxonomy mapping", unit="rows") as pbar:
            mapping_df = load_taxid_mapping(mapping_file, level)
            pbar.update(len(mapping_df))
        print(f"✅ Loaded {len(mapping_df)} mapping entries")

        # Ensure consistent data types and remove duplicates
        df['taxid'] = df['taxid'].astype(str)
        mapping_df['taxid'] = mapping_df['taxid'].astype(str)

        if mapping_df['taxid'].duplicated().any():
            print("⚠️  Found duplicated taxids in mapping file, removing duplicates...")
            mapping_df = mapping_df.drop_duplicates(subset=['taxid'], keep='first')
            print(f"✅ Cleaned mapping data: {len(mapping_df)} unique entries")

        # Merge assembly data with mapping data
        print("Merging assembly data with taxonomic mapping...")
        merged_df = pd.merge(df, mapping_df, on='taxid', how='inner')
        print(f"✅ Merged data: {len(merged_df)} entries")
        print(f"💾 Memory usage after merge: {get_memory_usage():.2f} GB")

        # Add genome source classification using the sophisticated classifier
        print("Classifying genome sources...")
        classifier = GenomeSourceClassifier()
        with tqdm(desc="Classifying genomes", total=len(merged_df), unit="rows") as pbar:
            merged_df['genome_source'] = classifier.classify_genomes_vectorized(merged_df)
            pbar.update(len(merged_df))

        # Show classification summary
        source_counts = merged_df['genome_source'].value_counts()
        print(f"📊 Genome source classification summary:")
        for source, count in source_counts.items():
            percentage = (count / len(merged_df)) * 100
            print(f"   {source}: {count:,} ({percentage:.1f}%)")

        # Note: Species subset is now handled in the counting logic above
        # The merged_df here contains all genomes for the comprehensive accession file

        # Generate comprehensive counts with both total genomes and species counts
        print(f"Generating comprehensive {level} counts...")

        # Total genome counts
        total_counts = merged_df.groupby([level, "domain"]).size().reset_index(name=f"{level}_genome_count")

        # Species subset counts (isolate-preferential)
        if 'species' not in merged_df.columns:
            # Extract species from organism_name (first two words)
            merged_df['species'] = merged_df['organism_name'].str.split().str[:2].str.join(' ')

        print("Creating isolate-preferential species subset for counting...")
        species_subset_df = create_isolate_preferential_species_subset(merged_df, 'species')
        species_counts = species_subset_df.groupby([level, "domain"]).size().reset_index(name=f"{level}_species_count")

        # Merge total and species counts
        counts = pd.merge(total_counts, species_counts, on=[level, "domain"], how="outer")
        counts = counts.fillna(0).astype({f"{level}_genome_count": int, f"{level}_species_count": int})
        counts = counts.sort_values(f"{level}_genome_count", ascending=False)

        # Print preview
        print(f"\n📈 Top 10 {level}s by genome count:")
        for _, row in counts.head(10).iterrows():
            print(f"  {row[level]} ({row['domain']}): {row[f'{level}_genome_count']} genomes, {row[f'{level}_species_count']} species")

        # Get taxids and lineages
        unique_names = counts[level].unique().tolist()
        print(f"\n🔍 Found {len(unique_names)} unique {level}s")

        print(f"Getting taxids for {level} names...")
        with tqdm(desc=f"Processing {level} names", total=len(unique_names), unit=f"{level}s") as pbar:
            name_to_taxid = get_taxids_from_names(unique_names)
            pbar.update(len(unique_names))

        counts['taxid'] = counts[level].map(name_to_taxid)

        # Get lineage information
        valid_taxids = [taxid for taxid in name_to_taxid.values() if taxid]
        print(f"Getting lineage information for {len(valid_taxids)} valid taxids...")
        with tqdm(desc="Processing lineages", total=len(valid_taxids), unit="taxids") as pbar:
            lineage_data = get_lineages_from_taxids(valid_taxids)
            pbar.update(len(valid_taxids))

        # Add lineage columns
        print("Adding lineage information to results...")
        lineages, lineage_ranks, lineage_taxids = [], [], []

        for _, row in tqdm(counts.iterrows(), desc="Processing lineage data", total=len(counts), unit="rows"):
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

        counts['lineage'] = lineages
        counts['lineage_ranks'] = lineage_ranks
        counts['lineage_taxids'] = lineage_taxids

        # Save summary counts file
        print(f"💾 Saving summary to {output_file}...")
        counts.to_csv(output_file, index=False)
        print(f"✅ Saved {len(counts)} {level}-domain combinations")

        # Generate comprehensive accessions file (all genomes with classification)
        comprehensive_accession_file = output_dir / f"ncbi_{level}_with_accessions_classified.csv"
        print(f"💾 Generating comprehensive classified accessions file: {comprehensive_accession_file}...")

        # Always include genome_source since it's the main value-add
        columns_to_include = ['assembly_accession', level, 'domain', 'taxid', 'genome_source']
        if 'excluded_from_refseq' in merged_df.columns:
            columns_to_include.append('excluded_from_refseq')
        if 'species' in merged_df.columns:
            columns_to_include.append('species')

        comprehensive_df = merged_df[columns_to_include].copy()
        comprehensive_df = comprehensive_df.rename(columns={'assembly_accession': 'accession_clean'})

        comprehensive_df.to_csv(comprehensive_accession_file, index=False)
        print(f"✅ Saved {len(comprehensive_df)} comprehensive classified accession mappings")

        # Generate species subset accessions file (isolate-preferential)
        if args.species_subset or True:  # Always generate species subset file for completeness
            species_subset_accession_file = output_dir / f"ncbi_{level}_species_subset_with_accessions_classified.csv"
            print(f"💾 Generating species subset classified accessions file: {species_subset_accession_file}...")

            # Prepare species subset DataFrame for accessions
            species_subset_accession_df = species_subset_df[columns_to_include].copy()
            species_subset_accession_df = species_subset_accession_df.rename(columns={'assembly_accession': 'accession_clean'})

            species_subset_accession_df.to_csv(species_subset_accession_file, index=False)
            print(f"✅ Saved {len(species_subset_accession_df)} species subset classified accession mappings")

        # Show genome source breakdown
        accession_source_counts = comprehensive_df['genome_source'].value_counts()
        print(f"📊 Comprehensive accession file genome source breakdown:")
        for source, count in accession_source_counts.items():
            percentage = (count / len(comprehensive_df)) * 100
            print(f"   {source}: {count:,} ({percentage:.1f}%)")

        if 'species_subset_accession_df' in locals():
            species_source_counts = species_subset_accession_df['genome_source'].value_counts()
            print(f"📊 Species subset accession file genome source breakdown:")
            for source, count in species_source_counts.items():
                percentage = (count / len(species_subset_accession_df)) * 100
                print(f"   {source}: {count:,} ({percentage:.1f}%)")

        logging.info(f"Unified NCBI parser completed successfully for {level} level")
        return 0

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        logging.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
