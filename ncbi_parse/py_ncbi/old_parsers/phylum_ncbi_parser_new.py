#!/usr/bin/env python3
"""
NCBI Phylum Parser

This script processes NCBI assembly data to extract phylum-level information.
Generates both total genome counts and species subset counts in a single run.
Uses species_taxid for true species-level subsetting.

Input:
- NCBI assembly summary file (00assembly_summary_genbank.txt)
- Taxonomic mapping file (taxid_to_phylum.csv)

Output:
- Total genome counts (ncbi_phylum_counts.csv)
- Species subset counts (ncbi_phylum_species_counts.csv)

Usage:
    python phylum_ncbi_parser.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
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

def classify_genome_source(organism_name: str, excluded_from_refseq: str) -> str:
    """
    Classify genome source as isolate or uncultured based on organism name and exclusion status.
    Uses comprehensive pattern matching from add_genome_source_classification.py
    """
    if pd.isna(organism_name):
        organism_name = ""
    if pd.isna(excluded_from_refseq):
        excluded_from_refseq = ""
    
    organism_lower = str(organism_name).lower()
    excluded_lower = str(excluded_from_refseq).lower()
    
    # Priority 1: Check excluded_from_refseq for metagenome indicators
    metagenome_indicators = [
        'derived from metagenome', 'metagenome', 'single cell', 
        'contaminated', 'derived from single cell'
    ]
    
    if any(indicator in excluded_lower for indicator in metagenome_indicators):
        return 'uncultured'
    
    # Priority 2: Check organism name for uncultured patterns
    uncultured_patterns = [
        'uncultured', 'unculture', 'environmental', 'metagenome', 
        'single cell', 'single-cell', 'mag', 'sag', 'candidatus',
        'marine group', 'deep sea', 'hydrothermal', 'hot spring',
        'bin ', 'contig ', 'scaffold ', 'assembly '
    ]
    
    if any(pattern in organism_lower for pattern in uncultured_patterns):
        return 'uncultured'
    
    return 'isolate'

def create_isolate_preferential_species_subset(df: pd.DataFrame, species_column: str) -> pd.DataFrame:
    """
    Create species subset using isolate-preferential strategy.
    For each species, prioritize isolate genomes over uncultured genomes.
    """
    if 'genome_source' not in df.columns:
        print("⚠️  No genome_source column found, using standard keep='first' strategy")
        return df.drop_duplicates(subset=species_column, keep='first')

    print("🧬 Applying isolate-preferential species subset strategy...")

    # Sort by species and genome_source (isolate comes before uncultured alphabetically)
    df_sorted = df.sort_values([species_column, 'genome_source'])

    # For each species, prioritize isolate over uncultured
    def select_best_genome(group):
        isolates = group[group['genome_source'] == 'isolate']
        if not isolates.empty:
            return isolates.iloc[0]
        else:
            return group.iloc[0]

    # Group by species and apply selection logic
    species_subset = df_sorted.groupby(species_column).apply(select_best_genome).reset_index(drop=True)

    # Count isolate vs uncultured selections
    isolate_count = (species_subset['genome_source'] == 'isolate').sum()
    uncultured_count = (species_subset['genome_source'] == 'uncultured').sum()

    print(f"✅ Species subset results:")
    print(f"   - Isolate genomes selected: {isolate_count:,}")
    print(f"   - Uncultured genomes selected: {uncultured_count:,}")
    print(f"   - Total species: {len(species_subset):,}")
    print(f"   - Isolate preference rate: {(isolate_count / len(species_subset) * 100):.1f}%")

    return species_subset

def get_taxids_from_names(phylum_names: List[str]) -> Dict[str, str]:
    """Get taxids for phylum names using taxonkit name2taxid."""
    if not phylum_names:
        return {}

    print(f"Getting taxids for {len(phylum_names)} phylum names...")
    name_to_taxid = {}

    # Create temporary file for phylum names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in phylum_names:
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

    print(f"Successfully processed {len(name_to_taxid)} phylum name->taxid mappings")
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

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse NCBI assembly data for phylum-level counts")
    parser.add_argument("--input-dir", help="Directory containing input files (default: ../metadata if exists, else script directory)")
    parser.add_argument("--output-dir", help="Directory for output files (default: ../ncbi_parse - safe for comparison)")
    parser.add_argument("--mapping-file", help="Path to taxid mapping file (default: ../taxonomic_mapping/taxid_to_phylum.csv)")
    return parser.parse_args()

def setup_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path]:
    """Set up input and output file paths."""
    script_dir = Path(__file__).resolve().parent
    
    # Default paths - prioritize metadata folder, fallback to script directory
    if args.input_dir:
        input_dir = Path(args.input_dir)
    else:
        metadata_dir = script_dir.parent / "metadata"
        input_dir = metadata_dir if metadata_dir.exists() else script_dir

    output_dir = Path(args.output_dir) if args.output_dir else script_dir.parent / "csv_ncbi"  # csv_ncbi directory
    output_dir.mkdir(exist_ok=True)

    # Set up paths for input files with fallback logic
    assembly_file = input_dir / "00assembly_summary_genbank.txt"
    if not assembly_file.exists() and input_dir != script_dir.parent / "metadata":
        metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
        if metadata_assembly.exists():
            assembly_file = metadata_assembly
            print(f"📁 Using assembly file from metadata folder: {assembly_file}")

    mapping_file = Path(args.mapping_file) if args.mapping_file else script_dir.parent / "taxonomic_mapping" / "taxid_to_phylum.csv"
    
    return assembly_file, mapping_file, output_dir, script_dir

def load_assembly_file(file_path: Path) -> pd.DataFrame:
    """Load NCBI assembly file with memory optimization."""
    try:
        dtype_dict = {
            'taxid': 'str',
            'assembly_accession': 'str',
            'organism_name': 'str',
            'excluded_from_refseq': 'str'
        }
        return pd.read_csv(file_path, sep="\t", skiprows=1, low_memory=False, dtype=dtype_dict)
    except Exception as e:
        print(f"❌ Error loading {file_path}: {str(e)}")
        raise

def load_taxid_mapping(file_path: Path) -> pd.DataFrame:
    """Load taxid mapping file with memory optimization."""
    try:
        dtype_dict = {
            'taxid': 'str',
            'phylum': 'str',
            'domain': 'str'
        }
        return pd.read_csv(file_path, dtype=dtype_dict)
    except Exception as e:
        print(f"❌ Error loading mapping file {file_path}: {str(e)}")
        raise

def main() -> int:
    """Main function to run the phylum parser."""
    args = parse_arguments()

    try:
        # Set up logging
        log_file = setup_logging("phylum_ncbi_parser")
        logging.info("Starting phylum NCBI parser")

        # Set up paths
        assembly_file, mapping_file, output_dir, script_dir = setup_paths(args)
        print(f"📂 Input assembly file: {assembly_file}")
        print(f"📂 Input mapping file: {mapping_file}")
        print(f"📂 Output directory: {output_dir}")
        print(f"📂 Log file: {log_file}")

        # Load assembly file
        print("Loading assembly file...")
        print(f"💾 Memory usage before loading: {get_memory_usage():.2f} GB")
        with tqdm(desc="Loading assembly data", unit="rows") as pbar:
            df = load_assembly_file(assembly_file)
            df = df.rename(columns={"#assembly_accession": "assembly_accession"})
            pbar.update(len(df))
        print(f"✅ Loaded {len(df)} assembly entries")
        print(f"💾 Memory usage after loading assembly: {get_memory_usage():.2f} GB")

        # Load taxid mapping file
        print("Loading taxid mapping file...")
        with tqdm(desc="Loading taxonomy mapping", unit="rows") as pbar:
            mapping_df = load_taxid_mapping(mapping_file)
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

        # Note: Genome source classification is now handled later by loading classified data

        # Debug: Check what columns are available
        print(f"📋 Available columns: {list(merged_df.columns)}")

        # Check for species_taxid column in the assembly file
        print("Checking for species_taxid column in the assembly file...")

        # If species_taxid is not in the data, we need to extract it from the taxid column
        if 'species_taxid' not in merged_df.columns:
            print("⚠️  species_taxid column not found in assembly file")
            print("🔍 Extracting species_taxid from taxid using taxonkit...")

            # Create a temporary file with taxids
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
                temp_filename = temp_file.name
                taxids = merged_df['taxid'].unique().tolist()
                for taxid in taxids:
                    temp_file.write(f"{taxid}\n")

            try:
                # Use taxonkit to get species taxids
                result = subprocess.run(
                    ["taxonkit", "reformat", "--taxid-field", "1", "--format", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}", temp_filename],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )

                if result.returncode == 0 and result.stdout.strip():
                    # Parse the output to get species taxids
                    taxid_to_species_taxid = {}
                    lines = result.stdout.strip().split('\n')

                    for line in lines:
                        parts = line.split('\t')
                        if len(parts) >= 8:
                            taxid = parts[0]
                            species_info = parts[6]  # Species field

                            # Extract species taxid from species info (format: "Species [taxid:123456]")
                            if "[taxid:" in species_info:
                                species_taxid = species_info.split("[taxid:")[1].split("]")[0]
                                taxid_to_species_taxid[taxid] = species_taxid

                    # Add species_taxid column to the dataframe
                    merged_df['species_taxid'] = merged_df['taxid'].map(taxid_to_species_taxid)
                    print(f"✅ Added species_taxid column to {len(merged_df)} entries")
                else:
                    print("⚠️  Failed to get species taxids from taxonkit")
                    print(f"Error: {result.stderr}")
            except Exception as e:
                print(f"⚠️  Error running taxonkit: {e}")
            finally:
                try:
                    os.unlink(temp_filename)
                except:
                    pass

        # Determine species column for subsetting - prioritize species_taxid
        species_column = None
        if 'species_taxid' in merged_df.columns and not merged_df['species_taxid'].isna().all():
            species_column = 'species_taxid'
            print(f"🧬 Using species_taxid column for true species-level subsetting")
            print(f"   This is the most biologically accurate approach")
        elif 'organism_name' in merged_df.columns:
            # Extract species from organism_name (first two words) for true species-level subsetting
            print(f"⚠️  No valid species_taxid found, extracting species from organism_name")
            merged_df['species'] = merged_df['organism_name'].str.split().str[:2].str.join(' ')
            species_column = 'species'
            print(f"🧬 Using extracted species column for species-level subsetting")
        elif 'scientific_name' in merged_df.columns:
            species_column = 'scientific_name'
            print(f"🧬 Using scientific_name column for species-level subsetting")

        if not species_column:
            print("❌ No suitable species column found!")
            return 1

        # Perform genome source classification directly from metadata (no external dependencies)
        print("🧬 Classifying genome sources (isolate vs uncultured) from metadata...")
        with tqdm(desc="Classifying genomes", total=len(merged_df), unit="genomes") as pbar:
            merged_df['genome_source'] = merged_df.apply(
                lambda row: classify_genome_source(row.get('organism_name'), row.get('excluded_from_refseq')),
                axis=1
            )
            pbar.update(len(merged_df))

        # Show genome source distribution
        source_counts = merged_df['genome_source'].value_counts()
        print("📊 Genome source distribution:")
        for source, count in source_counts.items():
            percentage = (count / len(merged_df)) * 100
            print(f"   - {source}: {count:,} ({percentage:.1f}%)")

        # Generate total genome counts
        print("Generating total genome counts...")
        total_counts = merged_df.groupby(["phylum", "domain"]).size().reset_index(name="phylum_count")
        total_counts = total_counts.sort_values("phylum_count", ascending=False)

        # Generate species subset counts using isolate-preferential strategy (matching old script logic)
        print(f"Subsetting to keep one genome per species using '{species_column}' column...")
        species_subset_df = create_isolate_preferential_species_subset(merged_df, species_column)
        print(f"✅ Reduced from {len(merged_df)} genomes to {len(species_subset_df)} unique species")

        # Group by phylum and domain with proper counting logic (matching old script)
        print("Grouping by phylum and domain...")
        if species_column == 'species_taxid':
            # For species_taxid, count unique taxids per phylum
            species_counts = species_subset_df.groupby(["phylum", "domain"]).agg({
                species_column: "nunique"  # Count unique species taxids
            }).reset_index()
        else:
            # For other species columns, count occurrences
            species_counts = species_subset_df.groupby(["phylum", "domain"]).agg({
                species_column: "count"
            }).reset_index()

        # Rename column for clarity
        species_counts = species_counts.rename(columns={species_column: "phylum_species_count"})
        species_counts = species_counts.sort_values("phylum_species_count", ascending=False)

        # Get unique phylum names for taxonkit processing
        unique_phyla = total_counts['phylum'].unique().tolist()
        print(f"\n🔍 Found {len(unique_phyla)} unique phyla")

        # Get taxids and lineages
        print("Getting taxids for phylum names...")
        with tqdm(desc="Processing phylum names", total=len(unique_phyla), unit="phyla") as pbar:
            phylum_to_taxid = get_taxids_from_names(unique_phyla)
            pbar.update(len(unique_phyla))

        # Add taxid column to both count files
        total_counts['taxid'] = total_counts['phylum'].map(phylum_to_taxid)
        species_counts['taxid'] = species_counts['phylum'].map(phylum_to_taxid)

        # Get lineage information
        valid_taxids = [taxid for taxid in phylum_to_taxid.values() if taxid]
        print(f"Getting lineage information for {len(valid_taxids)} valid taxids...")
        with tqdm(desc="Processing lineages", total=len(valid_taxids), unit="taxids") as pbar:
            lineage_data = get_lineages_from_taxids(valid_taxids)
            pbar.update(len(valid_taxids))

        # Add lineage columns to both count files
        for counts_df, file_type in [(total_counts, "total"), (species_counts, "species")]:
            print(f"Adding lineage information to {file_type} counts...")
            lineages, lineage_ranks, lineage_taxids = [], [], []

            for _, row in tqdm(counts_df.iterrows(), desc=f"Processing {file_type} lineage data", total=len(counts_df), unit="rows"):
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

        # Save count files
        total_output_file = output_dir / "ncbi_phylum_counts.csv"
        species_output_file = output_dir / "ncbi_phylum_species_counts.csv"

        print(f"💾 Saving total counts to {total_output_file}...")
        total_counts.to_csv(total_output_file, index=False)
        print(f"✅ Saved {len(total_counts)} phylum-domain combinations (total counts)")

        print(f"💾 Saving species counts to {species_output_file}...")
        species_counts.to_csv(species_output_file, index=False)
        print(f"✅ Saved {len(species_counts)} phylum-domain combinations (species counts)")

        # Print preview of results
        print(f"\n📈 Top 10 phyla by total genome count:")
        for _, row in total_counts.head(10).iterrows():
            print(f"  {row['phylum']} ({row['domain']}): {row['phylum_count']} genomes")

        print(f"\n🧬 Top 10 phyla by species count:")
        for _, row in species_counts.head(10).iterrows():
            print(f"  {row['phylum']} ({row['domain']}): {row['phylum_species_count']} species")

        logging.info("Phylum NCBI parser completed successfully")
        return 0

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        logging.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
