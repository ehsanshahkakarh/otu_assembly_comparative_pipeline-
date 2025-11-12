#!/usr/bin/env python3
"""
NCBI Family Parser (Species Subset Version)

This script processes NCBI assembly data to extract family-level information.
It uses taxonomic mapping and groups by family and domain, similar to the eukprot parser format.
This version subsets data to one genome per species before counting.

Input:
- NCBI assembly summary file (00assembly_summary_genbank.txt)
- Taxonomic mapping file (taxid_to_family.csv)

Output:
- Summary CSV with family species counts (ncbi_family_species_counts.csv)

Usage:
    python family_ncbi_parser_species_subset.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
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

def create_isolate_preferential_species_subset(df: pd.DataFrame, species_column: str) -> pd.DataFrame:
    """
    Create species subset using isolate-preferential strategy.

    Strategy:
    1. For each species, prioritize isolate genomes over uncultured genomes
    2. If a species has both isolate and uncultured genomes, keep the isolate
    3. If a species has only uncultured genomes, keep the first uncultured
    4. If genome_source column doesn't exist, fall back to keep='first'

    Args:
        df: DataFrame with species and genome_source information
        species_column: Column name containing species information

    Returns:
        DataFrame with one genome per species (isolate-preferential)
    """
    if 'genome_source' not in df.columns:
        print("⚠️  No genome_source column found, using standard keep='first' strategy")
        return df.drop_duplicates(subset=species_column, keep='first')

    print("🧬 Applying isolate-preferential species subset strategy...")

    # Sort by species and genome_source (isolate comes before uncultured alphabetically)
    # This ensures isolates are selected first when available
    df_sorted = df.sort_values([species_column, 'genome_source'])

    # For each species, prioritize isolate over uncultured
    def select_best_genome(group):
        # If any isolate genomes exist, select the first isolate
        isolates = group[group['genome_source'] == 'isolate']
        if not isolates.empty:
            return isolates.iloc[0]
        # Otherwise, select the first uncultured genome
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

def get_taxids_from_names(family_names: List[str]) -> Dict[str, str]:
    """
    Get taxids for family names using taxonkit name2taxid.

    Args:
        family_names: List of family names

    Returns:
        Dictionary mapping family names to their taxids
    """
    if not family_names:
        return {}

    print(f"Getting taxids for {len(family_names)} family names...")

    name_to_taxid = {}

    # Create temporary file for family names
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        for name in family_names:
            temp_file.write(f"{name}\n")

    try:
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit name2taxid {temp_filename}")

        # Run taxonkit name2taxid command
        result = subprocess.run(
            ["taxonkit", "name2taxid", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} name2taxid results")

            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    taxid = parts[1].strip()

                    if taxid and taxid != "":
                        name_to_taxid[name] = taxid

                        # Show first few examples
                        if i < 3:
                            print(f"  Example: {name} -> {taxid}")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No results from taxonkit name2taxid")

    except Exception as e:
        print(f"Error running taxonkit name2taxid: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(name_to_taxid)} family name->taxid mappings")
    return name_to_taxid

def get_lineages_from_taxids(taxids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Get lineages for taxids using taxonkit lineage -R -t.

    Args:
        taxids: List of taxids

    Returns:
        Dictionary mapping taxids to (lineage, lineage_ranks, lineage_taxids)
    """
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
        print(f"Created temporary file: {temp_filename}")
        print(f"Running: taxonkit lineage -R -t {temp_filename}")

        # Run taxonkit lineage command with -R flag for ranks and -t flag for taxids
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"Taxonkit return code: {result.returncode}")
        if result.stderr:
            print(f"Taxonkit stderr: {result.stderr}")

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            print(f"Got {len(lines)} lineage results")

            for i, line in enumerate(lines):
                if not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    taxid = parts[0].strip()
                    lineage = parts[1].strip()
                    lineage_taxids = parts[2].strip()
                    lineage_ranks = parts[3].strip()

                    # Store the lineage data
                    if lineage and lineage != taxid:
                        lineage_data[taxid] = (
                            lineage,
                            lineage_ranks,
                            lineage_taxids
                        )

                        # Show first few examples
                        if i < 3:
                            print(f"  Example: Taxid {taxid}:")
                            print(f"    Lineage: {lineage}")
                            print(f"    Ranks: {lineage_ranks}")
                            print(f"    Taxids: {lineage_taxids}")
                else:
                    print(f"Unexpected line format: {line}")
        else:
            print("No lineage results from taxonkit")

    except Exception as e:
        print(f"Error running taxonkit lineage: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

    print(f"Successfully processed {len(lineage_data)} lineages")
    return lineage_data

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse NCBI assembly data for family-level information")
    parser.add_argument("--input-dir", help="Directory containing input files (default: ../metadata if exists, else script directory)")
    parser.add_argument("--output-dir", help="Directory for output files (default: ../csv_ncbi)")
    parser.add_argument("--mapping-file", help="Path to taxid mapping file (default: ../taxonomic_mapping/taxid_to_family.csv)")
    return parser.parse_args()

def setup_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    """Set up input and output file paths."""
    script_dir = Path(__file__).resolve().parent

    # Default paths - prioritize metadata folder, fallback to script directory
    if args.input_dir:
        input_dir = Path(args.input_dir)
    else:
        # Try metadata folder first, then script directory
        metadata_dir = script_dir.parent / "metadata"
        input_dir = metadata_dir if metadata_dir.exists() else script_dir

    output_dir = Path(args.output_dir) if args.output_dir else script_dir.parent / "csv_ncbi"
    output_dir.mkdir(exist_ok=True)

    # Set up paths for input files with fallback logic
    assembly_file = input_dir / "00assembly_summary_genbank.txt"

    # If assembly file not found in input_dir, try metadata folder
    if not assembly_file.exists() and input_dir != script_dir.parent / "metadata":
        metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
        if metadata_assembly.exists():
            assembly_file = metadata_assembly
            print(f"📁 Using assembly file from metadata folder: {assembly_file}")

    mapping_file = Path(args.mapping_file) if args.mapping_file else script_dir.parent / "taxonomic_mapping" / "taxid_to_family.csv"
    output_file = output_dir / "ncbi_family_species_counts.csv"
    accession_output_file = output_dir / "ncbi_family_species_subset_with_accessions.csv"

    return assembly_file, mapping_file, output_file, accession_output_file

def load_assembly_file(file_path: Path) -> pd.DataFrame:
    """Load NCBI assembly file with memory optimization."""
    try:
        # Load with specific dtypes to reduce memory usage
        dtype_dict = {
            'taxid': 'str',  # Use string to avoid int overflow issues
            'assembly_accession': 'str',
            'organism_name': 'str'
        }
        return pd.read_csv(file_path, sep="\t", skiprows=1, low_memory=False, dtype=dtype_dict)
    except FileNotFoundError:
        print(f"❌ Error: File {file_path} not found")
        raise
    except pd.errors.EmptyDataError:
        print(f"❌ Error: File {file_path} is empty")
        raise
    except Exception as e:
        print(f"❌ Error loading {file_path}: {str(e)}")
        raise

def load_taxid_mapping(file_path: Path) -> pd.DataFrame:
    """Load taxid mapping file with memory optimization."""
    try:
        # Load with specific dtypes to reduce memory usage
        dtype_dict = {
            'taxid': 'str',  # Use string to match assembly file
            'family': 'str',
            'domain': 'str'
        }
        return pd.read_csv(file_path, dtype=dtype_dict)
    except FileNotFoundError:
        print(f"❌ Error: Mapping file {file_path} not found")
        raise
    except pd.errors.EmptyDataError:
        print(f"❌ Error: Mapping file {file_path} is empty")
        raise
    except Exception as e:
        print(f"❌ Error loading mapping file {file_path}: {str(e)}")
        raise

def main() -> int:
    """Main function to run the script."""
    args = parse_arguments()

    try:
        # Set up logging
        log_file = setup_logging("family_ncbi_parser_species_subset")
        logging.info("Starting family NCBI parser (species subset)")

        # Set up paths
        assembly_file, mapping_file, output_file, accession_output_file = setup_paths(args)
        print(f"📂 Input assembly file: {assembly_file}")
        print(f"📂 Input mapping file: {mapping_file}")
        print(f"📂 Output file: {output_file}")
        print(f"📂 Log file: {log_file}")

        logging.info(f"Input assembly file: {assembly_file}")
        logging.info(f"Input mapping file: {mapping_file}")
        logging.info(f"Output file: {output_file}")

        # Load assembly file
        print("Loading assembly file...")
        logging.info("Loading assembly file...")
        print(f"💾 Memory usage before loading: {get_memory_usage():.2f} GB")
        logging.info(f"Memory usage before loading: {get_memory_usage():.2f} GB")
        with tqdm(desc="Loading assembly data", unit="rows") as pbar:
            df = load_assembly_file(assembly_file)
            df = df.rename(columns={"#assembly_accession": "assembly_accession"})
            pbar.update(len(df))
        print(f"✅ Loaded {len(df)} assembly entries")
        logging.info(f"Loaded {len(df)} assembly entries")
        print(f"💾 Memory usage after loading assembly: {get_memory_usage():.2f} GB")
        logging.info(f"Memory usage after loading assembly: {get_memory_usage():.2f} GB")

        # Load taxid mapping file
        print("Loading taxid mapping file...")
        with tqdm(desc="Loading taxonomy mapping", unit="rows") as pbar:
            mapping_df = load_taxid_mapping(mapping_file)
            pbar.update(len(mapping_df))
        print(f"✅ Loaded {len(mapping_df)} mapping entries")
        print(f"💾 Memory usage after loading mapping: {get_memory_usage():.2f} GB")
        logging.info(f"Memory usage after loading mapping: {get_memory_usage():.2f} GB")

        # Debug: Check data types and duplicates before merge
        print("🔍 Checking data before merge...")
        print(f"Assembly df taxid dtype: {df['taxid'].dtype}")
        print(f"Mapping df taxid dtype: {mapping_df['taxid'].dtype}")
        print(f"Assembly df unique taxids: {df['taxid'].nunique()}")
        print(f"Mapping df unique taxids: {mapping_df['taxid'].nunique()}")
        print(f"Assembly df duplicated taxids: {df['taxid'].duplicated().sum()}")
        print(f"Mapping df duplicated taxids: {mapping_df['taxid'].duplicated().sum()}")

        # Ensure taxid columns have the same data type
        df['taxid'] = df['taxid'].astype(str)
        mapping_df['taxid'] = mapping_df['taxid'].astype(str)

        # Remove any duplicates in mapping_df to prevent Cartesian product
        if mapping_df['taxid'].duplicated().any():
            print("⚠️  Found duplicated taxids in mapping file, removing duplicates...")
            mapping_df = mapping_df.drop_duplicates(subset=['taxid'], keep='first')
            print(f"✅ Cleaned mapping data: {len(mapping_df)} unique entries")

        # Merge assembly data with mapping data
        print("Merging assembly data with taxonomic mapping...")
        merged_df = pd.merge(df, mapping_df, on='taxid', how='inner')
        print(f"✅ Merged data: {len(merged_df)} entries")
        print(f"💾 Memory usage after merge: {get_memory_usage():.2f} GB")
        logging.info(f"Memory usage after merge: {get_memory_usage():.2f} GB")

        # Debug: Check what columns are available
        print(f"📋 Available columns: {list(merged_df.columns)}")

        # Check if we have species-level columns for deduplication
        species_column = None
        if 'species_taxid' in merged_df.columns:
            species_column = 'species_taxid'
        elif 'organism_name' in merged_df.columns:
            species_column = 'organism_name'
        elif 'species' in merged_df.columns:
            species_column = 'species'
        elif 'scientific_name' in merged_df.columns:
            species_column = 'scientific_name'

        # Load classified accession data to get genome_source information
        print("Loading classified accession data for isolate/uncultured information...")
        output_dir_path = Path(args.output_dir) if args.output_dir else Path(__file__).resolve().parent.parent / "csv_ncbi"
        classified_file = output_dir_path / "ncbi_family_with_accessions_classified.csv"

        if classified_file.exists():
            print(f"📁 Loading classified data from: {classified_file}")
            classified_df = pd.read_csv(classified_file)

            # Merge with classified data to get genome_source information
            print("Merging with classified data to get genome_source information...")
            merged_df = pd.merge(merged_df, classified_df[['accession_clean', 'genome_source']],
                               left_on='assembly_accession', right_on='accession_clean', how='left')

            # Fill missing genome_source with 'isolate' (default assumption)
            merged_df['genome_source'] = merged_df['genome_source'].fillna('isolate')

            print(f"✅ Added genome_source information to {len(merged_df)} genomes")

            # Show genome source distribution
            source_counts = merged_df['genome_source'].value_counts()
            print("📊 Genome source distribution:")
            for source, count in source_counts.items():
                print(f"   - {source}: {count:,}")
        else:
            print(f"⚠️  Classified file not found: {classified_file}")
            print("⚠️  Will use standard keep='first' strategy without isolate preference")

        # Step 1: Subset to keep one genome per species using isolate-preferential strategy
        print(f"Subsetting to keep one genome per species using '{species_column}' column...")
        species_subset_df = create_isolate_preferential_species_subset(merged_df, species_column)
        print(f"✅ Reduced from {len(merged_df)} genomes to {len(species_subset_df)} unique species")

        # Step 2: Group by family and domain, and count how many unique species
        print("Grouping by family and domain...")
        if species_column == 'species_taxid':
            # For species_taxid, count unique taxids per family
            counts = species_subset_df.groupby(["family", "domain"]).agg({
                species_column: "nunique"  # Count unique species taxids
            }).reset_index()
        else:
            # For other species columns, count occurrences
            counts = species_subset_df.groupby(["family", "domain"]).agg({
                species_column: "count"
            }).reset_index()

        # Step 3: Rename column for clarity
        counts = counts.rename(columns={species_column: "family_species_count"})

        # Step 4: Sort descending by count
        counts = counts.sort_values("family_species_count", ascending=False)

        # Print preview
        print(f"\n📈 Top 10 families by species count:")
        for _, row in counts.head(10).iterrows():
            print(f"  {row['family']} ({row['domain']}): {row['family_species_count']} species")

        # Get unique family names for taxonkit processing
        unique_families = counts['family'].unique().tolist()
        print(f"\n🔍 Found {len(unique_families)} unique families")

        # Get taxids for family names using taxonkit name2taxid
        print("Getting taxids for family names...")
        with tqdm(desc="Processing family names", total=len(unique_families), unit="families") as pbar:
            family_to_taxid = get_taxids_from_names(unique_families)
            pbar.update(len(unique_families))

        # Add taxid column to counts
        print("Mapping taxids to family counts...")
        with tqdm(desc="Adding taxids", total=len(counts), unit="rows") as pbar:
            counts['taxid'] = counts['family'].map(family_to_taxid)
            pbar.update(len(counts))

        # Get lineage information for the taxids
        valid_taxids = [taxid for taxid in family_to_taxid.values() if taxid]
        print(f"Getting lineage information for {len(valid_taxids)} valid taxids...")
        with tqdm(desc="Processing lineages", total=len(valid_taxids), unit="taxids") as pbar:
            lineage_data = get_lineages_from_taxids(valid_taxids)
            pbar.update(len(valid_taxids))

        # Add lineage columns
        print("Adding lineage information to results...")
        lineages = []
        lineage_ranks = []
        lineage_taxids = []

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
        print(f"✅ Saved {len(counts)} family-domain combinations")

        # Save species subset accession file with isolate-preferential selection
        print(f"💾 Saving species subset accessions to {accession_output_file}...")

        # Select essential columns for the accession file
        accession_columns = ['assembly_accession', 'family', 'domain', 'taxid']
        if 'genome_source' in species_subset_df.columns:
            accession_columns.append('genome_source')
        if 'excluded_from_refseq' in species_subset_df.columns:
            accession_columns.append('excluded_from_refseq')

        # Filter to only include columns that exist
        available_columns = [col for col in accession_columns if col in species_subset_df.columns]
        accession_subset = species_subset_df[available_columns].copy()

        # Save the accession file
        accession_subset.to_csv(accession_output_file, index=False)
        print(f"✅ Saved {len(accession_subset)} species subset accessions")

        if 'genome_source' in accession_subset.columns:
            # Show final isolate/uncultured distribution in saved file
            final_source_counts = accession_subset['genome_source'].value_counts()
            print("📊 Final species subset genome source distribution:")
            for source, count in final_source_counts.items():
                percentage = (count / len(accession_subset)) * 100
                print(f"   - {source}: {count:,} ({percentage:.1f}%)")

        print("✅ Species subset version complete - generated both counts and accession files")

        logging.info("Family NCBI parser (species subset) completed successfully")
        return 0

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        logging.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())


