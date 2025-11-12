#!/usr/bin/env python3
"""
EukProt Source Metadata Validator and Processor

This script processes the source Eukprot_included_datasets.txt file directly,
validates organism names using Previous_Names column and known corrections,
and generates taxid mappings with comprehensive logging for manual review.

Key Features:
- Processes source metadata directly (no intermediate CSV needed)
- Uses Previous_Names column for validation
- Applies known spelling corrections
- Creates detailed logs for fishy matches
- Generates taxid mappings with match type tracking

Usage:
    python validate_eukprot_source.py [input_file] [output_file]

    input_file  - Source metadata file (default: Eukprot_included_datasets.txt)
    output_file - Output CSV file (default: eukprot_validated_with_taxids.csv)
"""

import pandas as pd
import os
import subprocess
import logging
import sys
import tempfile
from pathlib import Path
import time
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("eukprot_validation.log"),
        logging.StreamHandler()
    ]
)

# NCBI taxdump directory
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")

# Known name corrections for EukProt misspellings
NAME_CORRECTIONS = {
    "Luapelamoeba hula": "Luapeleamoeba hula",
    "Sappina": "Sappinia pedata",
    "Sapocribrum chincoteaguense": "Sapocribrum chincoteaguense",  # underscore fix
    # Add more corrections as needed
}

def process_source_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process the source EukProt metadata to create a working dataframe.

    Args:
        df: Raw metadata dataframe

    Returns:
        Processed dataframe with Name_to_use, Previous_Names, etc.
    """
    logging.info("🔧 Processing source metadata...")

    # Create working dataframe
    result_df = pd.DataFrame()

    # Extract and clean Name_to_Use (replace underscores with spaces)
    result_df['Name_to_use'] = df['Name_to_Use'].str.replace('_', ' ')

    # Extract Previous_Names for validation
    result_df['Previous_Names'] = df['Previous_Names'].fillna('')

    # Extract taxonomic components for reference
    result_df['Genus_UniEuk'] = df['Genus_UniEuk'].fillna('')
    result_df['Epithet_UniEuk'] = df['Epithet_UniEuk'].fillna('')
    result_df['Taxonomy_UniEuk'] = df['Taxonomy_UniEuk'].fillna('')  # For genus validation
    result_df['EukProt_ID'] = df['EukProt_ID']

    # Initialize taxid and match_type columns
    result_df['taxid'] = 'FAILED'
    result_df['match_type'] = 'none'

    logging.info(f"✅ Processed {len(result_df)} entries from source metadata")
    return result_df

def get_taxid_for_name(name: str, env: dict) -> str:
    """
    Get taxid for a single name using taxonkit.

    Args:
        name: Species name to look up
        env: Environment variables for subprocess

    Returns:
        Taxid string or empty string if not found
    """
    if not name or name.strip() == "":
        return ""

    # Create temporary file with the name
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_file.write(f"{name.strip()}\n")
        temp_file_path = temp_file.name

    try:
        # Create script to run taxonkit name2taxid
        script_path = temp_file_path + "_script.sh"
        with open(script_path, 'w') as script_file:
            script_file.write(f"""#!/bin/bash
cat {temp_file_path} | taxonkit name2taxid
""")

        os.chmod(script_path, 0o755)

        # Run the script
        result = subprocess.run(
            [script_path],
            capture_output=True,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 2 and parts[1].strip() != "":
                    return parts[1].strip()

        return ""

    except Exception as e:
        logging.warning(f"Error getting taxid for '{name}': {e}")
        return ""

    finally:
        # Clean up temporary files
        try:
            os.unlink(temp_file_path)
            if os.path.exists(script_path):
                os.unlink(script_path)
        except Exception as e:
            logging.warning(f"Error cleaning up temporary files: {e}")

def clean_name_for_taxonkit(name: str) -> str:
    """Clean organism name for taxonkit compatibility."""
    if not name:
        return ""

    # Remove common problematic patterns
    cleaned = name.strip()

    # Remove strain information (anything after 'strain')
    if 'strain' in cleaned.lower():
        cleaned = cleaned.split('strain')[0].strip()

    # Remove isolate information
    if 'isolate' in cleaned.lower():
        cleaned = cleaned.split('isolate')[0].strip()

    # Clean up extra spaces
    cleaned = ' '.join(cleaned.split())

    return cleaned

def extract_genus_from_name(name: str) -> str:
    """Extract genus from a species name."""
    if not name:
        return ""

    # Clean the name first
    cleaned_name = clean_name_for_taxonkit(name)
    parts = cleaned_name.split()

    if len(parts) >= 1:
        first_word = parts[0].strip()
        if first_word and first_word[0].isupper():
            return first_word

    return ""

def extract_genus_from_taxonomy(taxonomy: str) -> str:
    """
    Extract genus from Taxonomy_UniEuk lineage.

    Args:
        taxonomy: Taxonomy lineage string

    Returns:
        Genus name if found, empty string otherwise
    """
    if not taxonomy or pd.isna(taxonomy):
        return ""

    # Split taxonomy by semicolons
    parts = [part.strip() for part in str(taxonomy).split(';') if part.strip()]

    # Look for genus (usually second to last part before species)
    for part in reversed(parts):
        # Skip species names (containing spaces) and look for genus
        if ' ' not in part and part != 'Eukaryota' and len(part) > 2:
            # Check if it looks like a genus (starts with capital letter)
            if part[0].isupper():
                return part

    return ""

def validate_and_get_taxids(df: pd.DataFrame, env: dict) -> pd.DataFrame:
    """
    Validate organism names and get taxids using multiple strategies.

    Args:
        df: Dataframe with organism names
        env: Environment variables for subprocess

    Returns:
        Dataframe with validated taxids and match types
    """
    logging.info("🔍 Starting validation and taxid lookup...")

    df_result = df.copy()
    validation_log = []
    fishy_matches = []

    for idx, row in tqdm(df_result.iterrows(), total=len(df_result), desc="Validating entries"):
        original_name = row['Name_to_use']
        previous_names = row['Previous_Names']
        taxonomy = row['Taxonomy_UniEuk']

        taxid = ""
        match_type = "none"
        correction_applied = ""
        genus_used = ""
        species_to_append = ""

        # Clean the original name for taxonkit
        cleaned_name = clean_name_for_taxonkit(original_name)

        # Strategy 1: Try known corrections first
        corrected_name = NAME_CORRECTIONS.get(original_name, original_name)
        if corrected_name != original_name:
            correction_applied = f"Known correction: {original_name} → {corrected_name}"
            cleaned_corrected = clean_name_for_taxonkit(corrected_name)
            taxid = get_taxid_for_name(cleaned_corrected, env)
            if taxid:
                match_type = "corrected_full"
                species_to_append = corrected_name  # Use corrected name for lineage

        # Strategy 2: Try cleaned original name
        if not taxid:
            taxid = get_taxid_for_name(cleaned_name, env)
            if taxid:
                match_type = "full"
                species_to_append = original_name  # Use original name for lineage

        # Strategy 3: Try genus from Name_to_use
        if not taxid:
            genus = extract_genus_from_name(original_name)
            if genus:
                taxid = get_taxid_for_name(genus, env)
                if taxid:
                    match_type = "genus"
                    genus_used = genus
                    species_to_append = original_name  # Append full species name to lineage

        # Strategy 4: Try genus from Taxonomy_UniEuk (for misspelled genus validation)
        if not taxid and taxonomy:
            taxonomy_genus = extract_genus_from_taxonomy(taxonomy)
            if taxonomy_genus:
                # Check if this genus is different from Name_to_use genus (indicating misspelling)
                name_genus = extract_genus_from_name(original_name)
                if taxonomy_genus != name_genus:
                    correction_applied += f" | Genus from taxonomy: {name_genus} → {taxonomy_genus}"

                taxid = get_taxid_for_name(taxonomy_genus, env)
                if taxid:
                    match_type = "taxonomy_genus"
                    genus_used = taxonomy_genus
                    species_to_append = original_name  # Append original species name to lineage

        # Strategy 5: Try previous names (full species)
        if not taxid and previous_names:
            prev_names_list = [name.strip() for name in previous_names.split(',') if name.strip()]
            for prev_name in prev_names_list:
                cleaned_prev = clean_name_for_taxonkit(prev_name)
                taxid = get_taxid_for_name(cleaned_prev, env)
                if taxid:
                    match_type = "previous_name"
                    species_to_append = prev_name  # Use previous name for lineage
                    break

        # Strategy 6: Try genus from previous names
        if not taxid and previous_names:
            prev_names_list = [name.strip() for name in previous_names.split(',') if name.strip()]
            for prev_name in prev_names_list:
                prev_genus = extract_genus_from_name(prev_name)
                if prev_genus:
                    taxid = get_taxid_for_name(prev_genus, env)
                    if taxid:
                        match_type = "previous_name_genus"
                        genus_used = prev_genus
                        species_to_append = prev_name  # Append previous species name to lineage
                        break

        # Update the dataframe
        df_result.at[idx, 'taxid'] = taxid if taxid else 'FAILED'
        df_result.at[idx, 'match_type'] = match_type

        # Add species name for lineage appending (will be used later in lineage generation)
        if species_to_append:
            df_result.at[idx, 'species_for_lineage'] = species_to_append
        else:
            df_result.at[idx, 'species_for_lineage'] = original_name

        # Check for fishy matches
        is_fishy = False
        fishy_reasons = []

        if not taxid:
            is_fishy = True
            fishy_reasons.append("No taxid found with any strategy")

        if match_type in ['previous_name_genus', 'genus']:
            is_fishy = True
            fishy_reasons.append(f"Fallback match type: {match_type}")

        if correction_applied:
            fishy_reasons.append(correction_applied)

        # Log the result
        validation_log.append({
            'index': idx,
            'original_name': original_name,
            'previous_names': previous_names,
            'taxonomy': taxonomy,
            'final_taxid': taxid if taxid else 'FAILED',
            'match_type': match_type,
            'genus_used': genus_used,
            'species_for_lineage': species_to_append if species_to_append else original_name,
            'correction_applied': correction_applied,
            'is_fishy': is_fishy
        })

        if is_fishy:
            fishy_matches.append({
                'index': idx,
                'original_name': original_name,
                'previous_names': previous_names,
                'final_taxid': taxid if taxid else 'FAILED',
                'match_type': match_type,
                'reasons': fishy_reasons
            })

    # Write logs
    write_validation_logs(validation_log, fishy_matches)

    return df_result

def write_validation_logs(validation_log: list, fishy_matches: list):
    """Write validation and fishy match logs."""

    # Write comprehensive validation log
    with open("eukprot_validation_results.log", "w") as f:
        f.write("EukProt Validation Results\n")
        f.write("=" * 50 + "\n")

        success_count = sum(1 for entry in validation_log if entry['final_taxid'] != 'FAILED')
        total_count = len(validation_log)

        f.write(f"Total entries: {total_count}\n")
        f.write(f"Successful matches: {success_count} ({success_count/total_count*100:.1f}%)\n")
        f.write(f"Failed matches: {total_count - success_count} ({(total_count - success_count)/total_count*100:.1f}%)\n\n")

        for entry in validation_log:
            f.write(f"Row {entry['index']}: {entry['match_type'].upper()}\n")
            f.write(f"  Original: '{entry['original_name']}'\n")
            f.write(f"  Previous Names: '{entry['previous_names']}'\n")
            f.write(f"  Final Taxid: {entry['final_taxid']}\n")
            if entry['genus_used']:
                f.write(f"  Genus Used: {entry['genus_used']}\n")
            f.write(f"  Species for Lineage: {entry['species_for_lineage']}\n")
            if entry['correction_applied']:
                f.write(f"  Correction: {entry['correction_applied']}\n")
            if entry['is_fishy']:
                f.write(f"  ⚠️  FISHY - Needs Review\n")
            f.write("-" * 30 + "\n")

    # Write fishy matches log
    if fishy_matches:
        with open("fishy_matches_manual_review.log", "w") as f:
            f.write("🚨 FISHY MATCHES REQUIRING MANUAL REVIEW 🚨\n")
            f.write("=" * 60 + "\n")
            f.write(f"Found {len(fishy_matches)} entries that need attention:\n\n")

            for entry in fishy_matches:
                f.write(f"Row {entry['index']}: {entry['match_type'].upper()}\n")
                f.write(f"  Original Name: '{entry['original_name']}'\n")
                f.write(f"  Previous Names: '{entry['previous_names']}'\n")
                f.write(f"  Final Taxid: {entry['final_taxid']}\n")
                f.write(f"  Issues Found:\n")
                for reason in entry['reasons']:
                    f.write(f"    - {reason}\n")
                f.write("=" * 60 + "\n")

        logging.warning(f"⚠️  Found {len(fishy_matches)} fishy matches - see fishy_matches_manual_review.log")

def main():
    """Main function."""
    start_time = time.time()

    # Parse command line arguments
    if len(sys.argv) < 2:
        input_file = "Eukprot_included_datasets.txt"
        output_file = "eukprot_validated_with_taxids.csv"
    elif len(sys.argv) == 2:
        input_file = sys.argv[1]
        output_file = "eukprot_validated_with_taxids.csv"
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

    logging.info(f"Processing: {input_file} → {output_file}")

    # Load source metadata
    try:
        df = pd.read_csv(input_file, sep='\t')
        logging.info(f"✅ Loaded {len(df)} entries from source metadata")
    except Exception as e:
        logging.error(f"❌ Error loading {input_file}: {e}")
        sys.exit(1)

    # Process metadata
    df_processed = process_source_metadata(df)

    # Set up environment for taxonkit
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    # Validate and get taxids
    df_final = validate_and_get_taxids(df_processed, env)

    # Save results
    df_final.to_csv(output_file, index=False)

    # Calculate statistics
    elapsed_time = time.time() - start_time
    success_count = len(df_final[df_final['taxid'] != 'FAILED'])
    total_count = len(df_final)

    logging.info(f"✅ Processing complete in {elapsed_time:.2f} seconds")
    logging.info(f"✅ Found taxids for {success_count}/{total_count} entries ({success_count/total_count*100:.1f}%)")
    logging.info(f"✅ Results saved to {output_file}")

if __name__ == "__main__":
    main()
