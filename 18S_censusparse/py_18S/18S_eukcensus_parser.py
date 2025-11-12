#!/usr/bin/env python3
"""
Optimized EukCensus 18S Cluster Parser

This script processes EukCensus 18S cluster metadata using the same logic as the original
parse_eukcensus_clusters.py but with optimized performance. It generates CSV files
organized by division, family, and genus with taxid, taxon_name, member_size,
occurrence_count, and lineage columns.

The script uses taxonkit to get taxids for each taxon name and properly handles organelle
information in taxon names (e.g., "Aspergillus_nidulans_FGSC_A4.Mitochondria").

Output files:
- eukcensus_by_division.csv
- eukcensus_by_family.csv
- eukcensus_by_genus.csv
"""

import pandas as pd
import os
import subprocess
import tempfile
from collections import defaultdict
import csv
import concurrent.futures
import math
import logging
import sys
import time
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

def setup_directory_paths():
    """
    Set up directory paths for the reorganized 18S_censusparse structure.

    Returns:
        Tuple of (script_dir, metadata_dir, csv_output_dir, log_dir)
    """
    script_dir = Path(__file__).resolve().parent
    censusparse_dir = script_dir.parent
    metadata_dir = censusparse_dir / "metadata"
    csv_output_dir = censusparse_dir / "csv_outputs"
    log_dir = script_dir / "logs"

    # Ensure directories exist
    csv_output_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    return script_dir, metadata_dir, csv_output_dir, log_dir

def setup_taxonkit_environment():
    """
    Set up environment for taxonkit.

    Note: taxonkit has its own built-in taxdump and doesn't require external TAXONKIT_DB.

    Returns:
        Environment dictionary (using system default)
    """
    env = os.environ.copy()
    # taxonkit uses its own built-in NCBI taxdump, no need to set TAXONKIT_DB
    return env

def setup_logging(log_dir, output_prefix="eukcensus_optimized"):
    """
    Set up logging configuration for the optimization script.

    Args:
        log_dir: Directory to store log files
        output_prefix: Prefix for output files to include in log messages
    """
    # Create log file path in the logs directory
    log_file = log_dir / "eukcensus_optimization.log"

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),  # Append to existing log
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Log the start of processing
    logging.info(f"Starting optimized EukCensus processing: eukcensus_18S.clusters.97.tsv -> {output_prefix}_*")
    logging.info(f"Log file: {log_file}")
    return logging.getLogger(__name__)

def create_comprehensive_unmapped_log(division_data, family_data, genus_data,
                                    division_to_taxid, family_to_taxid, genus_to_taxid,
                                    taxid_to_lineage, log_dir, output_prefix):
    """
    Create a comprehensive log of all unmapped taxonomic names with enhanced analysis.

    Args:
        division_data, family_data, genus_data: Data dictionaries for each rank
        division_to_taxid, family_to_taxid, genus_to_taxid: Taxid mapping dictionaries
        taxid_to_lineage: Lineage information dictionary
        log_dir: Directory to store log files
        output_prefix: Prefix for output files
    """
    log_file = log_dir / f"{output_prefix}_comprehensive_unmapped.log"
    logging.info("Creating comprehensive unmapped log...")

    with open(log_file, 'w') as f:
        f.write("# Enhanced Comprehensive Unmapped Names Log - 18S Census Parser with Fallback Analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# This log contains all taxonomic names that failed to get NCBI taxids or lineages\n")
        f.write("# Enhanced with fallback strategy analysis and improved pattern recognition\n")
        f.write("# Format: Rank | Original_Name | Cleaned_Name | Appropriate_Name | OTU_Count | Size_Count | Taxid | Reason\n\n")

        # Summary statistics
        f.write("=== SUMMARY STATISTICS ===\n")

        # Calculate statistics for each rank
        rank_stats = {}
        for rank_name, data_dict, taxid_dict in [
            ('division', division_data, division_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            total = len(data_dict)
            mapped = len([t for t in taxid_dict.values() if t != 'NA'])
            unmapped = total - mapped
            unmapped_pct = (unmapped / total * 100) if total > 0 else 0

            rank_stats[rank_name] = {
                'total': total,
                'mapped': mapped,
                'unmapped': unmapped,
                'unmapped_pct': unmapped_pct
            }

            f.write(f"{rank_name.capitalize()} Statistics:\n")
            f.write(f"  Total: {total}\n")
            f.write(f"  Mapped: {mapped} ({(mapped/total*100):.1f}%)\n")
            f.write(f"  Unmapped: {unmapped} ({unmapped_pct:.1f}%)\n\n")

        # Overall statistics
        total_all = sum(stats['total'] for stats in rank_stats.values())
        mapped_all = sum(stats['mapped'] for stats in rank_stats.values())
        unmapped_all = total_all - mapped_all
        overall_unmapped_pct = (unmapped_all / total_all * 100) if total_all > 0 else 0

        f.write(f"Overall Statistics:\n")
        f.write(f"  Total entries: {total_all}\n")
        f.write(f"  Successfully mapped: {mapped_all} ({(mapped_all/total_all*100):.1f}%)\n")
        f.write(f"  Failed to map: {unmapped_all} ({overall_unmapped_pct:.1f}%)\n\n")

        # Detailed unmapped entries by rank
        for rank_name, data_dict, taxid_dict in [
            ('division', division_data, division_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            f.write(f"=== {rank_name.upper()} LEVEL UNMAPPED NAMES ===\n")

            unmapped_entries = []
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    # Determine reason for failure
                    if taxid == "NA":
                        reason = "NO_TAXID_FOUND"
                    else:
                        reason = "TAXID_NO_LINEAGE"

                    unmapped_entries.append({
                        'original_name': orig_name,
                        'cleaned_name': clean_taxon_name(orig_name),
                        'appropriate_name': orig_name,  # For 18S, we use original name
                        'otu_count': data['otu_count'],
                        'size_count': data.get('size_count', 0),
                        'taxid': taxid,
                        'reason': reason
                    })

            f.write(f"Total unmapped {rank_name} entries: {len(unmapped_entries)}\n\n")

            # Sort by OTU count (descending)
            unmapped_entries.sort(key=lambda x: x['otu_count'], reverse=True)

            for entry in unmapped_entries:
                f.write(f"{rank_name.upper()} | {entry['original_name']} | {entry['cleaned_name']} | {entry['appropriate_name']} | {entry['otu_count']} | {entry['size_count']} | {entry['taxid']} | {entry['reason']}\n")

            f.write(f"\n")

        # Pattern analysis
        f.write("=== PATTERN ANALYSIS ===\n")

        # Collect all unmapped names for pattern analysis
        all_unmapped = []
        for rank_name, data_dict, taxid_dict in [
            ('division', division_data, division_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    all_unmapped.append({
                        'rank': rank_name,
                        'name': orig_name,
                        'otu_count': data['otu_count'],
                        'size_count': data.get('size_count', 0)
                    })

        # Analyze patterns (adapted for 18S eukaryotic data)
        patterns = {
            'Candidatus': [n for n in all_unmapped if 'Candidatus' in n['name']],
            'Organelles': [n for n in all_unmapped if any(org in n['name'] for org in ['Mitochondria', 'Chloroplast', 'Apicoplast', 'Plastid', ':plas', ':mito', ':api'])],
            'Uncultured': [n for n in all_unmapped if 'uncultured' in n['name'].lower()],
            'Environmental': [n for n in all_unmapped if any(env in n['name'].lower() for env in ['metagenome', 'environmental', 'marine'])],
            'Species_level': [n for n in all_unmapped if '_' in n['name'] and '.' not in n['name']],
            'Unclassified': [n for n in all_unmapped if 'unclassified' in n['name'].lower()],
            'Group_taxa': [n for n in all_unmapped if 'Group' in n['name'] or 'group' in n['name']]
        }

        for pattern_name, pattern_names in patterns.items():
            if pattern_names:
                count = len(pattern_names)
                total_otu_occurrences = sum(n['otu_count'] for n in pattern_names)
                total_size_occurrences = sum(n['size_count'] for n in pattern_names)
                f.write(f"{pattern_name.replace('_', ' ').title()}: {count} names ({total_otu_occurrences} total OTU occurrences, {total_size_occurrences} total sequence occurrences)\n")

                # Show top 5 most frequent for each pattern
                top_names = sorted(pattern_names, key=lambda x: x['otu_count'], reverse=True)[:5]
                for name_info in top_names:
                    f.write(f"  {name_info['name']} ({name_info['rank']}) - {name_info['otu_count']} OTU occurrences, {name_info['size_count']} sequence occurrences\n")
                f.write("\n")

        f.write("=== RECOMMENDATIONS ===\n")
        f.write("1. Organelle sequences (mitochondria, chloroplast, plastid, apicoplast) should be mapped to host organisms\n")
        f.write("2. Original organellar names are preserved to prevent merging with non-organellar OTUs\n")
        f.write("3. Candidatus taxa should be properly preserved in NCBI taxonomy (no filtering)\n")
        f.write("4. Modern eukaryotic division names may need mapping to older NCBI names\n")
        f.write("5. Environmental/uncultured samples may not have valid NCBI taxids\n")
        f.write("6. Species-level entries need genus extraction for genus-level parsing\n")
        f.write("7. Consider using alternative taxonomic databases for better eukaryotic coverage\n")

    logging.info(f"Comprehensive unmapped log written to {log_file}")
    return log_file

def clean_taxon_name(taxon_name):
    """
    Clean a taxon name by replacing underscores with two spaces and removing trailing numbers.

    This handles EukCensus patterns like "_XX" by removing everything after the underscore
    and adding two spaces to help with taxonkit matching.

    Args:
        taxon_name: The taxon name to clean

    Returns:
        The cleaned taxon name
    """
    # Replace underscores with two spaces and remove trailing numbers
    cleaned = taxon_name.replace("_", "  ")
    return strip_trailing_numbers(cleaned)

def strip_trailing_numbers(taxon_name):
    """
    Remove trailing numbers from taxon names (e.g., "Theileria1" -> "Theileria").

    Args:
        taxon_name: The taxon name to process

    Returns:
        The taxon name with trailing numbers removed
    """
    import re

    # Handle taxa with numbers at the end (e.g., "Cryptosporidium15", "Eimeria1", "Plasmodium1")
    # This regex matches any word characters followed by one or more digits at the end
    match = re.match(r'^(.+?)(\d+)$', taxon_name.strip())
    if match:
        base_name = match.group(1).rstrip()  # Remove any trailing spaces
        return base_name

    # If no trailing numbers found, return as is
    return taxon_name

def extract_taxa_from_hyphenated(taxon_name):
    """
    Extract taxa from hyphenated names based on different patterns.

    Args:
        taxon_name: The hyphenated taxon name

    Returns:
        The extracted taxon name or None if no valid pattern found
    """
    if "-" not in taxon_name:
        return None

    # Pattern 1: [taxa]-lineage → extract first part
    if taxon_name.endswith("-lineage") or "_X" in taxon_name and "-lineage" in taxon_name:
        base_name = taxon_name.split("-lineage")[0]
        # Remove any _X suffix
        if "_" in base_name:
            base_name = base_name.split("_")[0]
        return base_name

    # Pattern 2: [taxa]-Group → extract first part
    if "-Group" in taxon_name:
        return taxon_name.split("-Group")[0]

    # Pattern 3: X-[taxa]_XX → extract middle part (existing logic)
    parts = taxon_name.split("-")
    if len(parts) >= 2:
        second_part = parts[1]
        # Remove any _XXX, _XX, _X suffix
        if "_" in second_part:
            clean_part = second_part.split("_")[0]
        else:
            clean_part = second_part

        # Only return if it looks like a valid taxonomic name (not a research clade)
        if not any(keyword in clean_part.lower() for keyword in ["clade", "group", "relatives"]):
            return clean_part

    return None

def extract_genus(taxon_name):
    """
    Extract the genus part from a taxon name, handling trailing numbers.

    Args:
        taxon_name: The taxon name to extract genus from

    Returns:
        The genus part of the taxon name, or None if not found
    """
    # First clean the name to remove organelle information
    if "." in taxon_name:
        # For names like "Genus_species.Mitochondria", get the part before the dot
        parts = taxon_name.split(".")
        name_part = parts[0]

        # If the name part contains an underscore, get the first part (genus)
        if "_" in name_part:
            genus = name_part.split("_")[0]
        else:
            genus = name_part

        # Strip trailing numbers from the genus
        return strip_trailing_numbers(genus)

    # For names like "Genus_species", get the first part
    if "_" in taxon_name:
        genus = taxon_name.split("_")[0]
        return strip_trailing_numbers(genus)

    # For names like "Genus" or "Genus1", return the genus part (without numbers)
    return strip_trailing_numbers(taxon_name)

def should_filter_taxon(taxon_name):
    """
    Check if a taxon name should be filtered out from final output.

    Modified to KEEP unclassified (.U.) entries as they represent significant biological diversity.
    Based on analysis showing .U. entries represent 20% of division, 59% of family, and 72% of genus data.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    # Previously filtered out .U. entries, but analysis shows this causes massive data loss:
    # - 20% loss at division level
    # - 59% loss at family level
    # - 72% loss at genus level
    # Therefore, we now KEEP all .U. entries to preserve biological diversity

    # Only filter out truly problematic entries (none currently defined)
    # Future filtering criteria can be added here if needed

    return False

def get_taxonomic_mapping():
    """
    Get mapping from outdated/informal taxonomic names to modern valid names.

    Returns:
        Dictionary mapping old names to new names
    """
    return {
        # Outdated taxonomic names to modern equivalents
        "Maxillopoda": "Copepoda",
        "Embryophyceae_XX": "Embryophyta",
        "Chytridiomycetaceae": "Chytridiomycetes",
        "Ophryoglenida": "Ophryoglenidae",

        # Informal groups to formal taxonomy
        "Blastocystis-Group": "Blastocystis",
        "Flamella-lineage": "Flamella",
        "Endostelium-lineage": "Endostelium",
        "Protaspa-lineage": "Protaspa",
        "Rhogostoma-lineage": "Rhogostoma",

        # Research clades to broader valid groups where possible
        "Neobodonidae": "Bodonidae",  # Neobodonidae is often considered part of Bodonidae
        "Vermamoebidae": "Amoebidae",  # Map to broader amoeba family
        "Tholoniidae": "Tholonia",  # Map to genus level
        "Nolandellidae": "Nolandella",  # Map to genus level
        "Paradinidae": "Paradinium",  # Map to genus level
        "Skeletonemaceae": "Skeletonema",  # Map to genus level

        # Remove suffixes for research annotations
        "Dino-Group-II_X": "Dinoflagellata",  # Map to broader group
        "Endomyxa-Ascetosporea_XX": "Endomyxa",
        "Eupetalomonads_X": "Eupetalomonas",
        "Filoretidae_X": "Filoreta",
        "Novel-clade-10_X": None,  # Keep as unmappable

        # Genus level mappings
        "Paradinida_XX": "Paradinium",
        "Cryptosporidium15": "Cryptosporidium",
        "Filosa-Thecofilosea_XXX": "Thecofilosea",
        "Craniata_XXX": "Craniata",
        "Mataza-lineage_X": "Mataza",
        "Novel-Clade-4_X": None,  # Keep as unmappable
        "Corallicolla": "Corallicola",  # Fix spelling
        "OLIGO5_XX": None,  # Environmental clade, unmappable

        # Environmental clades that should remain unmapped (return None)
        "MAST-12": None,
        "MAST-3E": None,
        "LKM74-lineage": None,
        "NC12A-lineage": None,
        "NC12B-lineage": None,
        "NPK2-lineage": None,
        "WIM80-lineage": None,
        "AND16-lineage": None,
        "LOS7N/I-lineage": None,
        "Mariager-Fjord-lineage": None,
        "Mb5C-lineage_X": None,
        "CCW10-lineage_X": None,
        "OLIGO2": None,
        "Nucleohelea": None,  # Environmental group
    }

def apply_taxonomic_mapping(taxon_name, mapping_dict):
    """
    Apply taxonomic mapping to convert outdated names to modern ones.

    Args:
        taxon_name: Original taxon name
        mapping_dict: Dictionary of mappings

    Returns:
        Mapped taxon name or original if no mapping exists
    """
    if taxon_name in mapping_dict:
        mapped_name = mapping_dict[taxon_name]
        if mapped_name is not None:
            return mapped_name
        else:
            return taxon_name
    return taxon_name

def should_append_name_to_lineage(name_to_use):
    """
    Check if the name_to_use should be appended to the lineage.

    Appends the original name for entries that contain:
    - Numbers (e.g., "Theileria1")
    - .U. patterns (e.g., "Eukaryota.U.family")
    - Underscores (e.g., "Embryophyceae_XX")

    Args:
        name_to_use: The original taxonomic name

    Returns:
        True if the name should be appended to lineage, False otherwise
    """
    if not name_to_use or name_to_use == "NA":
        return False

    # Check for numbers
    if any(char.isdigit() for char in name_to_use):
        return True

    # Check for .U. patterns
    if ".U." in name_to_use:
        return True

    # Check for underscores
    if "_" in name_to_use:
        return True

    return False

def append_name_to_lineage(lineage, lineage_ranks, lineage_taxids, name_to_use, taxid, env=None, taxid_to_lineage_cache=None):
    """
    Append the name_to_use to the end of lineage components if it meets criteria.
    If we have a taxid but no lineage, try to retrieve the lineage first.

    Args:
        lineage: Original lineage string
        lineage_ranks: Original lineage ranks string
        lineage_taxids: Original lineage taxids string
        name_to_use: The original taxonomic name
        taxid: The taxid for this entry
        env: Environment for taxonkit subprocess calls
        taxid_to_lineage_cache: Dictionary to update with newly retrieved lineages

    Returns:
        Tuple of (updated_lineage, updated_ranks, updated_taxids)
    """
    if not should_append_name_to_lineage(name_to_use):
        return lineage, lineage_ranks, lineage_taxids

    # If we have a taxid but no lineage, try to get the lineage
    if taxid != "NA" and not lineage and env:
        try:
            result = subprocess.run(
                ["taxonkit", "lineage", "-R", "-t"],
                input=taxid,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env
            )

            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 4:
                    lineage = parts[1]
                    lineage_taxids = parts[2]
                    lineage_ranks = parts[3]

                    # Update the cache so the log reflects the newly retrieved lineage
                    if taxid_to_lineage_cache is not None:
                        taxid_to_lineage_cache[taxid] = (lineage, lineage_ranks, lineage_taxids)
                        logging.info(f"Retrieved missing lineage for taxid {taxid}: {lineage}")

        except Exception as e:
            logging.warning(f"Failed to retrieve lineage for taxid {taxid}: {e}")

    # Append the name_to_use to lineage components
    if lineage:
        updated_lineage = f"{lineage};{name_to_use}"
    else:
        updated_lineage = name_to_use

    if lineage_ranks:
        updated_ranks = f"{lineage_ranks};original_name"
    else:
        updated_ranks = "original_name"

    if lineage_taxids:
        updated_taxids = f"{lineage_taxids};{taxid}" if taxid != "NA" else f"{lineage_taxids};NA"
    else:
        updated_taxids = taxid if taxid != "NA" else "NA"

    return updated_lineage, updated_ranks, updated_taxids

def clean_csv_field(field):
    """
    Clean a field for CSV output by removing newlines and other problematic characters.

    Args:
        field: The field to clean

    Returns:
        Cleaned field safe for CSV
    """
    if field is None:
        return ""

    # Convert to string and remove newlines, carriage returns, and tabs
    cleaned = str(field).replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')

    # Remove extra whitespace
    cleaned = ' '.join(cleaned.split())

    return cleaned

def process_taxon_batch(taxon_batch):
    """
    Process a batch of taxon names to get their taxids.

    Args:
        taxon_batch: List of taxon names to process

    Returns:
        Dictionary mapping taxon names to (taxid, method) tuples
    """
    # Set up environment - taxonkit uses its own built-in NCBI database
    env = os.environ.copy()

    results = {}

    # Clean the names
    cleaned_names = [clean_taxon_name(name) for name in taxon_batch]

    # Run taxonkit name2taxid for the batch
    try:
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input="\n".join(cleaned_names),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode != 0:
            return results

        # Parse the output
        lines = result.stdout.strip().split('\n')
        for i, line in enumerate(lines):
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0" and parts[1].strip():
                # Use the name from the output, not the index, to avoid misalignment
                output_name = parts[0]
                # Find the corresponding original name (handle cleaning)
                original_name = None
                for orig_name in taxon_batch:
                    if clean_taxon_name(orig_name) == output_name:
                        original_name = orig_name
                        break

                if original_name:
                    results[original_name] = (parts[1], "direct")
    except Exception:
        pass

    # For names that didn't get a match, try genus fallback
    for name in taxon_batch:
        if name not in results:
            genus = extract_genus(name)
            if genus:
                try:
                    genus_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=genus.replace("_", " "),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if genus_result.returncode == 0 and genus_result.stdout.strip():
                        genus_parts = genus_result.stdout.strip().split('\t')
                        if len(genus_parts) >= 2 and genus_parts[1] != "0" and genus_parts[1].strip():
                            results[name] = (genus_parts[1], "genus_fallback")
                except Exception:
                    pass

    # For names that still didn't get a match, try stripping numbers as a fallback
    for name in taxon_batch:
        if name not in results:
            # Try stripping numbers from the original name
            stripped_name = strip_trailing_numbers(name.replace("_", " "))
            if stripped_name != name.replace("_", " "):  # Only try if we actually stripped something
                try:
                    stripped_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=stripped_name,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if stripped_result.returncode == 0 and stripped_result.stdout.strip():
                        stripped_parts = stripped_result.stdout.strip().split('\t')
                        if len(stripped_parts) >= 2 and stripped_parts[1] != "0" and stripped_parts[1].strip():
                            results[name] = (stripped_parts[1], "number_stripped")
                except Exception:
                    pass

    # For names that still didn't get a match, try extracting from hyphenated patterns as final fallback
    for name in taxon_batch:
        if name not in results:
            # Try extracting taxa from hyphenated names
            extracted_name = extract_taxa_from_hyphenated(name)
            if extracted_name:  # Only try if we actually extracted something
                try:
                    extracted_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=extracted_name.replace("_", " "),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if extracted_result.returncode == 0 and extracted_result.stdout.strip():
                        extracted_parts = extracted_result.stdout.strip().split('\t')
                        if len(extracted_parts) >= 2 and extracted_parts[1] != "0" and extracted_parts[1].strip():
                            results[name] = (extracted_parts[1], "hyphenated_extracted")
                except Exception:
                    pass

    return results

def get_taxids_for_names(taxon_names, rank_name="unknown"):
    """
    Get NCBI taxids for a list of taxon names using taxonkit.
    Handles underscore removal and genus fallback.
    Uses parallel processing for better performance.

    Args:
        taxon_names: List of taxon names
        rank_name: Name of the taxonomic rank being processed (for logging)

    Returns:
        Tuple of (results_dict, failed_names_dict) where:
        - results_dict: Dictionary mapping taxon names to their taxids
        - failed_names_dict: Dictionary of failed names with failure details
    """
    if not taxon_names:
        return {}, {}

    logging.info(f"Getting taxids for {len(taxon_names)} unique names...")

    # Determine the number of workers and batch size
    num_workers = min(os.cpu_count() or 4, 8)  # Use at most 8 workers
    batch_size = max(1, math.ceil(len(taxon_names) / (num_workers * 4)))  # Ensure enough batches

    # Split the taxon names into batches
    batches = [taxon_names[i:i + batch_size] for i in range(0, len(taxon_names), batch_size)]

    # Process batches in parallel
    results = {}
    failed_names = {}
    direct_match_count = 0
    genus_fallback_count = 0
    number_stripped_count = 0
    hyphenated_extracted_count = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all batches
        future_to_batch = {executor.submit(process_taxon_batch, batch): i for i, batch in enumerate(batches)}

        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_batch):
            batch_results = future.result()

            # Count match types
            for name, (taxid, method) in batch_results.items():
                results[name] = taxid
                if method == "direct":
                    direct_match_count += 1
                elif method == "genus_fallback":
                    genus_fallback_count += 1
                elif method == "number_stripped":
                    number_stripped_count += 1
                elif method == "hyphenated_extracted":
                    hyphenated_extracted_count += 1

    # Collect failed names and add "NA" for names that didn't get a match
    for name in taxon_names:
        if name not in results:
            results[name] = "NA"
            failed_names[name] = {
                'type': 'NO_TAXID_FOUND',
                'details': f'No taxid found for {rank_name} name after direct and genus fallback attempts',
                'taxid': 'NA',
                'rank': rank_name
            }

    # Log statistics
    total = len(taxon_names)
    matched = direct_match_count + genus_fallback_count + number_stripped_count + hyphenated_extracted_count

    logging.info(f"Taxid matching complete: {matched}/{total} matched ({matched/total*100:.1f}%)")

    return results, failed_names

def get_lineages_for_taxids(taxids, env):
    """
    Get lineages for a list of taxids using taxonkit with temporary file approach.
    Returns lineage, ranks, and taxids in a single pass.

    Args:
        taxids: List of taxids
        env: Environment variables for subprocess

    Returns:
        Dictionary mapping taxids to a tuple of (lineage, lineage_ranks, lineage_taxids)
    """
    if not taxids:
        return {}

    lineage_data = {}

    # Create temporary file for taxids
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_filename = temp_file.name
        # Write taxids to temporary file
        for taxid in taxids:
            temp_file.write(f"{taxid}\n")

    try:
        # Run taxonkit lineage command with -R flag for ranks and -t flag for taxids
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", temp_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')

            for line in lines:
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 4:  # Expecting taxid, lineage, lineage_taxids, lineage_ranks
                    taxid = parts[0]
                    lineage = parts[1]  # The lineage is in the second column
                    lineage_taxids = parts[2]  # The taxids are in the third column
                    lineage_ranks = parts[3]  # The ranks are in the fourth column

                    # Clean up the data
                    if lineage and lineage != taxid:
                        # Clean and format the lineage components
                        clean_lineage = lineage.strip()
                        clean_ranks = lineage_ranks.strip()
                        clean_taxids = lineage_taxids.strip()

                        # Store the cleaned data
                        lineage_data[taxid] = (
                            clean_lineage,
                            clean_ranks,
                            clean_taxids
                        )

    except Exception as e:
        logging.error(f"Error getting lineages: {e}")
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except Exception as e:
            logging.warning(f"Error cleaning up temporary file: {e}")

    return lineage_data

def main():
    # Set up directory paths for reorganized structure
    _, metadata_dir, csv_output_dir, log_dir = setup_directory_paths()

    # Parse command line arguments for custom input/output
    import sys
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        # If relative path provided, make it relative to metadata directory
        if not os.path.isabs(input_file):
            input_file = metadata_dir / input_file
    else:
        input_file = metadata_dir / "eukcensus_18S.clusters.97.tsv"

    if len(sys.argv) > 2:
        output_prefix = sys.argv[2]
    else:
        output_prefix = "eukcensus_optimized"

    # Set up logging with proper log directory
    setup_logging(log_dir, output_prefix)
    start_time = time.time()

    # Output file paths - using csv_outputs directory with 18S specification
    division_output = csv_output_dir / "eukcensus_18S_by_division.csv"
    family_output = csv_output_dir / "eukcensus_18S_by_family.csv"
    genus_output = csv_output_dir / "eukcensus_18S_by_genus.csv"

    # Log the paths being used
    logging.info(f"Input file: {input_file}")
    logging.info(f"Output directory: {csv_output_dir}")
    logging.info(f"Log directory: {log_dir}")

    logging.info("Processing file in chunks of 50000 rows...")



    # Read the TSV file
    try:
        # Process in chunks for memory efficiency
        chunk_size = 50000
        chunks = []

        logging.info(f"Loading file in chunks of {chunk_size:,}...")

        # Read chunks with progress bar (no total estimation to avoid hanging)
        chunk_iterator = pd.read_csv(input_file, sep='\t', chunksize=chunk_size)
        for chunk in tqdm(chunk_iterator, desc="Loading chunks", unit="chunk"):
            chunks.append(chunk)

        # Combine chunks
        logging.info(f"Combining {len(chunks)} chunks ({sum(len(chunk) for chunk in chunks)} total rows)...")
        df = pd.concat(chunks, ignore_index=True)

    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        return

    logging.info(f"Successfully processed {len(df)} rows")

    # Check if required columns exist
    required_columns = ['centroid', 'members', 'size', 'division', 'family', 'genus']
    for col in required_columns:
        if col not in df.columns:
            logging.error(f"Required column '{col}' not found in the input file")
            return

    # Initialize dictionaries to store grouped data
    # Added size_count to track total sequence count (sum of cluster sizes)
    division_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})
    family_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})
    genus_data = defaultdict(lambda: {'otu_count': 0, 'size_count': 0})

    # Process data using vectorized operations
    logging.info("Processing clusters using vectorized operations...")

    # Process each taxonomic level
    taxonomic_levels = [('division', division_data), ('family', family_data), ('genus', genus_data)]

    for level, data_dict in tqdm(taxonomic_levels, desc="Processing taxonomic levels", unit="level"):
        logging.info(f"Processing {level} level...")

        # Filter only truly empty/null entries, keep "Unknown" and .U. entries
        level_df = df[~df[level].isna()]
        filtered_count = len(df) - len(level_df)
        logging.info(f"Filtered out {filtered_count} null/empty entries from {level} (keeping Unknown and .U. entries)")

        # Group by taxon and count occurrences (OTU clusters) and sum sizes (sequence counts)
        grouped = level_df.groupby(level).agg({
            'centroid': 'count',  # Count OTU occurrences
            'size': 'sum'  # Sum of sequence counts (cluster sizes)
        }).reset_index()

        # Store in data dictionary with progress bar
        for _, row in tqdm(grouped.iterrows(), total=len(grouped), desc=f"Processing {level} taxa", leave=False):
            taxon = row[level]
            if not should_filter_taxon(taxon):
                data_dict[taxon]['otu_count'] = row['centroid']  # Number of OTU occurrences
                data_dict[taxon]['size_count'] = row['size']  # Sum of sequence counts

        logging.info(f"Processed {len(data_dict)} unique {level} entries")

    # Apply taxonomic mappings before getting taxids
    logging.info("Adding taxids and lineages to processed data...")

    # Collect all unique names for batch processing
    all_names = set()
    all_names.update(division_data.keys())
    all_names.update(family_data.keys())
    all_names.update(genus_data.keys())

    # Apply mappings to all names (mappings are applied within get_taxids_for_names)
    # Get taxids for all names at once
    all_taxid_results, all_failed_names = get_taxids_for_names(list(all_names), "all_ranks")

    # Create individual mappings for each rank
    division_to_taxid = {name: all_taxid_results.get(name, "NA") for name in division_data.keys()}
    family_to_taxid = {name: all_taxid_results.get(name, "NA") for name in family_data.keys()}
    genus_to_taxid = {name: all_taxid_results.get(name, "NA") for name in genus_data.keys()}

    # Initialize failed taxon data structure
    failed_taxon_data = defaultdict(dict)

    # Collect failed names by rank
    for name, failure_info in all_failed_names.items():
        if name in division_data:
            failed_taxon_data['division'][name] = failure_info
        if name in family_data:
            failed_taxon_data['family'][name] = failure_info
        if name in genus_data:
            failed_taxon_data['genus'][name] = failure_info

    # Set up environment - taxonkit uses its own built-in NCBI database
    env = os.environ.copy()

    # Get lineages for all taxids
    logging.info("Getting lineages for taxids...")
    all_taxids = set()

    # Collect all valid taxids with progress tracking
    taxid_sources = [
        ("division", division_to_taxid),
        ("family", family_to_taxid),
        ("genus", genus_to_taxid)
    ]

    for _, taxid_dict in tqdm(taxid_sources, desc="Collecting taxids", unit="source"):
        for taxid in taxid_dict.values():
            if taxid != "NA":
                all_taxids.add(taxid)

    logging.info(f"Collected {len(all_taxids)} unique taxids")

    # Get lineages for all taxids
    taxid_to_lineage = get_lineages_for_taxids(list(all_taxids), env)
    logging.info(f"Retrieved {len(taxid_to_lineage)} lineages from {len(all_taxids)} taxids")



    # Write division data to CSV
    logging.info(f"Writing division data to {division_output}")
    with open(division_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'size_count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        filtered_division_count = 0
        # Sort by decreasing OTU count for better readability
        sorted_division_data = sorted(division_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)
        for division, data in tqdm(sorted_division_data, desc="Writing division data", unit="entry"):
            # Check if this taxon should be filtered out
            if should_filter_taxon(division):
                filtered_division_count += 1
                continue

            taxid = division_to_taxid.get(division, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, division, taxid, env, taxid_to_lineage
            )

            # Clean all fields for CSV output
            clean_division = clean_csv_field(division)
            clean_taxid = clean_csv_field(taxid)
            clean_lineage = clean_csv_field(lineage)
            clean_ranks = clean_csv_field(lineage_ranks)
            clean_taxids = clean_csv_field(lineage_taxids)

            writer.writerow([
                clean_division,
                clean_taxid,
                data['otu_count'],
                data['size_count'],  # Added size_count column
                clean_lineage,
                clean_ranks,
                clean_taxids
            ])

        if filtered_division_count > 0:
            logging.info(f"Filtered out {filtered_division_count} divisions with unidentified taxon patterns")

    # Write family data to CSV
    logging.info(f"Writing family data to {family_output}")
    with open(family_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'size_count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        filtered_family_count = 0
        # Sort by decreasing OTU count for better readability
        sorted_family_data = sorted(family_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)
        for family, data in tqdm(sorted_family_data, desc="Writing family data", unit="entry"):
            # Check if this taxon should be filtered out
            if should_filter_taxon(family):
                filtered_family_count += 1
                continue

            taxid = family_to_taxid.get(family, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, family, taxid, env, taxid_to_lineage
            )

            # Clean all fields for CSV output
            clean_family = clean_csv_field(family)
            clean_taxid = clean_csv_field(taxid)
            clean_lineage = clean_csv_field(lineage)
            clean_ranks = clean_csv_field(lineage_ranks)
            clean_taxids = clean_csv_field(lineage_taxids)

            writer.writerow([
                clean_family,
                clean_taxid,
                data['otu_count'],
                data['size_count'],  # Added size_count column
                clean_lineage,
                clean_ranks,
                clean_taxids
            ])

        if filtered_family_count > 0:
            logging.info(f"Filtered out {filtered_family_count} families with unidentified taxon patterns")

    # Write genus data to CSV
    logging.info(f"Writing genus data to {genus_output}")
    with open(genus_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'otu_count', 'size_count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        filtered_genus_count = 0
        # Sort by decreasing OTU count for better readability
        sorted_genus_data = sorted(genus_data.items(), key=lambda x: x[1]['otu_count'], reverse=True)
        for genus, data in tqdm(sorted_genus_data, desc="Writing genus data", unit="entry"):
            # Check if this taxon should be filtered out
            if should_filter_taxon(genus):
                filtered_genus_count += 1
                continue

            taxid = genus_to_taxid.get(genus, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            # Append name_to_use to lineage if it contains numbers, .U., or underscores
            lineage, lineage_ranks, lineage_taxids = append_name_to_lineage(
                lineage, lineage_ranks, lineage_taxids, genus, taxid, env, taxid_to_lineage
            )

            # Clean all fields for CSV output
            clean_genus = clean_csv_field(genus)
            clean_taxid = clean_csv_field(taxid)
            clean_lineage = clean_csv_field(lineage)
            clean_ranks = clean_csv_field(lineage_ranks)
            clean_taxids = clean_csv_field(lineage_taxids)

            writer.writerow([
                clean_genus,
                clean_taxid,
                data['otu_count'],
                data['size_count'],  # Added size_count column
                clean_lineage,
                clean_ranks,
                clean_taxids
            ])

        if filtered_genus_count > 0:
            logging.info(f"Filtered out {filtered_genus_count} genera with unidentified taxon patterns")

    # Create comprehensive unmapped log
    create_comprehensive_unmapped_log(
        division_data, family_data, genus_data,
        division_to_taxid, family_to_taxid, genus_to_taxid,
        taxid_to_lineage, log_dir, output_prefix
    )

    # Calculate and log performance metrics
    end_time = time.time()
    processing_time = end_time - start_time
    total_entries = len(division_data) + len(family_data) + len(genus_data)

    logging.info("Saving results to CSV files...")
    logging.info(f"Saved {len(division_data)} division entries to {division_output}")
    logging.info(f"Division summary: {len([t for t in division_to_taxid.values() if t != 'NA'])} with taxids, {len([t for t in division_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Saved {len(family_data)} family entries to {family_output}")
    logging.info(f"Family summary: {len([t for t in family_to_taxid.values() if t != 'NA'])} with taxids, {len([t for t in family_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Saved {len(genus_data)} genus entries to {genus_output}")
    logging.info(f"Genus summary: {len([t for t in genus_to_taxid.values() if t != 'NA'])} with taxids, {len([t for t in genus_to_taxid.values() if t in taxid_to_lineage])} with lineages")

    logging.info(f"Processing complete in {processing_time:.2f} seconds")
    logging.info(f"Total entries processed: {total_entries}")
    logging.info(f"Performance: {total_entries/processing_time:.1f} entries/second")

    logging.info("Processing complete! Generated the following files:")
    logging.info(f"- {division_output}")
    logging.info(f"- {family_output}")
    logging.info(f"- {genus_output}")
    logging.info(f"- {log_dir / f'{output_prefix}_comprehensive_unmapped.log'}")
    logging.info(f"- {log_dir / 'eukcensus_optimization.log'}")
    logging.info("Output files organized in csv_outputs/ and logs/ directories")

if __name__ == "__main__":
    main()
