#!/usr/bin/env python3
"""
Enhanced EukCensus 16S Cluster Parser with Organelle Handling and Taxonomic Rank Filtering

This script processes EukCensus 16S cluster metadata with enhanced features:
1. Better organelle entry handling (e.g., Vitis_vinifera:plas.Chloroplast)
2. Taxonomic rank filtering to remove inappropriate rank entries
3. Improved species name cleaning for organellar sequences

Key improvements:
- Enhanced organelle detection and cleaning
- Taxonomic rank validation against lineage information
- Better handling of species-level entries in genus/family parsing
- Improved logging and verification

Output files:
- eukcensus_by_phylum_enhanced.csv
- eukcensus_by_family_enhanced.csv
- eukcensus_by_genus_enhanced.csv
"""

import pandas as pd
import os
import subprocess
import tempfile
from collections import defaultdict
import csv
import logging
import sys
import time
from datetime import datetime
import re
from tqdm import tqdm

def setup_logging(output_prefix="eukcensus_enhanced"):
    """
    Set up logging configuration for the enhanced script.

    Args:
        output_prefix: Prefix for output files to include in log messages
    """
    # Create log file path
    log_file = "eukcensus_enhanced_processing.log"

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
    logging.info(f"🚀 Starting enhanced EukCensus processing with organelle handling: {output_prefix}_*")
    return logging.getLogger(__name__)

def detect_organelle_type(taxon_name):
    """
    Detect the type of organelle from taxon name.
    
    Args:
        taxon_name: The taxon name to analyze
        
    Returns:
        Tuple of (is_organelle, organelle_type, cleaned_name)
    """
    organelle_patterns = {
        'chloroplast': ['.Chloroplast', ':plas.Chloroplast', '.plas.Chloroplast', 'chloroplast'],
        'mitochondria': ['.Mitochondria', ':mito.Mitochondria', '.mito.Mitochondria', 'mitochondria'],
        'plastid': ['.Plastid', ':plas.Plastid', '.plas.Plastid', 'plastid'],
        'apicoplast': ['.Apicoplast', ':api.Apicoplast', '.api.Apicoplast', 'apicoplast']
    }
    
    taxon_lower = taxon_name.lower()
    
    for organelle_type, patterns in organelle_patterns.items():
        for pattern in patterns:
            if pattern.lower() in taxon_lower:
                # Extract the base name before the organelle indicator
                if '.' in taxon_name:
                    base_name = taxon_name.split('.')[0]
                elif ':' in taxon_name:
                    # Handle cases like "Vitis_vinifera:plas.Chloroplast"
                    base_name = taxon_name.split(':')[0]
                else:
                    base_name = taxon_name
                
                # Further clean any remaining organelle indicators
                if ':plas' in base_name:
                    base_name = base_name.replace(':plas', '')
                if ':mito' in base_name:
                    base_name = base_name.replace(':mito', '')
                if ':api' in base_name:
                    base_name = base_name.replace(':api', '')
                
                return True, organelle_type, base_name
    
    return False, None, taxon_name

def extract_meaningful_taxonomic_part(taxon_name):
    """
    Extract meaningful taxonomic information from complex names.

    Handles cases like:
    - uncultured_Alphaproteobacteria_bacterium -> Alphaproteobacteria
    - uncultured_Rickettsia_sp -> Rickettsia
    - marine_metagenome -> None (not taxonomically meaningful)

    Args:
        taxon_name: The taxon name to analyze

    Returns:
        Meaningful taxonomic part or None if no meaningful part found
    """
    if not taxon_name or pd.isna(taxon_name):
        return None

    # Skip purely environmental/technical terms
    environmental_terms = ['metagenome', 'environmental', 'sample', 'clone', 'specimen']
    if any(term in taxon_name.lower() for term in environmental_terms):
        return None

    # For uncultured names, try to extract the taxonomic part
    if 'uncultured' in taxon_name.lower():
        parts = taxon_name.replace('_', ' ').split()

        # Look for meaningful taxonomic terms (not generic terms)
        generic_terms = ['uncultured', 'bacterium', 'organism', 'eukaryote', 'sp', 'species']
        meaningful_parts = []

        for part in parts:
            if part.lower() not in generic_terms and len(part) > 2:
                # Check if this looks like a taxonomic name
                if part[0].isupper() and part[1:].islower():
                    meaningful_parts.append(part)
                elif part.endswith('bacteria') or part.endswith('proteobacteria'):
                    meaningful_parts.append(part)

        if meaningful_parts:
            # Return the most specific meaningful part
            return meaningful_parts[-1] if len(meaningful_parts) == 1 else ' '.join(meaningful_parts[:2])

    return None

def clean_organelle_taxon_name(taxon_name):
    """
    Enhanced cleaning for organelle-containing taxon names with improved Candidatus handling.

    Handles cases like:
    - Vitis_vinifera:plas.Chloroplast -> Vitis vinifera
    - uncultured_bacterium.Mitochondria -> uncultured bacterium
    - Genus_species.Plastid -> Genus species
    - Candidatus names -> preserve (no stripping needed with updated NCBI taxonomy)

    Args:
        taxon_name: The taxon name to clean

    Returns:
        Cleaned taxon name suitable for taxonomic lookup
    """
    if pd.isna(taxon_name) or taxon_name == "Unknown":
        return taxon_name

    # Detect and handle organelles
    is_organelle, organelle_type, base_name = detect_organelle_type(taxon_name)

    if is_organelle:
        # Quietly handle organelle detection
        taxon_name = base_name

    # Candidatus taxa are now preserved in NCBI taxonomy - no stripping needed
    # The NCBI taxonomic mapping scripts have been updated to handle Candidatus taxa properly

    # Replace underscores with spaces for species-level names
    if '_' in taxon_name:
        # Check if this looks like a binomial species name
        parts = taxon_name.split('_')
        if len(parts) == 2 and not any(char.isdigit() for char in parts[1]):
            # Likely a species name like Genus_species
            cleaned = taxon_name.replace('_', ' ')
            return strip_trailing_numbers(cleaned)
        else:
            # Handle other underscore cases
            cleaned = taxon_name.replace('_', ' ')
            return strip_trailing_numbers(cleaned)

    # Handle trailing numbers
    return strip_trailing_numbers(taxon_name)

def extract_genus_from_species(species_name):
    """
    Extract genus name from a species name.
    
    Args:
        species_name: Full species name (e.g., "Vitis vinifera")
        
    Returns:
        Genus name (e.g., "Vitis")
    """
    if not species_name or pd.isna(species_name):
        return species_name
    
    parts = species_name.strip().split()
    if len(parts) >= 1:
        return parts[0]
    
    return species_name

def recover_organelle_taxonomy(taxon_name, target_rank):
    """
    Recover taxonomic information from organelle sequences by extracting host organism info.

    NOTE: This function is currently DISABLED in the main processing workflow to prevent
    performance bottlenecks. It's kept for potential future use or debugging.
    Organellar sequences are now handled by vectorized_organelle_detection() instead.

    Args:
        taxon_name: Original taxonomic name (e.g., "Prasiola_crispa.Chloroplast")
        target_rank: Target rank (phylum, family, genus)

    Returns:
        Tuple of (host_species_name, appropriate_rank_name) or (None, None) if recovery fails
    """
    # Quick check if this is an organelle sequence - return early if not
    organelle_patterns = ['.Chloroplast', '.Mitochondria', '.Apicoplast', '.Plastid', ':Chloroplast', ':Mitochondria', ':Apicoplast', ':Plastid']

    host_species = None
    for pattern in organelle_patterns:
        if pattern in taxon_name:
            # Extract host species name
            host_species = taxon_name.split(pattern)[0]
            break

    # Early return if no organelle pattern found - this avoids expensive subprocess calls
    if not host_species:
        return None, None

    # Clean the host species name
    host_species = host_species.replace('_', ' ').strip()

    # Get taxid and lineage for host species
    try:
        # Use taxonkit to get taxid for host species
        env = os.environ.copy()
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=host_species,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
            cwd="."
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0":
                host_taxid = parts[1]

                # Get lineage for host taxid
                lineage_result = subprocess.run(
                    ["taxonkit", "lineage", "-R"],
                    input=host_taxid,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                    cwd="."
                )

                if lineage_result.returncode == 0 and lineage_result.stdout.strip():
                    lineage_parts = lineage_result.stdout.strip().split('\t')
                    if len(lineage_parts) >= 4:  # taxid, lineage, lineage_taxids, lineage_ranks
                        lineage = lineage_parts[1]
                        lineage_ranks = lineage_parts[3]

                        # Extract appropriate rank from lineage
                        lineage_list = lineage.split(';')
                        ranks_list = lineage_ranks.split(';')

                        # Find the target rank in the lineage
                        for i, rank in enumerate(ranks_list):
                            if rank.lower() == target_rank.lower():
                                if i < len(lineage_list):
                                    return host_species, lineage_list[i]

                        # If exact rank not found, use fallback logic
                        if target_rank == "genus" and len(lineage_list) > 0:
                            # For genus, try to extract genus from species name
                            genus = host_species.split()[0] if ' ' in host_species else host_species
                            return host_species, genus
                        elif target_rank == "family" and "family" in [r.lower() for r in ranks_list]:
                            family_idx = [r.lower() for r in ranks_list].index("family")
                            if family_idx < len(lineage_list):
                                return host_species, lineage_list[family_idx]
                        elif target_rank == "phylum" and "phylum" in [r.lower() for r in ranks_list]:
                            phylum_idx = [r.lower() for r in ranks_list].index("phylum")
                            if phylum_idx < len(lineage_list):
                                return host_species, lineage_list[phylum_idx]

    except Exception as e:
        # If organelle recovery fails, return None
        pass

    return None, None

def extract_appropriate_rank_name(taxon_name, target_rank, lineage_info=None):
    """
    Extract the appropriate taxonomic rank name from a taxon entry with DISABLED organelle recovery.

    For example, if parsing at genus level but entry is species-level,
    extract the genus portion.

    NOTE: Organelle recovery is now DISABLED here to prevent performance bottlenecks.
    Organellar sequences are handled by the vectorized approach in get_taxids_using_taxonkit_optimized().

    Args:
        taxon_name: Original taxon name
        target_rank: Target taxonomic rank ('genus', 'family', 'phylum')
        lineage_info: Optional lineage information tuple (lineage, ranks, taxids)

    Returns:
        Appropriate rank name or None if should be filtered
    """
    # PERFORMANCE FIX: Disable individual organelle recovery here
    # Organellar sequences will be handled by vectorized approach later
    # This prevents thousands of individual subprocess calls during initial processing

    cleaned_name = clean_organelle_taxon_name(taxon_name)
    
    # If we have lineage information, use it to validate rank appropriateness
    if lineage_info and lineage_info[1]:  # lineage_ranks exists
        lineage, ranks, taxids = lineage_info
        rank_list = ranks.split(';') if ranks else []
        
        # Check if the entry matches our target rank
        if target_rank in rank_list:
            # Find the position of our target rank
            try:
                target_idx = rank_list.index(target_rank)
                lineage_parts = lineage.split(';') if lineage else []
                if target_idx < len(lineage_parts):
                    return lineage_parts[target_idx]
            except (ValueError, IndexError):
                pass
    
    # Fallback to name-based extraction
    if target_rank == 'genus':
        # For genus level, extract genus from species names
        if ' ' in cleaned_name and len(cleaned_name.split()) >= 2:
            # Likely a species name, extract genus
            return extract_genus_from_species(cleaned_name)
        else:
            # Already genus-level or higher
            return cleaned_name  
    
    elif target_rank == 'family':
        # For family level, we generally keep the name as-is unless it's clearly species
        if ' ' in cleaned_name and len(cleaned_name.split()) >= 2:
            # This is a species name, which is inappropriate for family-level parsing
            return None
        return cleaned_name
    
    elif target_rank == 'phylum':
        # For phylum level, filter out species and genus names
        if ' ' in cleaned_name:
            # This is likely a species name
            return None
        # Could still be genus or family, but we'll let taxonkit determine appropriateness
        return cleaned_name
    
    return cleaned_name

def strip_trailing_numbers(taxon_name):
    """
    Remove trailing numbers from taxon names (e.g., "Theileria1" -> "Theileria").

    Args:
        taxon_name: The taxon name to process

    Returns:
        The taxon name with trailing numbers removed
    """
    if not taxon_name:
        return taxon_name
        
    # Handle taxa with numbers at the end (e.g., "Cryptosporidium15", "Eimeria1", "Plasmodium1")
    match = re.match(r'^(.+?)(\d+)$', taxon_name.strip())
    if match:
        base_name = match.group(1).rstrip()  # Remove any trailing spaces
        return base_name

    # If no trailing numbers found, return as is
    return taxon_name

def should_filter_taxon(taxon_name):
    """
    Check if a taxon name should be filtered out.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    if not taxon_name or pd.isna(taxon_name):
        return True
        
    # Filter out entries with ".U.phylum", ".U.genus", etc.
    if any(pattern in taxon_name for pattern in [".U.phylum", ".U.genus", ".U.family", ".U.order", ".U.class", ".U.species", ".U.division"]):
        return True
    
    # Filter out clearly unidentified entries
    unidentified_patterns = ['unidentified', 'unknown', 'uncultured', 'environmental']
    taxon_lower = taxon_name.lower()
    
    # Only filter if the entire name is just an unidentified pattern
    if taxon_lower in unidentified_patterns:
        return True

    return False

def validate_rank_appropriateness(taxon_name, target_rank, lineage_info):
    """
    Validate if a taxon entry is appropriate for the target taxonomic rank.
    
    Args:
        taxon_name: Original taxon name
        target_rank: Target rank ('genus', 'family', 'phylum')
        lineage_info: Tuple of (lineage, lineage_ranks, lineage_taxids)
        
    Returns:
        True if appropriate, False if should be filtered
    """
    if not lineage_info or not lineage_info[1]:  # No rank information
        return True  # Let it through, will be handled by name-based logic
    
    lineage, ranks, taxids = lineage_info
    rank_list = ranks.split(';') if ranks else []
    
    # Check if the lineage contains our target rank
    if target_rank not in rank_list:
        logging.info(f"⚠️ Taxon '{taxon_name}' lacks {target_rank} rank in lineage, filtering out")
        return False
    
    # For genus-level parsing, filter out entries that are clearly family or higher
    if target_rank == 'genus':
        if 'species' in rank_list:
            # Species-level entry is OK for genus parsing (we extract genus)
            return True
        elif 'genus' in rank_list:
            # Genus-level entry is perfect
            return True
        else:
            # Higher-level entry (family, order, etc.) - might not be appropriate
            logging.info(f"⚠️ Taxon '{taxon_name}' appears to be higher than genus level, filtering out")
            return False
    
    # For family-level parsing
    elif target_rank == 'family':
        if 'species' in rank_list or 'genus' in rank_list:
            # Too specific for family-level parsing
            logging.info(f"⚠️ Taxon '{taxon_name}' is too specific for family-level parsing, filtering out")
            return False
        elif 'family' in rank_list:
            return True
        else:
            # Higher level might be OK
            return True
    
    # For phylum-level parsing
    elif target_rank == 'phylum':
        if any(rank in rank_list for rank in ['species', 'genus', 'family']):
            # Too specific for phylum-level parsing
            logging.info(f"⚠️ Taxon '{taxon_name}' is too specific for phylum-level parsing, filtering out")
            return False
        elif 'phylum' in rank_list:
            return True
        else:
            return True
    
    return True

def get_single_taxid_with_fallbacks(taxon_name, env):
    """
    Get taxid for a single taxon name with multiple fallback strategies.

    Args:
        taxon_name: The taxon name to lookup
        env: Environment variables for subprocess

    Returns:
        Tuple of (taxid, method_used) where method_used indicates which approach worked
    """
    fallback_strategies = []

    # Strategy 1: Original name (especially important for Candidatus)
    fallback_strategies.append(("original", taxon_name))

    # Strategy 2: Cleaned name (organelle removal, underscore handling)
    cleaned_name = clean_organelle_taxon_name(taxon_name)
    if cleaned_name != taxon_name:
        fallback_strategies.append(("cleaned", cleaned_name))

    # Strategy 3: Candidatus names are now handled properly by NCBI taxonomy
    # No need to strip Candidatus prefix - NCBI scripts preserve these taxa

    # Strategy 4: Extract meaningful taxonomic part for complex uncultured names
    meaningful_part = extract_meaningful_taxonomic_part(taxon_name)
    if meaningful_part and meaningful_part not in [s[1] for s in fallback_strategies]:
        fallback_strategies.append(("meaningful_part", meaningful_part))

    # Strategy 5: For names with trailing numbers, try without numbers
    no_numbers = strip_trailing_numbers(taxon_name)
    if no_numbers != taxon_name and no_numbers not in [s[1] for s in fallback_strategies]:
        fallback_strategies.append(("no_numbers", no_numbers))

    # Try each strategy
    for method, name_to_try in fallback_strategies:
        if not name_to_try or name_to_try.lower() in ['uncultured', 'unknown', 'environmental']:
            continue

        try:
            result = subprocess.run(
                ["taxonkit", "name2taxid"],
                input=name_to_try,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                cwd=".",
                timeout=5
            )

            if result.returncode == 0 and result.stdout.strip():
                # Handle multiple results by taking the first line
                first_line = result.stdout.strip().split('\n')[0]
                parts = first_line.strip().split('\t')
                if len(parts) >= 2 and parts[1] != "0" and parts[1].strip():
                    return parts[1].strip(), method

        except Exception:
            continue

    return "NA", "failed_all_strategies"

def vectorized_organelle_detection(taxon_names):
    """
    Vectorized organelle detection and host extraction using pandas operations.
    Much faster than individual name processing.

    Args:
        taxon_names: List of taxon names

    Returns:
        tuple: (organellar_names_dict, all_lookup_names)
            - organellar_names_dict: {original_name: host_organism_name}
            - all_lookup_names: list of names to use for taxid lookup
    """
    import pandas as pd

    # Convert to pandas Series for vectorized operations
    names_series = pd.Series(taxon_names)

    # Define organelle patterns
    organelle_patterns = ['.Chloroplast', '.Mitochondria', '.Apicoplast', '.Plastid',
                         ':Chloroplast', ':Mitochondria', ':Apicoplast', ':Plastid']

    # Vectorized organelle detection
    is_organellar = names_series.str.contains('|'.join([p.replace('.', r'\.').replace(':', r':') for p in organelle_patterns]), na=False)

    organellar_names = {}
    lookup_names = []

    if is_organellar.any():
        # Process organellar sequences
        organellar_series = names_series[is_organellar]

        # Vectorized host extraction
        for pattern in organelle_patterns:
            mask = organellar_series.str.contains(pattern.replace('.', r'\.').replace(':', r':'), na=False)
            if mask.any():
                # Extract host organisms for this pattern
                hosts = organellar_series[mask].str.split(pattern, expand=True)[0].str.replace('_', ' ').str.strip()

                # Map original names to host names
                for orig_name, host_name in zip(organellar_series[mask].index, hosts):
                    original_name = names_series.iloc[orig_name]
                    if host_name:
                        organellar_names[original_name] = host_name
                        lookup_names.append(host_name)
                    else:
                        lookup_names.append(original_name)

    # Add non-organellar names
    non_organellar_names = names_series[~is_organellar].tolist()
    lookup_names.extend(non_organellar_names)

    print(f"🧬 Vectorized organelle detection: {len(organellar_names)} organellar sequences found")
    if len(organellar_names) > 0:
        print(f"⚡ PERFORMANCE BOOST: Avoided {len(organellar_names) * 2} individual subprocess calls!")

    return organellar_names, lookup_names


def get_taxids_using_taxonkit_optimized(taxon_names, rank_name):
    """
    Optimized taxid lookup with batch processing and smart fallback strategies.
    Enhanced with organelle handling to infer taxonomy from host organisms.

    Speed optimizations:
    1. Batch all fallback strategies together instead of individual calls
    2. Pre-generate all possible name variants
    3. Single large taxonkit call with all variants
    4. Smart mapping back to original names
    5. Special handling for organellar sequences

    Args:
        taxon_names: List of taxon names
        rank_name: Name of the taxonomic rank (for logging)

    Returns:
        Dictionary mapping taxon names to taxids
    """
    if not taxon_names:
        return {}

    print(f"🔍 Getting taxids for {len(taxon_names)} {rank_name} names (optimized with organelle handling)...")

    # Use default environment
    env = os.environ.copy()

    # STEP 1: Vectorized organelle detection and variant generation
    print(f"📝 Step 1: Vectorized organelle detection and variant generation...")

    # Use vectorized organelle detection for much better performance
    organellar_names, primary_lookup_names = vectorized_organelle_detection(taxon_names)

    name_variants = {}  # original_name -> [list of variants to try]
    all_variants = []   # flat list of all variants
    variant_to_original = {}  # variant -> original_name

    for original_name in taxon_names:
        variants = []

        # Use the result from vectorized organelle detection
        if original_name in organellar_names:
            # For organellar sequences, use host organism name as primary variant
            host_organism = organellar_names[original_name]
            variants.append(host_organism)
        else:
            # Variant 1: Original name (for non-organellar sequences)
            variants.append(original_name)

        # Variant 2: Cleaned name
        cleaned = clean_organelle_taxon_name(original_name)
        if cleaned != original_name and cleaned not in variants:
            variants.append(cleaned)

        # Variant 3: Candidatus names are now preserved in NCBI taxonomy
        # No need to strip Candidatus prefix - NCBI scripts handle these properly

        # Variant 4: Meaningful part for uncultured
        meaningful = extract_meaningful_taxonomic_part(original_name)
        if meaningful and meaningful not in variants:
            variants.append(meaningful)

        # Variant 5: Numbers stripped
        no_numbers = strip_trailing_numbers(original_name)
        if no_numbers != original_name and no_numbers not in variants:
            variants.append(no_numbers)

        # Store variants for this name
        name_variants[original_name] = variants

        # Add to flat list and create reverse mapping
        for variant in variants:
            if variant and variant.lower() not in ['uncultured', 'unknown', 'environmental']:
                all_variants.append(variant)
                if variant not in variant_to_original:
                    variant_to_original[variant] = []
                variant_to_original[variant].append(original_name)

    print(f"📊 Generated {len(all_variants)} total variants for {len(taxon_names)} names")

    # STEP 2: Single batch taxonkit call for ALL variants
    print(f"📝 Step 2: Running single batch taxonkit call...")

    variant_to_taxid = {}

    if all_variants:
        variants_input = "\n".join(all_variants)

        try:
            with tqdm(total=1, desc=f"Running optimized batch taxonkit", leave=False, ncols=80) as pbar:
                result = subprocess.run(
                    ["taxonkit", "name2taxid"],
                    input=variants_input,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                    cwd="."
                )
                pbar.update(1)

            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                valid_lines = [line for line in lines if line.strip()]

                # Parse results
                line_idx = 0
                for variant in all_variants:
                    if line_idx < len(valid_lines):
                        line = valid_lines[line_idx]
                        parts = line.strip().split('\t')

                        if len(parts) >= 2 and parts[1] != "0" and parts[1].strip():
                            variant_to_taxid[variant] = parts[1].strip()
                            line_idx += 1

                            # Skip multiple results for same variant
                            while (line_idx < len(valid_lines) and
                                   valid_lines[line_idx].startswith(parts[0] + '\t')):
                                line_idx += 1
                        else:
                            line_idx += 1

        except Exception as e:
            logging.error(f"❌ Error in optimized taxonkit call: {e}")

    # STEP 3: Map results back to original names using priority order
    print(f"📝 Step 3: Mapping results back to original names...")

    name_to_taxid = {}
    success_count = 0

    for original_name in taxon_names:
        found_taxid = None

        # Try variants in priority order
        for variant in name_variants[original_name]:
            if variant in variant_to_taxid:
                found_taxid = variant_to_taxid[variant]
                break

        if found_taxid:
            name_to_taxid[original_name] = found_taxid
            success_count += 1
        else:
            name_to_taxid[original_name] = "NA"

    # Count organellar sequences processed
    organellar_count = len(organellar_names)
    if organellar_count > 0:
        print(f"🧬 Processed {organellar_count} organellar sequences (taxonomy inferred from host organisms)")

    print(f"✅ Optimized processing: Successfully mapped {success_count}/{len(taxon_names)} {rank_name} names")
    print(f"⚡ Speed improvement: Single batch call vs {len(taxon_names)} individual calls")

    return name_to_taxid

# Keep the old function as fallback
def get_taxids_using_taxonkit(taxon_names, rank_name):
    """Wrapper to use optimized version by default."""
    return get_taxids_using_taxonkit_optimized(taxon_names, rank_name)

def get_lineages_using_taxonkit(taxids):
    """
    Get lineages for taxids using taxonkit.

    Args:
        taxids: List of taxids

    Returns:
        Dictionary mapping taxids to (lineage, lineage_ranks, lineage_taxids) tuples
    """
    if not taxids:
        return {}

    # Filter out "NA" taxids
    valid_taxids = [tid for tid in taxids if tid != "NA"]
    if not valid_taxids:
        return {}

    print(f"🧬 Getting lineages for {len(valid_taxids)} taxids...")

    # Use default environment - let taxonkit find its own database
    env = os.environ.copy()

    taxid_to_lineage = {}

    try:
        # Run taxonkit lineage
        taxids_input = "\n".join(valid_taxids)

        with tqdm(total=1, desc="Running taxonkit lineage", leave=False, ncols=80) as pbar:
            result = subprocess.run(
                ["taxonkit", "lineage", "-R", "-t"],
                input=taxids_input,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                cwd="."
            )
            pbar.update(1)

        if result.returncode == 0 and result.stdout.strip():
            # Parse results
            lines = result.stdout.strip().split('\n')
            for line in tqdm(lines, desc="Parsing lineage results", leave=False, ncols=80):
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 4:  # taxid, lineage, lineage_taxids, lineage_ranks
                    taxid = parts[0]
                    lineage = parts[1]
                    lineage_taxids = parts[2]
                    lineage_ranks = parts[3]

                    taxid_to_lineage[taxid] = (lineage, lineage_ranks, lineage_taxids)
        else:
            print(f"❌ taxonkit lineage failed: {result.stderr}")

        print(f"✅ Successfully retrieved lineages for {len(taxid_to_lineage)} taxids")

    except Exception as e:
        print(f"❌ Error in taxonkit lineage: {e}")

    return taxid_to_lineage

def create_comprehensive_unmapped_log(phylum_data, family_data, genus_data,
                                    phylum_to_taxid, family_to_taxid, genus_to_taxid,
                                    taxid_to_lineage, output_prefix):
    """
    Create a comprehensive log of all unmapped taxonomic names with enhanced analysis.

    Args:
        phylum_data, family_data, genus_data: Data dictionaries for each rank
        phylum_to_taxid, family_to_taxid, genus_to_taxid: Taxid mapping dictionaries
        taxid_to_lineage: Lineage information dictionary
        output_prefix: Prefix for output files
    """
    log_file = f"{output_prefix}_comprehensive_unmapped.log"
    print(f"📝 Creating enhanced comprehensive unmapped log...")

    with open(log_file, 'w') as f:
        f.write("# Enhanced Comprehensive Unmapped Names Log - 16S Census Parser with Fallback Analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# This log contains all taxonomic names that failed to get NCBI taxids or lineages\n")
        f.write("# Enhanced with fallback strategy analysis and improved pattern recognition\n")
        f.write("# Format: Rank | Original_Name | Cleaned_Name | Appropriate_Name | Size_Count | Occurrence_Count | Taxid | Reason\n\n")

        # Summary statistics
        f.write("=== SUMMARY STATISTICS ===\n")

        # Calculate statistics for each rank
        rank_stats = {}
        for rank_name, data_dict, taxid_dict in [
            ('phylum', phylum_data, phylum_to_taxid),
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
            ('phylum', phylum_data, phylum_to_taxid),
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
                        'cleaned_name': data.get('cleaned_name', ''),
                        'appropriate_name': data.get('appropriate_name', ''),
                        'size_count': data['size_count'],
                        'count': data['count'],
                        'taxid': taxid,
                        'reason': reason
                    })

            f.write(f"Total unmapped {rank_name} entries: {len(unmapped_entries)}\n\n")

            # Sort by occurrence count (descending)
            unmapped_entries.sort(key=lambda x: x['count'], reverse=True)

            for entry in unmapped_entries:
                f.write(f"{rank_name.upper()} | {entry['original_name']} | {entry['cleaned_name']} | {entry['appropriate_name']} | {entry['size_count']} | {entry['count']} | {entry['taxid']} | {entry['reason']}\n")

            f.write(f"\n")

        # Pattern analysis
        f.write("=== PATTERN ANALYSIS ===\n")

        # Collect all unmapped names for pattern analysis
        all_unmapped = []
        for rank_name, data_dict, taxid_dict in [
            ('phylum', phylum_data, phylum_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    all_unmapped.append({
                        'rank': rank_name,
                        'name': orig_name,
                        'count': data['count']
                    })

        # Analyze patterns
        patterns = {
            'Candidatus': [n for n in all_unmapped if 'Candidatus' in n['name']],
            'Organelles': [n for n in all_unmapped if any(org in n['name'] for org in ['Mitochondria', 'Chloroplast', 'Apicoplast', 'Plastid', ':plas', ':mito', ':api'])],
            'Uncultured': [n for n in all_unmapped if 'uncultured' in n['name'].lower()],
            'Environmental': [n for n in all_unmapped if any(env in n['name'].lower() for env in ['metagenome', 'environmental', 'marine'])],
            'Species_level': [n for n in all_unmapped if '_' in n['name'] and '.' not in n['name']],
            'Modern_phyla': [n for n in all_unmapped if n['name'].endswith('ota') or n['name'].endswith('eia')],
            'Candidate_division': [n for n in all_unmapped if 'candidate division' in n['name'].lower()]
        }

        for pattern_name, pattern_names in patterns.items():
            if pattern_names:
                count = len(pattern_names)
                total_occurrences = sum(n['count'] for n in pattern_names)
                f.write(f"{pattern_name.replace('_', ' ').title()}: {count} names ({total_occurrences} total occurrences)\n")

                # Show top 5 most frequent for each pattern
                top_names = sorted(pattern_names, key=lambda x: x['count'], reverse=True)[:5]
                for name_info in top_names:
                    f.write(f"  {name_info['name']} ({name_info['rank']}) - {name_info['count']} occurrences\n")
                f.write("\n")

        f.write("=== RECOMMENDATIONS ===\n")
        f.write("1. Organelle sequences (mitochondria, chloroplast, plastid, apicoplast) are now mapped to host organisms\n")
        f.write("2. Original organellar names are preserved to prevent merging with non-organellar OTUs\n")
        f.write("3. Candidatus taxa are now properly preserved in NCBI taxonomy (no filtering)\n")
        f.write("4. Modern bacterial phylum names may need mapping to older NCBI names\n")
        f.write("5. Environmental/uncultured samples may not have valid NCBI taxids\n")
        f.write("6. Species-level entries need genus extraction for genus-level parsing\n")
        f.write("7. Consider using GTDB taxonomy database for better bacterial coverage\n")

    print(f"✅ Comprehensive unmapped log written to {log_file}")
    return log_file

def process_taxonomic_level(df, col_name, target_rank):
    """
    Process a specific taxonomic level with enhanced filtering and rank validation.

    Args:
        df: Input dataframe
        col_name: Column name for this taxonomic level
        target_rank: Target taxonomic rank

    Returns:
        Dictionary of processed data for this level, keyed by ORIGINAL names
    """
    print(f"📊 Processing {target_rank} level...")

    # Filter out unidentified taxa and NaN values
    level_df = df[~df[col_name].isna() & (df[col_name] != "Unknown")]
    filtered_count = len(df) - len(level_df)

    # Initialize data dictionary - KEY CHANGE: Use original names as keys
    data_dict = defaultdict(lambda: {'size_count': 0, 'count': 0, 'cleaned_name': '', 'appropriate_name': ''})

    # Process each entry with progress bar
    for _, row in tqdm(level_df.iterrows(), total=len(level_df), desc=f"Processing {target_rank} entries", leave=False, ncols=80):
        original_taxon = row[col_name]

        if should_filter_taxon(original_taxon):
            continue

        # Extract appropriate rank name for taxonkit lookup
        appropriate_name = extract_appropriate_rank_name(original_taxon, target_rank)

        if appropriate_name is None:
            continue  # Filtered out

        # Store data using ORIGINAL name as key
        data_dict[original_taxon]['size_count'] += row['size']
        data_dict[original_taxon]['count'] += 1
        data_dict[original_taxon]['cleaned_name'] = clean_organelle_taxon_name(original_taxon)
        data_dict[original_taxon]['appropriate_name'] = appropriate_name

    print(f"✅ Processed {len(data_dict)} unique {target_rank} entries (filtered {filtered_count} unidentified)")

    return dict(data_dict)

def main():
    """
    PERFORMANCE OPTIMIZED VERSION:
    - Disabled individual organelle recovery during initial processing
    - All organellar sequences now handled by vectorized approach
    - This prevents thousands of individual subprocess calls
    - Expected speedup: 10-100x for datasets with many organellar sequences
    """
    # Parse command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "eukcensus_16S.clusters.97.tsv"

    if len(sys.argv) > 2:
        output_prefix = sys.argv[2]
    else:
        output_prefix = "eukcensus_enhanced"

    # Set up logging
    setup_logging(output_prefix)
    start_time = time.time()

    # Output file names
    output_dir = "."
    phylum_output = os.path.join(output_dir, f"{output_prefix}_by_phylum.csv")
    family_output = os.path.join(output_dir, f"{output_prefix}_by_family.csv")
    genus_output = os.path.join(output_dir, f"{output_prefix}_by_genus.csv")

    print(f"📁 Processing file: {input_file}")

    # Read the TSV file
    try:
        print(f"📖 Reading input file...")
        df = pd.read_csv(input_file, sep='\t')
        print(f"✅ Successfully loaded {len(df)} rows")
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        return

    # Check if required columns exist
    required_columns = ['centroid', 'members', 'size', 'phylum', 'familiy', 'genus']
    for col in required_columns:
        if col not in df.columns:
            print(f"❌ Required column '{col}' not found in the input file")
            return

    # Process each taxonomic level with enhanced logic
    phylum_data = process_taxonomic_level(df, 'phylum', 'phylum')
    family_data = process_taxonomic_level(df, 'familiy', 'family')
    genus_data = process_taxonomic_level(df, 'genus', 'genus')

    # Get taxids for each taxonomic level using appropriate names for lookup
    print("🔬 Getting taxids using enhanced taxonkit processing...")

    # OPTIMIZATION: Vectorized collection and processing of unique names
    print("🚀 Vectorized optimization: Collecting unique names across all ranks...")

    all_unique_names = set()
    all_unique_names.update(phylum_data.keys())
    all_unique_names.update(family_data.keys())
    all_unique_names.update(genus_data.keys())

    print(f"📊 Total unique names across all ranks: {len(all_unique_names)}")

    # Single vectorized lookup for all unique names (with optimized organelle handling)
    all_names_to_taxid = get_taxids_using_taxonkit_optimized(list(all_unique_names), "all_ranks")

    # Map results to individual ranks
    phylum_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in phylum_data.keys()}
    family_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in family_data.keys()}
    genus_to_taxid = {name: all_names_to_taxid.get(name, "NA") for name in genus_data.keys()}

    # Collect all valid taxids
    all_taxids = set()
    for taxid_dict in [phylum_to_taxid, family_to_taxid, genus_to_taxid]:
        for taxid in taxid_dict.values():
            if taxid != "NA":
                all_taxids.add(taxid)

    print(f"🧬 Collected {len(all_taxids)} unique taxids for lineage retrieval")

    # Get lineages for all taxids
    taxid_to_lineage = get_lineages_using_taxonkit(list(all_taxids))

    # Write results to CSV files with enhanced information
    print("💾 Writing results to CSV files...")

    # Write phylum data
    with open(phylum_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'size_count', 'count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_phyla = sorted(phylum_data.items(), key=lambda x: x[1]['count'], reverse=True)

        for phylum, data in tqdm(sorted_phyla, desc="Writing phylum data", leave=False, ncols=80):
            taxid = phylum_to_taxid.get(phylum, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            writer.writerow([
                phylum,  # Original name preserved
                taxid,
                data['size_count'],
                data['count'],
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Write family data
    with open(family_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'size_count', 'count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_families = sorted(family_data.items(), key=lambda x: x[1]['count'], reverse=True)

        for family, data in tqdm(sorted_families, desc="Writing family data", leave=False, ncols=80):
            taxid = family_to_taxid.get(family, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            writer.writerow([
                family,  # Original name preserved
                taxid,
                data['size_count'],
                data['count'],
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Write genus data
    with open(genus_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name_to_use', 'taxid', 'size_count', 'count', 'lineage', 'lineage_ranks', 'lineage_taxids'])

        sorted_genera = sorted(genus_data.items(), key=lambda x: x[1]['count'], reverse=True)

        for genus, data in tqdm(sorted_genera, desc="Writing genus data", leave=False, ncols=80):
            taxid = genus_to_taxid.get(genus, "NA")
            lineage_info = taxid_to_lineage.get(taxid, ("", "", "")) if taxid != "NA" else ("", "", "")
            lineage, lineage_ranks, lineage_taxids = lineage_info

            writer.writerow([
                genus,  # Original name preserved
                taxid,
                data['size_count'],
                data['count'],
                lineage,
                lineage_ranks,
                lineage_taxids
            ])

    # Create comprehensive unmapped log
    unmapped_log = create_comprehensive_unmapped_log(
        phylum_data, family_data, genus_data,
        phylum_to_taxid, family_to_taxid, genus_to_taxid,
        taxid_to_lineage, output_prefix
    )

    # Calculate and log performance metrics
    end_time = time.time()
    processing_time = end_time - start_time
    total_entries = len(phylum_data) + len(family_data) + len(genus_data)

    # Calculate success rates
    phylum_success = len([t for t in phylum_to_taxid.values() if t != 'NA'])
    family_success = len([t for t in family_to_taxid.values() if t != 'NA'])
    genus_success = len([t for t in genus_to_taxid.values() if t != 'NA'])

    phylum_success_rate = (phylum_success / len(phylum_data) * 100) if len(phylum_data) > 0 else 0
    family_success_rate = (family_success / len(family_data) * 100) if len(family_data) > 0 else 0
    genus_success_rate = (genus_success / len(genus_data) * 100) if len(genus_data) > 0 else 0

    print(f"✅ Saved {len(phylum_data)} phylum entries to {phylum_output}")
    print(f"📊 Phylum success rate: {phylum_success}/{len(phylum_data)} ({phylum_success_rate:.1f}%)")

    print(f"✅ Saved {len(family_data)} family entries to {family_output}")
    print(f"📊 Family success rate: {family_success}/{len(family_data)} ({family_success_rate:.1f}%)")

    print(f"✅ Saved {len(genus_data)} genus entries to {genus_output}")
    print(f"📊 Genus success rate: {genus_success}/{len(genus_data)} ({genus_success_rate:.1f}%)")

    print(f"🎉 Vectorized processing complete in {processing_time:.2f} seconds")
    print(f"📊 Total entries processed: {total_entries}")
    print(f"⚡ Performance: {total_entries/processing_time:.1f} entries/second")
    print(f"🚀 Vectorized organelle handling and optimized batch processing enabled")

    print("Enhanced processing complete! Generated the following files:")
    print(f"- {phylum_output}")
    print(f"- {family_output}")
    print(f"- {genus_output}")
    print(f"- {output_prefix}_comprehensive_unmapped.log")
    print(f"- {output_prefix}_processing.log")

if __name__ == "__main__":
    main()