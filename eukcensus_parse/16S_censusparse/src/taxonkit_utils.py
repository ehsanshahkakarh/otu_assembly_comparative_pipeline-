"""
Taxonkit Integration for 16S Census Parser
==========================================

Handles NCBI taxonomy ID lookup and lineage retrieval using taxonkit
with optimized batch processing and smart fallback strategies.
"""

import subprocess
import os
import logging
from tqdm import tqdm

from .taxon_cleaner import (
    clean_organelle_taxon_name,
    extract_meaningful_taxonomic_part,
    strip_trailing_numbers
)
from .organelle_handler import vectorized_organelle_detection


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


def get_taxids_using_taxonkit(taxon_names, rank_name):
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
    organellar_names, _ = vectorized_organelle_detection(taxon_names)

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

                # Parse results - more robust parsing that matches names to taxids
                for line in lines:
                    if not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        variant_name = parts[0]
                        taxid = parts[1].strip()

                        # Only store valid taxids (not "0" or empty)
                        if taxid and taxid != "0":
                            # Only store if this variant was actually requested
                            if variant_name in all_variants:
                                variant_to_taxid[variant_name] = taxid

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
            failed_taxids = []

            for line in tqdm(lines, desc="Parsing lineage results", leave=False, ncols=80):
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 4:  # taxid, lineage, lineage_taxids, lineage_ranks
                    taxid = parts[0]
                    lineage = parts[1]
                    lineage_taxids = parts[2]
                    lineage_ranks = parts[3]

                    # Only store if we have actual lineage data
                    if lineage and lineage.strip():
                        taxid_to_lineage[taxid] = (lineage, lineage_ranks, lineage_taxids)
                    else:
                        failed_taxids.append(taxid)
                else:
                    # Malformed line - try to extract taxid for error reporting
                    if parts:
                        failed_taxids.append(parts[0])

            # Report any failures
            if failed_taxids:
                print(f"⚠️  Failed to get lineages for {len(failed_taxids)} taxids")
                if len(failed_taxids) <= 10:  # Show details for small numbers
                    print(f"Failed taxids: {', '.join(failed_taxids)}")
        else:
            print(f"❌ taxonkit lineage failed: {result.stderr}")

        success_count = len(taxid_to_lineage)
        total_requested = len(valid_taxids)
        print(f"✅ Successfully retrieved lineages for {success_count}/{total_requested} taxids ({success_count/total_requested*100:.1f}%)")

    except Exception as e:
        print(f"❌ Error in taxonkit lineage: {e}")

    return taxid_to_lineage

