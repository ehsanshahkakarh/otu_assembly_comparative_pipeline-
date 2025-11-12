#!/usr/bin/env python3
"""
Parse EukCensus 18S Clusters TSV File

This script parses the eukcensus_18S.clusters.97.tsv file (metadata) and generates three CSV files
organized by division, family, and genus. Each file contains columns for taxid, taxon_name,
domain, and member_size.

The script uses taxonkit to get taxids for each taxon name and properly handles organelle information
in taxon names (e.g., "Aspergillus_nidulans_FGSC_A4.Mitochondria").

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
from tqdm import tqdm

def clean_taxon_name(taxon_name):
    """
    Clean a taxon name by removing organelle information, unclassified indicators, and underscores.

    Args:
        taxon_name: The taxon name to clean

    Returns:
        The cleaned taxon name
    """
    # Handle special cases with organelle information or unclassified indicators
    if "." in taxon_name:
        # For names like "Genus_species.Mitochondria", get the part before the dot
        parts = taxon_name.split(".")
        name_part = parts[0]

        # Replace underscores with spaces
        return name_part.replace("_", " ")

    # For regular names, just replace underscores with spaces
    return taxon_name.replace("_", " ")

def extract_genus(taxon_name):
    """
    Extract the genus part from a taxon name.

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
            return name_part.split("_")[0]
        return name_part

    # For names like "Genus_species", get the first part
    if "_" in taxon_name:
        return taxon_name.split("_")[0]

    # For names like "Genus", return as is
    return taxon_name

def should_filter_taxon(taxon_name):
    """
    Check if a taxon name should be filtered out.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    # Filter out entries with ".U.phylum", ".U.genus", etc.
    if any(pattern in taxon_name for pattern in [".U.phylum", ".U.genus", ".U.family", ".U.order", ".U.class", ".U.species", ".U.division"]):
        return True

    return False

def process_taxon_batch(taxon_batch):
    """
    Process a batch of taxon names to get their taxids.

    Args:
        taxon_batch: List of taxon names to process

    Returns:
        Dictionary mapping taxon names to (taxid, method) tuples
    """
    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir

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
            if len(parts) >= 2 and parts[1] != "0":
                original_name = taxon_batch[i]
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
                        if len(genus_parts) >= 2 and genus_parts[1] != "0":
                            results[name] = (genus_parts[1], "genus_fallback")
                except Exception:
                    pass

    return results

def get_taxids_for_names(taxon_names):
    """
    Get NCBI taxids for a list of taxon names using taxonkit.
    Handles underscore removal and genus fallback.
    Uses parallel processing for better performance.

    Args:
        taxon_names: List of taxon names

    Returns:
        Dictionary mapping taxon names to their taxids
    """
    if not taxon_names:
        return {}

    print(f"Getting taxids for {len(taxon_names)} taxon names...")

    # Determine the number of workers and batch size
    num_workers = min(os.cpu_count() or 4, 8)  # Use at most 8 workers
    batch_size = max(1, math.ceil(len(taxon_names) / (num_workers * 4)))  # Ensure enough batches

    # Split the taxon names into batches
    batches = [taxon_names[i:i + batch_size] for i in range(0, len(taxon_names), batch_size)]

    # Process batches in parallel
    results = {}
    direct_match_count = 0
    genus_fallback_count = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all batches
        future_to_batch = {executor.submit(process_taxon_batch, batch): i for i, batch in enumerate(batches)}

        # Process results as they complete
        for future in tqdm(concurrent.futures.as_completed(future_to_batch), total=len(batches), desc="Processing taxon batches"):
            batch_results = future.result()

            # Count match types
            for name, (taxid, method) in batch_results.items():
                results[name] = taxid
                if method == "direct":
                    direct_match_count += 1
                elif method == "genus_fallback":
                    genus_fallback_count += 1

    # Add "NA" for names that didn't get a match
    for name in taxon_names:
        if name not in results:
            results[name] = "NA"

    # Print statistics
    total = len(taxon_names)
    matched = direct_match_count + genus_fallback_count
    not_matched = total - matched

    print("\nTaxid Matching Statistics:")
    print(f"Total taxon names: {total}")
    print(f"Direct matches: {direct_match_count} ({direct_match_count/total*100:.1f}%)")
    print(f"Genus fallback matches: {genus_fallback_count} ({genus_fallback_count/total*100:.1f}%)")
    print(f"Total matched: {matched} ({matched/total*100:.1f}%)")
    print(f"Not matched: {not_matched} ({not_matched/total*100:.1f}%)")

    return results

def get_domains_batch(taxids, env):
    """
    Get domains for a batch of taxids using taxonkit.

    Args:
        taxids: List of taxids to check
        env: Environment variables for subprocess

    Returns:
        Dictionary mapping taxids to domain names
    """
    # Filter out "NA" taxids
    valid_taxids = [taxid for taxid in taxids if taxid != "NA"]

    if not valid_taxids:
        return {}

    domains = {}

    try:
        # Run taxonkit lineage to get the lineage with names for all taxids at once
        result = subprocess.run(
            ["taxonkit", "lineage", "--show-name"],
            input="\n".join(valid_taxids),
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
                if len(parts) >= 3:
                    taxid = parts[0]
                    lineage = parts[2]  # The lineage with names is in the third column

                    # Determine the domain
                    if "Eukaryota" in lineage:
                        domains[taxid] = "Eukaryota"
                    elif "Bacteria" in lineage:
                        domains[taxid] = "Bacteria"
                    elif "Archaea" in lineage:
                        domains[taxid] = "Archaea"
                    elif "Viruses" in lineage:
                        domains[taxid] = "Virus"
                    else:
                        domains[taxid] = "Unknown"

    except Exception as e:
        print(f"Error getting domains for taxids: {e}")

    return domains

def get_domain(taxid, env, domain_cache=None):
    """
    Get the domain (Bacteria, Archaea, Virus, or Eukaryota) for a taxid using taxonkit.
    Uses a cache to avoid redundant lookups.

    Args:
        taxid: The taxid to check
        env: Environment variables for subprocess
        domain_cache: Optional cache of taxid to domain mappings

    Returns:
        The domain name as a string, or "Unknown" if it can't be determined
    """
    if taxid == "NA":
        return "Unknown"

    # Check cache first if provided
    if domain_cache is not None and taxid in domain_cache:
        return domain_cache[taxid]

    try:
        # Run taxonkit lineage to get the lineage with names
        result = subprocess.run(
            ["taxonkit", "lineage", "--show-name"],
            input=taxid,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            lineage = result.stdout.strip()

            # Check for each domain
            if "Eukaryota" in lineage:
                return "Eukaryota"
            elif "Bacteria" in lineage:
                return "Bacteria"
            elif "Archaea" in lineage:
                return "Archaea"
            elif "Viruses" in lineage:
                return "Virus"

            # If we can't determine the domain from the lineage
            return "Unknown"

    except Exception as e:
        print(f"Error getting domain for taxid {taxid}: {e}")
        pass

    return "Unknown"

def main():
    # Input and output file paths
    input_file = "eukcensus_18S.clusters.97.tsv"
    output_dir = "."

    # Output file names
    division_output = os.path.join(output_dir, "eukcensus_by_division.csv")
    family_output = os.path.join(output_dir, "eukcensus_by_family.csv")
    genus_output = os.path.join(output_dir, "eukcensus_by_genus.csv")

    print(f"Reading input file: {input_file}")

    # Read the TSV file
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    print(f"Successfully read {len(df)} rows from the input file")

    # Check if required columns exist
    required_columns = ['centroid', 'members', 'size', 'division', 'family', 'genus']
    for col in required_columns:
        if col not in df.columns:
            print(f"Error: Required column '{col}' not found in the input file")
            return

    # Initialize dictionaries to store grouped data
    division_data = defaultdict(lambda: {'size': 0, 'centroids': [], 'count': 0})
    family_data = defaultdict(lambda: {'size': 0, 'centroids': [], 'count': 0})
    genus_data = defaultdict(lambda: {'size': 0, 'centroids': [], 'count': 0})

    # Process each row in the dataframe
    print("Processing data and grouping by taxonomic ranks...")
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Grouping by taxonomic ranks"):
        # Extract values
        centroid = row['centroid']

        # Convert size to integer if it's not already
        try:
            size = int(row['size'])
        except (ValueError, TypeError):
            size = 0
            # Don't print warnings for every conversion error to reduce verbosity

        division = row['division'] if not pd.isna(row['division']) else "Unknown"
        family = row['family'] if not pd.isna(row['family']) else "Unknown"
        genus = row['genus'] if not pd.isna(row['genus']) else "Unknown"

        # Update division data
        division_data[division]['size'] += size
        division_data[division]['centroids'].append(centroid)
        division_data[division]['count'] += 1  # Increment count for each occurrence

        # Update family data
        family_data[family]['size'] += size
        family_data[family]['centroids'].append(centroid)
        family_data[family]['count'] += 1  # Increment count for each occurrence

        # Update genus data
        genus_data[genus]['size'] += size
        genus_data[genus]['centroids'].append(centroid)
        genus_data[genus]['count'] += 1  # Increment count for each occurrence

    # Get taxids for divisions
    print("Getting taxids for divisions...")
    division_names = list(division_data.keys())
    division_to_taxid = get_taxids_for_names(division_names)

    # Get taxids for families
    print("Getting taxids for families...")
    family_names = list(family_data.keys())
    family_to_taxid = get_taxids_for_names(family_names)

    # Get taxids for genera
    print("Getting taxids for genera...")
    genus_names = list(genus_data.keys())
    genus_to_taxid = get_taxids_for_names(genus_names)

    # Set up environment with TAXONKIT_DB for domain check
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir

    # Create a domain cache to avoid redundant lookups
    domain_cache = {}

    # Collect all taxids that need domain lookup
    all_taxids = set()

    # Process division taxids
    division_entries = []
    filtered_division_count = 0
    for division, data in division_data.items():
        # Check if this taxon should be filtered out
        if should_filter_taxon(division):
            filtered_division_count += 1
            continue

        taxid = division_to_taxid.get(division, "NA")

        # Special handling for organelle information
        if "." in division:
            organism_part = division.split(".")[0]
            if organism_part in division_to_taxid:
                organism_taxid = division_to_taxid[organism_part]
                if organism_taxid != "NA":
                    all_taxids.add(organism_taxid)

        if taxid != "NA":
            all_taxids.add(taxid)

        division_entries.append((division, data, taxid))

    # Process family taxids
    family_entries = []
    filtered_family_count = 0
    for family, data in family_data.items():
        # Check if this taxon should be filtered out
        if should_filter_taxon(family):
            filtered_family_count += 1
            continue

        taxid = family_to_taxid.get(family, "NA")

        # Special handling for organelle information
        if "." in family:
            organism_part = family.split(".")[0]
            if organism_part in family_to_taxid:
                organism_taxid = family_to_taxid[organism_part]
                if organism_taxid != "NA":
                    all_taxids.add(organism_taxid)

        if taxid != "NA":
            all_taxids.add(taxid)

        family_entries.append((family, data, taxid))

    # Process genus taxids
    genus_entries = []
    filtered_genus_count = 0
    for genus, data in genus_data.items():
        # Check if this taxon should be filtered out
        if should_filter_taxon(genus):
            filtered_genus_count += 1
            continue

        taxid = genus_to_taxid.get(genus, "NA")

        # Special handling for organelle information
        if "." in genus:
            organism_part = genus.split(".")[0]
            if organism_part in genus_to_taxid:
                organism_taxid = genus_to_taxid[organism_part]
                if organism_taxid != "NA":
                    all_taxids.add(organism_taxid)

        if taxid != "NA":
            all_taxids.add(taxid)

        genus_entries.append((genus, data, taxid))

    # Batch lookup domains for all taxids
    print(f"Looking up domains for {len(all_taxids)} unique taxids...")
    all_taxids_list = list(all_taxids)

    # Process in batches of 1000 to avoid command line length limits
    batch_size = 1000
    for i in range(0, len(all_taxids_list), batch_size):
        batch = all_taxids_list[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(all_taxids_list) + batch_size - 1)//batch_size}...")
        batch_domains = get_domains_batch(batch, env)
        domain_cache.update(batch_domains)

    # Write division data to CSV
    print(f"Writing division data to {division_output}")
    with open(division_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['taxid', 'taxon_name', 'domain', 'member_size', 'occurrence_count'])

        for division, data, taxid in division_entries:
            # Since this is 18S data, default to Eukaryota
            domain = "Eukaryota"

            # Special handling for organelle information
            if "." in division:
                organism_part = division.split(".")[0]
                if organism_part in division_to_taxid:
                    organism_taxid = division_to_taxid[organism_part]
                    if organism_taxid in domain_cache:
                        organism_domain = domain_cache[organism_taxid]
                        if organism_domain == "Eukaryota":
                            print(f"✅ Keeping organelle entry in division data: {division} (organism part: {organism_part} is eukaryotic)")
                            taxid = organism_taxid

            # Check if we have a valid taxid with a different domain in the cache
            if taxid != "NA" and taxid in domain_cache and domain_cache[taxid] != "Unknown":
                # Only override the default if we have a definitive non-Unknown domain
                domain = domain_cache[taxid]

            writer.writerow([taxid, division, domain, data['size'], data['count']])

        if filtered_division_count > 0:
            print(f"Filtered out {filtered_division_count} divisions with unidentified taxon patterns")

    # Write family data to CSV
    print(f"Writing family data to {family_output}")
    with open(family_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['taxid', 'taxon_name', 'domain', 'member_size', 'occurrence_count'])

        for family, data, taxid in family_entries:
            # Since this is 18S data, default to Eukaryota
            domain = "Eukaryota"

            # Special handling for organelle information
            if "." in family:
                organism_part = family.split(".")[0]
                if organism_part in family_to_taxid:
                    organism_taxid = family_to_taxid[organism_part]
                    if organism_taxid in domain_cache:
                        organism_domain = domain_cache[organism_taxid]
                        if organism_domain == "Eukaryota":
                            print(f"✅ Keeping organelle entry in family data: {family} (organism part: {organism_part} is eukaryotic)")
                            taxid = organism_taxid

            # Check if we have a valid taxid with a different domain in the cache
            if taxid != "NA" and taxid in domain_cache and domain_cache[taxid] != "Unknown":
                # Only override the default if we have a definitive non-Unknown domain
                domain = domain_cache[taxid]

            writer.writerow([taxid, family, domain, data['size'], data['count']])

        if filtered_family_count > 0:
            print(f"Filtered out {filtered_family_count} families with unidentified taxon patterns")

    # Write genus data to CSV
    print(f"Writing genus data to {genus_output}")
    with open(genus_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['taxid', 'taxon_name', 'domain', 'member_size', 'occurrence_count'])

        for genus, data, taxid in genus_entries:
            # Since this is 18S data, default to Eukaryota
            domain = "Eukaryota"

            # Special handling for organelle information
            if "." in genus:
                organism_part = genus.split(".")[0]
                if organism_part in genus_to_taxid:
                    organism_taxid = genus_to_taxid[organism_part]
                    if organism_taxid in domain_cache:
                        organism_domain = domain_cache[organism_taxid]
                        if organism_domain == "Eukaryota":
                            print(f"✅ Keeping organelle entry in genus data: {genus} (organism part: {organism_part} is eukaryotic)")
                            taxid = organism_taxid

            # Check if we have a valid taxid with a different domain in the cache
            if taxid != "NA" and taxid in domain_cache and domain_cache[taxid] != "Unknown":
                # Only override the default if we have a definitive non-Unknown domain
                domain = domain_cache[taxid]

            writer.writerow([taxid, genus, domain, data['size'], data['count']])

        if filtered_genus_count > 0:
            print(f"Filtered out {filtered_genus_count} genera with unidentified taxon patterns")

    print("Done! Generated the following files:")
    print(f"- {division_output}")
    print(f"- {family_output}")
    print(f"- {genus_output}")

if __name__ == "__main__":
    main()
