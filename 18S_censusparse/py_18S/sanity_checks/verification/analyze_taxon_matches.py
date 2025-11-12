#!/usr/bin/env python3
"""
Analyze Taxon Name Matches Across Taxonomic Levels

This script analyzes the taxon names (division/phylum, family, and genus) directly from
the original 18S clusters TSV file, identifies exact matches, merges matching entries,
and provides occurrence counts.

The script reads the original TSV file:
- eukcensus_18S.clusters.97.tsv

It then:
1. Extracts taxon names at each taxonomic level
2. Identifies exact matches across taxonomic levels
3. Merges matching entries and counts occurrences
4. Generates reports showing which taxa were successfully merged

Output files:
- merged_taxa_report.csv: Contains all merged taxa with occurrence counts
- division_matches.csv: Division-specific matches
- family_matches.csv: Family-specific matches
- genus_matches.csv: Genus-specific matches
- multi_level_taxa.csv: Taxa present at multiple taxonomic levels
"""

import pandas as pd
import os
import subprocess
import tempfile
from collections import defaultdict
# Try to import tqdm, fall back to simple progress if not available
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    def tqdm(iterable, desc="Processing", total=None):
        """Simple fallback progress reporter"""
        if total is None:
            total = len(iterable) if hasattr(iterable, '__len__') else 0

        for i, item in enumerate(iterable):
            if total > 0 and (i % max(1, total // 20) == 0 or i == total - 1):
                percent = (i + 1) / total * 100
                print(f"\r{desc}: {i+1}/{total} ({percent:.1f}%)", end="", flush=True)
            yield item
        print()  # New line at the end

def clean_taxon_name(taxon_name):
    """
    Clean a taxon name by removing organelle information, unclassified indicators, and underscores.
    Extract only the important taxonomic portion.

    Args:
        taxon_name: The taxon name to clean

    Returns:
        The cleaned taxon name
    """
    if pd.isna(taxon_name):
        return None

    import re

    # Convert to string and strip whitespace
    name = str(taxon_name).strip()

    # Handle special cases with organelle information or unclassified indicators
    if "." in name:
        # For names like "Genus_species.Mitochondria", get the part before the dot
        parts = name.split(".")
        name = parts[0]

    # Remove underscores and take only the part before the first underscore
    if "_" in name:
        name = name.split("_")[0]

    # Remove any trailing numbers, hyphens, or special characters
    # Keep only alphabetic characters at the start
    name = re.sub(r'[^a-zA-Z\s].*$', '', name).strip()

    # Remove extra whitespace
    name = re.sub(r'\s+', ' ', name).strip()

    # If the name is too short or contains only numbers/special chars, return None
    if len(name) < 2 or not re.match(r'^[a-zA-Z]', name):
        return None

    return name

def extract_clean_taxid(taxid_value):
    """
    Extract a clean numeric taxid from a potentially messy string.

    Args:
        taxid_value: The taxid value to clean

    Returns:
        Clean numeric taxid or "NA" if not found
    """
    if pd.isna(taxid_value) or taxid_value == "NA":
        return "NA"

    import re

    # Convert to string
    taxid_str = str(taxid_value).strip()

    # Extract all numeric sequences
    numbers = re.findall(r'\d+', taxid_str)

    if numbers:
        # Take the first (and usually longest) numeric sequence
        # Filter out very short numbers (likely not taxids)
        valid_numbers = [num for num in numbers if len(num) >= 3]
        if valid_numbers:
            return valid_numbers[0]
        elif numbers:
            return numbers[0]

    return "NA"

def get_taxid_for_name(taxon_name, failed_log=None):
    """
    Get NCBI taxid for a taxon name using taxonkit.

    Args:
        taxon_name: The taxon name to look up
        failed_log: List to append failed lookups to

    Returns:
        The taxid as a string, or "NA" if not found
    """
    if pd.isna(taxon_name) or not taxon_name:
        if failed_log is not None:
            failed_log.append(f"Empty/NaN name: {taxon_name}")
        return "NA"

    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir

    try:
        # Clean the name for taxonkit
        cleaned_name = clean_taxon_name(taxon_name)
        if not cleaned_name:
            if failed_log is not None:
                failed_log.append(f"Failed cleaning: '{taxon_name}' -> None")
            return "NA"

        # Run taxonkit name2taxid
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=cleaned_name,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0":
                # Extract clean taxid from the result
                clean_taxid = extract_clean_taxid(parts[1])
                if clean_taxid != "NA":
                    return clean_taxid

        # If we get here, the lookup failed
        if failed_log is not None:
            failed_log.append(f"Taxonkit lookup failed: '{taxon_name}' -> '{cleaned_name}'")

    except Exception as e:
        if failed_log is not None:
            failed_log.append(f"Error with '{taxon_name}': {str(e)}")

    return "NA"

def should_filter_taxon(taxon_name):
    """
    Check if a taxon name should be filtered out.

    Args:
        taxon_name: The taxon name to check

    Returns:
        True if the taxon should be filtered out, False otherwise
    """
    if pd.isna(taxon_name):
        return True

    # Filter out entries with ".U.phylum", ".U.genus", etc.
    if any(pattern in taxon_name for pattern in [".U.phylum", ".U.genus", ".U.family", ".U.order", ".U.class", ".U.species", ".U.division"]):
        return True

    return False

def get_ranks_for_taxids(taxids):
    """
    Get taxonomic ranks for a list of taxids using taxonkit.

    Args:
        taxids: List of taxids to look up

    Returns:
        Dictionary mapping taxids to their ranks
    """
    if not taxids:
        return {}

    # Set up environment with TAXONKIT_DB
    env = os.environ.copy()
    taxdump_dir = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"
    env["TAXONKIT_DB"] = taxdump_dir

    # Create a temporary file with taxids
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_path = temp_file.name
        for taxid in taxids:
            if str(taxid) != "NA":
                temp_file.write(f"{taxid}\n")

    try:
        # Run taxonkit lineage with rank information
        result = subprocess.run(
            ["taxonkit", "lineage", "-r", temp_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        taxid_to_rank = {}

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            for line in lines:
                if not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    taxid = parts[0]
                    rank = parts[2]
                    taxid_to_rank[taxid] = rank
                elif len(parts) == 2:
                    taxid = parts[0]
                    taxid_to_rank[taxid] = "unknown"

        # Add "unknown" for any taxids that weren't found
        for taxid in taxids:
            if str(taxid) != "NA" and str(taxid) not in taxid_to_rank:
                taxid_to_rank[str(taxid)] = "not_found"

        return taxid_to_rank

    except Exception as e:
        print(f"Error getting ranks: {e}")
        return {}

    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_path)
        except Exception:
            pass

def analyze_taxon_matches():
    """
    Analyze taxon name matches across division, family, and genus directly from the TSV file.
    """
    # Input file path
    input_file = "eukcensus_18S.clusters.97.tsv"

    # Initialize failed names log
    failed_names_log = []

    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        return

    print(f"Reading input file: {input_file}")

    # Read the TSV file
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully read {len(df)} rows from {input_file}")
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Check if required columns exist
    required_columns = ['centroid', 'members', 'size', 'division', 'family', 'genus']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Input file is missing columns: {missing_columns}")
        return

    # Initialize dictionaries to store grouped data
    division_data = defaultdict(lambda: {'size': 0, 'centroids': set(), 'count': 0})
    family_data = defaultdict(lambda: {'size': 0, 'centroids': set(), 'count': 0})
    genus_data = defaultdict(lambda: {'size': 0, 'centroids': set(), 'count': 0})

    # Process each row in the dataframe
    print("Processing data and grouping by taxonomic ranks...")
    if TQDM_AVAILABLE:
        iterator = tqdm(df.iterrows(), total=len(df), desc="Grouping by taxonomic ranks",
                       ncols=80, leave=True, position=0)
    else:
        iterator = tqdm(df.iterrows(), desc="Grouping by taxonomic ranks", total=len(df))

    for _, row in iterator:
        # Extract values
        centroid = row['centroid']

        # Convert size to integer if it's not already
        try:
            size = int(row['size'])
        except (ValueError, TypeError):
            size = 0

        # Process division
        division = row['division'] if not pd.isna(row['division']) else None
        if division and not should_filter_taxon(division):
            division_data[division]['size'] += size
            division_data[division]['centroids'].add(centroid)
            division_data[division]['count'] += 1

        # Process family
        family = row['family'] if not pd.isna(row['family']) else None
        if family and not should_filter_taxon(family):
            family_data[family]['size'] += size
            family_data[family]['centroids'].add(centroid)
            family_data[family]['count'] += 1

        # Process genus
        genus = row['genus'] if not pd.isna(row['genus']) else None
        if genus and not should_filter_taxon(genus):
            genus_data[genus]['size'] += size
            genus_data[genus]['centroids'].add(centroid)
            genus_data[genus]['count'] += 1

    # Create dictionaries to store taxon information by name
    division_taxa = {}
    family_taxa = {}
    genus_taxa = {}

    # Process division data
    print("Processing division data...")
    if TQDM_AVAILABLE:
        div_iterator = tqdm(division_data.items(), desc="Processing divisions",
                           ncols=80, leave=True, position=0)
    else:
        div_iterator = tqdm(division_data.items(), desc="Processing divisions",
                           total=len(division_data))

    for name, data in div_iterator:
        taxid = get_taxid_for_name(name, failed_names_log)
        division_taxa[name] = {
            'taxon_name': name,
            'taxid': taxid,
            'domain': 'Eukaryota',  # Default for 18S data
            'member_size': data['size'],
            'occurrence_count': data['count']
        }

    # Process family data
    print("Processing family data...")
    if TQDM_AVAILABLE:
        fam_iterator = tqdm(family_data.items(), desc="Processing families",
                           ncols=80, leave=True, position=0)
    else:
        fam_iterator = tqdm(family_data.items(), desc="Processing families",
                           total=len(family_data))

    for name, data in fam_iterator:
        taxid = get_taxid_for_name(name, failed_names_log)
        family_taxa[name] = {
            'taxon_name': name,
            'taxid': taxid,
            'domain': 'Eukaryota',  # Default for 18S data
            'member_size': data['size'],
            'occurrence_count': data['count']
        }

    # Process genus data
    print("Processing genus data...")
    if TQDM_AVAILABLE:
        gen_iterator = tqdm(genus_data.items(), desc="Processing genera",
                           ncols=80, leave=True, position=0)
    else:
        gen_iterator = tqdm(genus_data.items(), desc="Processing genera",
                           total=len(genus_data))

    for name, data in gen_iterator:
        taxid = get_taxid_for_name(name, failed_names_log)
        genus_taxa[name] = {
            'taxon_name': name,
            'taxid': taxid,
            'domain': 'Eukaryota',  # Default for 18S data
            'member_size': data['size'],
            'occurrence_count': data['count']
        }

    # Find exact matches across all taxonomic levels
    all_taxon_names = set(division_taxa.keys()) | set(family_taxa.keys()) | set(genus_taxa.keys())
    print(f"Found {len(all_taxon_names)} unique taxon names across all taxonomic levels")

    # Create a dictionary to store occurrence counts and merged data
    taxon_occurrences = defaultdict(lambda: {'division': False, 'family': False, 'genus': False,
                                            'occurrence_count': 0, 'member_size': 0,
                                            'taxid': 'NA', 'domain': 'Eukaryota'})

    # Process division taxa
    for name, data in division_taxa.items():
        taxon_occurrences[name]['division'] = True
        taxon_occurrences[name]['occurrence_count'] += data['occurrence_count']
        taxon_occurrences[name]['member_size'] += data['member_size']
        taxon_occurrences[name]['taxid'] = data['taxid']

    # Process family taxa
    for name, data in family_taxa.items():
        taxon_occurrences[name]['family'] = True
        # Only add to counts if this is a new occurrence (not already counted in division)
        if not taxon_occurrences[name]['division']:
            taxon_occurrences[name]['occurrence_count'] += data['occurrence_count']
            taxon_occurrences[name]['member_size'] += data['member_size']
            taxon_occurrences[name]['taxid'] = data['taxid']

    # Process genus taxa
    for name, data in genus_taxa.items():
        taxon_occurrences[name]['genus'] = True
        # Only add to counts if this is a new occurrence (not already counted in division or family)
        if not taxon_occurrences[name]['division'] and not taxon_occurrences[name]['family']:
            taxon_occurrences[name]['occurrence_count'] += data['occurrence_count']
            taxon_occurrences[name]['member_size'] += data['member_size']
            taxon_occurrences[name]['taxid'] = data['taxid']

    # Get ranks for all taxids
    print("Getting taxonomic ranks using taxonkit...")
    all_taxids = []
    taxid_to_name = {}

    for name, data in taxon_occurrences.items():
        taxid = data['taxid']
        if taxid != 'NA':
            all_taxids.append(taxid)
            taxid_to_name[taxid] = name

    # Get ranks for all taxids at once
    taxid_to_rank = get_ranks_for_taxids(all_taxids)

    # Create a DataFrame from the merged data
    merged_data = []
    for name, data in taxon_occurrences.items():
        taxid = data['taxid']
        rank = taxid_to_rank.get(str(taxid), "NA") if taxid != 'NA' else "NA"

        merged_data.append({
            'taxon_name': name,
            'taxid': taxid,
            'occurrence_count': data['occurrence_count'],
            'rank': rank,
            'division': data['division'],
            'family': data['family'],
            'genus': data['genus']
        })

    # Convert to DataFrame and sort by occurrence count
    merged_df = pd.DataFrame(merged_data)
    merged_df = merged_df.sort_values(by=['occurrence_count'], ascending=[False])

    # Save the complete merged data to a CSV file
    complete_df = merged_df[['taxon_name', 'taxid', 'occurrence_count', 'rank']].copy()
    complete_df.to_csv('taxa_with_ranks.csv', index=False)
    print(f"Saved complete taxa report with ranks to taxa_with_ranks.csv with {len(complete_df)} entries")

    # Create separate reports for each taxonomic level with ranks
    # Division matches
    division_matches = merged_df[merged_df['division'] == True][['taxon_name', 'taxid', 'occurrence_count', 'rank']].copy()
    division_matches.to_csv('division_matches.csv', index=False)
    print(f"Saved division matches to division_matches.csv with {len(division_matches)} entries")

    # Family matches
    family_matches = merged_df[merged_df['family'] == True][['taxon_name', 'taxid', 'occurrence_count', 'rank']].copy()
    family_matches.to_csv('family_matches.csv', index=False)
    print(f"Saved family matches to family_matches.csv with {len(family_matches)} entries")

    # Genus matches
    genus_matches = merged_df[merged_df['genus'] == True][['taxon_name', 'taxid', 'occurrence_count', 'rank']].copy()
    genus_matches.to_csv('genus_matches.csv', index=False)
    print(f"Saved genus matches to genus_matches.csv with {len(genus_matches)} entries")

    # Generate summary statistics
    print("\nSummary Statistics:")
    print(f"Total unique taxa: {len(merged_df)}")
    print(f"Taxa present in division: {len(division_matches)}")
    print(f"Taxa present in family: {len(family_matches)}")
    print(f"Taxa present in genus: {len(genus_matches)}")

    # Show rank distribution
    rank_counts = merged_df['rank'].value_counts()
    print(f"\nRank distribution:")
    for rank, count in rank_counts.items():
        print(f"  {rank}: {count}")

    # Save failed names log
    if failed_names_log:
        failed_log_file = "failed_taxon_names.log"
        with open(failed_log_file, 'w') as f:
            f.write("# Failed Taxon Name Lookups\n")
            f.write(f"# Total failed: {len(failed_names_log)}\n\n")
            for i, failure in enumerate(failed_names_log, 1):
                f.write(f"{i}. {failure}\n")
        print(f"\nSaved {len(failed_names_log)} failed taxon name lookups to {failed_log_file}")
    else:
        print("\nNo failed taxon name lookups to report!")

    # Create a visualization of the overlap
    try:
        create_overlap_visualization(merged_df)
    except ImportError:
        print("Warning: Visualization libraries (matplotlib, seaborn) not available. Skipping visualizations.")

    return merged_df

def create_overlap_visualization(merged_df):
    """
    Create a simple bar chart of the top taxa by occurrence count.

    Args:
        merged_df: DataFrame containing the merged taxa data
    """
    try:
        import matplotlib.pyplot as plt

        # Get the top 20 taxa by occurrence count
        top_taxa = merged_df.nlargest(20, 'occurrence_count')

        # Create a simple bar chart
        plt.figure(figsize=(12, 8))
        plt.bar(top_taxa['taxon_name'], top_taxa['occurrence_count'], color='skyblue')
        plt.title('Top 20 Taxa by Occurrence Count')
        plt.xlabel('Taxon Name')
        plt.ylabel('Occurrence Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('top_taxa_occurrence_counts.png')
        print("Created visualization: top_taxa_occurrence_counts.png")

    except Exception as e:
        print(f"Error creating visualization: {e}")

if __name__ == "__main__":
    analyze_taxon_matches()
