#!/usr/bin/env python3
"""
Clean 18S EukCensus Parser with Comprehensive Taxa Preservation

Simple, streamlined parser for 18S EukCensus data that:
1. Processes division, family, and genus levels
2. PRESERVES .U. entries and unidentified taxa for downstream visualization
3. Filters out only _X entries (technical artifacts)
4. Gets taxids and lineages using taxonkit
5. Outputs sorted CSV files by decreasing OTU count

UPDATED: Now preserves .U. entries and unidentified taxa for comprehensive visualization

Output files:
- eukcensus_18S_by_division.csv
- eukcensus_18S_by_family.csv
- eukcensus_18S_by_genus.csv
"""

import pandas as pd
import subprocess
import csv
import tempfile
import os
from pathlib import Path

def setup_paths():
    """Setup input/output paths"""
    script_dir = Path(__file__).resolve().parent
    censusparse_dir = script_dir.parent
    metadata_dir = censusparse_dir / "metadata"
    csv_output_dir = censusparse_dir / "csv_outputs"
    
    # Ensure output directory exists
    csv_output_dir.mkdir(exist_ok=True)
    
    return metadata_dir, csv_output_dir

def get_taxids_batch(taxon_names):
    """Get taxids for a list of taxon names using taxonkit with temporary files"""
    if not taxon_names:
        return {}

    # Clean all names and maintain mapping
    clean_names = []
    name_mapping = {}  # Maps clean_name back to original_name

    for name in taxon_names:
        clean_name = name.replace("_", " ")
        # Only split on dots if it looks like a taxonomic suffix (e.g., "sp.", "cf.")
        # Don't split family names like "Lipomycetaceae"
        if "." in clean_name and any(suffix in clean_name.lower() for suffix in [" sp.", " cf.", " aff.", " gen."]):
            clean_name = clean_name.split(".")[0]
        clean_names.append(clean_name)
        name_mapping[clean_name] = name

    print(f"    Querying taxonkit for {len(clean_names)} names...")

    # Create temporary files for input and output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as input_file:
        input_filename = input_file.name
        # Write clean names to temporary file
        for clean_name in clean_names:
            input_file.write(f"{clean_name}\n")

    output_filename = input_filename.replace('.txt', '_output.txt')

    try:
        # Run taxonkit with file input/output
        result = subprocess.run(
            ["taxonkit", "name2taxid", input_filename, "-o", output_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        taxid_map = {}

        # Read results from output file
        if result.returncode == 0 and os.path.exists(output_filename):
            with open(output_filename, 'r') as f:
                lines = f.readlines()

            for i, line in enumerate(lines):
                if i < len(clean_names) and line.strip():
                    parts = line.strip().split('\t')
                    clean_name = clean_names[i]
                    original_name = name_mapping[clean_name]

                    if len(parts) >= 2 and parts[1] != "0":
                        taxid_map[original_name] = parts[1]
                    else:
                        taxid_map[original_name] = "NA"
                        print(f"      No taxid found for: {original_name} (cleaned: {clean_name})")
        else:
            print(f"      taxonkit command failed with return code: {result.returncode}")
            if result.stderr:
                print(f"      Error: {result.stderr}")

        # Fill in any missing entries
        for name in taxon_names:
            if name not in taxid_map:
                taxid_map[name] = "NA"
                print(f"      Missing from results: {name}")

        successful_lookups = sum(1 for v in taxid_map.values() if v != "NA")
        print(f"    Successfully found taxids for {successful_lookups}/{len(taxon_names)} names")

        return taxid_map

    except Exception as e:
        print(f"    Error in taxonkit name2taxid: {e}")
        return {name: "NA" for name in taxon_names}
    finally:
        # Clean up temporary files
        try:
            if os.path.exists(input_filename):
                os.unlink(input_filename)
            if os.path.exists(output_filename):
                os.unlink(output_filename)
        except:
            pass

def get_lineages_batch(taxids):
    """Get lineages for a list of taxids using taxonkit with temporary files"""
    if not taxids:
        return {}

    # Filter out "NA" taxids
    valid_taxids = [tid for tid in taxids if tid != "NA"]
    if not valid_taxids:
        return {tid: ("", "", "") for tid in taxids}

    print(f"    Querying lineages for {len(valid_taxids)} valid taxids...")

    # Create temporary files for input and output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as input_file:
        input_filename = input_file.name
        # Write taxids to temporary file
        for taxid in valid_taxids:
            input_file.write(f"{taxid}\n")

    output_filename = input_filename.replace('.txt', '_lineage_output.txt')

    try:
        # Run taxonkit with file input/output
        result = subprocess.run(
            ["taxonkit", "lineage", "-R", "-t", input_filename, "-o", output_filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        lineage_map = {}

        # Read results from output file
        if result.returncode == 0 and os.path.exists(output_filename):
            with open(output_filename, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        taxid = parts[0]
                        lineage = parts[1]
                        lineage_taxids = parts[2]
                        lineage_ranks = parts[3]
                        lineage_map[taxid] = (lineage, lineage_ranks, lineage_taxids)
                    else:
                        print(f"      Incomplete lineage data for line: {line}")
        else:
            print(f"      taxonkit lineage command failed with return code: {result.returncode}")
            if result.stderr:
                print(f"      Error: {result.stderr}")

        # Fill in missing entries and NA taxids
        for tid in taxids:
            if tid not in lineage_map:
                lineage_map[tid] = ("", "", "")
                if tid != "NA":
                    print(f"      No lineage found for taxid: {tid}")

        successful_lineages = sum(1 for tid, (lineage, _, _) in lineage_map.items() if lineage and tid != "NA")
        print(f"    Successfully found lineages for {successful_lineages}/{len(valid_taxids)} taxids")

        return lineage_map

    except Exception as e:
        print(f"    Error in taxonkit lineage: {e}")
        return {tid: ("", "", "") for tid in taxids}
    finally:
        # Clean up temporary files
        try:
            if os.path.exists(input_filename):
                os.unlink(input_filename)
            if os.path.exists(output_filename):
                os.unlink(output_filename)
        except:
            pass

def process_taxonomic_level(df, level_name):
    """Process a single taxonomic level - UPDATED: preserves .U. entries for visualization"""
    print(f"Processing {level_name} level...")

    # Filter out null entries and _X entries only - preserve .U. entries for visualization
    level_df = df[~df[level_name].isna()]
    print(f"  Found {len(level_df):,} entries after removing null values")

    # REMOVED: .U. filtering - these entries are now preserved for downstream visualization
    # Filter out entries with _X (technical artifacts)
    level_df = level_df[~level_df[level_name].str.contains('_X', na=False)]
    print(f"  Found {len(level_df):,} entries after removing _X entries (preserved .U. entries)")

    # Group by taxon and aggregate
    grouped = level_df.groupby(level_name).agg({
        'centroid': 'count',  # OTU count
        'size': 'sum'         # Sequence count
    }).reset_index()

    # Sort by OTU count (descending)
    grouped = grouped.sort_values('centroid', ascending=False)

    print(f"  Found {len(grouped)} unique {level_name} taxa")

    return grouped

def add_taxonomy_info(grouped_df, level_name):
    """Add taxid and lineage information (vectorized)"""
    print(f"  Getting taxonomy info for {level_name}...")

    # Get all unique taxon names
    taxon_names = grouped_df[level_name].tolist()

    # Get taxids for all names at once
    print(f"    Getting taxids for {len(taxon_names)} taxa...")
    taxid_map = get_taxids_batch(taxon_names)

    # Get all unique taxids
    unique_taxids = list(set(taxid_map.values()))

    # Get lineages for all taxids at once
    print(f"    Getting lineages for {len(unique_taxids)} taxids...")
    lineage_map = get_lineages_batch(unique_taxids)

    # Build results
    results = []
    for _, row in grouped_df.iterrows():
        taxon_name = row[level_name]
        otu_count = row['centroid']
        size_count = row['size']

        # Get taxid and lineage from maps
        taxid = taxid_map[taxon_name]
        lineage, lineage_ranks, lineage_taxids = lineage_map[taxid]

        results.append({
            'Name_to_use': taxon_name,
            'taxid': taxid,
            'otu_count': otu_count,
            'size_count': size_count,
            'lineage': lineage,
            'lineage_ranks': lineage_ranks,
            'lineage_taxids': lineage_taxids
        })

    return results

def write_csv(results, output_file, total_otu_count, total_size_count):
    """Write results to CSV file with percentage columns"""
    print(f"  Writing {len(results)} entries to {output_file}")

    # Add percentage calculations to each result
    for result in results:
        # Calculate percentages
        otu_percentage = round((result['otu_count'] / total_otu_count * 100), 2) if total_otu_count > 0 else 0
        size_percentage = round((result['size_count'] / total_size_count * 100), 2) if total_size_count > 0 else 0

        # Add percentage columns
        result['otu_percentage'] = otu_percentage
        result['size_percentage'] = size_percentage

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'Name_to_use', 'taxid', 'otu_count', 'otu_percentage', 'size_count', 'size_percentage',
            'lineage', 'lineage_ranks', 'lineage_taxids'
        ])
        writer.writeheader()
        writer.writerows(results)

def main():
    """Main processing function"""
    print("🚀 Starting clean 18S EukCensus parsing...")
    
    # Setup paths
    metadata_dir, csv_output_dir = setup_paths()
    input_file = metadata_dir / "eukcensus_18S.clusters.97.tsv"
    
    # Check input file
    if not input_file.exists():
        print(f"❌ Error: Input file not found: {input_file}")
        return
    
    # Read data
    print(f"📊 Reading {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    print(f"✅ Loaded {len(df):,} records")

    # Calculate database totals for percentage calculations
    print("📊 Calculating database totals for percentage calculations...")
    total_otu_count = len(df)  # Each row represents one OTU
    total_size_count = df['size'].sum()  # Sum of all sequence counts

    print(f"📊 Total OTUs in database: {total_otu_count:,}")
    print(f"📊 Total size count in database: {total_size_count:,}")

    # Process each taxonomic level
    levels = ['division', 'family', 'genus']
    
    for level in levels:
        print(f"\n🔍 Processing {level}...")
        
        # Process the level
        grouped_df = process_taxonomic_level(df, level)
        
        # Add taxonomy information
        results = add_taxonomy_info(grouped_df, level)
        
        # Write output
        output_file = csv_output_dir / f"eukcensus_18S_by_{level}.csv"
        write_csv(results, output_file, total_otu_count, total_size_count)
        
        print(f"✅ Completed {level}: {len(results)} entries")
    
    print(f"\n🎉 Processing complete!")
    print(f"📁 Output files in: {csv_output_dir}")
    for level in levels:
        print(f"  • eukcensus_18S_by_{level}.csv")

if __name__ == "__main__":
    main()
