#!/usr/bin/env python3
"""
Verify taxonomy classifications using taxonKit with GTDB taxonomy.
This script checks if taxa in CSV files match their expected domain and rank
using the GTDB taxdump.
"""

import subprocess
import csv
from pathlib import Path
import sys
import os
import re
import argparse
import shutil

def verify_taxa(input_file, expected_domain, expected_rank, taxdump_dir=None, output_dir=None):
    """
    Verify taxonomy classifications using taxonKit with GTDB taxonomy.

    Args:
        input_file (str): Path to input CSV file
        expected_domain (str): Expected domain (e.g., 'Bacteria', 'Archaea')
        expected_rank (str): Expected taxonomic rank (e.g., 'phylum', 'family', 'genus')
        taxdump_dir (str, optional): Path to GTDB taxdump directory. If None, uses default.
        output_dir (str, optional): Directory for output log files. Defaults to same as input.

    Returns:
        tuple: (total_taxa, mismatches_count)
    """
    input_path = Path(input_file).resolve()
    assert input_path.exists(), f"File not found: {input_path}"

    # Set output directory to input directory if not specified
    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir).resolve()
        os.makedirs(output_dir, exist_ok=True)

    # Create output log filename based on input filename
    base_name = input_path.stem
    output_log = output_dir / f"{base_name}_verification.log"

    print(f"🔍 Verifying {input_path.name} for domain={expected_domain}, rank={expected_rank}")

    # Read taxa from CSV file (assuming format with header and taxon_name column)
    taxa = []
    try:
        with open(input_path, 'r') as f:
            reader = csv.DictReader(f)
            if 'taxon_name' not in reader.fieldnames:
                print(f"❌ Error: CSV file {input_path.name} does not have 'taxon_name' column")
                return 0, 0
            taxa = [row['taxon_name'] for row in reader if row['taxon_name'].strip()]
    except Exception as e:
        print(f"❌ Error reading CSV file: {e}")
        return 0, 0

    if not taxa:
        print(f"⚠️ No taxa found in {input_path.name}")
        return 0, 0

    print(f"📊 Found {len(taxa)} taxa to verify")

    mismatches = []
    not_found = []

    # Process taxa in batches to avoid command line length limits
    batch_size = 100
    for i in range(0, len(taxa), batch_size):
        batch = taxa[i:i+batch_size]

        # Create temporary file with batch taxa
        temp_file = output_dir / f"temp_taxa_batch_{i}.txt"
        with open(temp_file, 'w') as f:
            f.write('\n'.join(batch))

        try:
            # Set up taxonkit command with GTDB taxdump if provided
            taxonkit_cmd = ["taxonkit"]
            if taxdump_dir:
                # Ensure the taxdump directory exists
                taxdump_path = Path(taxdump_dir).resolve()
                if not taxdump_path.exists():
                    print(f"❌ GTDB taxdump directory not found: {taxdump_path}")
                    for taxon in batch:
                        mismatches.append((taxon, "❌ GTDB taxdump directory not found"))
                    continue

                # Check if required files exist
                required_files = ["nodes.dmp", "names.dmp"]
                missing_files = [f for f in required_files if not (taxdump_path / f).exists()]
                if missing_files:
                    print(f"❌ Missing required files in taxdump: {', '.join(missing_files)}")
                    for taxon in batch:
                        mismatches.append((taxon, f"❌ Missing taxdump files: {', '.join(missing_files)}"))
                    continue

                # Set TAXONKIT_DB environment variable
                os.environ["TAXONKIT_DB"] = str(taxdump_path)
                print(f"🔧 Using GTDB taxdump at: {taxdump_path}")

            # Run taxonkit lineage
            lineage_result = subprocess.run(
                taxonkit_cmd + ["lineage", str(temp_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ
            )

            if lineage_result.returncode != 0:
                print(f"❌ taxonkit lineage error: {lineage_result.stderr}")
                for taxon in batch:
                    mismatches.append((taxon, "❌ taxonkit lineage error"))
                continue

            # Run taxonkit reformat to get ranks
            reformat_result = subprocess.run(
                taxonkit_cmd + ["reformat", "-f", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}"],
                input=lineage_result.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ
            )

            if reformat_result.returncode != 0:
                print(f"❌ taxonkit reformat error: {reformat_result.stderr}")
                for taxon in batch:
                    mismatches.append((taxon, "❌ taxonkit reformat error"))
                continue

            # Process results
            lineage_lines = lineage_result.stdout.strip().split('\n')
            reformat_lines = reformat_result.stdout.strip().split('\n')

            for line_idx, (lineage_line, reformat_line) in enumerate(zip(lineage_lines, reformat_lines)):
                lineage_parts = lineage_line.strip().split('\t')
                reformat_parts = reformat_line.strip().split('\t')

                if len(lineage_parts) < 2 or not lineage_parts[1]:
                    taxon = lineage_parts[0] if lineage_parts else batch[line_idx]
                    not_found.append(taxon)
                    mismatches.append((taxon, "❌ Not found in NCBI taxonomy"))
                    continue

                taxon = lineage_parts[0]
                lineage = lineage_parts[1] if len(lineage_parts) > 1 else ""
                ranks = '\t'.join(reformat_parts[1:]) if len(reformat_parts) > 1 else ""

                issues = []

                # Check domain
                if expected_domain.lower() == "bacteria":
                    if not re.search(r'\bbacteria\b', lineage.lower()):
                        issues.append(f"⚠️ Domain mismatch: Not in Bacteria")
                elif expected_domain.lower() == "archaea":
                    if not re.search(r'\barchaea\b', lineage.lower()):
                        issues.append(f"⚠️ Domain mismatch: Not in Archaea")

                # Check for Candidatus status (which should be excluded)
                if "candidatus" in taxon.lower():
                    issues.append(f"⚠️ Candidatus taxon should be excluded")

                # Check rank based on filename pattern
                rank_index = {"phylum": 1, "family": 4, "genus": 5}
                if expected_rank.lower() in rank_index:
                    idx = rank_index[expected_rank.lower()]
                    if len(reformat_parts) <= idx or not reformat_parts[idx]:
                        issues.append(f"⚠️ Rank missing: No {expected_rank} classification")

                if issues:
                    mismatches.append((taxon, "; ".join(issues)))

        except Exception as e:
            print(f"❌ Error processing batch: {e}")
            for taxon in batch:
                mismatches.append((taxon, f"❌ Error: {e}"))

        # Clean up temp file
        if temp_file.exists():
            temp_file.unlink()

    # Write results to log file
    with open(output_log, "w") as log:
        log.write(f"Taxon\tIssue\n")
        for taxon, issue in mismatches:
            log.write(f"{taxon}\t{issue}\n")

    if not_found:
        print(f"⚠️ {len(not_found)} taxa not found in NCBI taxonomy")

    if mismatches:
        print(f"⚠️ {len(mismatches)} mismatches found. See {output_log}")
    else:
        print(f"✅ All {len(taxa)} taxa passed verification")

    return len(taxa), len(mismatches)

def find_gtdb_taxdump():
    """Find the GTDB taxdump directory in the repository."""
    # Start with known locations
    potential_paths = [
        # GTDB taxdump
        Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/gtdb_parse/taxdump_gtdp/gtdb-taxdump-R226"),
        # Relative path from script location
        Path(__file__).resolve().parents[3] / "taxdump_gtdp" / "gtdb-taxdump-R226"
    ]

    for path in potential_paths:
        if path.exists() and (path / "nodes.dmp").exists() and (path / "names.dmp").exists():
            return path

    # If not found in known locations, search in repository
    repo_root = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes")
    if repo_root.exists():
        for path in repo_root.glob("**/gtdb-taxdump-*/"):
            if (path / "nodes.dmp").exists() and (path / "names.dmp").exists():
                return path

    return None

def main():
    """Process all CSV files in the directory."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Verify taxonomy classifications using taxonKit with GTDB taxonomy.")
    parser.add_argument("--taxdump", help="Path to GTDB taxdump directory")
    parser.add_argument("--output-dir", help="Directory for output files")
    args = parser.parse_args()

    # Get current directory
    current_dir = Path(__file__).parent.resolve()

    # Find GTDB taxdump directory
    taxdump_dir = None
    if args.taxdump:
        taxdump_dir = Path(args.taxdump).resolve()
    else:
        taxdump_dir = find_gtdb_taxdump()

    if taxdump_dir and taxdump_dir.exists():
        print(f"🔧 Using GTDB taxdump at: {taxdump_dir}")
    else:
        print("⚠️ GTDB taxdump directory not found. Will use default NCBI taxonomy.")
        taxdump_dir = None

    # Set output directory
    output_dir = current_dir
    if args.output_dir:
        output_dir = Path(args.output_dir).resolve()
        os.makedirs(output_dir, exist_ok=True)

    # Define file patterns and their expected domains and ranks
    file_patterns = [
        {"pattern": "archaeal_phyla.csv", "domain": "Archaea", "rank": "phylum"},
        {"pattern": "archaeal_families.csv", "domain": "Archaea", "rank": "family"},
        {"pattern": "archaeal_genera.csv", "domain": "Archaea", "rank": "genus"},
        {"pattern": "bacterial_phyla.csv", "domain": "Bacteria", "rank": "phylum"},
        {"pattern": "bacterial_families.csv", "domain": "Bacteria", "rank": "family"},
        {"pattern": "bacterial_genera.csv", "domain": "Bacteria", "rank": "genus"}
    ]

    # Create a summary file
    summary_file = output_dir / "verification_summary.txt"
    with open(summary_file, "w") as summary:
        summary.write("File\tTotal Taxa\tMismatches\tPercentage Correct\n")

    # Process each file
    for pattern in file_patterns:
        file_path = current_dir / pattern["pattern"]
        if file_path.exists():
            total, mismatches = verify_taxa(
                file_path,
                pattern["domain"],
                pattern["rank"],
                taxdump_dir,
                output_dir
            )

            # Update summary
            if total > 0:
                percentage_correct = ((total - mismatches) / total) * 100
            else:
                percentage_correct = 0

            with open(summary_file, "a") as summary:
                summary.write(f"{pattern['pattern']}\t{total}\t{mismatches}\t{percentage_correct:.2f}%\n")
        else:
            print(f"⚠️ File not found: {file_path}")

    print(f"✅ Verification complete. Summary written to {summary_file}")

if __name__ == "__main__":
    main()
