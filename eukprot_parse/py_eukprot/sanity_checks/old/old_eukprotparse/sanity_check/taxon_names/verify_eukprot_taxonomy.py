#!/usr/bin/env python3
"""
Verify EukProt taxonomy classifications using taxonKit.
This script checks if taxa in CSV files match their expected domain and rank.
"""

import subprocess
import csv
from pathlib import Path
import sys
import os
import re
import argparse
import logging
from datetime import datetime
import difflib  # For fuzzy matching

def setup_logging(log_dir):
    """Set up logging configuration"""
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / f"eukprot_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file

def verify_taxa(input_file, expected_domain, expected_rank, taxdump_dir=None, output_dir=None):
    """
    Verify taxonomy classifications using taxonKit.

    Args:
        input_file (str): Path to input CSV file
        expected_domain (str): Expected domain (e.g., 'Eukaryota')
        expected_rank (str): Expected taxonomic rank (e.g., 'phylum', 'family', 'genus')
        taxdump_dir (str, optional): Path to NCBI taxdump directory. If None, uses default.
        output_dir (str, optional): Directory for output log files. Defaults to same as input.

    Returns:
        tuple: (total_taxa, mismatches_count)
    """
    input_path = Path(input_file).resolve()
    if not input_path.exists():
        logging.error(f"File not found: {input_path}")
        return 0, 0

    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir).resolve()
        os.makedirs(output_dir, exist_ok=True)

    base_name = input_path.stem
    output_log = output_dir / f"{base_name}_verification.log"

    logging.info(f"Verifying {input_path.name} for domain={expected_domain}, rank={expected_rank}")

    # Read taxa from CSV file
    taxa = []
    taxa_with_taxids = {}
    try:
        with open(input_path, 'r') as f:
            reader = csv.DictReader(f)
            if 'taxon_name' not in reader.fieldnames:
                logging.error(f"CSV file {input_path.name} does not have 'taxon_name' column")
                return 0, 0

            # Check if taxid column exists
            has_taxid = 'taxid' in reader.fieldnames

            if has_taxid:
                logging.info(f"Found taxid column in {input_path.name}, will use taxids when available")
                for row in reader:
                    if row['taxon_name'].strip():
                        taxa.append(row['taxon_name'])
                        if row['taxid'] and row['taxid'].strip():
                            taxa_with_taxids[row['taxon_name']] = row['taxid']
            else:
                taxa = [row['taxon_name'] for row in reader if row['taxon_name'].strip()]

            logging.info(f"Found {len(taxa_with_taxids)} taxa with taxids out of {len(taxa)} total taxa")
    except Exception as e:
        logging.error(f"Error reading CSV file: {e}")
        return 0, 0

    if not taxa:
        logging.warning(f"No taxa found in {input_path.name}")
        return 0, 0

    logging.info(f"Found {len(taxa)} taxa to verify")

    mismatches = []
    not_found = []

    batch_size = 100
    for i in range(0, len(taxa), batch_size):
        batch = taxa[i:i+batch_size]

        # Prepare input - use taxids when available
        batch_with_taxids = []
        for taxon in batch:
            if taxon in taxa_with_taxids:
                # Use taxid directly if available
                batch_with_taxids.append(taxa_with_taxids[taxon])
            else:
                batch_with_taxids.append(taxon)

        # Create input text - try both taxon names and taxids
        input_text = "\n".join(batch)
        input_text_with_taxids = "\n".join(batch_with_taxids)

        env = os.environ.copy()
        if taxdump_dir:
            env["TAXONKIT_DB"] = str(taxdump_dir)
            logging.info(f"Using taxdump directory: {taxdump_dir}")

        try:
            # Log the first few taxa we're checking
            if i == 0:
                sample_taxa = batch[:5]
                sample_taxids = batch_with_taxids[:5]
                logging.info(f"Sample taxa being checked: {sample_taxa}")
                logging.info(f"Sample taxids being checked: {sample_taxids}")

            # Run taxonkit lineage with verbose output
            logging.info(f"Running taxonkit lineage on batch {i//batch_size + 1} ({len(batch)} taxa)")

            # First try with taxids if available
            if any(taxon in taxa_with_taxids for taxon in batch):
                logging.info(f"Trying with taxids first for batch {i//batch_size + 1}")
                try:
                    lineage_result = subprocess.run(
                        ["taxonkit", "lineage"],
                        input=input_text_with_taxids,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if i == 0:
                        # Log sample results from lineage with taxids
                        sample_lines = lineage_result.stdout.strip().split('\n')[:5]
                        logging.info(f"Sample lineage results with taxids: {sample_lines}")

                    # If successful, use this result
                    if lineage_result.returncode == 0 and lineage_result.stdout.strip():
                        logging.info(f"Successfully used taxids for batch {i//batch_size + 1}")
                        # Continue with this result
                    else:
                        # Fall back to using names
                        logging.warning(f"Taxid lookup failed, falling back to names for batch {i//batch_size + 1}")
                        lineage_result = None
                except Exception as e:
                    logging.error(f"Error running lineage with taxids: {e}")
                    lineage_result = None
            else:
                lineage_result = None

            # If taxid approach didn't work, try with name2taxid
            if lineage_result is None:
                try:
                    name2taxid_result = subprocess.run(
                        ["taxonkit", "name2taxid"],
                        input=input_text,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env
                    )

                    if i == 0:
                        # Log sample results from name2taxid
                        sample_lines = name2taxid_result.stdout.strip().split('\n')[:5]
                        logging.info(f"Sample name2taxid results: {sample_lines}")

                        if name2taxid_result.stderr:
                            logging.warning(f"name2taxid warnings: {name2taxid_result.stderr}")
                except Exception as e:
                    logging.error(f"Error running name2taxid: {e}")

                # Now run lineage with names
                lineage_result = subprocess.run(
                    ["taxonkit", "lineage"],
                    input=input_text,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                    check=True
                )

            if i == 0:
                # Log sample results from lineage
                sample_lines = lineage_result.stdout.strip().split('\n')[:5]
                logging.info(f"Sample lineage results: {sample_lines}")

            if lineage_result.stderr:
                logging.warning(f"taxonkit lineage warnings: {lineage_result.stderr}")

            # Run taxonkit reformat
            reformat_result = subprocess.run(
                ["taxonkit", "reformat", "-f", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}"],
                input=lineage_result.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                check=True
            )
            if reformat_result.stderr:
                logging.error(f"taxonkit reformat error: {reformat_result.stderr}")
                for taxon in batch:
                    mismatches.append((taxon, "taxonkit reformat error"))
                continue

            lineage_lines = lineage_result.stdout.strip().split('\n')
            reformat_lines = reformat_result.stdout.strip().split('\n')

            for line_idx, (lineage_line, reformat_line) in enumerate(zip(lineage_lines, reformat_lines)):
                lineage_parts = lineage_line.strip().split('\t')
                reformat_parts = reformat_line.strip().split('\t')

                if len(lineage_parts) < 2 or not lineage_parts[1]:
                    taxon = lineage_parts[0] if lineage_parts else batch[line_idx]
                    not_found.append(taxon)
                    mismatches.append((taxon, "Not found in NCBI taxonomy"))
                    continue

                taxon = lineage_parts[0]
                lineage = lineage_parts[1]

                issues = []

                if expected_domain.lower() == "eukaryota":
                    if not re.search(r'\beukaryota\b', lineage.lower()):
                        issues.append("Domain mismatch: Not in Eukaryota")

                if "candidatus" in taxon.lower():
                    issues.append("Candidatus taxon should be excluded")

                rank_index = {"phylum": 1, "family": 4, "genus": 5}
                if expected_rank.lower() in rank_index:
                    idx = rank_index[expected_rank.lower()]
                    if len(reformat_parts) <= idx or not reformat_parts[idx]:
                        issues.append(f"Rank missing: No {expected_rank} classification")

                if issues:
                    mismatches.append((taxon, "; ".join(issues)))

        except subprocess.CalledProcessError as e:
            logging.error(f"TaxonKit error: {e}")
            for taxon in batch:
                mismatches.append((taxon, f"TaxonKit error: {e}"))
        except Exception as e:
            logging.error(f"Unexpected error: {e}")
            for taxon in batch:
                mismatches.append((taxon, f"Unexpected error: {e}"))

    # Write results
    with open(output_log, "w") as log:
        log.write("Taxon\tIssue\n")
        for taxon, issue in mismatches:
            log.write(f"{taxon}\t{issue}\n")

    if not_found:
        logging.warning(f"{len(not_found)} taxa not found in NCBI taxonomy")
    if mismatches:
        logging.warning(f"{len(mismatches)} mismatches found. See {output_log}")
    else:
        logging.info(f"All {len(taxa)} taxa passed verification")

    return len(taxa), len(mismatches)

def find_ncbi_taxdump():
    """Find the NCBI taxdump directory in the repository."""
    # First check the default taxonkit directory
    home_dir = Path.home()
    taxonkit_dir = home_dir /".taxonkit"

    if taxonkit_dir.exists() and (taxonkit_dir / "nodes.dmp").exists() and (taxonkit_dir / "names.dmp").exists():
        return taxonkit_dir

    # Check other known locations
    potential_paths = [
        Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi"),
        Path(__file__).resolve().parents[3] / "ncbi_parse_scripts" / "taxonomic_mapping" / "taxdump_ncbi"
    ]

    for path in potential_paths:
        if path.exists() and (path / "nodes.dmp").exists() and (path / "names.dmp").exists():
            return path

    # If not found in known locations, search in repository
    repo_root = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes")
    if repo_root.exists():
        for path in repo_root.glob("**/taxdump_ncbi/"):
            if (path / "nodes.dmp").exists() and (path / "names.dmp").exists():
                return path

    return None

def main():
    """Process all CSV files in the directory."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Verify EukProt taxonomy classifications using taxonKit.")
    parser.add_argument("--taxdump", help="Path to NCBI taxdump directory")
    parser.add_argument("--output-dir", help="Directory for output files")
    args = parser.parse_args()

    # Get current directory
    script_dir = Path(__file__).parent.resolve()
    log_dir = script_dir / "logs"
    log_file = setup_logging(log_dir)

    # Find NCBI taxdump directory
    taxdump_dir = None
    if args.taxdump:
        taxdump_dir = Path(args.taxdump).resolve()
    else:
        taxdump_dir = find_ncbi_taxdump()

    if taxdump_dir and taxdump_dir.exists():
        logging.info(f"Using NCBI taxdump at: {taxdump_dir}")
    else:
        logging.warning("NCBI taxdump directory not found. Will use default NCBI taxonomy.")
        taxdump_dir = None

    # Set output directory
    output_dir = script_dir
    if args.output_dir:
        output_dir = Path(args.output_dir).resolve()
        os.makedirs(output_dir, exist_ok=True)

    # Define file patterns and their expected domains and ranks
    file_patterns = [
        {"pattern": "eukaryotic_phyla.csv", "domain": "Eukaryota", "rank": "phylum"},
        {"pattern": "eukaryotic_families.csv", "domain": "Eukaryota", "rank": "family"},
        {"pattern": "eukaryotic_genera.csv", "domain": "Eukaryota", "rank": "genus"}
    ]

    # Create a summary file
    summary_file = output_dir / "verification_summary.txt"
    with open(summary_file, "w") as summary:
        summary.write("File\tTotal Taxa\tMismatches\tPercentage Correct\n")

    # Process each file
    for pattern in file_patterns:
        file_path = script_dir / pattern["pattern"]
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
            logging.warning(f"File not found: {file_path}")

    logging.info(f"Verification complete. Summary written to {summary_file}")
    logging.info(f"Log file: {log_file}")

if __name__ == "__main__":
    main()
