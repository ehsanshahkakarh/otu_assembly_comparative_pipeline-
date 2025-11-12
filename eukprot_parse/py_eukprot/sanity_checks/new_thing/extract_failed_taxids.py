#!/usr/bin/env python3
"""
Extract Failed Taxid Entries

This script reads the eukprot_names_with_taxids.csv file and extracts
only the entries that failed to get a taxid (marked as "FAILED").
It outputs these entries to a new CSV file for further analysis.

Usage:
    python extract_failed_taxids.py

Output:
    A CSV file named 'eukprot_failed_taxids.csv' containing only the
    entries that failed to get a taxid.
"""

import pandas as pd
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("extract_failed_taxids.log"),
        logging.StreamHandler()
    ]
)

def main():
    """Extract entries that failed to get a taxid."""
    # Get the current directory
    script_dir = Path(__file__).resolve().parent

    # Input file path (now in the current directory)
    input_file = script_dir / "eukprot_names_with_taxids.csv"

    # Output file path (also in the current directory)
    output_file = script_dir / "eukprot_failed_taxids.csv"

    try:
        # Check if input file exists
        if not input_file.exists():
            logging.error(f"❌ Input file not found: {input_file}")
            return

        # Load the data
        logging.info(f"Loading data from {input_file}...")
        df = pd.read_csv(input_file)

        # Get column names
        if len(df.columns) < 2:
            logging.error(f"❌ Input file does not have the expected format. Found columns: {df.columns}")
            return

        name_column = df.columns[0]  # First column should be the name
        taxid_column = df.columns[1]  # Second column should be the taxid

        # Filter for failed entries
        failed_df = df[df[taxid_column] == "FAILED"].copy()

        # Count failed entries
        total_entries = len(df)
        failed_entries = len(failed_df)
        success_entries = total_entries - failed_entries

        logging.info(f"Total entries: {total_entries}")
        logging.info(f"Successful entries: {success_entries} ({success_entries/total_entries*100:.1f}%)")
        logging.info(f"Failed entries: {failed_entries} ({failed_entries/total_entries*100:.1f}%)")

        # Save failed entries to CSV
        failed_df.to_csv(output_file, index=False)
        logging.info(f"✅ Saved {failed_entries} failed entries to {output_file}")

        # Print preview of failed entries
        if not failed_df.empty:
            logging.info("\nPreview of failed entries:")
            for i, (_, row) in enumerate(failed_df.head(10).iterrows()):
                logging.info(f"{i+1}. {row[name_column]}")

            if len(failed_df) > 10:
                logging.info(f"... and {len(failed_df) - 10} more")

            # Analyze failed entries
            analyze_failed_entries(failed_df, name_column)

    except Exception as e:
        logging.error(f"❌ Error: {e}")

def analyze_failed_entries(failed_df, name_column):
    """Analyze failed entries to identify potential patterns."""
    logging.info("\nAnalyzing failed entries...")

    # Check for patterns in names
    has_sp = 0
    has_underscore = 0
    has_special_chars = 0

    for name in failed_df[name_column]:
        if " sp " in name.lower() or name.lower().endswith(" sp"):
            has_sp += 1
        if "_" in name:
            has_underscore += 1
        if any(c in name for c in "()[]{}.,;:!?@#$%^&*+="):
            has_special_chars += 1

    logging.info(f"Entries containing 'sp' (species): {has_sp} ({has_sp/len(failed_df)*100:.1f}%)")
    logging.info(f"Entries containing underscores: {has_underscore} ({has_underscore/len(failed_df)*100:.1f}%)")
    logging.info(f"Entries containing special characters: {has_special_chars} ({has_special_chars/len(failed_df)*100:.1f}%)")

    # Suggest potential fixes
    logging.info("\nPotential fixes to try:")
    logging.info("1. For entries with 'sp', try removing the 'sp' part and search for the genus only")
    logging.info("2. For entries with underscores, try replacing them with spaces")
    logging.info("3. For entries with special characters, try removing them")
    logging.info("4. Try using fuzzy matching with a higher similarity threshold")
    logging.info("5. Consider manual curation for the most important entries")

if __name__ == "__main__":
    main()
