#!/usr/bin/env python3
"""
Extract Name_to_Use Column with Previous Names

This script extracts all values from the Name_to_Use column in the
Eukprot_included_datasets.txt file and saves them to a CSV file with
a column named "Name_to_use". It also replaces underscores with spaces
in the organism names for better readability.

Additionally, it adds a "Previous_Names" column for specific organisms
that have alternative names that might be useful for taxonomic mapping.

Usage:
    python extract_name_to_use.py

Output:
    A CSV file named 'eukprot_names_to_use.csv' with columns 'Name_to_use'
    and 'Previous_Names' containing the organism names with underscores
    replaced by spaces.
"""

import pandas as pd
import os
from pathlib import Path

def main():
    """Extract and save the Name_to_Use column."""
    # Get the current directory
    script_dir = Path(__file__).resolve().parent

    # Input file path
    input_file = script_dir / "Eukprot_included_datasets.txt"

    # Output file path
    output_file = script_dir / "eukprot_names_to_use.csv"

    try:
        # Load the data
        print(f"Loading data from {input_file}...")
        df = pd.read_csv(input_file, sep='\t')

        # Extract the Name_to_Use column, replace underscores with spaces, and create a new DataFrame
        names = df["Name_to_Use"].dropna()
        # Replace underscores with spaces for better readability
        print(f"Replacing underscores with spaces in {len(names)} organism names...")
        names = names.str.replace('_', ' ')

        # Create DataFrame with Name_to_use column
        names_df = pd.DataFrame({"Name_to_use": names})

        # Add Previous_Names column (initialized as empty)
        names_df["Previous_Names"] = ""

        # Define mapping of specific organisms to their previous names
        previous_names_map = {
            "Manchomonas bermudensis": "Acanthamoeba bermudensis",
            "Chromosphaera perkinsii": "Mayorella perkinsii",
            "Luapelamoeba hula": "Korotnevella hula",
            "Sappina pedata": "Thecamoeba pedata",
            "Chlorochromonas danica": "Heterochromonas danica, Ochromonas danica",  # Added Ochromonas danica as alternative
            "Eupelagonema oceanica": "Mayorella oceanica",
            "Sapocribrum chincoteaguense": "Sapocribrum chincoteaguensis",
            "Vannellida sp DIVA3-517-6-12": "Vannella sp DIVA3-517-6-12"
        }

        # Update Previous_Names for specific organisms
        print("Adding previous names for specific organisms...")
        for i, row in names_df.iterrows():
            name = row["Name_to_use"]
            if name in previous_names_map:
                names_df.at[i, "Previous_Names"] = previous_names_map[name]
                print(f"  ✅ Added previous name for {name}: {previous_names_map[name]}")

        # Save to CSV file
        names_df.to_csv(output_file, index=False)

        print(f"✅ Successfully extracted {len(names_df)} names to {output_file}")
        print(f"✅ All underscores have been replaced with spaces in the output file")

        # Count how many entries have previous names
        previous_names_count = (names_df["Previous_Names"] != "").sum()
        print(f"✅ Added previous names for {previous_names_count} organisms")

        # Print the first 10 names as a preview
        print("\nPreview of extracted names:")
        for i, (_, row) in enumerate(names_df.head(10).iterrows()):
            name = row["Name_to_use"]
            prev_name = row["Previous_Names"]
            if prev_name:
                print(f"{i+1}. {name} (Previous name: {prev_name})")
            else:
                print(f"{i+1}. {name}")

        if len(names_df) > 10:
            print(f"... and {len(names_df) - 10} more")

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
