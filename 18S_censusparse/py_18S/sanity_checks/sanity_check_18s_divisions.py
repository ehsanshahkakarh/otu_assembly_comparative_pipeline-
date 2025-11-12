#!/usr/bin/env python3
"""
Extract 18S Division Names

This script parses the eukcensus_by_division.csv file and extracts only the division names.

Usage:
    python sanity_check_18s_divisions.py
"""

import csv
import pandas as pd

def main():
    # Input file
    input_file = "eukcensus_by_division.csv"

    try:
        # Read the CSV file using pandas
        df = pd.read_csv(input_file)

        # Extract division names
        division_names = df['taxon_name'].tolist()

        # Sort alphabetically
        division_names.sort()

        # Print division names
        print("18S Division Names:")
        for i, name in enumerate(division_names):
            print(f"{i+1}. {name}")

        # Write division names to a text file
        with open("18s_division_names.txt", "w") as f:
            for name in division_names:
                f.write(f"{name}\n")

        print(f"\nTotal divisions: {len(division_names)}")
        print("Division names written to 18s_division_names.txt")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
