import pandas as pd
import argparse

# Set up CLI args
parser = argparse.ArgumentParser(description="Extract Taxogroup2_UniEuk from a tab-delimited .txt metadata file.")
parser.add_argument("input_file", help="Path to the metadata .txt file (tab-delimited).")
parser.add_argument("output_file", help="Path to save unique phylum names (one per line).")
args = parser.parse_args()

# Load tab-delimited text file
df = pd.read_csv(args.input_file, sep="\t")

# Ensure the expected column exists
if 'Taxogroup2_UniEuk' not in df.columns:
    raise ValueError("Column 'Taxogroup2_UniEuk' not found in the input file.")

# Drop duplicates and NA, sort, and save
unique_phyla = df['Taxogroup2_UniEuk'].dropna().drop_duplicates().sort_values()
unique_phyla.to_csv(args.output_file, index=False, header=False)

print(f"✅ Unique phylum names saved to: {args.output_file}")

