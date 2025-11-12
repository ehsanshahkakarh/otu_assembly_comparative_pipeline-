import pandas as pd
import argparse

# Set up CLI args
parser = argparse.ArgumentParser(description="Extract taxon_name column from metadata file.")
parser.add_argument("input_file", help="Path to the metadata file (CSV or TSV).")
parser.add_argument("output_file", help="Path to save extracted taxon names.")
parser.add_argument("--sep", default="\t", help="Column separator (default: tab). Use ',' for CSV.")
args = parser.parse_args()

# Load file
df = pd.read_csv(args.input_file, sep=args.sep)

# Check if column exists
if 'family' not in df.columns:
    raise ValueError("Column 'genus' not found in the input file.")


# Extract and save both taxon_name and data_source to the same file
df[['family']].drop_duplicates().to_csv(args.output_file, index=False, header=True)
print(f"✅ Extracted genus columns saved to: {args.output_file}")

# If you still want separate files for each column, uncomment these lines:
# df[['taxon_name']].drop_duplicates().to_csv(args.output_file.replace(".txt", "_taxon_only.txt"), index=False, header=True)
# print(f"✅ Extracted taxon_name column saved to: {args.output_file.replace('.txt', '_taxon_only.txt')}")
# df[['data_source']].drop_duplicates().to_csv(args.output_file.replace(".txt", "_source_only.txt"), index=False, header=True)
# print(f"✅ Extracted data_source column saved to: {args.output_file.replace('.txt', '_source_only.txt')}")
