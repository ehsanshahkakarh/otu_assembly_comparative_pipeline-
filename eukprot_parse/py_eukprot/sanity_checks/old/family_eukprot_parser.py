import pandas as pd
import os
from pathlib import Path

# Set up paths
script_dir = Path(__file__).resolve().parent
csv_output_dir = script_dir.parent / "csv_eukprot"
input_file = script_dir / "Eukprot_included_datasets.txt"
output_file = csv_output_dir / "eukprot_family_counts.csv"
detailed_output_file = csv_output_dir / "eukprot_family_with_accessions.csv"

# Ensure output directory exists
csv_output_dir.mkdir(exist_ok=True)

def process_eukprot_data():
    """Process the EukProt dataset to extract accurate family-level genome counts."""
    try:
        # Load EukProt dataset
        df = pd.read_csv(input_file, sep='\t')
        print(f"✅ Loaded {len(df)} rows from EukProt metadata")

        # Parse family from Taxonomy_UniEuk (index 4 if present)
        df["taxonomy_parts"] = df["Taxonomy_UniEuk"].str.split(';')
        df["family"] = df["taxonomy_parts"].apply(
            lambda x: x[4].replace("'", "").strip() if isinstance(x, list) and len(x) > 4 else None
        )

        # Clean up
        df["domain"] = "Eukaryota"
        df["accession_clean"] = df["EukProt_ID"]
        df = df.dropna(subset=["family"])

        # Save detailed file
        df[["accession_clean", "family", "domain"]].to_csv(detailed_output_file, index=False)
        print(f"✅ Saved detailed accession-family file with {len(df)} rows")

        # Group to get accurate genome counts per family
        counts = (
            df.groupby(["family", "domain"])
              .agg(
                  eukprot_genome_count=("accession_clean", "nunique"),
                  accessions=("accession_clean", lambda x: list(sorted(x.unique())))
              )
              .reset_index()
              .sort_values("eukprot_genome_count", ascending=False)
        )

        # Save summary table
        counts.to_csv(output_file, index=False)
        print(f"✅ Saved family summary to: {output_file}")

        # Preview
        print("\n🔝 Top 10 families:")
        print(counts[["family", "eukprot_genome_count"]].head(10))

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    process_eukprot_data()
