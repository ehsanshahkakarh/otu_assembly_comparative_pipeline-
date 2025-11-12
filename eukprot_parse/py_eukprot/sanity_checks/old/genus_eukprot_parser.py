import pandas as pd
import os
from pathlib import Path

# Paths
script_dir = Path(__file__).resolve().parent
csv_output_dir = script_dir.parent / "csv_eukprot"
input_file = script_dir / "Eukprot_included_datasets.txt"
output_file = csv_output_dir / "eukprot_genus_counts.csv"
detailed_output_file = csv_output_dir / "eukprot_genus_with_accessions.csv"

# Ensure output directory exists
csv_output_dir.mkdir(exist_ok=True)

def process_eukprot_data():
    """Extract genus information from EukProt metadata and summarize counts."""
    try:
        # Load data
        df = pd.read_csv(input_file, sep='\t')
        print(f"✅ Loaded {len(df)} rows from EukProt metadata.")

        if "Genus_UniEuk" not in df.columns:
            raise ValueError("Missing required column 'Genus_UniEuk' in the dataset.")

        # Rename for consistency
        df = df.rename(columns={"Genus_UniEuk": "genus"})

        # Add domain and accession fields
        df["domain"] = "Eukaryota"
        df["accession_clean"] = df["EukProt_ID"]

        # Drop missing genus
        df = df.dropna(subset=["genus"])
        print(f"✅ Rows with valid genus: {len(df)}")

        # Save detailed file
        df[["accession_clean", "genus", "domain"]].to_csv(detailed_output_file, index=False)
        print(f"✅ Saved detailed accession-genus file to: {detailed_output_file}")

        # Group and count unique genomes per genus
        counts = (
            df.groupby(["genus", "domain"])
              .agg(
                  eukprot_genome_count=('accession_clean', 'nunique'),
                  accessions=('accession_clean', lambda x: list(sorted(x.unique())))
              )
              .reset_index()
              .sort_values("eukprot_genome_count", ascending=False)
        )

        # Save summary
        counts.to_csv(output_file, index=False)
        print(f"✅ Saved genus summary to: {output_file}")

        # Preview
        print("\n🔝 Top 10 genera:")
        print(counts[["genus", "eukprot_genome_count"]].head(10))

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    process_eukprot_data()
