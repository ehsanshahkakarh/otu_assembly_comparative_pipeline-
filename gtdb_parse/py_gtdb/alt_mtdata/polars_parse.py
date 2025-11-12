import polars as pl
import os

def check_file_exists(filepath):
    if not os.path.exists(filepath):
        print(f"Error: File not found: {filepath}")
        return False
    return True

# Check if input files exist
bac_file = "00bac120_taxonomy.tsv"
arc_file = "00ar53_taxonomy.tsv"

if not (check_file_exists(bac_file) and check_file_exists(arc_file)):
    print("Please check your input file paths")
    exit(1)

try:
    # Load GTDB taxonomy files
    print(f"Loading bacterial taxonomy from {bac_file}")
    df_bac = pl.read_csv(
        bac_file,
        separator="\t",
        has_header=False,
        new_columns=["accession", "taxonomy"]
    )
    print(f"Loaded {len(df_bac)} bacterial records")

    print(f"Loading archaeal taxonomy from {arc_file}")
    df_arc = pl.read_csv(
        arc_file,
        separator="\t",
        has_header=False,
        new_columns=["accession", "taxonomy"]
    )
    print(f"Loaded {len(df_arc)} archaeal records")

    # Combine bacteria and archaea
    print("Concatenating dataframes")
    df = pl.concat([df_bac, df_arc])
    print(f"Combined dataframe has {len(df)} records")

    # Clean accession and extract phylum
    print("Cleaning and processing data")
    df_cleaned = df.with_columns([
        pl.col("accession").str.replace(r"^[A-Z]{2}_", "").alias("accession"),
        pl.col("taxonomy").str.extract(r"p__([A-Za-z0-9_\-]+)", 1).alias("phylum")
    ]).drop_nulls(subset=["phylum"])
    print(f"After cleaning, {len(df_cleaned)} records remain")

    # Count genomes per phylum and sort by count (descending)
    print("Counting and sorting phyla")
    phylum_counts = (
        df_cleaned.group_by("phylum")
                 .agg(pl.count().alias("gtdb_genome_count"))
                 .sort("gtdb_genome_count", descending=True)
    )
    print(f"Found {len(phylum_counts)} unique phyla")

    # Save the output CSV
    output_file = "gtdb_phylum_counts_polars_fixed.csv"
    print(f"Saving results to {output_file}")
    phylum_counts.write_csv(output_file)

    # Print preview
    print("\nTop 10 phyla by genome count:")
    print(phylum_counts.head(10))

except Exception as e:
    print(f"An error occurred: {str(e)}")

