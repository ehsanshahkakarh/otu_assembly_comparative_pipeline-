import polars as pl
from pathlib import Path

# Get absolute path to the current script
script_dir = Path(__file__).resolve().parent

# Define sibling directory paths
csv_output_dir = script_dir.parent / "csv_ncbi"
mapping_dir = script_dir.parent / "taxonomic_mapping"

# Define file paths
assembly_file = script_dir / "00assembly_summary_genbank.txt"
taxmap_file = mapping_dir / "taxid_to_phylum.csv"
output_file = csv_output_dir / "ncbi_phylum_counts_polars.csv"

# Ensure output directory exists
csv_output_dir.mkdir(exist_ok=True)

# Print paths for debugging
print(f"Assembly file path: {assembly_file}")
print(f"Taxonomy mapping file path: {taxmap_file}")
print(f"Output file path: {output_file}")

def process_with_domain():
    """Process the assembly file with domain information included"""
    try:
        # Load assembly file
        print("Loading assembly summary file...")
        df = pl.read_csv(
            assembly_file,
            separator="\t",
            skip_rows=1,
            infer_schema_length=10000,
            null_values=["na", "NA", ""],
            ignore_errors=True,
            quote_char=None
        )
        print(f"Successfully loaded assembly file with {df.height} rows")

        # Rename column
        df = df.rename({"#assembly_accession": "assembly_accession"})

        # Load taxonomy mapping with domain information
        print("Loading taxonomy mapping file...")
        taxmap = pl.read_csv(taxmap_file)
        print(f"Successfully loaded taxonomy mapping with {taxmap.height} rows")
        
        # Check if domain column exists
        if 'domain' not in taxmap.columns:
            print("Warning: 'domain' column not found in taxonomy mapping file")
            return
        
        # Merge and drop null phyla
        print("Merging datasets...")
        df = df.join(taxmap, on="taxid", how="left")
        print(f"After merge: {df.height} rows")
        
        df = df.drop_nulls("phylum")
        print(f"After dropping nulls: {df.height} rows")
        
        # Group by phylum and domain, then count
        print("Counting genomes per phylum and domain...")
        counts = (
            df.group_by(["phylum", "domain"])
              .agg(pl.count().alias("ncbi_genome_count"))
              .sort("ncbi_genome_count", descending=True)
        )
        print(f"Found {counts.height} unique phylum-domain combinations")
        
        # Save output
        counts.write_csv(output_file)
        print(f"✅ Saved to {output_file} with domain information")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"Error: {e}")

# Run the processing function
process_with_domain()
