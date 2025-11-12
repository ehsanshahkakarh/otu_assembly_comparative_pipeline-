import pandas as pd
import os
import subprocess
import logging
from pathlib import Path
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("phylum_eukprot_parser.log"),
        logging.StreamHandler()
    ]
)

# Get absolute path to the current script
script_dir = Path(__file__).resolve().parent

# Define output directory paths
csv_output_dir = script_dir.parent / "csv_eukprot"
input_file = script_dir / "Eukprot_included_datasets.txt"
output_file = csv_output_dir / "eukprot_phylum_counts.csv"
detailed_output_file = csv_output_dir / "eukprot_phylum_with_accessions.csv"
verification_dir = csv_output_dir / "sanity_check" / "taxon_names"

# NCBI taxdump directory
TAXDUMP_DIR = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi")

# Ensure output directories exist
csv_output_dir.mkdir(exist_ok=True)
verification_dir.mkdir(parents=True, exist_ok=True)

def try_get_taxid(name, env=None):
    """
    Try to get NCBI taxid for a name using taxonkit.
    Returns the taxid if found, None otherwise.
    """
    if not name or pd.isna(name) or name == 'nan':
        return None

    if env is None:
        env = os.environ.copy()
        env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    try:
        # Run taxonkit name2taxid
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=name,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1].strip() and parts[1] != "0":
                return parts[1].strip()
    except Exception as e:
        logging.debug(f"Error getting taxid for {name}: {e}")

    return None

def map_phyla_to_taxids(phyla):
    """Map phylum names to NCBI taxids."""
    # Set up environment for taxonkit
    env = os.environ.copy()
    env["TAXONKIT_DB"] = str(TAXDUMP_DIR)

    # Check if taxonkit and database are available
    try:
        result = subprocess.run(
            ["taxonkit", "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        if result.returncode != 0:
            logging.warning("⚠️ taxonkit not available or database not found. Skipping taxid mapping.")
            return {}
        logging.info(f"✅ Using taxonkit {result.stdout.strip()}")
    except Exception as e:
        logging.warning(f"⚠️ taxonkit not available: {e}. Skipping taxid mapping.")
        return {}

    # Map phyla to taxids
    taxid_map = {}
    for phylum in phyla:
        if phylum and phylum != 'nan':
            taxid = try_get_taxid(phylum, env)
            if taxid:
                taxid_map[phylum] = taxid

    logging.info(f"✅ Mapped {len(taxid_map)} out of {len(phyla)} phyla to NCBI taxids")
    return taxid_map

def generate_verification_file(phyla_with_taxids):
    """Generate a file for taxonomic verification."""
    try:
        # Create the verification file
        verify_file = verification_dir / "eukaryotic_phyla.csv"

        # Create DataFrame with taxon_name and taxid columns
        verify_df = pd.DataFrame({
            "taxon_name": list(phyla_with_taxids.keys()),
            "taxid": list(phyla_with_taxids.values())
        })

        # Save to CSV
        verify_df.to_csv(verify_file, index=False)
        logging.info(f"✅ Generated verification file with {len(verify_df)} phyla at {verify_file}")

    except Exception as e:
        logging.error(f"❌ Error generating verification file: {e}")

def process_eukprot_data():
    """Process the EukProt dataset file to extract phylum information with taxids"""
    try:
        # Load EukProt dataset file
        df = pd.read_csv(input_file, sep='\t')
        logging.info(f"✅ Successfully loaded EukProt file with {len(df)} rows")

        # Extract phylum directly from Taxogroup2_UniEuk
        df["phylum"] = df["Taxogroup2_UniEuk"].astype(str).str.strip().str.replace("'", "")

        # Add domain column (all are Eukaryota)
        df["domain"] = "Eukaryota"

        # Create accession column from EukProt_ID
        df["accession_clean"] = df["EukProt_ID"]

        # Get unique phyla for taxid mapping
        unique_phyla = df["phylum"].dropna().unique()
        unique_phyla = [p for p in unique_phyla if p != 'nan']
        logging.info(f"✅ Found {len(unique_phyla)} unique phyla")

        # Map phyla to taxids
        phyla_to_taxids = map_phyla_to_taxids(unique_phyla)

        # Add taxid column if we have mappings
        if phyla_to_taxids:
            df["taxid"] = df["phylum"].map(phyla_to_taxids)
            logging.info(f"✅ Added taxids to {df['taxid'].notna().sum()} rows")

            # Generate verification file
            generate_verification_file(phyla_to_taxids)

        # Save detailed file with all accessions
        detailed_df = df.dropna(subset=["phylum"])

        # Select columns for detailed output
        if phyla_to_taxids:
            detailed_df[["accession_clean", "phylum", "domain", "taxid"]].to_csv(detailed_output_file, index=False)
        else:
            detailed_df[["accession_clean", "phylum", "domain"]].to_csv(detailed_output_file, index=False)

        logging.info(f"✅ Saved detailed file with {len(detailed_df)} entries to {detailed_output_file}")

        # Drop rows with null phylum values
        df = df.dropna(subset=["phylum"])
        logging.info(f"✅ After dropping null phylum values: {len(df)} rows")

        # Group by phylum and domain (and taxid if available), then count
        if phyla_to_taxids:
            counts = df.groupby(['phylum', 'domain', 'taxid']).size().reset_index(name='eukprot_genome_count')
        else:
            counts = df.groupby(['phylum', 'domain']).size().reset_index(name='eukprot_genome_count')

        logging.info(f"✅ Found {len(counts)} unique phylum combinations")

        # Add accession lists to the counts dataframe
        if phyla_to_taxids:
            accession_lists = df.groupby(['phylum', 'domain', 'taxid'])["accession_clean"].apply(list).reset_index()
            counts = pd.merge(counts, accession_lists, on=['phylum', 'domain', 'taxid'])
        else:
            accession_lists = df.groupby(['phylum', 'domain'])["accession_clean"].apply(list).reset_index()
            counts = pd.merge(counts, accession_lists, on=['phylum', 'domain'])

        # Sort by count (descending)
        counts = counts.sort_values('eukprot_genome_count', ascending=False)

        # Save
        counts.to_csv(output_file, index=False)
        logging.info(f"✅ Saved phylum counts to {output_file}")

        # Print preview
        logging.info("\n🔝 Top 10 phyla:")
        if phyla_to_taxids:
            logging.info(counts[["phylum", "domain", "taxid", "eukprot_genome_count"]].head(10).to_string())
        else:
            logging.info(counts[["phylum", "domain", "eukprot_genome_count"]].head(10).to_string())

        # Validate the data
        assert df[['accession_clean']].isna().sum().sum() == 0, "Missing accessions found"

    except Exception as e:
        logging.error(f"❌ Error: {e}")

if __name__ == "__main__":
    process_eukprot_data()
