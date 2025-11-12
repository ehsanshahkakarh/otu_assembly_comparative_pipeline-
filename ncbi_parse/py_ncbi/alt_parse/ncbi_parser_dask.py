import dask.dataframe as dd
import pandas as pd
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Get absolute path to the current script
script_dir = Path(__file__).resolve().parent

# Define sibling directory paths
csv_output_dir = script_dir.parent / "csv_ncbi"
mapping_dir = script_dir.parent / "taxonomic_mapping"

# Define file paths
assembly_file = script_dir / "00assembly_summary_genbank.txt"
taxmap_file = mapping_dir / "taxid_to_phylum.csv"
output_file = csv_output_dir / "ncbi_phylum_counts_dask.csv"

# Ensure output directory exists
csv_output_dir.mkdir(exist_ok=True)

# Print paths for debugging
logger.info(f"Assembly file path: {assembly_file}")
logger.info(f"Taxonomy mapping file path: {taxmap_file}")
logger.info(f"Output file path: {output_file}")

def check_input_files():
    """Verify input files exist before processing"""
    required_files = [assembly_file, taxmap_file]
    for file in required_files:
        if not file.exists():
            raise FileNotFoundError(f"Required file {file} not found")

def process_with_domain():
    """Process the assembly file with domain information included"""
    try:
        # Check files exist
        check_input_files()
        
        logger.info("Loading NCBI assembly summary file...")
        # Load the NCBI assembly summary file with proper data types
        df = dd.read_csv(
            assembly_file, 
            sep="\t", 
            skiprows=1, 
            assume_missing=True,
            dtype={
                'taxid': 'Int64',  # Handle potential NaN in taxid
                'assembly_level': 'object',
                'version_status': 'object',
                'relation_to_type_material': 'object',
                'non_coding_gene_count': 'object',
                'protein_coding_gene_count': 'object',
                'total_gene_count': 'object',
                'genome_rep': 'object',
                'seq_rel_date': 'object',
                'submitter': 'object'
            },
            na_values=['na', 'NA', '', 'null'],
            low_memory=False
        )
        logger.info(f"Successfully loaded assembly file")

        # Handle potential column name variations
        if "#assembly_accession" in df.columns:
            df = df.rename(columns={"#assembly_accession": "assembly_accession"})
        elif "# assembly_accession" in df.columns:
            df = df.rename(columns={"# assembly_accession": "assembly_accession"})

        logger.info("Loading taxid to phylum mapping...")
        # Load taxid → phylum mapping (with Pandas since it's expected to be smaller)
        taxmap = pd.read_csv(
            taxmap_file,
            dtype={'taxid': 'Int64'}  # Ensure consistent type with df
        )
        logger.info(f"Successfully loaded taxonomy mapping with {len(taxmap)} rows")
        
        # Check if domain column exists
        if 'domain' not in taxmap.columns:
            logger.warning("'domain' column not found in taxonomy mapping file")
            return

        logger.info("Merging datasets...")
        # Merge and handle missing data
        df = df.merge(taxmap, on="taxid", how="left")
        df = df.dropna(subset=["phylum"])
        logger.info("Successfully merged and filtered datasets")

        logger.info("Counting genomes per phylum and domain...")
        # Group by phylum and domain, then count
        # Note: For Dask, we need to compute() to materialize the result
        counts = (df.groupby(['phylum', 'domain'])
                 .size()
                 .compute()  # Trigger actual computation
                 .reset_index(name='ncbi_genome_count'))
        
        # Sort by count in descending order
        counts = counts.sort_values("ncbi_genome_count", ascending=False)
        logger.info(f"Found {len(counts)} unique phylum-domain combinations")

        logger.info("Saving results...")
        # Save output with proper encoding
        counts.to_csv(output_file, index=False, encoding='utf-8')
        logger.info(f"✅ Dask output saved successfully to {output_file}")

        # Print summary
        print("\nSummary:")
        print(f"Total phyla processed: {len(counts['phylum'].unique())}")
        print(f"Total genomes counted: {counts['ncbi_genome_count'].sum():,}")
        print("\nTop 5 phyla by genome count:")
        print(counts.head().to_string(index=False))

    except FileNotFoundError as e:
        logger.error(f"File not found error: {e}")
        raise
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    try:
        process_with_domain()
    except Exception as e:
        logger.error(f"Script failed: {e}")
        raise

