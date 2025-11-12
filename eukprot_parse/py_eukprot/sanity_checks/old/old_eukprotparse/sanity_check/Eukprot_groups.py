import pandas as pd
from pathlib import Path
import re
from tqdm import tqdm
import os
import subprocess
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('taxid_lookup.log'),
        logging.StreamHandler()
    ]
)

def standardize_taxon_name(taxon: str) -> str:
    """
    Standardize taxon names by removing suffixes like _A, _B, etc.

    Args:
        taxon (str): Original taxon name

    Returns:
        str: Standardized taxon name
    """
    # Check if the taxon name ends with _[A-Z]
    if re.search(r'_[A-Z]$', taxon):
        # Remove the suffix
        return taxon.rsplit('_', 1)[0]
    return taxon

def get_taxid_from_taxonkit(taxon_name: str, taxdump_dir: str) -> str:
    """
    Get taxid for a taxon name using taxonkit.

    Args:
        taxon_name (str): Taxon name to look up
        taxdump_dir (str): Path to taxdump directory

    Returns:
        str: Taxid if found, empty string otherwise
    """
    # Try different formats of the taxon name
    name_variants = [
        taxon_name,  # Original name
        taxon_name.lower(),  # Lowercase
        taxon_name.replace('ae', 'a'),  # Remove 'ae' suffix
        taxon_name.replace('ceae', 'a'),  # Remove 'ceae' suffix
        taxon_name.replace('ales', 'a'),  # Remove 'ales' suffix
    ]
    
    for name in name_variants:
        try:
            # Run taxonkit name2taxid command
            cmd = ['taxonkit', 'name2taxid', '-d', taxdump_dir, '-i', name]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse output
            if result.stdout.strip():
                # Output format is typically: name<TAB>taxid
                taxid = result.stdout.strip().split('\t')[1]
                if taxid:
                    logging.info(f"Found taxid {taxid} for {taxon_name} (using variant: {name})")
                    return taxid
        except subprocess.CalledProcessError:
            continue
        except Exception as e:
            logging.warning(f"Error getting taxid for {taxon_name} (variant: {name}): {e}")
            continue
    
    logging.warning(f"No taxid found for {taxon_name} after trying all variants")
    return ''

def merge_related_taxa(df: pd.DataFrame, level: str) -> pd.DataFrame:
    """
    Merge counts for related taxa (e.g., Bacillota and Bacillota_A).

    Args:
        df (pd.DataFrame): DataFrame with taxon counts
        level (str): Taxonomic level (genus, family, or phylum)

    Returns:
        pd.DataFrame: DataFrame with merged taxon counts
    """
    print(f"Merging related {level}...")

    # Create a copy to avoid modifying the original
    df = df.copy()

    # Create a new column with standardized taxon names
    df[f'{level}_base'] = df[level].apply(standardize_taxon_name)

    # Find taxa that need to be merged (those where taxon_base != taxon)
    taxa_to_merge = df[df[f'{level}_base'] != df[level]][f'{level}_base'].unique()

    if len(taxa_to_merge) == 0:
        print(f"No related {level} found to merge")
        return df

    # For each base taxon that has variants
    merged_rows = []
    for base_taxon in tqdm(taxa_to_merge, desc=f"Merging {level}", unit="taxon"):
        # Get all rows for this base taxon and its variants
        related_rows = df[df[f'{level}_base'] == base_taxon]

        # Group by domain and sum the genome counts
        for domain, domain_group in related_rows.groupby('domain'):
            # Sum the genome counts
            total_count = domain_group['eukprot_genome_count'].sum()

            # Create a new row with the merged data
            merged_row = {
                level: base_taxon,
                'domain': domain,
                'eukprot_genome_count': total_count,
                f'{level}_base': base_taxon  # Keep this for filtering later
            }
            merged_rows.append(merged_row)

            # Log the merge
            variants = domain_group[level].tolist()
            print(f"Merged {len(variants)} variants of {base_taxon} ({domain}): {', '.join(variants)}")

    # Create a DataFrame from the merged rows
    merged_df = pd.DataFrame(merged_rows)

    # Remove the original rows for the merged taxa (both variants and base taxa)
    df = df[~(df[f'{level}_base'].isin(taxa_to_merge))]

    # Add the merged rows
    df = pd.concat([df, merged_df], ignore_index=True)

    # Remove the temporary column and sort by genome count
    df = df.drop(columns=[f'{level}_base'])
    df = df.sort_values('eukprot_genome_count', ascending=False)

    print(f"Merged {len(taxa_to_merge)} {level} groups")
    return df

def process_eukprot_groups():
    """Process CSV files to extract unique taxonomic groups for Eukaryotes with their taxids."""
    # Setup paths using pathlib
    script_dir = Path(__file__).resolve().parent
    parent_dir = script_dir.parent  # csv_eukprot directory
    output_dir = script_dir / "taxon_names"
    taxdump_dir = "/clusterfs/jgi/groups/science/homes/ehsankakar/.taxonkit"
    
    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)
    
    # Define input and output files
    counts_files = {
        'genus': parent_dir / "eukprot_genus_counts.csv",
        'family': parent_dir / "eukprot_family_counts.csv",
        'phylum': parent_dir / "eukprot_phylum_counts.csv"
    }
    
    output_files = {
        'genus': output_dir / "eukprot_genera.csv",
        'family': output_dir / "eukprot_families.csv",
        'phylum': output_dir / "eukprot_phyla.csv"
    }
    
    # Create a summary file for taxid lookup results
    summary_file = output_dir / "taxid_lookup_summary.csv"
    summary_data = []
    
    print("📥 Processing eukaryotic taxonomic groups...")
    
    for level in counts_files.keys():
        try:
            print(f"\n🔍 Processing {level} level...")
            print(f"  Input counts: {counts_files[level].name}")
            print(f"  Output: {output_files[level].name}")
            
            # Load counts CSV
            counts_df = pd.read_csv(counts_files[level])
            
            # Check required columns
            required_columns = ['domain', level, 'eukprot_genome_count']
            missing_columns = [col for col in required_columns if col not in counts_df.columns]
            if missing_columns:
                raise ValueError(f"Counts file missing required columns: {', '.join(missing_columns)}")
            
            # Filter for Eukaryota (case insensitive, no leading/trailing spaces)
            filtered_counts = counts_df[counts_df['domain'].str.strip().str.lower() == 'eukaryota']
            
            # Merge related taxa and their counts
            merged_df = merge_related_taxa(filtered_counts, level)
            
            # Extract unique taxa and create result dataframe
            result_data = []
            for _, row in tqdm(merged_df.iterrows(), desc="Getting taxids", total=len(merged_df)):
                taxon_name = row[level]
                # Get taxid using taxonkit
                taxid = get_taxid_from_taxonkit(taxon_name, taxdump_dir)
                result_data.append({
                    'taxon_name': taxon_name,
                    'taxid': taxid
                })
                
                # Add to summary
                summary_data.append({
                    'level': level,
                    'taxon_name': taxon_name,
                    'taxid': taxid,
                    'found': bool(taxid)
                })
            
            # Convert to DataFrame, deduplicate, and sort
            result_df = pd.DataFrame(result_data)
            result_df = result_df.drop_duplicates(subset=['taxon_name']).sort_values('taxon_name')
            
            # Save to output CSV with headers
            result_df.to_csv(output_files[level], index=False)
            
            # Log statistics
            total_taxa = len(result_df)
            found_taxids = result_df['taxid'].notna().sum()
            logging.info(f"Level {level}: Found {found_taxids}/{total_taxa} taxids ({found_taxids/total_taxa*100:.1f}%)")
            
            print(f"✅ Saved {len(result_df)} unique Eukaryotic {level} with taxids to: {output_files[level].name}")
            
        except FileNotFoundError as e:
            logging.error(f"File not found: {e}")
        except ValueError as e:
            logging.error(f"Value error: {e}")
        except Exception as e:
            logging.error(f"Error processing {level}: {e}")
    
    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_file, index=False)
    logging.info(f"Saved taxid lookup summary to {summary_file}")

if __name__ == "__main__":
    process_eukprot_groups() 