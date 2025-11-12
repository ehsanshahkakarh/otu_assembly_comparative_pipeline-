import pandas as pd
from pathlib import Path
import re
from tqdm import tqdm

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
            total_count = domain_group['gtdb_genome_count'].sum()

            # Create a new row with the merged data
            merged_row = {
                level: base_taxon,
                'domain': domain,
                'gtdb_genome_count': total_count,
                f'{level}_base': base_taxon  # Keep this for filtering later
            }

            # Combine taxids if available
            if 'taxid' in domain_group.columns:
                all_taxids = []
                for taxid in domain_group['taxid']:
                    if pd.notna(taxid):
                        all_taxids.append(taxid)
                merged_row['taxid'] = all_taxids[0] if all_taxids else ''
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
    df = df.sort_values('gtdb_genome_count', ascending=False)

    print(f"Merged {len(taxa_to_merge)} {level} groups")
    return df

def process_bacterial_groups():
    """Process CSV files to extract unique taxonomic groups for Bacteria with their taxids."""
    # Setup paths using pathlib
    script_dir = Path(__file__).resolve().parent
    parent_dir = script_dir.parent  # csv_gtdb directory
    output_dir = script_dir / "taxon_names"
    
    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)
    
    # Define input and output files
    counts_files = {
        'genus': parent_dir / "gtdb_genus_counts.csv",
        'family': parent_dir / "gtdb_family_counts.csv",
        'phylum': parent_dir / "gtdb_phylum_counts.csv"
    }
    
    accessions_files = {
        'genus': parent_dir / "gtdb_genus_with_accessions.csv",
        'family': parent_dir / "gtdb_family_with_accessions.csv",
        'phylum': parent_dir / "gtdb_phylum_with_accessions.csv"
    }
    
    output_files = {
        'genus': output_dir / "bacterial_genera.csv",
        'family': output_dir / "bacterial_families.csv",
        'phylum': output_dir / "bacterial_phyla.csv"
    }
    
    print("📥 Processing bacterial taxonomic groups...")
    
    for level in counts_files.keys():
        try:
            print(f"\n🔍 Processing {level} level...")
            print(f"  Input counts: {counts_files[level].name}")
            print(f"  Input taxids: {accessions_files[level].name}")
            print(f"  Output: {output_files[level].name}")
            
            # Load counts CSV
            counts_df = pd.read_csv(counts_files[level])
            
            # Load accessions CSV for taxids
            accessions_df = pd.read_csv(accessions_files[level])
            
            # Check required columns
            required_columns = ['domain', level, 'gtdb_genome_count']
            missing_columns = [col for col in required_columns if col not in counts_df.columns]
            if missing_columns:
                raise ValueError(f"Counts file missing required columns: {', '.join(missing_columns)}")
            
            # Filter for Bacteria (case insensitive, no leading/trailing spaces)
            filtered_counts = counts_df[counts_df['domain'].str.strip().str.lower() == 'bacteria']
            filtered_accessions = accessions_df[accessions_df['domain'].str.strip().str.lower() == 'bacteria']
            
            # Get taxids from accessions file
            taxid_map = filtered_accessions.groupby([level, 'domain'])['taxid'].first().reset_index()
            
            # Merge counts with taxids
            filtered_counts = pd.merge(
                filtered_counts,
                taxid_map,
                on=[level, 'domain'],
                how='left'
            )
            
            # Merge related taxa and their counts
            merged_df = merge_related_taxa(filtered_counts, level)
            
            # Extract unique taxa and create result dataframe
            result_data = []
            for _, row in merged_df.iterrows():
                taxon_name = row[level]
                taxid = row['taxid']
                result_data.append({
                    'taxon_name': taxon_name,
                    'taxid': taxid
                })
            
            # Convert to DataFrame, deduplicate, and sort
            result_df = pd.DataFrame(result_data)
            result_df = result_df.drop_duplicates(subset=['taxon_name']).sort_values('taxon_name')
            
            # Save to output CSV with headers
            result_df.to_csv(output_files[level], index=False)
            
            print(f"✅ Saved {len(result_df)} unique Bacterial {level} with taxids to: {output_files[level].name}")
            
        except FileNotFoundError as e:
            print(f"❌ Error: File not found - {e}")
        except ValueError as e:
            print(f"❌ Error: {e}")
        except Exception as e:
            print(f"❌ Error processing {level}: {e}")

if __name__ == "__main__":
    process_bacterial_groups()
