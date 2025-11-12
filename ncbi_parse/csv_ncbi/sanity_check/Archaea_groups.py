import pandas as pd
from pathlib import Path

def process_archaeal_groups():
    """Process CSV files to extract unique taxonomic groups for Archaea with their taxids."""
    # Setup paths using pathlib
    script_dir = Path(__file__).resolve().parent
    parent_dir = script_dir.parent  # csv_ncbi directory
    output_dir = script_dir / "taxon_names"

    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)

    # Define input and output files
    input_files = {
        'genus': parent_dir / "ncbi_genus_counts.csv",
        'family': parent_dir / "ncbi_family_counts.csv",
        'phylum': parent_dir / "ncbi_phylum_counts.csv"
    }

    output_files = {
        'genus': output_dir / "archaeal_genera.csv",
        'family': output_dir / "archaeal_families.csv",
        'phylum': output_dir / "archaeal_phyla.csv"
    }

    print("📥 Processing archaeal taxonomic groups...")

    for level, input_file in input_files.items():
        try:
            print(f"\n🔍 Processing {level} level...")
            print(f"  Input: {input_file.name}")
            print(f"  Output: {output_files[level].name}")

            # Load CSV
            df = pd.read_csv(input_file)

            # Check required columns
            required_columns = ['domain', level, 'taxid']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(f"Input file missing required columns: {', '.join(missing_columns)}")

            # Filter for Archaea (case insensitive, no leading/trailing spaces)
            filtered_df = df[df['domain'].str.strip().str.lower() == 'archaea']

            # For each unique taxon, get the first taxid (since taxid is stored as a list in the CSV)
            result_data = []
            for _, row in filtered_df.iterrows():
                taxon_name = row[level]
                # Handle the case where taxid might be a string representation of a list
                taxid_value = row['taxid']
                if isinstance(taxid_value, str) and taxid_value.startswith('[') and taxid_value.endswith(']'):
                    # Extract the first taxid from the list representation
                    try:
                        # Convert string representation of list to actual list
                        taxid_list = eval(taxid_value)
                        taxid = taxid_list[0] if taxid_list else None
                    except:
                        taxid = None
                else:
                    taxid = taxid_value

                result_data.append({
                    'taxon_name': taxon_name,
                    'taxid': taxid
                })

            # Convert to DataFrame, deduplicate, and sort
            result_df = pd.DataFrame(result_data)
            result_df = result_df.drop_duplicates(subset=['taxon_name']).sort_values('taxon_name')

            # Save to output CSV with headers
            result_df.to_csv(output_files[level], index=False)

            print(f"✅ Saved {len(result_df)} unique Archaeal {level} with taxids to: {output_files[level].name}")

        except FileNotFoundError as e:
            print(f"❌ Error: File not found - {e}")
        except ValueError as e:
            print(f"❌ Error: {e}")
        except Exception as e:
            print(f"❌ Error processing {level}: {e}")

if __name__ == "__main__":
    process_archaeal_groups()