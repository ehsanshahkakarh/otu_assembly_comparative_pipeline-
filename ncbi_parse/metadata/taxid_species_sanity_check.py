#!/usr/bin/env python3
"""
NCBI Taxid and Species_taxid Sanity Check Script
Created: 2025-08-18

This script parses the NCBI assembly file and extracts taxid, species_taxid, and organism_name
columns to create two output files for sanity checking:

1. No grouping: All records with taxid, species_taxid, organism_name
2. Grouped: Group by species_taxid, remove taxid column, show unique organism names per species

Usage:
    python taxid_species_sanity_check.py [--sample-size N]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys
from tqdm import tqdm

class TaxidSpeciesSanityCheck:
    """Sanity check for taxid and species_taxid relationships."""
    
    def __init__(self, input_dir=None, output_dir=None):
        self.script_dir = Path(__file__).resolve().parent
        self.input_dir = Path(input_dir) if input_dir else self.script_dir
        self.output_dir = Path(output_dir) if output_dir else self.script_dir
        
        # Find assembly file
        self.assembly_file = self._find_assembly_file()
        if not self.assembly_file:
            raise FileNotFoundError("Could not find 00assembly_summary_genbank.txt")
            
        print(f"📁 Input file: {self.assembly_file}")
        print(f"📁 Output directory: {self.output_dir}")
    
    def _find_assembly_file(self):
        """Find the NCBI assembly summary file."""
        possible_files = [
            self.input_dir / "00assembly_summary_genbank.txt",
            self.script_dir / "00assembly_summary_genbank.txt"
        ]
        
        for file_path in possible_files:
            if file_path.exists():
                return file_path
        return None
    
    def load_assembly_data(self, sample_size=None):
        """Load assembly data with focus on taxid columns."""
        print("📂 Loading NCBI assembly data...")
        
        # Read header first
        with open(self.assembly_file, 'r') as f:
            lines = f.readlines()
        
        # Find header line
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#assembly_accession'):
                header_line = line.strip().lstrip('#')
                data_start = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find header line")
        
        # Load data with sampling if requested
        if sample_size:
            print(f"📊 Loading sample of {sample_size:,} records...")
            # Read sample of lines
            data_lines = [header_line]
            step = max(1, (len(lines) - data_start) // sample_size)
            for i in range(data_start, len(lines), step):
                if len(data_lines) >= sample_size + 1:  # +1 for header
                    break
                if not lines[i].startswith('#'):
                    data_lines.append(lines[i])
            
            # Create temporary dataframe
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp_file:
                tmp_file.writelines(data_lines)
                tmp_filename = tmp_file.name
            
            df = pd.read_csv(tmp_filename, sep='\t', low_memory=False)
            Path(tmp_filename).unlink()
        else:
            # Load full dataset
            print("⚠️  Loading full dataset - this may take several minutes...")
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp_file:
                tmp_file.write(header_line + '\n')
                for line in tqdm(lines[data_start:], desc="Processing lines"):
                    if not line.startswith('#'):
                        tmp_file.write(line)
                tmp_filename = tmp_file.name
            
            df = pd.read_csv(tmp_filename, sep='\t', low_memory=False)
            Path(tmp_filename).unlink()
        
        print(f"✅ Loaded {len(df):,} assembly records")
        
        # Extract only the columns we need
        required_columns = ['taxid', 'species_taxid', 'organism_name']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Select and clean the required columns
        df_clean = df[required_columns].copy()
        
        # Convert taxid columns to numeric, keeping as strings for now to preserve original values
        print("🔧 Processing taxid columns...")
        
        # Basic statistics
        print(f"📊 Data quality check:")
        print(f"   Total records: {len(df_clean):,}")
        print(f"   Records with taxid: {df_clean['taxid'].notna().sum():,}")
        print(f"   Records with species_taxid: {df_clean['species_taxid'].notna().sum():,}")
        print(f"   Records with organism_name: {df_clean['organism_name'].notna().sum():,}")
        
        return df_clean
    
    def create_no_grouping_output(self, df):
        """Create output file with no grouping - all records as-is."""
        print("📝 Creating no-grouping output file...")
        
        # Sort by species_taxid, then taxid for better organization
        df_sorted = df.sort_values(['species_taxid', 'taxid'], na_position='last')
        
        # Output file
        output_file = self.output_dir / "taxid_species_no_grouping.csv"
        
        # Save with all records
        df_sorted.to_csv(output_file, index=False)
        
        print(f"💾 Saved no-grouping file: {output_file}")
        print(f"   Records: {len(df_sorted):,}")
        
        # Show sample
        print(f"\n📋 Sample of no-grouping output (first 10 rows):")
        print(df_sorted.head(10).to_string(index=False))
        
        return output_file
    
    def get_species_names_from_taxonkit(self, species_taxids):
        """Get species names using taxonkit lineage."""
        print("🔍 Getting species names from taxonkit...")

        import tempfile
        import subprocess

        # Create temporary file with species taxids
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
            temp_filename = temp_file.name
            for taxid in species_taxids:
                temp_file.write(f"{taxid}\n")

        try:
            # Run taxonkit lineage to get species names
            result = subprocess.run([
                "taxonkit", "lineage", "--show-lineage-taxids", temp_filename
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                print(f"⚠️  Taxonkit warning: {result.stderr}")
                # Fallback: create mapping with taxid as name
                return {taxid: f"species_{taxid}" for taxid in species_taxids}

            # Parse results to extract species names
            taxid_to_species_name = {}
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        taxid = int(parts[0])
                        lineage = parts[1]

                        # Extract species name from lineage (last part before strain info)
                        lineage_parts = lineage.split(';')
                        species_name = lineage_parts[-1].strip() if lineage_parts else f"species_{taxid}"

                        # Clean up species name
                        if species_name and species_name != 'unclassified':
                            taxid_to_species_name[taxid] = species_name
                        else:
                            taxid_to_species_name[taxid] = f"species_{taxid}"
                    else:
                        taxid_to_species_name[int(parts[0])] = f"species_{parts[0]}"

            # Clean up temp file
            Path(temp_filename).unlink()

            print(f"✅ Retrieved names for {len(taxid_to_species_name):,} species")
            return taxid_to_species_name

        except Exception as e:
            print(f"❌ Error running taxonkit: {e}")
            # Clean up temp file
            Path(temp_filename).unlink()
            # Fallback: create mapping with taxid as name
            return {taxid: f"species_{taxid}" for taxid in species_taxids}

    def create_grouped_output(self, df):
        """Create output file grouped by species_taxid, using taxonkit for species names."""
        print("📝 Creating grouped output file...")

        # Group by species_taxid only - one row per species
        print("🔄 Grouping by species_taxid...")

        # Group by species_taxid and count records
        grouped_stats = df.groupby('species_taxid').agg({
            'taxid': 'count'  # Total number of records (assemblies) for this species
        }).reset_index()

        # Rename columns
        grouped_stats.columns = ['species_taxid', 'total_records']

        # Get species names using taxonkit
        unique_species_taxids = grouped_stats['species_taxid'].unique()
        species_name_mapping = self.get_species_names_from_taxonkit(unique_species_taxids)

        # Add species names to the dataframe
        grouped_stats['species_name'] = grouped_stats['species_taxid'].map(species_name_mapping)

        # Reorder columns
        grouped_stats = grouped_stats[['species_taxid', 'species_name', 'total_records']]

        # Sort by total_records (descending) to show species with most assemblies first
        grouped_stats = grouped_stats.sort_values('total_records', ascending=False)

        # Output file
        output_file = self.output_dir / "taxid_species_grouped.csv"

        # Save grouped data
        grouped_stats.to_csv(output_file, index=False)

        print(f"💾 Saved grouped file: {output_file}")
        print(f"   Unique species: {len(grouped_stats):,}")
        print(f"   Total records processed: {grouped_stats['total_records'].sum():,}")

        # Show sample
        print(f"\n📋 Sample of grouped output (first 10 rows):")
        print(grouped_stats.head(10).to_string(index=False))

        # Show some interesting statistics
        print(f"\n📊 Grouping statistics:")
        print(f"   Average records per species: {grouped_stats['total_records'].mean():.2f}")
        print(f"   Species with most records: {grouped_stats['total_records'].max():,}")
        print(f"   Species with only 1 record: {(grouped_stats['total_records'] == 1).sum():,}")

        # Show top species by record count
        print(f"\n🔝 Top 10 species by record count:")
        for i, row in grouped_stats.head(10).iterrows():
            print(f"   {row['species_taxid']}: {row['total_records']:,} records - {row['species_name']}")

        return output_file
    
    def run_sanity_check(self, sample_size=None):
        """Run the complete sanity check process."""
        try:
            print("🚀 Starting taxid/species_taxid sanity check...")
            
            # Load data
            df = self.load_assembly_data(sample_size)
            
            # Create both output files
            no_grouping_file = self.create_no_grouping_output(df)
            grouped_file = self.create_grouped_output(df)
            
            print(f"\n" + "="*60)
            print("SANITY CHECK COMPLETE")
            print("="*60)
            print(f"📁 Output files created:")
            print(f"   1. No grouping: {no_grouping_file}")
            print(f"   2. Grouped: {grouped_file}")
            
            if sample_size:
                print(f"📝 Note: Results based on sample of {sample_size:,} records")
                print(f"   For complete analysis, run without --sample-size")
            
            print(f"\n🔍 Files ready for manual inspection in: {self.output_dir}")
            
        except Exception as e:
            print(f"❌ Error during sanity check: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description='NCBI taxid/species_taxid sanity check')
    parser.add_argument('--sample-size', type=int, help='Sample size for testing (default: full dataset)')
    parser.add_argument('--input-dir', help='Directory containing assembly file')
    parser.add_argument('--output-dir', help='Directory for output files')
    
    args = parser.parse_args()
    
    checker = TaxidSpeciesSanityCheck(
        input_dir=args.input_dir,
        output_dir=args.output_dir
    )
    checker.run_sanity_check(sample_size=args.sample_size)

if __name__ == "__main__":
    main()
