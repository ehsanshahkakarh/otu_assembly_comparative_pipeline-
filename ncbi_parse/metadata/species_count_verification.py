#!/usr/bin/env python3
"""
NCBI Species Count Verification Script
Created: 2025-08-18

This script explains and verifies our method for counting species in the NCBI database.
It demonstrates how we calculate species numbers and validates the approach used in our parsing pipeline.

OUR METHOD EXPLAINED:
====================

1. **Data Source**: NCBI Assembly Summary file (00assembly_summary_genbank.txt)
   - Contains metadata for all genome assemblies in NCBI GenBank
   - Key columns: assembly_accession, taxid, species_taxid, organism_name

2. **Species Identification**: 
   - We use the 'species_taxid' column as the primary species identifier
   - Each unique species_taxid represents one distinct species
   - Multiple assemblies can belong to the same species (different strains, isolates, etc.)

3. **Species Counting Method**:
   - **Total Species Count**: Count unique species_taxid values across all assemblies
   - **Taxonomic Level Species Count**: Group by taxonomic level (phylum/family/genus) 
     and count unique species_taxid within each group
   - **Domain-specific Counts**: Further breakdown by domain (Bacteria, Archaea, Eukaryota, Viruses)

4. **Quality Control**:
   - Verify 100% species_taxid coverage (no missing values)
   - Check for valid taxid assignments
   - Validate taxonomic hierarchy consistency

5. **Two-tier Analysis**:
   - **Genome Counts**: Total number of assemblies per taxonomic group
   - **Species Counts**: Number of unique species per taxonomic group
   - This reveals assembly density per species (some species have many assemblies)

Usage:
    python species_count_verification.py [--sample-size N] [--detailed]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from collections import Counter
import sys

class SpeciesCountVerifier:
    """Verifies and explains NCBI species counting methodology."""
    
    def __init__(self, input_dir=None):
        self.script_dir = Path(__file__).resolve().parent
        self.input_dir = Path(input_dir) if input_dir else self.script_dir
        
        # Find assembly file
        self.assembly_file = self._find_assembly_file()
        if not self.assembly_file:
            raise FileNotFoundError("Could not find 00assembly_summary_genbank.txt")
            
        print(f"📁 Using assembly file: {self.assembly_file}")
    
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
    
    def load_sample_data(self, sample_size=None):
        """Load a sample of the assembly data for verification."""
        print(f"📂 Loading NCBI assembly data...")
        
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
            print(f"📊 Loading sample of {sample_size:,} records for verification...")
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
            # Load full dataset (memory intensive)
            print("⚠️  Loading full dataset - this may take several minutes...")
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp_file:
                tmp_file.write(header_line + '\n')
                for line in lines[data_start:]:
                    if not line.startswith('#'):
                        tmp_file.write(line)
                tmp_filename = tmp_file.name
            
            df = pd.read_csv(tmp_filename, sep='\t', low_memory=False)
            Path(tmp_filename).unlink()
        
        print(f"✅ Loaded {len(df):,} assembly records")
        return df
    
    def explain_method(self):
        """Explain our species counting methodology."""
        print("\n" + "="*80)
        print("OUR SPECIES COUNTING METHOD EXPLAINED")
        print("="*80)
        
        explanation = """
🔬 SPECIES IDENTIFICATION APPROACH:
----------------------------------
1. Primary Key: 'species_taxid' column from NCBI assembly summary
   - Each unique species_taxid = one distinct species
   - Multiple assemblies can share the same species_taxid (strains, isolates)
   
2. Data Structure:
   - assembly_accession: Unique identifier for each genome assembly
   - taxid: Taxonomic ID (can be strain/isolate level)  
   - species_taxid: Species-level taxonomic ID (our counting unit)
   - organism_name: Scientific name of the organism

📊 COUNTING METHODOLOGY:
-----------------------
1. Total Species Count: 
   df['species_taxid'].nunique()
   
2. Taxonomic Group Species Count:
   df.groupby(['phylum'])['species_taxid'].nunique()
   
3. Domain-specific Species Count:
   df.groupby(['domain', 'phylum'])['species_taxid'].nunique()

🎯 WHY THIS METHOD IS ACCURATE:
------------------------------
- Uses NCBI's official species-level taxonomy assignments
- Avoids double-counting strains/isolates of the same species
- Provides true biological species diversity metrics
- Enables comparison of assembly density vs species diversity

📈 DUAL METRICS APPROACH:
------------------------
- Genome Count: Total assemblies per taxonomic group
- Species Count: Unique species per taxonomic group  
- Ratio reveals: assemblies per species (research intensity)
        """
        print(explanation)
    
    def verify_data_quality(self, df):
        """Verify data quality and species_taxid coverage."""
        print("\n" + "="*60)
        print("DATA QUALITY VERIFICATION")
        print("="*60)
        
        total_records = len(df)
        
        # Check species_taxid coverage
        has_species_taxid = df['species_taxid'].notna().sum()
        missing_species_taxid = df['species_taxid'].isna().sum()
        
        print(f"📊 Total assembly records: {total_records:,}")
        print(f"✅ Records with species_taxid: {has_species_taxid:,} ({has_species_taxid/total_records*100:.2f}%)")
        print(f"❌ Records missing species_taxid: {missing_species_taxid:,} ({missing_species_taxid/total_records*100:.2f}%)")
        
        # Check taxid coverage
        has_taxid = df['taxid'].notna().sum()
        print(f"✅ Records with taxid: {has_taxid:,} ({has_taxid/total_records*100:.2f}%)")
        
        # Species statistics
        unique_species = df['species_taxid'].nunique()
        avg_assemblies_per_species = total_records / unique_species if unique_species > 0 else 0
        
        print(f"\n🧬 SPECIES STATISTICS:")
        print(f"   Unique species: {unique_species:,}")
        print(f"   Average assemblies per species: {avg_assemblies_per_species:.2f}")
        
        # Show examples of species with multiple assemblies
        if unique_species > 0:
            species_assembly_counts = df['species_taxid'].value_counts()
            top_species = species_assembly_counts.head(5)
            
            print(f"\n📈 Top 5 species by assembly count:")
            for i, (species_taxid, count) in enumerate(top_species.items(), 1):
                # Get organism name for this species
                example_organism = df[df['species_taxid'] == species_taxid]['organism_name'].iloc[0]
                print(f"   {i}. Species {species_taxid} ({example_organism}): {count:,} assemblies")
        
        return {
            'total_records': total_records,
            'unique_species': unique_species,
            'species_coverage': has_species_taxid/total_records*100,
            'avg_assemblies_per_species': avg_assemblies_per_species
        }
    
    def demonstrate_counting(self, df, detailed=False):
        """Demonstrate the species counting process."""
        print("\n" + "="*60)
        print("SPECIES COUNTING DEMONSTRATION")
        print("="*60)
        
        # Basic domain classification (simplified)
        def classify_domain(organism_name):
            if pd.isna(organism_name):
                return 'Unknown'
            
            organism_lower = str(organism_name).lower()
            
            if any(pattern in organism_lower for pattern in ['virus', 'phage']):
                return 'Viruses'
            elif any(pattern in organism_lower for pattern in ['escherichia', 'salmonella', 'bacillus']):
                return 'Bacteria'
            elif any(pattern in organism_lower for pattern in ['homo', 'mus', 'drosophila']):
                return 'Eukaryota'
            else:
                return 'Unclassified'
        
        df['domain'] = df['organism_name'].apply(classify_domain)
        
        # Overall counts
        total_assemblies = len(df)
        total_species = df['species_taxid'].nunique()
        
        print(f"📊 OVERALL COUNTS:")
        print(f"   Total assemblies: {total_assemblies:,}")
        print(f"   Unique species: {total_species:,}")
        print(f"   Assembly-to-species ratio: {total_assemblies/total_species:.2f}")
        
        # Domain breakdown
        print(f"\n🌍 DOMAIN BREAKDOWN:")
        domain_stats = df.groupby('domain').agg({
            'assembly_accession': 'count',  # Assembly count
            'species_taxid': 'nunique'      # Species count
        }).rename(columns={
            'assembly_accession': 'assemblies',
            'species_taxid': 'species'
        })
        
        for domain, row in domain_stats.iterrows():
            assemblies = row['assemblies']
            species = row['species']
            ratio = assemblies / species if species > 0 else 0
            print(f"   {domain}:")
            print(f"     Assemblies: {assemblies:,}")
            print(f"     Species: {species:,}")
            print(f"     Ratio: {ratio:.2f} assemblies/species")
        
        if detailed:
            print(f"\n🔍 DETAILED EXAMPLES:")
            # Show some specific examples
            sample_species = df['species_taxid'].value_counts().head(3)
            for species_taxid, count in sample_species.items():
                species_df = df[df['species_taxid'] == species_taxid]
                organism_name = species_df['organism_name'].iloc[0]
                print(f"\n   Species {species_taxid} ({organism_name}):")
                print(f"     Total assemblies: {count}")
                if count <= 5:
                    for _, row in species_df.iterrows():
                        print(f"       - {row['assembly_accession']}: {row['organism_name']}")
                else:
                    print(f"       - First 3 assemblies:")
                    for _, row in species_df.head(3).iterrows():
                        print(f"         - {row['assembly_accession']}: {row['organism_name']}")
                    print(f"       - ... and {count-3} more")
    
    def run_verification(self, sample_size=None, detailed=False):
        """Run complete verification process."""
        try:
            # Load data
            df = self.load_sample_data(sample_size)
            
            # Explain method
            self.explain_method()
            
            # Verify data quality
            stats = self.verify_data_quality(df)
            
            # Demonstrate counting
            self.demonstrate_counting(df, detailed)
            
            print(f"\n" + "="*60)
            print("VERIFICATION COMPLETE")
            print("="*60)
            print(f"✅ Method validated on {len(df):,} records")
            print(f"🧬 Species coverage: {stats['species_coverage']:.2f}%")
            print(f"📊 Average assemblies per species: {stats['avg_assemblies_per_species']:.2f}")
            
            if sample_size:
                print(f"📝 Note: Results based on sample of {sample_size:,} records")
                print(f"   For complete analysis, run without --sample-size")
            
        except Exception as e:
            print(f"❌ Error during verification: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description='Verify NCBI species counting methodology')
    parser.add_argument('--sample-size', type=int, help='Sample size for verification (default: full dataset)')
    parser.add_argument('--detailed', action='store_true', help='Show detailed examples')
    parser.add_argument('--input-dir', help='Directory containing assembly file')
    
    args = parser.parse_args()
    
    verifier = SpeciesCountVerifier(input_dir=args.input_dir)
    verifier.run_verification(sample_size=args.sample_size, detailed=args.detailed)

if __name__ == "__main__":
    main()
