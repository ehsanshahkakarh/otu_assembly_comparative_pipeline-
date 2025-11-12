#!/usr/bin/env python3
"""
Test Comprehensive Output for Unified NCBI Parser

This script tests the new comprehensive output that includes both genome counts
and species subset counts in a single file.

Usage:
    python test_comprehensive_output.py --test-genus
    python test_comprehensive_output.py --compare-with-original
"""

import subprocess
import pandas as pd
import argparse
from pathlib import Path

def run_unified_parser_test(level: str = "genus") -> bool:
    """Test the unified parser with comprehensive output."""
    print(f"🧪 Testing unified parser comprehensive output for {level} level")
    print("=" * 60)
    
    try:
        # Run unified parser
        result = subprocess.run(
            ["python", "unified_ncbi_parser.py", "--level", level],
            check=True,
            capture_output=True,
            text=True
        )
        
        print("✅ Unified parser completed successfully")
        
        # Check output files
        script_dir = Path(__file__).resolve().parent
        output_dir = script_dir.parent  # ncbi_parse directory
        
        expected_files = [
            f"ncbi_{level}_comprehensive_counts.csv",
            f"ncbi_{level}_with_accessions_classified.csv",
            f"ncbi_{level}_species_subset_with_accessions_classified.csv"
        ]
        
        print(f"\n📁 Checking output files in {output_dir}:")
        all_files_exist = True
        
        for filename in expected_files:
            file_path = output_dir / filename
            if file_path.exists():
                size_info = f"({file_path.stat().st_size / 1024:.1f} KB)" if file_path.stat().st_size < 1024*1024 else f"({file_path.stat().st_size / (1024*1024):.1f} MB)"
                print(f"   ✅ {filename} {size_info}")
            else:
                print(f"   ❌ {filename} - NOT FOUND")
                all_files_exist = False
        
        if all_files_exist:
            # Analyze the comprehensive counts file
            analyze_comprehensive_counts(output_dir / f"ncbi_{level}_comprehensive_counts.csv")
            return True
        else:
            print(f"❌ Some expected files are missing")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Unified parser failed: {e}")
        if e.stderr:
            print(f"Error output: {e.stderr}")
        return False

def analyze_comprehensive_counts(file_path: Path):
    """Analyze the comprehensive counts file."""
    print(f"\n📊 Analyzing comprehensive counts file: {file_path.name}")
    print("-" * 50)
    
    try:
        df = pd.read_csv(file_path)
        
        print(f"📋 File structure:")
        print(f"   Rows: {len(df):,}")
        print(f"   Columns: {list(df.columns)}")
        
        # Check for expected columns
        level = file_path.name.split('_')[1]  # Extract level from filename
        expected_cols = [level, 'domain', f'{level}_genome_count', f'{level}_species_count']
        
        missing_cols = [col for col in expected_cols if col not in df.columns]
        if missing_cols:
            print(f"   ⚠️  Missing expected columns: {missing_cols}")
        else:
            print(f"   ✅ All expected columns present")
        
        # Show top entries
        print(f"\n📈 Top 5 entries by genome count:")
        if f'{level}_genome_count' in df.columns:
            top_entries = df.nlargest(5, f'{level}_genome_count')
            for _, row in top_entries.iterrows():
                genome_count = row[f'{level}_genome_count']
                species_count = row[f'{level}_species_count'] if f'{level}_species_count' in row else 'N/A'
                print(f"   {row[level]} ({row['domain']}): {genome_count:,} genomes, {species_count} species")
        
        # Summary statistics
        if f'{level}_genome_count' in df.columns and f'{level}_species_count' in df.columns:
            total_genomes = df[f'{level}_genome_count'].sum()
            total_species = df[f'{level}_species_count'].sum()
            print(f"\n📊 Summary statistics:")
            print(f"   Total genomes: {total_genomes:,}")
            print(f"   Total species: {total_species:,}")
            print(f"   Average genomes per {level}: {total_genomes / len(df):.1f}")
            print(f"   Average species per {level}: {total_species / len(df):.1f}")
            
            # Genome to species ratio
            df_with_species = df[df[f'{level}_species_count'] > 0]
            if len(df_with_species) > 0:
                avg_ratio = (df_with_species[f'{level}_genome_count'] / df_with_species[f'{level}_species_count']).mean()
                print(f"   Average genomes per species: {avg_ratio:.1f}")
        
    except Exception as e:
        print(f"❌ Error analyzing file: {e}")

def compare_with_original():
    """Compare comprehensive output with original separate files."""
    print(f"🔍 Comparing comprehensive output with original files")
    print("=" * 60)
    
    script_dir = Path(__file__).resolve().parent
    new_output_dir = script_dir.parent  # ncbi_parse directory
    original_dir = script_dir.parent / "csv_ncbi"
    
    level = "genus"  # Test with genus
    
    # File paths
    comprehensive_file = new_output_dir / f"ncbi_{level}_comprehensive_counts.csv"
    original_genome_file = original_dir / f"ncbi_{level}_counts.csv"
    original_species_file = original_dir / f"ncbi_{level}_species_counts.csv"
    
    print(f"📁 File locations:")
    print(f"   Comprehensive: {comprehensive_file}")
    print(f"   Original genome counts: {original_genome_file}")
    print(f"   Original species counts: {original_species_file}")
    
    # Check if files exist
    files_exist = {
        'comprehensive': comprehensive_file.exists(),
        'original_genome': original_genome_file.exists(),
        'original_species': original_species_file.exists()
    }
    
    print(f"\n📋 File existence:")
    for file_type, exists in files_exist.items():
        status = "✅ Exists" if exists else "❌ Missing"
        print(f"   {file_type}: {status}")
    
    if not files_exist['comprehensive']:
        print(f"\n⚠️  Run unified parser first: python unified_ncbi_parser.py --level {level}")
        return False
    
    if not (files_exist['original_genome'] and files_exist['original_species']):
        print(f"\n⚠️  Original files missing - cannot compare")
        return False
    
    try:
        # Load files
        comprehensive_df = pd.read_csv(comprehensive_file)
        original_genome_df = pd.read_csv(original_genome_file)
        original_species_df = pd.read_csv(original_species_file)
        
        print(f"\n📊 Comparison results:")
        print(f"   Comprehensive file: {len(comprehensive_df):,} rows")
        print(f"   Original genome file: {len(original_genome_df):,} rows")
        print(f"   Original species file: {len(original_species_df):,} rows")
        
        # Compare totals if possible
        if f'{level}_genome_count' in comprehensive_df.columns:
            comp_total_genomes = comprehensive_df[f'{level}_genome_count'].sum()
            orig_total_genomes = original_genome_df[f'{level}_count'].sum() if f'{level}_count' in original_genome_df.columns else 0
            
            print(f"\n🔢 Total genome counts:")
            print(f"   Comprehensive: {comp_total_genomes:,}")
            print(f"   Original: {orig_total_genomes:,}")
            
            if comp_total_genomes == orig_total_genomes:
                print(f"   ✅ Genome counts match!")
            else:
                diff = abs(comp_total_genomes - orig_total_genomes)
                print(f"   ⚠️  Difference: {diff:,} genomes")
        
        if f'{level}_species_count' in comprehensive_df.columns:
            comp_total_species = comprehensive_df[f'{level}_species_count'].sum()
            orig_total_species = original_species_df[f'{level}_species_count'].sum() if f'{level}_species_count' in original_species_df.columns else 0
            
            print(f"\n🧬 Total species counts:")
            print(f"   Comprehensive: {comp_total_species:,}")
            print(f"   Original: {orig_total_species:,}")
            
            if comp_total_species == orig_total_species:
                print(f"   ✅ Species counts match!")
            else:
                diff = abs(comp_total_species - orig_total_species)
                print(f"   ⚠️  Difference: {diff:,} species")
        
        return True
        
    except Exception as e:
        print(f"❌ Error comparing files: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Test comprehensive output of unified NCBI parser")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--test-genus", action="store_true", help="Test unified parser with genus level")
    group.add_argument("--test-family", action="store_true", help="Test unified parser with family level")
    group.add_argument("--compare-with-original", action="store_true", help="Compare with original separate files")
    
    args = parser.parse_args()
    
    if args.test_genus:
        success = run_unified_parser_test("genus")
        exit(0 if success else 1)
    elif args.test_family:
        success = run_unified_parser_test("family")
        exit(0 if success else 1)
    elif args.compare_with_original:
        success = compare_with_original()
        exit(0 if success else 1)

if __name__ == "__main__":
    main()
