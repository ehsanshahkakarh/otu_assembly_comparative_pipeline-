#!/usr/bin/env python3
"""
Example Usage Script for Unified NCBI Parser

This script demonstrates how to use the unified parser with genome source classification
to replace all the redundant individual parsers.

Usage:
    python run_unified_parser_examples.py --help
    python run_unified_parser_examples.py --run-all
    python run_unified_parser_examples.py --level genus --species-subset
"""

import subprocess
import sys
import argparse
from pathlib import Path
from typing import List, Optional

def run_command(cmd: List[str], description: str) -> bool:
    """Run a command and return success status."""
    print(f"\n🚀 {description}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with return code {e.returncode}")
        return False
    except FileNotFoundError:
        print(f"❌ Command not found: {cmd[0]}")
        return False

def run_unified_parser(level: str, species_subset: bool = False, 
                      input_dir: Optional[str] = None, 
                      output_dir: Optional[str] = None) -> bool:
    """Run the unified parser for a specific taxonomic level."""
    
    cmd = ["python", "unified_ncbi_parser.py", "--level", level]
    
    if species_subset:
        cmd.append("--species-subset")
    
    if input_dir:
        cmd.extend(["--input-dir", input_dir])
    
    if output_dir:
        cmd.extend(["--output-dir", output_dir])
    
    subset_text = " (species subset)" if species_subset else ""
    description = f"Running unified parser for {level} level{subset_text}"
    
    return run_command(cmd, description)

def compare_with_old_parser(level: str, species_subset: bool = False) -> bool:
    """Compare unified parser output with old parser output (if available)."""
    
    # Determine old parser script name
    if species_subset:
        old_script = f"{level}_ncbi_parser_species_subset.py"
        output_suffix = "_species_counts"
    else:
        old_script = f"{level}_ncbi_parser.py"
        output_suffix = "_counts"
    
    old_script_path = Path(old_script)
    
    if not old_script_path.exists():
        print(f"⚠️  Old parser {old_script} not found - skipping comparison")
        return True
    
    print(f"\n🔍 Comparing unified parser with old parser for {level} level")
    
    # Run old parser
    old_success = run_command(
        ["python", old_script], 
        f"Running old {level} parser for comparison"
    )
    
    if not old_success:
        print(f"❌ Old parser failed - cannot compare")
        return False
    
    # Compare output files
    script_dir = Path(__file__).resolve().parent
    csv_dir = script_dir.parent / "csv_ncbi"
    
    old_file = csv_dir / f"ncbi_{level}{output_suffix}.csv"
    new_file = csv_dir / f"ncbi_{level}{output_suffix}.csv"  # Same name
    
    if old_file.exists() and new_file.exists():
        print(f"📊 Output files exist - manual comparison recommended:")
        print(f"   Old: {old_file}")
        print(f"   New: {new_file}")
        return True
    else:
        print(f"⚠️  Output files not found for comparison")
        return False

def run_all_levels(species_subset: bool = False) -> None:
    """Run unified parser for all taxonomic levels."""
    
    levels = ["phylum", "family", "genus"]
    results = {}
    
    print(f"🎯 Running unified parser for all levels{'(species subset)' if species_subset else ''}")
    
    for level in levels:
        success = run_unified_parser(level, species_subset)
        results[level] = success
    
    # Summary
    print(f"\n📋 Summary of unified parser runs:")
    for level, success in results.items():
        status = "✅ Success" if success else "❌ Failed"
        subset_text = " (species subset)" if species_subset else ""
        print(f"   {level}{subset_text}: {status}")
    
    successful_runs = sum(results.values())
    print(f"\n🎉 Completed {successful_runs}/{len(levels)} parser runs successfully")

def show_output_files() -> None:
    """Show generated output files."""

    script_dir = Path(__file__).resolve().parent
    # Check both locations: new default (ncbi_parse) and old location (csv_ncbi)
    ncbi_parse_dir = script_dir.parent  # ncbi_parse directory
    csv_ncbi_dir = script_dir.parent / "csv_ncbi"  # old csv_ncbi directory

    # Determine which directory to show
    if ncbi_parse_dir.exists() and list(ncbi_parse_dir.glob("ncbi_*_counts.csv")):
        output_dir = ncbi_parse_dir
        location_note = "(unified parser output - safe for comparison)"
    elif csv_ncbi_dir.exists():
        output_dir = csv_ncbi_dir
        location_note = "(original csv_ncbi directory)"
    else:
        print(f"❌ No output directories found")
        print(f"   Checked: {ncbi_parse_dir}")
        print(f"   Checked: {csv_ncbi_dir}")
        return

    print(f"\n📁 Output files in {output_dir} {location_note}:")
    
    # Group files by type
    count_files = list(output_dir.glob("*_counts.csv"))
    species_count_files = list(output_dir.glob("*_species_counts.csv"))
    accession_files = list(output_dir.glob("*_with_accessions_classified.csv"))
    species_accession_files = list(output_dir.glob("*_species_subset_with_accessions_classified.csv"))
    
    if count_files:
        print(f"   📊 Count files ({len(count_files)}):")
        for f in sorted(count_files):
            size_kb = f.stat().st_size / 1024
            print(f"      {f.name} ({size_kb:.1f} KB)")
    
    if species_count_files:
        print(f"   🧬 Species count files ({len(species_count_files)}):")
        for f in sorted(species_count_files):
            size_kb = f.stat().st_size / 1024
            print(f"      {f.name} ({size_kb:.1f} KB)")
    
    if accession_files:
        print(f"   🔗 Classified accession files ({len(accession_files)}):")
        for f in sorted(accession_files):
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"      {f.name} ({size_mb:.1f} MB)")

    if species_accession_files:
        print(f"   🧬 Classified species subset accession files ({len(species_accession_files)}):")
        for f in sorted(species_accession_files):
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"      {f.name} ({size_mb:.1f} MB)")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Example usage of unified NCBI parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run for single level
  python run_unified_parser_examples.py --level genus
  
  # Run with species subset
  python run_unified_parser_examples.py --level family --species-subset
  
  # Run all levels
  python run_unified_parser_examples.py --run-all
  
  # Run all levels with species subset
  python run_unified_parser_examples.py --run-all --species-subset
  
  # Show output files
  python run_unified_parser_examples.py --show-files
  
  # Compare with old parser
  python run_unified_parser_examples.py --level genus --compare
        """
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--level", choices=["phylum", "family", "genus"], 
                      help="Run parser for specific taxonomic level")
    group.add_argument("--run-all", action="store_true", 
                      help="Run parser for all taxonomic levels")
    group.add_argument("--show-files", action="store_true", 
                      help="Show generated output files")
    
    parser.add_argument("--species-subset", action="store_true", 
                       help="Use species subset mode (one genome per species)")
    parser.add_argument("--compare", action="store_true", 
                       help="Compare with old parser output (requires --level)")
    parser.add_argument("--input-dir", help="Input directory for assembly files")
    parser.add_argument("--output-dir", help="Output directory for results")
    
    args = parser.parse_args()
    
    # Check if unified parser exists
    unified_parser = Path("unified_ncbi_parser.py")
    if not unified_parser.exists():
        print(f"❌ Unified parser not found: {unified_parser}")
        print("Make sure you're running this script from the correct directory")
        sys.exit(1)
    
    if args.show_files:
        show_output_files()
    elif args.run_all:
        run_all_levels(args.species_subset)
        show_output_files()
    elif args.level:
        if args.compare:
            # Run unified parser first
            success = run_unified_parser(args.level, args.species_subset, 
                                       args.input_dir, args.output_dir)
            if success:
                compare_with_old_parser(args.level, args.species_subset)
        else:
            run_unified_parser(args.level, args.species_subset, 
                             args.input_dir, args.output_dir)
        show_output_files()

if __name__ == "__main__":
    main()
