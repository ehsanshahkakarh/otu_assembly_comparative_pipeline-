#!/usr/bin/env python3
"""
Compare Output Locations

This script shows the difference between old and new output locations
for the genus parser and accessions scripts.

Usage:
    python compare_output_locations.py --show-locations
    python compare_output_locations.py --show-files
"""

import argparse
from pathlib import Path

def show_directory_structure():
    """Show the directory structure and output locations."""
    script_dir = Path(__file__).resolve().parent
    ncbi_parse_dir = script_dir.parent
    csv_ncbi_dir = ncbi_parse_dir / "csv_ncbi"
    
    print("📁 DIRECTORY STRUCTURE & OUTPUT LOCATIONS")
    print("=" * 60)
    
    print(f"\n🏠 Project structure:")
    print(f"   ncbi_parse/                           # Main directory")
    print(f"   ├── py_ncbi/                          # Scripts location")
    print(f"   │   ├── genus_ncbi_parser_new.py      # New parser script")
    print(f"   │   ├── genus_ncbi_accessions_new.py  # New accessions script")
    print(f"   │   └── ...                           # Other scripts")
    print(f"   ├── csv_ncbi/                         # Old output location")
    print(f"   │   ├── ncbi_genus_counts.csv         # Old outputs")
    print(f"   │   └── ...                           # (preserved)")
    print(f"   ├── ncbi_genus_counts.csv             # 🆕 New outputs (safe location)")
    print(f"   └── ncbi_genus_with_accessions_classified.csv")
    
    print(f"\n📂 Output location comparison:")
    print(f"   Old scripts output to: {csv_ncbi_dir}")
    print(f"   New scripts output to: {ncbi_parse_dir}")
    print(f"   Benefit: Safe comparison without overriding existing files")

def show_file_comparison():
    """Show comparison of files in both locations."""
    script_dir = Path(__file__).resolve().parent
    ncbi_parse_dir = script_dir.parent
    csv_ncbi_dir = ncbi_parse_dir / "csv_ncbi"
    
    print("📊 FILE COMPARISON")
    print("=" * 60)
    
    # Expected file patterns
    genus_files = [
        "ncbi_genus_counts.csv",
        "ncbi_genus_species_counts.csv",
        "ncbi_genus_with_accessions_classified.csv",
        "ncbi_genus_species_subset_with_accessions_classified.csv"
    ]
    
    print(f"\n📁 Old location ({csv_ncbi_dir}):")
    if csv_ncbi_dir.exists():
        old_files_found = 0
        for filename in genus_files:
            file_path = csv_ncbi_dir / filename
            if file_path.exists():
                size_info = f"({file_path.stat().st_size / (1024*1024):.1f} MB)" if file_path.stat().st_size > 1024*1024 else f"({file_path.stat().st_size / 1024:.1f} KB)"
                print(f"   ✅ {filename} {size_info}")
                old_files_found += 1
            else:
                print(f"   ❌ {filename} - not found")
        
        if old_files_found == 0:
            print(f"   📝 No genus files found in old location")
    else:
        print(f"   📝 Directory does not exist")
    
    print(f"\n📁 New location ({ncbi_parse_dir}):")
    new_files_found = 0
    for filename in genus_files:
        file_path = ncbi_parse_dir / filename
        if file_path.exists():
            size_info = f"({file_path.stat().st_size / (1024*1024):.1f} MB)" if file_path.stat().st_size > 1024*1024 else f"({file_path.stat().st_size / 1024:.1f} KB)"
            print(f"   ✅ {filename} {size_info}")
            new_files_found += 1
        else:
            print(f"   ❌ {filename} - not found")
    
    if new_files_found == 0:
        print(f"   📝 No genus files found in new location")
        print(f"   💡 Run new scripts to generate files:")
        print(f"      python genus_ncbi_parser_new.py")
        print(f"      python genus_ncbi_accessions_new.py")
    
    # Summary
    print(f"\n📊 Summary:")
    print(f"   Old location files: {old_files_found if csv_ncbi_dir.exists() else 0}/{len(genus_files)}")
    print(f"   New location files: {new_files_found}/{len(genus_files)}")

def show_usage_examples():
    """Show usage examples for the new scripts."""
    print("🚀 USAGE EXAMPLES")
    print("=" * 60)
    
    print(f"\n📝 Running new genus scripts:")
    print(f"   cd py_ncbi/")
    print(f"   ")
    print(f"   # Generate count files (outputs to ncbi_parse/)")
    print(f"   python genus_ncbi_parser_new.py")
    print(f"   ")
    print(f"   # Generate accession files (outputs to ncbi_parse/)")
    print(f"   python genus_ncbi_accessions_new.py")
    print(f"   ")
    print(f"   # Test both scripts")
    print(f"   python test_new_genus_scripts.py --test-both")
    
    print(f"\n📊 Comparing outputs:")
    print(f"   # Check file sizes")
    print(f"   ls -la ../csv_ncbi/ncbi_genus_*")
    print(f"   ls -la ../ncbi_genus_*")
    print(f"   ")
    print(f"   # Compare row counts")
    print(f"   wc -l ../csv_ncbi/ncbi_genus_counts.csv ../ncbi_genus_counts.csv")
    print(f"   ")
    print(f"   # Compare headers")
    print(f"   head -1 ../csv_ncbi/ncbi_genus_counts.csv")
    print(f"   head -1 ../ncbi_genus_counts.csv")
    
    print(f"\n🔧 Custom output directory:")
    print(f"   # Output to specific directory")
    print(f"   python genus_ncbi_parser_new.py --output-dir /path/to/custom/dir")
    print(f"   python genus_ncbi_accessions_new.py --output-dir /path/to/custom/dir")

def show_benefits():
    """Show benefits of the new output location."""
    print("✅ BENEFITS OF NEW OUTPUT LOCATION")
    print("=" * 60)
    
    benefits = [
        "Safe comparison - won't override existing csv_ncbi/ files",
        "Side-by-side analysis - compare old vs new outputs easily",
        "Clean separation - original files preserved, new files clearly identified",
        "Easy testing - test new scripts without risk of data loss",
        "Flexible - can still specify custom output directory if needed",
        "Organized - keeps new experimental outputs separate from production files"
    ]
    
    for i, benefit in enumerate(benefits, 1):
        print(f"   {i}. {benefit}")
    
    print(f"\n💡 Workflow recommendation:")
    print(f"   1. Run new scripts (outputs to ncbi_parse/)")
    print(f"   2. Compare with existing files in csv_ncbi/")
    print(f"   3. Validate results and performance")
    print(f"   4. Once satisfied, replace old scripts")
    print(f"   5. Optionally move new outputs to csv_ncbi/ or update downstream scripts")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Compare old and new output locations")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--show-locations", action="store_true", help="Show directory structure and locations")
    group.add_argument("--show-files", action="store_true", help="Show file comparison between locations")
    group.add_argument("--show-usage", action="store_true", help="Show usage examples")
    group.add_argument("--show-benefits", action="store_true", help="Show benefits of new location")
    group.add_argument("--show-all", action="store_true", help="Show all information")
    
    args = parser.parse_args()
    
    if args.show_all:
        show_directory_structure()
        print("\n")
        show_file_comparison()
        print("\n")
        show_usage_examples()
        print("\n")
        show_benefits()
    elif args.show_locations:
        show_directory_structure()
    elif args.show_files:
        show_file_comparison()
    elif args.show_usage:
        show_usage_examples()
    elif args.show_benefits:
        show_benefits()

if __name__ == "__main__":
    main()
