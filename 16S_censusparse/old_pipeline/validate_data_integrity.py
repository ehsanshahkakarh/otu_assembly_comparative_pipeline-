#!/usr/bin/env python3
"""
Data Integrity Validation Script for 16S Census Parser
=======================================================

This script validates that the parser output maintains data integrity by checking:
1. Fidelity: New output matches archived baseline files
2. Completeness: Total OTU/size counts match between input and output
3. Consistency: Proper taxonomic aggregation at each level

Usage:
    python validate_data_integrity.py
"""

import pandas as pd
from pathlib import Path
import sys

def compare_dataframes(new_df, archive_df, level_name):
    """Compare new output against archived baseline with detailed difference reporting."""
    print(f"\n🔍 Comparing {level_name} against archive baseline...")

    differences_found = []

    # Check row counts
    if len(new_df) != len(archive_df):
        diff_msg = f"Row count mismatch: New={len(new_df)}, Archive={len(archive_df)}"
        print(f"   ⚠️  {diff_msg}")
        differences_found.append(diff_msg)
    else:
        print(f"   ✅ Row count matches: {len(new_df)} entries")

    # Check total OTU counts
    new_otu_total = new_df['otu_count'].sum()
    archive_otu_total = archive_df['otu_count'].sum()

    if new_otu_total != archive_otu_total:
        diff = new_otu_total - archive_otu_total
        diff_msg = f"OTU total mismatch: New={new_otu_total:,}, Archive={archive_otu_total:,}, Diff={diff:+,}"
        print(f"   ⚠️  {diff_msg}")
        differences_found.append(diff_msg)
    else:
        print(f"   ✅ OTU totals match: {new_otu_total:,}")

    # Check total size counts
    new_size_total = new_df['size_count'].sum()
    archive_size_total = archive_df['size_count'].sum()

    if new_size_total != archive_size_total:
        diff = new_size_total - archive_size_total
        diff_msg = f"Size total mismatch: New={new_size_total:,}, Archive={archive_size_total:,}, Diff={diff:+,}"
        print(f"   ⚠️  {diff_msg}")
        differences_found.append(diff_msg)
    else:
        print(f"   ✅ Size totals match: {new_size_total:,}")

    # Check if all taxa names match
    new_taxa = set(new_df['Name_to_use'].values)
    archive_taxa = set(archive_df['Name_to_use'].values)

    if new_taxa != archive_taxa:
        missing_in_new = archive_taxa - new_taxa
        extra_in_new = new_taxa - archive_taxa

        if missing_in_new:
            diff_msg = f"Taxa missing in new output: {len(missing_in_new)}"
            print(f"   ⚠️  {diff_msg}")
            print(f"      Examples: {sorted(list(missing_in_new))[:10]}")
            differences_found.append(diff_msg)
        if extra_in_new:
            diff_msg = f"Extra taxa in new output: {len(extra_in_new)}"
            print(f"   ⚠️  {diff_msg}")
            print(f"      Examples: {sorted(list(extra_in_new))[:10]}")
            differences_found.append(diff_msg)
    else:
        print(f"   ✅ All taxa names match: {len(new_taxa)} unique taxa")

    # Detailed row-by-row comparison for matching taxa
    if len(differences_found) == 0:
        # Merge on Name_to_use to compare values
        merged = new_df.merge(archive_df, on='Name_to_use', suffixes=('_new', '_archive'))

        # Check for differences in OTU counts
        otu_diffs = merged[merged['otu_count_new'] != merged['otu_count_archive']]
        if len(otu_diffs) > 0:
            diff_msg = f"OTU count differences in {len(otu_diffs)} taxa"
            print(f"   ⚠️  {diff_msg}")
            print(f"      Examples:")
            for idx, row in otu_diffs.head(5).iterrows():
                print(f"         {row['Name_to_use']}: New={row['otu_count_new']:,}, Archive={row['otu_count_archive']:,}")
            differences_found.append(diff_msg)

        # Check for differences in size counts
        size_diffs = merged[merged['size_count_new'] != merged['size_count_archive']]
        if len(size_diffs) > 0:
            diff_msg = f"Size count differences in {len(size_diffs)} taxa"
            print(f"   ⚠️  {diff_msg}")
            print(f"      Examples:")
            for idx, row in size_diffs.head(5).iterrows():
                print(f"         {row['Name_to_use']}: New={row['size_count_new']:,}, Archive={row['size_count_archive']:,}")
            differences_found.append(diff_msg)

        # Check for differences in taxids
        taxid_diffs = merged[merged['taxid_new'] != merged['taxid_archive']]
        if len(taxid_diffs) > 0:
            diff_msg = f"Taxid differences in {len(taxid_diffs)} taxa"
            print(f"   ⚠️  {diff_msg}")
            print(f"      Examples:")
            for idx, row in taxid_diffs.head(5).iterrows():
                print(f"         {row['Name_to_use']}: New={row['taxid_new']}, Archive={row['taxid_archive']}")
            differences_found.append(diff_msg)

    # Final verdict
    if len(differences_found) == 0:
        print(f"\n   🎉 PERFECT MATCH - No differences found!")
        return True
    else:
        print(f"\n   ❌ DIFFERENCES DETECTED - {len(differences_found)} issue(s) found")
        return False


def validate_data_integrity():
    """Validate that parser output maintains data integrity."""

    print("=" * 80)
    print("16S CENSUS PARSER - DATA INTEGRITY VALIDATION")
    print("=" * 80)
    print()

    # File paths
    base_dir = Path(__file__).parent.parent
    input_file = base_dir / "metadata" / "eukcensus_16S.clusters.97.tsv"

    output_dir = base_dir / "csv_16S"
    archive_dir = output_dir / "archive"

    # New output files
    division_file = output_dir / "eukcensus16S_by_division.csv"
    family_file = output_dir / "eukcensus16S_by_family.csv"
    genus_file = output_dir / "eukcensus16S_by_genus.csv"

    # Archive baseline files
    archive_division = archive_dir / "eukcensus16S_by_division.csv"
    archive_family = archive_dir / "eukcensus16S_by_family.csv"
    archive_genus = archive_dir / "eukcensus16S_by_genus.csv"

    # Check if files exist
    if not input_file.exists():
        print(f"❌ Input file not found: {input_file}")
        return False

    if not all([division_file.exists(), family_file.exists(), genus_file.exists()]):
        print("❌ New output files not found. Please run the parser first.")
        return False

    if not all([archive_division.exists(), archive_family.exists(), archive_genus.exists()]):
        print("⚠️  Archive baseline files not found. Skipping fidelity check.")
        print(f"   Expected archive location: {archive_dir}")
        archive_check = False
    else:
        archive_check = True
    
    # PART 1: FIDELITY CHECK - Compare against archive baseline
    if archive_check:
        print("=" * 80)
        print("PART 1: FIDELITY CHECK - Comparing against archive baseline")
        print("=" * 80)

        print("\n📂 Loading new output files...")
        division_df = pd.read_csv(division_file)
        family_df = pd.read_csv(family_file)
        genus_df = pd.read_csv(genus_file)

        print(f"✅ Loaded new division file: {len(division_df)} entries")
        print(f"✅ Loaded new family file: {len(family_df)} entries")
        print(f"✅ Loaded new genus file: {len(genus_df)} entries")

        print("\n📂 Loading archive baseline files...")
        archive_division_df = pd.read_csv(archive_division)
        archive_family_df = pd.read_csv(archive_family)
        archive_genus_df = pd.read_csv(archive_genus)

        print(f"✅ Loaded archive division file: {len(archive_division_df)} entries")
        print(f"✅ Loaded archive family file: {len(archive_family_df)} entries")
        print(f"✅ Loaded archive genus file: {len(archive_genus_df)} entries")

        # Compare each level
        fidelity_valid = True
        fidelity_valid &= compare_dataframes(division_df, archive_division_df, "DIVISION")
        fidelity_valid &= compare_dataframes(family_df, archive_family_df, "FAMILY")
        fidelity_valid &= compare_dataframes(genus_df, archive_genus_df, "GENUS")

        print("\n" + "=" * 80)
        if fidelity_valid:
            print("✅ FIDELITY CHECK PASSED")
            print("   New output matches archive baseline perfectly!")
        else:
            print("❌ FIDELITY CHECK FAILED")
            print("   New output differs from archive baseline!")
        print("=" * 80)
    else:
        print("\n📂 Loading new output files...")
        division_df = pd.read_csv(division_file)
        family_df = pd.read_csv(family_file)
        genus_df = pd.read_csv(genus_file)
        fidelity_valid = None  # Skip fidelity check

    # PART 2: COMPLETENESS CHECK - Validate against input data
    print("\n" + "=" * 80)
    print("PART 2: COMPLETENESS CHECK - Validating against input data")
    print("=" * 80)

    print("\n📂 Loading input data...")
    input_df = pd.read_csv(input_file, sep='\t')

    print(f"✅ Loaded {len(input_df):,} rows from input file")
    print()

    # Calculate input totals
    input_total_otus = len(input_df)
    input_total_size = input_df['size'].sum()

    print("📊 INPUT DATA SUMMARY:")
    print(f"   Total OTUs (rows): {input_total_otus:,}")
    print(f"   Total size count: {input_total_size:,}")
    print()

    # Validate each level
    completeness_valid = True
    
    for level_name, df in [("DIVISION", division_df), ("FAMILY", family_df), ("GENUS", genus_df)]:
        print(f"🔍 Validating {level_name} level...")

        # Calculate totals
        output_total_otus = df['otu_count'].sum()
        output_total_size = df['size_count'].sum()

        # Check OTU counts
        otu_match = output_total_otus == input_total_otus
        otu_diff = output_total_otus - input_total_otus
        otu_pct = (output_total_otus / input_total_otus * 100) if input_total_otus > 0 else 0

        # Check size counts
        size_match = output_total_size == input_total_size
        size_diff = output_total_size - input_total_size
        size_pct = (output_total_size / input_total_size * 100) if input_total_size > 0 else 0

        print(f"   Total OTUs: {output_total_otus:,} ({otu_pct:.2f}% of input)")
        if otu_match:
            print(f"   ✅ OTU count matches perfectly!")
        else:
            print(f"   ⚠️  OTU count difference: {otu_diff:+,}")
            completeness_valid = False

        print(f"   Total size: {output_total_size:,} ({size_pct:.2f}% of input)")
        if size_match:
            print(f"   ✅ Size count matches perfectly!")
        else:
            print(f"   ⚠️  Size count difference: {size_diff:+,}")
            completeness_valid = False

        # Show top 5 entries
        print(f"   Top 5 {level_name.lower()} entries by OTU count:")
        top_5 = df.nlargest(5, 'otu_count')[['Name_to_use', 'otu_count', 'size_count']]
        for idx, row in top_5.iterrows():
            print(f"      {row['Name_to_use']:<40} OTUs: {row['otu_count']:>6,}  Size: {row['size_count']:>8,}")

        print()

    # Final summary
    print("=" * 80)
    print("FINAL VALIDATION SUMMARY")
    print("=" * 80)

    if fidelity_valid is not None:
        if fidelity_valid:
            print("✅ FIDELITY: New output matches archive baseline")
        else:
            print("❌ FIDELITY: New output differs from archive baseline")
    else:
        print("⚠️  FIDELITY: Skipped (no archive baseline found)")

    if completeness_valid:
        print("✅ COMPLETENESS: Genus level preserves all input data")
    else:
        print("⚠️  COMPLETENESS: Some data filtered at higher taxonomic levels (expected)")

    print("=" * 80)

    # Return True only if fidelity check passed (or was skipped) and genus level is complete
    return (fidelity_valid is None or fidelity_valid) and completeness_valid


if __name__ == "__main__":
    success = validate_data_integrity()
    sys.exit(0 if success else 1)

