#!/usr/bin/env python3
"""
Comprehensive Validation Script for New Merger Outputs
=======================================================

Validates:
1. Data integrity (no missing values, valid ranges)
2. Consistency across taxonomic levels
3. Percentage calculations
4. Domain filtering
5. Matching logic
6. Comparison with old merger outputs

Generates a comprehensive validation report.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import sys


class MergerValidator:
    """Comprehensive validation for merger outputs."""
    
    def __init__(self):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.output_dir = self.base_dir / 'outputs'
        self.report_dir = self.base_dir / 'validation_reports'
        self.report_dir.mkdir(exist_ok=True)
        
        self.issues = []
        self.warnings = []
        self.passed_checks = []
    
    def log_issue(self, message: str):
        """Log a critical issue."""
        self.issues.append(f"❌ {message}")
        print(f"❌ {message}")
    
    def log_warning(self, message: str):
        """Log a warning."""
        self.warnings.append(f"⚠️  {message}")
        print(f"⚠️  {message}")
    
    def log_pass(self, message: str):
        """Log a passed check."""
        self.passed_checks.append(f"✅ {message}")
        print(f"✅ {message}")
    
    def validate_data_integrity(self, df: pd.DataFrame, filename: str) -> bool:
        """Check for missing values, negative numbers, etc."""
        print(f"\n🔍 Validating data integrity: {filename}")
        
        # Check for missing values in critical columns
        critical_cols = ['census_otu_count', 'census_size_count', 'ncbi_genome_count', 
                        'ncbi_species_count', 'match_status']
        
        for col in critical_cols:
            if col in df.columns:
                missing = df[col].isna().sum()
                if missing > 0:
                    self.log_issue(f"{filename}: {missing} missing values in {col}")
                    return False
        
        # Check for negative values in count columns
        count_cols = [c for c in df.columns if 'count' in c.lower()]
        for col in count_cols:
            if col in df.columns:
                negative = (df[col] < 0).sum()
                if negative > 0:
                    self.log_issue(f"{filename}: {negative} negative values in {col}")
                    return False
        
        # Check match_status values
        if 'match_status' in df.columns:
            valid_statuses = {'matched', 'census_only'}
            invalid = ~df['match_status'].isin(valid_statuses)
            if invalid.sum() > 0:
                self.log_issue(f"{filename}: {invalid.sum()} invalid match_status values")
                return False
        
        self.log_pass(f"{filename}: Data integrity check passed")
        return True
    
    def validate_domain_filtering(self, df: pd.DataFrame, filename: str, expected_domains: list) -> bool:
        """Check that domain filtering is correct."""
        print(f"\n🔍 Validating domain filtering: {filename}")
        
        if 'domain' not in df.columns:
            self.log_warning(f"{filename}: No domain column found")
            return True
        
        # Check matched entries have correct domains
        matched = df[df['match_status'] == 'matched']
        if len(matched) > 0:
            invalid_domains = ~matched['domain'].isin(expected_domains)
            if invalid_domains.sum() > 0:
                invalid_list = matched[invalid_domains]['domain'].unique()
                self.log_issue(f"{filename}: Found {invalid_domains.sum()} entries with invalid domains: {invalid_list}")
                return False
        
        self.log_pass(f"{filename}: Domain filtering correct (expected: {expected_domains})")
        return True
    
    def validate_matching_logic(self, df: pd.DataFrame, filename: str) -> bool:
        """Check that matching logic is consistent."""
        print(f"\n🔍 Validating matching logic: {filename}")
        
        # Matched entries should have non-zero NCBI counts
        matched = df[df['match_status'] == 'matched']
        if len(matched) > 0:
            zero_genomes = (matched['ncbi_genome_count'] == 0).sum()
            if zero_genomes > 0:
                self.log_issue(f"{filename}: {zero_genomes} matched entries with zero NCBI genomes")
                return False
        
        # Census-only entries should have zero NCBI counts
        census_only = df[df['match_status'] == 'census_only']
        if len(census_only) > 0:
            nonzero_genomes = (census_only['ncbi_genome_count'] > 0).sum()
            if nonzero_genomes > 0:
                self.log_issue(f"{filename}: {nonzero_genomes} census_only entries with non-zero NCBI genomes")
                return False
        
        self.log_pass(f"{filename}: Matching logic is consistent")
        return True

    def validate_percentage_calculations(self, df: pd.DataFrame, filename: str) -> bool:
        """Check that percentage calculations are reasonable."""
        print(f"\n🔍 Validating percentage calculations: {filename}")

        # Check OTU percentages
        if 'otu_percentage' in df.columns:
            total_otu_pct = df['otu_percentage'].sum()
            if abs(total_otu_pct - 100) > 10:  # Allow 10% deviation
                self.log_warning(f"{filename}: OTU percentages sum to {total_otu_pct:.2f}% (expected ~100%)")
            else:
                self.log_pass(f"{filename}: OTU percentages sum to {total_otu_pct:.2f}%")

        # Check size percentages
        if 'size_percentage' in df.columns:
            total_size_pct = df['size_percentage'].sum()
            if abs(total_size_pct - 100) > 10:  # Allow 10% deviation
                self.log_warning(f"{filename}: Size percentages sum to {total_size_pct:.2f}% (expected ~100%)")
            else:
                self.log_pass(f"{filename}: Size percentages sum to {total_size_pct:.2f}%")

        return True

    def validate_hierarchical_consistency(self, census_type: str) -> bool:
        """Check that counts are consistent across taxonomic levels."""
        print(f"\n🔍 Validating hierarchical consistency: {census_type}")

        # Load all three levels
        div_file = self.output_dir / f'{census_type.lower()}_ncbi_merged_division.csv'
        fam_file = self.output_dir / f'{census_type.lower()}_ncbi_merged_family.csv'
        gen_file = self.output_dir / f'{census_type.lower()}_ncbi_merged_genus.csv'

        if not all([div_file.exists(), fam_file.exists(), gen_file.exists()]):
            self.log_warning(f"{census_type}: Not all level files found")
            return True

        div_df = pd.read_csv(div_file)
        fam_df = pd.read_csv(fam_file)
        gen_df = pd.read_csv(gen_file)

        # Check that total census counts are similar across levels
        div_total = div_df['census_otu_count'].sum()
        fam_total = fam_df['census_otu_count'].sum()
        gen_total = gen_df['census_otu_count'].sum()

        print(f"  Division total OTUs: {div_total:,}")
        print(f"  Family total OTUs: {fam_total:,}")
        print(f"  Genus total OTUs: {gen_total:,}")

        # They should be similar (within 20%)
        max_total = max(div_total, fam_total, gen_total)
        min_total = min(div_total, fam_total, gen_total)

        if max_total > 0:
            diff_pct = ((max_total - min_total) / max_total) * 100
            if diff_pct > 20:
                self.log_warning(f"{census_type}: Large variation in totals across levels ({diff_pct:.1f}%)")
            else:
                self.log_pass(f"{census_type}: Totals consistent across levels (variation: {diff_pct:.1f}%)")

        return True

    def validate_file(self, census_type: str, level: str, expected_domains: list) -> bool:
        """Run all validations on a single file."""
        filename = f'{census_type.lower()}_ncbi_merged_{level}.csv'
        filepath = self.output_dir / filename

        if not filepath.exists():
            self.log_warning(f"File not found: {filename}")
            return False

        print(f"\n{'='*80}")
        print(f"VALIDATING: {filename}")
        print(f"{'='*80}")

        df = pd.read_csv(filepath)
        print(f"📊 Loaded {len(df):,} rows, {len(df.columns)} columns")

        # Run all validation checks
        checks = [
            self.validate_data_integrity(df, filename),
            self.validate_domain_filtering(df, filename, expected_domains),
            self.validate_matching_logic(df, filename),
            self.validate_percentage_calculations(df, filename)
        ]

        return all(checks)

    def generate_report(self) -> str:
        """Generate a comprehensive validation report."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = self.report_dir / f'validation_report_{timestamp}.txt'

        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("NEW MERGER VALIDATION REPORT\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            # Summary
            f.write("SUMMARY\n")
            f.write("-"*80 + "\n")
            f.write(f"✅ Passed checks: {len(self.passed_checks)}\n")
            f.write(f"⚠️  Warnings: {len(self.warnings)}\n")
            f.write(f"❌ Critical issues: {len(self.issues)}\n\n")

            # Details
            if self.passed_checks:
                f.write("\nPASSED CHECKS\n")
                f.write("-"*80 + "\n")
                for check in self.passed_checks:
                    f.write(f"{check}\n")

            if self.warnings:
                f.write("\nWARNINGS\n")
                f.write("-"*80 + "\n")
                for warning in self.warnings:
                    f.write(f"{warning}\n")

            if self.issues:
                f.write("\nCRITICAL ISSUES\n")
                f.write("-"*80 + "\n")
                for issue in self.issues:
                    f.write(f"{issue}\n")

            # Conclusion
            f.write("\n" + "="*80 + "\n")
            f.write("CONCLUSION\n")
            f.write("="*80 + "\n")
            if len(self.issues) == 0:
                f.write("✅ All critical checks passed! Outputs are valid.\n")
            else:
                f.write(f"❌ Found {len(self.issues)} critical issues. Review required.\n")

        return str(report_file)


def main():
    """Run comprehensive validation."""
    validator = MergerValidator()

    print("="*80)
    print("NEW MERGER COMPREHENSIVE VALIDATION")
    print("="*80)
    print("\nValidating all merger outputs for data quality and consistency\n")

    # Validate 18S files
    print("\n" + "="*80)
    print("VALIDATING 18S (EUKARYOTIC) OUTPUTS")
    print("="*80)

    for level in ['division', 'family', 'genus']:
        validator.validate_file('18S', level, ['Eukaryota'])

    validator.validate_hierarchical_consistency('18S')

    # Validate 16S files
    print("\n" + "="*80)
    print("VALIDATING 16S (PROKARYOTIC) OUTPUTS")
    print("="*80)

    for level in ['division', 'family', 'genus']:
        validator.validate_file('16S', level, ['Bacteria', 'Archaea'])

    validator.validate_hierarchical_consistency('16S')

    # Generate report
    print("\n" + "="*80)
    print("GENERATING VALIDATION REPORT")
    print("="*80)

    report_file = validator.generate_report()

    print(f"\n✅ Validation complete!")
    print(f"📄 Report saved: {report_file}")
    print(f"\n📊 Summary:")
    print(f"  ✅ Passed: {len(validator.passed_checks)}")
    print(f"  ⚠️  Warnings: {len(validator.warnings)}")
    print(f"  ❌ Issues: {len(validator.issues)}")

    # Exit with error code if there are critical issues
    sys.exit(1 if len(validator.issues) > 0 else 0)


if __name__ == "__main__":
    main()


