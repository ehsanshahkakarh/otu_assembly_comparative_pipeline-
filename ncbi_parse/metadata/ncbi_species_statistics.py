#!/usr/bin/env python3
"""
NCBI Species Assignment Statistics Script (Updated: 2025-01-13)
================================================================

Analyzes the NCBI assembly summary file to provide comprehensive statistics on:
1. Total species count in the database
2. Species assignment success/failure rates
3. Taxonomic distribution of unassigned species
4. Domain-specific species statistics
5. Assembly level vs species assignment correlation

This script helps identify data quality issues and understand the scope
of species-level taxonomic coverage in the NCBI database.

Input: 00assembly_summary_genbank.txt (NCBI assembly summary file)
Output: Detailed statistics report and CSV files with breakdowns

Author: Species Statistics Team
Date: 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from collections import Counter
from tqdm import tqdm

class NCBISpeciesStatistics:
    """Comprehensive NCBI species assignment statistics analyzer."""
    
    def __init__(self, input_dir=None, output_dir=None):
        self.script_dir = Path(__file__).resolve().parent
        
        # Set up directories
        self.input_dir = Path(input_dir) if input_dir else self.script_dir
        self.output_dir = Path(output_dir) if output_dir else self.script_dir
        self.output_dir.mkdir(exist_ok=True)
        
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
            self.input_dir / "assembly_summary_genbank.txt",
            self.script_dir / "00assembly_summary_genbank.txt",
            self.script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
        ]
        
        for file_path in possible_files:
            if file_path.exists():
                return file_path
        return None
    
    def load_assembly_data(self):
        """Load and preprocess the assembly data."""
        print("📂 Loading NCBI assembly data...")

        # Read the file properly, skipping comment lines
        with open(self.assembly_file, 'r') as f:
            lines = f.readlines()

        # Find the header line (starts with #assembly_accession)
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#assembly_accession'):
                header_line = line.strip().lstrip('#')
                data_start = i + 1
                break

        if header_line is None:
            raise ValueError("Could not find header line in assembly file")

        # Write clean data to temporary file for pandas
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp_file:
            tmp_file.write(header_line + '\n')
            for line in lines[data_start:]:
                if not line.startswith('#'):
                    tmp_file.write(line)
            tmp_filename = tmp_file.name

        # Read with pandas
        df = pd.read_csv(tmp_filename, sep='\t', low_memory=False)

        # Clean up temp file
        Path(tmp_filename).unlink()

        print(f"📊 Loaded {len(df):,} total assembly records")
        print(f"📋 Key columns: taxid, species_taxid, organism_name, assembly_level")

        # Basic data cleaning
        df['taxid'] = pd.to_numeric(df['taxid'], errors='coerce')
        df['species_taxid'] = pd.to_numeric(df['species_taxid'], errors='coerce')
        
        # Add derived columns for analysis
        df['has_taxid'] = df['taxid'].notna()
        df['has_species_taxid'] = df['species_taxid'].notna()
        df['species_assigned'] = df['has_species_taxid'] & (df['species_taxid'] > 0)
        
        # Classify assignment status
        df['assignment_status'] = 'unknown'
        df.loc[df['species_assigned'], 'assignment_status'] = 'species_assigned'
        df.loc[df['has_taxid'] & ~df['species_assigned'], 'assignment_status'] = 'taxid_only'
        df.loc[~df['has_taxid'], 'assignment_status'] = 'no_taxid'
        
        return df
    
    def analyze_overall_statistics(self, df):
        """Generate overall species assignment statistics."""
        print("📊 Analyzing overall statistics...")
        
        total_assemblies = len(df)
        
        # Basic counts
        has_taxid = df['has_taxid'].sum()
        has_species_taxid = df['has_species_taxid'].sum()
        species_assigned = df['species_assigned'].sum()
        
        # Unique species count
        unique_species = df[df['species_assigned']]['species_taxid'].nunique()
        
        # Assignment status breakdown
        status_counts = df['assignment_status'].value_counts()
        
        stats = {
            'total_assemblies': total_assemblies,
            'has_taxid': has_taxid,
            'has_species_taxid': has_species_taxid,
            'species_assigned': species_assigned,
            'unique_species': unique_species,
            'no_taxid': status_counts.get('no_taxid', 0),
            'taxid_only': status_counts.get('taxid_only', 0),
            'species_assigned_count': status_counts.get('species_assigned', 0)
        }
        
        # Calculate percentages
        stats['has_taxid_pct'] = (has_taxid / total_assemblies) * 100
        stats['has_species_taxid_pct'] = (has_species_taxid / total_assemblies) * 100
        stats['species_assigned_pct'] = (species_assigned / total_assemblies) * 100
        stats['no_taxid_pct'] = (stats['no_taxid'] / total_assemblies) * 100
        stats['taxid_only_pct'] = (stats['taxid_only'] / total_assemblies) * 100
        
        return stats
    
    def analyze_domain_statistics(self, df):
        """Analyze species assignment by domain (inferred from organism name)."""
        print("📊 Analyzing domain-specific statistics...")
        
        # Simple domain classification based on organism patterns
        def classify_domain(organism_name):
            if pd.isna(organism_name):
                return 'Unknown'
            
            organism_lower = str(organism_name).lower()
            
            # Virus patterns
            if any(pattern in organism_lower for pattern in ['virus', 'phage', 'viroid']):
                return 'Viruses'
            
            # Bacteria patterns (most common)
            bacteria_indicators = ['escherichia', 'salmonella', 'staphylococcus', 'streptococcus', 
                                 'bacillus', 'clostridium', 'pseudomonas', 'mycobacterium']
            if any(indicator in organism_lower for indicator in bacteria_indicators):
                return 'Bacteria'
            
            # Archaea patterns
            archaea_indicators = ['methanococcus', 'thermococcus', 'pyrococcus', 'sulfolobus']
            if any(indicator in organism_lower for indicator in archaea_indicators):
                return 'Archaea'
            
            # Eukaryote patterns
            eukaryote_indicators = ['homo', 'mus', 'drosophila', 'caenorhabditis', 'arabidopsis',
                                  'saccharomyces', 'candida', 'aspergillus', 'plasmodium']
            if any(indicator in organism_lower for indicator in eukaryote_indicators):
                return 'Eukaryota'
            
            return 'Unclassified'
        
        # Apply domain classification
        tqdm.pandas(desc="Classifying domains")
        df['inferred_domain'] = df['organism_name'].progress_apply(classify_domain)
        
        # Calculate domain statistics
        domain_stats = []
        for domain in df['inferred_domain'].unique():
            domain_df = df[df['inferred_domain'] == domain]
            total = len(domain_df)
            species_assigned = domain_df['species_assigned'].sum()
            unique_species = domain_df[domain_df['species_assigned']]['species_taxid'].nunique()
            
            domain_stats.append({
                'domain': domain,
                'total_assemblies': total,
                'species_assigned': species_assigned,
                'species_assigned_pct': (species_assigned / total) * 100 if total > 0 else 0,
                'unique_species': unique_species,
                'avg_assemblies_per_species': species_assigned / unique_species if unique_species > 0 else 0
            })
        
        return pd.DataFrame(domain_stats).sort_values('total_assemblies', ascending=False)
    
    def analyze_assembly_level_correlation(self, df):
        """Analyze correlation between assembly level and species assignment."""
        print("📊 Analyzing assembly level vs species assignment...")
        
        level_stats = []
        for level in df['assembly_level'].unique():
            if pd.isna(level):
                level = 'Unknown'
            
            level_df = df[df['assembly_level'] == level]
            total = len(level_df)
            species_assigned = level_df['species_assigned'].sum()
            
            level_stats.append({
                'assembly_level': level,
                'total_assemblies': total,
                'species_assigned': species_assigned,
                'species_assigned_pct': (species_assigned / total) * 100 if total > 0 else 0,
                'no_assignment': total - species_assigned,
                'no_assignment_pct': ((total - species_assigned) / total) * 100 if total > 0 else 0
            })
        
        return pd.DataFrame(level_stats).sort_values('total_assemblies', ascending=False)
    
    def analyze_problematic_entries(self, df):
        """Identify and analyze entries with assignment problems."""
        print("📊 Analyzing problematic entries...")
        
        # Get entries with problems
        no_taxid = df[~df['has_taxid']]
        taxid_only = df[df['has_taxid'] & ~df['species_assigned']]
        
        problems = {
            'no_taxid_count': len(no_taxid),
            'taxid_only_count': len(taxid_only),
            'no_taxid_top_organisms': no_taxid['organism_name'].value_counts().head(10).to_dict(),
            'taxid_only_top_organisms': taxid_only['organism_name'].value_counts().head(10).to_dict()
        }
        
        return problems

    def generate_report(self, df, overall_stats, domain_stats, level_stats, problems):
        """Generate comprehensive statistics report."""
        print("📝 Generating comprehensive report...")

        report_lines = []
        report_lines.append("NCBI SPECIES ASSIGNMENT STATISTICS REPORT")
        report_lines.append("=" * 60)
        report_lines.append(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Input file: {self.assembly_file.name}")
        report_lines.append("")

        # Overall Statistics
        report_lines.append("OVERALL STATISTICS")
        report_lines.append("-" * 30)
        report_lines.append(f"Total assembly records: {overall_stats['total_assemblies']:,}")
        report_lines.append(f"Records with taxid: {overall_stats['has_taxid']:,} ({overall_stats['has_taxid_pct']:.2f}%)")
        report_lines.append(f"Records with species_taxid: {overall_stats['has_species_taxid']:,} ({overall_stats['has_species_taxid_pct']:.2f}%)")
        report_lines.append(f"Records with valid species assignment: {overall_stats['species_assigned']:,} ({overall_stats['species_assigned_pct']:.2f}%)")
        report_lines.append(f"Unique species represented: {overall_stats['unique_species']:,}")
        report_lines.append("")

        # Assignment Status Breakdown
        report_lines.append("ASSIGNMENT STATUS BREAKDOWN")
        report_lines.append("-" * 35)
        report_lines.append(f"✅ Species assigned: {overall_stats['species_assigned_count']:,} ({overall_stats['species_assigned_pct']:.2f}%)")
        report_lines.append(f"⚠️  Taxid only (no species): {overall_stats['taxid_only']:,} ({overall_stats['taxid_only_pct']:.2f}%)")
        report_lines.append(f"❌ No taxid: {overall_stats['no_taxid']:,} ({overall_stats['no_taxid_pct']:.2f}%)")
        report_lines.append("")

        # Species Coverage Metrics
        avg_assemblies_per_species = overall_stats['species_assigned'] / overall_stats['unique_species'] if overall_stats['unique_species'] > 0 else 0
        report_lines.append("SPECIES COVERAGE METRICS")
        report_lines.append("-" * 30)
        report_lines.append(f"Average assemblies per species: {avg_assemblies_per_species:.2f}")
        report_lines.append(f"Species assignment success rate: {overall_stats['species_assigned_pct']:.2f}%")
        report_lines.append("")

        # Domain Statistics
        report_lines.append("DOMAIN-SPECIFIC STATISTICS")
        report_lines.append("-" * 35)
        for _, row in domain_stats.iterrows():
            report_lines.append(f"{row['domain']}:")
            report_lines.append(f"  Total assemblies: {row['total_assemblies']:,}")
            report_lines.append(f"  Species assigned: {row['species_assigned']:,} ({row['species_assigned_pct']:.2f}%)")
            report_lines.append(f"  Unique species: {row['unique_species']:,}")
            report_lines.append(f"  Avg assemblies/species: {row['avg_assemblies_per_species']:.2f}")
            report_lines.append("")

        # Assembly Level Statistics
        report_lines.append("ASSEMBLY LEVEL vs SPECIES ASSIGNMENT")
        report_lines.append("-" * 45)
        for _, row in level_stats.iterrows():
            report_lines.append(f"{row['assembly_level']}:")
            report_lines.append(f"  Total: {row['total_assemblies']:,}")
            report_lines.append(f"  Species assigned: {row['species_assigned']:,} ({row['species_assigned_pct']:.2f}%)")
            report_lines.append(f"  No assignment: {row['no_assignment']:,} ({row['no_assignment_pct']:.2f}%)")
            report_lines.append("")

        # Problematic Entries
        report_lines.append("PROBLEMATIC ENTRIES ANALYSIS")
        report_lines.append("-" * 40)
        report_lines.append(f"Entries with no taxid: {problems['no_taxid_count']:,}")
        report_lines.append("Top organisms with no taxid:")
        for org, count in list(problems['no_taxid_top_organisms'].items())[:5]:
            report_lines.append(f"  {org}: {count:,}")
        report_lines.append("")

        report_lines.append(f"Entries with taxid but no species assignment: {problems['taxid_only_count']:,}")
        report_lines.append("Top organisms with taxid only:")
        for org, count in list(problems['taxid_only_top_organisms'].items())[:5]:
            report_lines.append(f"  {org}: {count:,}")
        report_lines.append("")

        return "\n".join(report_lines)

    def save_detailed_csvs(self, df, domain_stats, level_stats):
        """Save detailed CSV files for further analysis."""
        print("💾 Saving detailed CSV files...")

        # 1. Domain statistics
        domain_file = self.output_dir / "ncbi_species_stats_by_domain.csv"
        domain_stats.to_csv(domain_file, index=False)

        # 2. Assembly level statistics
        level_file = self.output_dir / "ncbi_species_stats_by_assembly_level.csv"
        level_stats.to_csv(level_file, index=False)

        # 3. Problematic entries summary
        problematic_summary = []

        # No taxid entries
        no_taxid_df = df[~df['has_taxid']]
        if len(no_taxid_df) > 0:
            no_taxid_summary = no_taxid_df.groupby(['organism_name', 'assembly_level']).size().reset_index(name='count')
            no_taxid_summary['problem_type'] = 'no_taxid'
            problematic_summary.append(no_taxid_summary)

        # Taxid only entries
        taxid_only_df = df[df['has_taxid'] & ~df['species_assigned']]
        if len(taxid_only_df) > 0:
            taxid_only_summary = taxid_only_df.groupby(['organism_name', 'assembly_level']).size().reset_index(name='count')
            taxid_only_summary['problem_type'] = 'taxid_only'
            problematic_summary.append(taxid_only_summary)

        if problematic_summary:
            problematic_df = pd.concat(problematic_summary, ignore_index=True)
            problematic_file = self.output_dir / "ncbi_species_problematic_entries.csv"
            problematic_df.to_csv(problematic_file, index=False)

        print(f"📁 Saved CSV files to: {self.output_dir}")

    def run_analysis(self):
        """Run the complete species statistics analysis."""
        try:
            # Load data
            df = self.load_assembly_data()

            # Run analyses
            overall_stats = self.analyze_overall_statistics(df)
            domain_stats = self.analyze_domain_statistics(df)
            level_stats = self.analyze_assembly_level_correlation(df)
            problems = self.analyze_problematic_entries(df)

            # Generate report
            report = self.generate_report(df, overall_stats, domain_stats, level_stats, problems)

            # Save report
            report_file = self.output_dir / "ncbi_species_statistics_report.txt"
            with open(report_file, 'w') as f:
                f.write(report)

            # Save detailed CSVs
            self.save_detailed_csvs(df, domain_stats, level_stats)

            # Print summary to console
            print("\n" + "="*60)
            print("ANALYSIS COMPLETE!")
            print("="*60)
            print(f"📊 Total assemblies: {overall_stats['total_assemblies']:,}")
            print(f"🧬 Unique species: {overall_stats['unique_species']:,}")
            print(f"✅ Species assignment rate: {overall_stats['species_assigned_pct']:.2f}%")
            print(f"⚠️  Problematic entries: {problems['no_taxid_count'] + problems['taxid_only_count']:,}")
            print(f"📁 Report saved: {report_file}")
            print("="*60)

        except Exception as e:
            print(f"❌ Error during analysis: {e}")
            raise

def main():
    """Main function with command line argument support."""
    import argparse

    parser = argparse.ArgumentParser(description="Analyze NCBI species assignment statistics")
    parser.add_argument("--input-dir", help="Directory containing assembly summary file")
    parser.add_argument("--output-dir", help="Directory for output files")

    args = parser.parse_args()

    # Run analysis
    analyzer = NCBISpeciesStatistics(
        input_dir=args.input_dir,
        output_dir=args.output_dir
    )
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
