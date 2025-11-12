#!/usr/bin/env python3
"""
Species Analysis Tools - Direct Assembly Data Analysis
Created: 2025-01-27
Updated: 2025-01-27

This script provides comprehensive analysis of NCBI assembly data directly from the
assembly_summary_genbank.txt file. It generates species-level statistics, taxonomic
distribution analysis, and quality checks without relying on intermediate files.

Usage:
    python species_analysis_tools.py [--assembly-file FILE] [--analysis TYPE]
"""

import pandas as pd
from pathlib import Path
import argparse
import sys
from datetime import datetime

def load_assembly_data(assembly_file):
    """Load and process the NCBI assembly data directly."""
    print("📂 Loading NCBI assembly data...")

    # Read the assembly file, handling the header correctly
    with open(assembly_file, 'r') as f:
        lines = f.readlines()

    # Find the header line (starts with #assembly_accession)
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#assembly_accession'):
            header_line = line.strip().lstrip('#').split('\t')
            data_start = i + 1
            break

    if header_line is None:
        raise ValueError("Could not find header line in assembly file")

    # Read the data with the correct header
    assembly_data = pd.read_csv(assembly_file, sep='\t', skiprows=data_start,
                               names=header_line, low_memory=False)

    print(f"✅ Loaded assembly data: {len(assembly_data):,} records")

    # Clean and prepare the data
    assembly_data = assembly_data.copy()

    # Fill NaN values in isolate column
    assembly_data['isolate'] = assembly_data['isolate'].fillna('')

    # Determine if genome is isolate or uncultured
    assembly_data['is_isolate'] = assembly_data['isolate'] != ''
    assembly_data['is_uncultured'] = ~assembly_data['is_isolate']

    print(f"📊 Isolate genomes: {assembly_data['is_isolate'].sum():,}")
    print(f"📊 Uncultured genomes: {assembly_data['is_uncultured'].sum():,}")

    return assembly_data

def load_taxonomic_mappings(base_dir):
    """Load taxonomic mapping files."""
    print("🏷️  Loading taxonomic mappings...")

    mapping_dir = Path(base_dir) / "taxonomic_mapping"
    mappings = {}

    # Load available mapping files
    mapping_files = {
        'phylum': mapping_dir / "taxid_to_phylum.csv",
        'family': mapping_dir / "taxid_to_family.csv",
        'genus': mapping_dir / "taxid_to_genus.csv"
    }

    for rank, file_path in mapping_files.items():
        if file_path.exists():
            mappings[rank] = pd.read_csv(file_path)
            print(f"✅ Loaded {rank} mappings: {len(mappings[rank]):,} records")
        else:
            print(f"⚠️  {rank} mapping file not found: {file_path}")

    return mappings

def generate_species_statistics(assembly_data, taxonomic_mappings):
    """Generate species-level statistics from assembly data."""
    print("\n📊 Generating species-level statistics...")

    # Group by species_taxid
    species_groups = assembly_data.groupby('species_taxid').agg({
        'organism_name': 'first',  # Representative organism name
        'is_isolate': 'sum',       # Count of isolate genomes
        'is_uncultured': 'sum',    # Count of uncultured genomes
        'assembly_accession': 'count',  # Total genome count
        'genome_size': ['mean', 'std', 'min', 'max'],  # Genome size statistics
        'gc_percent': ['mean', 'std'],  # GC content statistics
        'seq_rel_date': ['min', 'max']  # Date range
    }).reset_index()

    # Flatten column names
    species_groups.columns = [
        'species_taxid', 'representative_organism_name', 'isolate_count',
        'uncultured_count', 'total_genomes', 'mean_genome_size', 'std_genome_size',
        'min_genome_size', 'max_genome_size', 'mean_gc_percent', 'std_gc_percent',
        'earliest_date', 'latest_date'
    ]

    # Add taxonomic information
    for rank, mapping in taxonomic_mappings.items():
        if not mapping.empty:
            # Ensure data types match for merging
            mapping = mapping.copy()
            mapping['taxid'] = mapping['taxid'].astype(str)
            species_groups['species_taxid_str'] = species_groups['species_taxid'].astype(str)

            # Merge taxonomic information
            merge_columns = ['taxid', rank]
            if 'domain' in mapping.columns:
                merge_columns.append('domain')

            species_groups = species_groups.merge(
                mapping[merge_columns],
                left_on='species_taxid_str',
                right_on='taxid',
                how='left'
            ).drop(['taxid', 'species_taxid_str'], axis=1, errors='ignore')

    print(f"✅ Generated statistics for {len(species_groups):,} species")
    return species_groups

def analyze_taxonomic_distribution(species_stats):
    """Analyze taxonomic distribution across different ranks."""
    print("\n🏷️  Taxonomic Distribution Analysis")
    print("=" * 50)

    ranks = ['domain', 'phylum', 'family', 'genus']

    for rank in ranks:
        if rank in species_stats.columns:
            # Count non-empty entries
            non_empty = species_stats[species_stats[rank].notna() & (species_stats[rank] != '')]
            unique_count = non_empty[rank].nunique()
            coverage = len(non_empty) / len(species_stats) * 100

            print(f"\n{rank.capitalize()}:")
            print(f"  Unique {rank}s: {unique_count:,}")
            print(f"  Coverage: {coverage:.1f}% ({len(non_empty):,}/{len(species_stats):,})")

            # Show top 10
            if unique_count > 0:
                top_counts = non_empty.groupby(rank)['total_genomes'].sum().sort_values(ascending=False).head(10)
                print(f"  Top 10 by genome count:")
                for taxa, count in top_counts.items():
                    print(f"    {taxa}: {count:,} genomes")

def analyze_isolate_vs_uncultured(species_stats):
    """Analyze isolate vs uncultured distribution."""
    print("\n🧪 Isolate vs Uncultured Analysis")
    print("=" * 40)

    total_species = len(species_stats)
    isolate_only = (species_stats['isolate_count'] > 0) & (species_stats['uncultured_count'] == 0)
    uncultured_only = (species_stats['isolate_count'] == 0) & (species_stats['uncultured_count'] > 0)
    mixed = (species_stats['isolate_count'] > 0) & (species_stats['uncultured_count'] > 0)

    print(f"Total species: {total_species:,}")
    print(f"Isolate only: {isolate_only.sum():,} ({isolate_only.sum()/total_species*100:.1f}%)")
    print(f"Uncultured only: {uncultured_only.sum():,} ({uncultured_only.sum()/total_species*100:.1f}%)")
    print(f"Mixed (both): {mixed.sum():,} ({mixed.sum()/total_species*100:.1f}%)")

    # Genome distribution
    total_genomes = species_stats['total_genomes'].sum()
    isolate_genomes = species_stats['isolate_count'].sum()
    uncultured_genomes = species_stats['uncultured_count'].sum()

    print(f"\nGenome distribution:")
    print(f"Total genomes: {total_genomes:,}")
    print(f"Isolate genomes: {isolate_genomes:,} ({isolate_genomes/total_genomes*100:.1f}%)")
    print(f"Uncultured genomes: {uncultured_genomes:,} ({uncultured_genomes/total_genomes*100:.1f}%)")

def analyze_species_diversity(species_stats):
    """Analyze species diversity patterns."""
    print("\n🌿 Species Diversity Analysis")
    print("=" * 35)

    # Genome count distribution
    genome_counts = species_stats['total_genomes']

    print(f"Genome count per species:")
    print(f"  Mean: {genome_counts.mean():.1f}")
    print(f"  Median: {genome_counts.median():.1f}")
    print(f"  Min: {genome_counts.min()}")
    print(f"  Max: {genome_counts.max():,}")
    print(f"  Std: {genome_counts.std():.1f}")

    # Distribution bins
    bins = [1, 2, 5, 10, 50, 100, 500, 1000, float('inf')]
    labels = ['1', '2-4', '5-9', '10-49', '50-99', '100-499', '500-999', '1000+']

    binned = pd.cut(genome_counts, bins=bins, labels=labels, right=False)
    distribution = binned.value_counts().sort_index()

    print(f"\nGenome count distribution:")
    for label, count in distribution.items():
        pct = count / len(species_stats) * 100
        print(f"  {label} genomes: {count:,} species ({pct:.1f}%)")

def analyze_genome_quality_metrics(species_stats):
    """Analyze genome quality and assembly metrics."""
    print("\n🔬 Genome Quality Analysis")
    print("=" * 35)

    # Genome size analysis
    if 'mean_genome_size' in species_stats.columns:
        genome_sizes = species_stats['mean_genome_size'].dropna()
        if len(genome_sizes) > 0:
            print(f"Genome size statistics (mean per species):")
            print(f"  Mean: {genome_sizes.mean()/1e6:.2f} Mb")
            print(f"  Median: {genome_sizes.median()/1e6:.2f} Mb")
            print(f"  Min: {genome_sizes.min()/1e6:.2f} Mb")
            print(f"  Max: {genome_sizes.max()/1e6:.2f} Mb")

    # GC content analysis
    if 'mean_gc_percent' in species_stats.columns:
        gc_content = species_stats['mean_gc_percent'].dropna()
        if len(gc_content) > 0:
            print(f"\nGC content statistics (mean per species):")
            print(f"  Mean: {gc_content.mean():.1f}%")
            print(f"  Median: {gc_content.median():.1f}%")
            print(f"  Min: {gc_content.min():.1f}%")
            print(f"  Max: {gc_content.max():.1f}%")

def find_potential_issues(species_stats):
    """Identify potential data quality issues."""
    print("\n⚠️  Potential Data Quality Issues")
    print("=" * 40)

    issues_found = 0

    # Missing taxonomic information
    ranks = ['domain', 'phylum', 'family', 'genus']
    for rank in ranks:
        if rank in species_stats.columns:
            missing = (species_stats[rank] == '') | species_stats[rank].isna()
            if missing.any():
                print(f"Missing {rank}: {missing.sum():,} species ({missing.sum()/len(species_stats)*100:.1f}%)")
                issues_found += 1

    # Unusual genome counts
    high_genome_species = species_stats[species_stats['total_genomes'] > 1000]
    if len(high_genome_species) > 0:
        print(f"\nSpecies with >1000 genomes: {len(high_genome_species)}")
        for _, row in high_genome_species.head(5).iterrows():
            print(f"  {row['representative_organism_name']}: {row['total_genomes']:,} genomes")
        issues_found += 1

    # Species with only uncultured genomes but high counts
    uncultured_high = species_stats[
        (species_stats['isolate_count'] == 0) &
        (species_stats['total_genomes'] > 100)
    ]
    if len(uncultured_high) > 0:
        print(f"\nSpecies with >100 uncultured genomes (no isolates): {len(uncultured_high)}")
        for _, row in uncultured_high.head(3).iterrows():
            print(f"  {row['representative_organism_name']}: {row['total_genomes']:,} genomes")
        issues_found += 1

    # Check for extreme genome sizes
    if 'mean_genome_size' in species_stats.columns:
        genome_sizes = species_stats['mean_genome_size'].dropna()
        if len(genome_sizes) > 0:
            very_small = genome_sizes < 100000  # < 100kb
            very_large = genome_sizes > 50000000  # > 50Mb

            if very_small.any():
                print(f"\nSpecies with very small genomes (<100kb): {very_small.sum()}")
                issues_found += 1

            if very_large.any():
                print(f"Species with very large genomes (>50Mb): {very_large.sum()}")
                issues_found += 1

    if issues_found == 0:
        print("✅ No obvious data quality issues detected")

def save_species_statistics(species_stats, output_dir):
    """Save the generated species statistics to CSV file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save main statistics file
    stats_file = output_path / f"species_statistics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    species_stats.to_csv(stats_file, index=False)
    print(f"💾 Species statistics saved: {stats_file}")

    # Save isolate-preferred subset
    isolate_preferred = species_stats[species_stats['isolate_count'] > 0].copy()
    if len(isolate_preferred) > 0:
        isolate_file = output_path / f"species_isolate_preferred_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        isolate_preferred.to_csv(isolate_file, index=False)
        print(f"💾 Isolate-preferred subset saved: {isolate_file}")

    return stats_file

def generate_summary_report(species_stats, output_dir):
    """Generate a comprehensive summary report."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    report_file = output_path / f"species_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    with open(report_file, 'w') as f:
        f.write("NCBI Assembly Data - Species Analysis Report\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write(f"Dataset Overview:\n")
        f.write(f"  Total species: {len(species_stats):,}\n")
        f.write(f"  Total genomes: {species_stats['total_genomes'].sum():,}\n")
        f.write(f"  Isolate genomes: {species_stats['isolate_count'].sum():,}\n")
        f.write(f"  Uncultured genomes: {species_stats['uncultured_count'].sum():,}\n\n")

        # Top species by genome count
        f.write("Top 20 Species by Genome Count:\n")
        top_species = species_stats.nlargest(20, 'total_genomes')
        for _, row in top_species.iterrows():
            f.write(f"  {row['representative_organism_name']}: {row['total_genomes']:,} genomes\n")

        # Taxonomic distribution summary
        f.write(f"\nTaxonomic Coverage:\n")
        ranks = ['domain', 'phylum', 'family', 'genus']
        for rank in ranks:
            if rank in species_stats.columns:
                non_empty = species_stats[species_stats[rank].notna() & (species_stats[rank] != '')]
                coverage = len(non_empty) / len(species_stats) * 100
                unique_count = non_empty[rank].nunique()
                f.write(f"  {rank.capitalize()}: {unique_count:,} unique ({coverage:.1f}% coverage)\n")

    print(f"📄 Summary report saved: {report_file}")
    return report_file

def main():
    parser = argparse.ArgumentParser(description="Analyze NCBI assembly data directly")
    parser.add_argument("--assembly-file", type=str,
                       default="../00assembly_summary_genbank.txt",
                       help="Path to NCBI assembly summary file")
    parser.add_argument("--base-dir", type=str, default="..",
                       help="Base directory containing taxonomic_mapping folder")
    parser.add_argument("--output-dir", type=str, default=".",
                       help="Output directory for results")
    parser.add_argument("--analysis", type=str,
                       choices=['all', 'taxonomy', 'isolates', 'diversity', 'quality', 'issues'],
                       default='all', help="Type of analysis to perform")
    parser.add_argument("--save-data", action='store_true',
                       help="Save generated species statistics to CSV files")
    args = parser.parse_args()

    print("🔬 NCBI Assembly Data Analysis Tools")
    print("=" * 40)
    print(f"📅 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Load assembly data
    if not Path(args.assembly_file).exists():
        print(f"❌ Assembly file not found: {args.assembly_file}")
        sys.exit(1)

    assembly_data = load_assembly_data(args.assembly_file)

    # Load taxonomic mappings
    taxonomic_mappings = load_taxonomic_mappings(args.base_dir)

    # Generate species statistics
    species_stats = generate_species_statistics(assembly_data, taxonomic_mappings)

    # Perform requested analysis
    if args.analysis in ['all', 'taxonomy']:
        analyze_taxonomic_distribution(species_stats)

    if args.analysis in ['all', 'isolates']:
        analyze_isolate_vs_uncultured(species_stats)

    if args.analysis in ['all', 'diversity']:
        analyze_species_diversity(species_stats)

    if args.analysis in ['all', 'quality']:
        analyze_genome_quality_metrics(species_stats)

    if args.analysis in ['all', 'issues']:
        find_potential_issues(species_stats)

    # Save data and generate report if requested
    if args.save_data or args.analysis == 'all':
        save_species_statistics(species_stats, args.output_dir)
        generate_summary_report(species_stats, args.output_dir)

    print(f"\n✅ Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📊 Analyzed {len(species_stats):,} species from {len(assembly_data):,} genomes")

if __name__ == "__main__":
    main()
