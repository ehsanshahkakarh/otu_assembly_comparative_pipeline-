#!/usr/bin/env python3
"""
Enhanced Data Quality Checker for Taxonomic Databases

This script validates the completeness and quality of taxonomic data across
GTDB and NCBI databases with comprehensive progress tracking and failure logging.

INPUT FILES:
- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/gtdb_parse/csv_gtdb/gtdb_phylum_with_accessions.csv
  Format: phylum,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: GTDB phylum-level taxonomic data with individual genome accessions for cross-database comparison

- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/gtdb_parse/csv_gtdb/gtdb_family_with_accessions.csv
  Format: family,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: GTDB family-level taxonomic data with individual genome accessions for cross-database comparison

- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/gtdb_parse/csv_gtdb/gtdb_genus_with_accessions.csv
  Format: genus,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: GTDB genus-level taxonomic data with individual genome accessions for cross-database comparison

- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/ncbi_parse/csv_ncbi/ncbi_phylum_with_accessions.csv
  Format: phylum,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: NCBI phylum-level taxonomic data with individual genome accessions for cross-database comparison

- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/ncbi_parse/csv_ncbi/ncbi_family_with_accessions.csv
  Format: family,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: NCBI family-level taxonomic data with individual genome accessions for cross-database comparison

- /clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/ncbi_parse/csv_ncbi/ncbi_genus_with_accessions.csv
  Format: genus,domain,accession_clean,taxid,lineage,lineage_ranks,lineage_taxids
  Description: NCBI genus-level taxonomic data with individual genome accessions for cross-database comparison

OUTPUT FILES:
- ./logs/data_quality_check_YYYYMMDD_HHMMSS.log
  Description: Main log file with timestamped quality check results and statistics

- ./logs/failed_mappings/gtdb_[rank]_failed_mappings_YYYYMMDD_HHMMSS.csv
  Description: GTDB entries with missing taxonomic classifications at specified rank

- ./logs/failed_mappings/ncbi_[rank]_failed_mappings_YYYYMMDD_HHMMSS.csv
  Description: NCBI entries with missing taxonomic classifications at specified rank

- ./logs/taxonomic_inconsistencies_[rank]_YYYYMMDD_HHMMSS.csv
  Description: Accessions with different taxonomic assignments between GTDB and NCBI

- [optional] --output specified file (e.g., quality_report.txt)
  Description: Comprehensive data quality summary report with statistics and recommendations

DIRECTORY STRUCTURE:
- gtdb_parse/csv_gtdb/: GTDB taxonomic count files
- ncbi_parse/csv_ncbi/: NCBI taxonomic count files
- merged_tables/: Location of this script
- merged_tables/logs/: Generated log files and error reports
- merged_tables/logs/failed_mappings/: CSV files with problematic entries

Usage:
    python enhanced_checker.py [--rank RANK] [--detailed] [--output OUTPUT_FILE] [--log-level LEVEL]

Arguments:
    --rank: Taxonomic rank to check (phylum, family, genus, or all) [default: all]
    --detailed: Show detailed statistics and sample problematic entries
    --output: Save comprehensive report to specified file
    --log-level: Logging level (DEBUG, INFO, WARNING, ERROR) [default: INFO]

Examples:
    python enhanced_checker.py --rank all --detailed
    python enhanced_checker.py --rank genus --output genus_quality_report.txt
    python enhanced_checker.py --detailed --log-level DEBUG
"""

import pandas as pd
import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import time

# Get absolute path to the current script
script_dir = Path(__file__).resolve().parent

# Define paths to source directories
ncbi_csv_dir = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/ncbi_parse/csv_ncbi")
gtdb_csv_dir = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/gtdb_parse/csv_gtdb")
eukprot_csv_dir = Path("/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/metadata_proj/parse_repaa_table/eukprot_parse/py_eukprot")

# Define taxonomic ranks to process
RANKS = ["phylum", "family", "genus"]

# Setup logging
def setup_logging(log_level="INFO", log_file=None):
    """Setup logging configuration with optional file output"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    # Create logs directory if it doesn't exist
    log_dir = script_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    
    # Default log file with timestamp
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"data_quality_check_{timestamp}.log"
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format=log_format,
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return log_file

def log_and_print(message, level="info"):
    """
    Print to terminal and log to file simultaneously

    Args:
        message (str): Message to print and log
        level (str): Logging level (info, warning, error, debug)
    """
    print(message)
    getattr(logging, level.lower())(message)

def log_failed_mappings(failed_data, source_name, rank, log_file_path):
    """
    Log detailed information about failed mappings to a separate file

    Args:
        failed_data (pd.DataFrame): DataFrame containing failed entries
        source_name (str): Name of the data source (GTDB, NCBI, EukProt)
        rank (str): Taxonomic rank being processed
        log_file_path (Path): Path to the main log file
    """
    if failed_data.empty:
        return

    # Create failed mappings log file
    failed_log_dir = log_file_path.parent / "failed_mappings"
    failed_log_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    failed_log_file = failed_log_dir / f"{source_name.lower()}_{rank}_failed_mappings_{timestamp}.csv"

    try:
        # Save failed entries to CSV
        failed_data.to_csv(failed_log_file, index=False)

        # Log summary statistics
        total_failed = len(failed_data)
        log_and_print(f"💾 Saved {total_failed} failed {source_name} {rank} entries to: {failed_log_file}", "warning")

        # Log sample of failed entries for immediate review
        if total_failed > 0:
            log_and_print(f"📋 Sample failed {source_name} {rank} entries:", "warning")
            sample_size = min(3, total_failed)
            for idx, row in failed_data.head(sample_size).iterrows():
                log_and_print(f"   Row {idx}: {dict(row)}", "warning")

    except Exception as e:
        log_and_print(f"❌ Error saving failed mappings for {source_name} {rank}: {e}", "error")

def load_source_files(rank, detailed=False, log_file_path=None):
    """
    Load and validate source files for a specific taxonomic rank with progress tracking
    
    Args:
        rank (str): Taxonomic rank to process (phylum, family, genus)
        detailed (bool): Whether to show detailed statistics and samples
        log_file_path (Path): Path to log file for failed mappings
        
    Returns:
        tuple: (sources dict, stats dict) containing loaded dataframes and statistics
    """
    log_and_print(f"\n{'='*70}")
    log_and_print(f"🔍 QUALITY CHECK: {rank.upper()} LEVEL")
    log_and_print(f"{'='*70}")
    log_and_print(f"Starting quality check for {rank} level", "info")
    
    sources = {}
    stats = {}
    
    # Define source configurations (GTDB and NCBI only - comprehensive prokaryotic and eukaryotic coverage)
    # Using with_accessions files to enable cross-database accession comparison
    source_configs = [
        {
            "name": "gtdb",
            "file_path": gtdb_csv_dir / f"gtdb_{rank}_with_accessions.csv",
            "display_name": "GTDB",
            "expected_columns": [rank, "domain", "accession_clean", "taxid"]
        },
        {
            "name": "ncbi",
            "file_path": ncbi_csv_dir / f"ncbi_{rank}_with_accessions.csv",
            "display_name": "NCBI",
            "expected_columns": [rank, "domain", "accession_clean", "taxid"]
        }
    ]
    
    # Process each source with progress bar
    for config in tqdm(source_configs, desc="Loading data sources", unit="source"):
        source_name = config["name"]
        file_path = config["file_path"]
        display_name = config["display_name"]
        
        if file_path.exists():
            try:
                # Load file with progress indication
                load_msg = f"📂 Loading {display_name} data from {file_path.name}..."
                tqdm.write(load_msg)
                logging.info(load_msg)

                # Simulate loading progress for large files
                with tqdm(total=100, desc=f"Reading {display_name}", leave=False, unit="%") as pbar:
                    sources[source_name] = pd.read_csv(file_path)
                    pbar.update(50)

                    # Validate columns
                    missing_cols = [col for col in config["expected_columns"]
                                  if col not in sources[source_name].columns]
                    pbar.update(25)

                    if missing_cols:
                        warning_msg = f"⚠️  {display_name}: Missing expected columns: {missing_cols}"
                        tqdm.write(warning_msg)
                        logging.warning(warning_msg)

                    pbar.update(25)

                df = sources[source_name]
                success_msg = f"✅ Loaded {display_name} with {len(df):,} rows and {len(df.columns)} columns"
                tqdm.write(success_msg)
                logging.info(success_msg)
                
                # Analyze data quality with progress (GTDB and NCBI have same structure)
                missing_taxa = df[df[rank].isna()]
                missing_pct = (len(missing_taxa) / len(df)) * 100

                stats[source_name] = {
                    "total_records": len(df),
                    "missing_taxa": len(missing_taxa),
                    "missing_percentage": missing_pct,
                    "complete_records": len(df) - len(missing_taxa),
                    "file_path": str(file_path)
                }

                # Log failed mappings
                if len(missing_taxa) > 0:
                    log_failed_mappings(missing_taxa, display_name, rank, log_file_path)

                stats_msg = f"   📊 {display_name}: {len(missing_taxa):,} missing {rank} ({missing_pct:.1f}%)"
                tqdm.write(stats_msg)
                logging.info(stats_msg)

                # Show detailed sample if requested
                if detailed:
                    missing_data = df[df[rank].isna()]
                    if len(missing_data) > 0:
                        sample_msg = f"   📋 Sample {display_name} entries with missing {rank}:"
                        tqdm.write(sample_msg)
                        logging.info(sample_msg)

                        sample_cols = [col for col in ['accession_clean', 'domain', rank, 'taxid']
                                     if col in missing_data.columns]
                        if sample_cols:
                            sample_df = missing_data[sample_cols].head(3)
                            sample_table = f"   {sample_df.to_string(index=False)}"
                            tqdm.write(sample_table)
                            logging.info(sample_table)
                
            except Exception as e:
                error_msg = f"Error loading {display_name} {rank} file: {e}"
                tqdm.write(f"❌ {error_msg}")
                logging.error(error_msg)
                stats[source_name] = {"error": str(e), "file_path": str(file_path)}
        else:
            warning_msg = f"{display_name} {rank} file not found at {file_path}"
            tqdm.write(f"⚠️  {warning_msg}")
            logging.warning(warning_msg)
            stats[source_name] = {"error": "File not found", "file_path": str(file_path)}
    
    logging.info(f"Completed quality check for {rank} level")
    return sources, stats

def extract_rank_from_lineage(df, rank):
    """
    Extract rank-specific data from EukProt lineage information
    
    Args:
        df (pd.DataFrame): EukProt dataframe with lineage column
        rank (str): Target taxonomic rank
        
    Returns:
        pd.DataFrame: Filtered dataframe with rank-specific data
    """
    rank_mapping = {
        "phylum": 1,  # Usually second level after domain
        "family": -2, # Usually second to last
        "genus": -1   # Usually last before species
    }
    
    if rank not in rank_mapping:
        return pd.DataFrame()
    
    rank_data = []
    for idx, row in df.iterrows():
        if pd.notna(row["lineage"]):
            lineage_parts = row["lineage"].split(";")
            rank_idx = rank_mapping[rank]
            
            if abs(rank_idx) < len(lineage_parts):
                rank_taxon = lineage_parts[rank_idx].strip()
                if rank_taxon and rank_taxon != "nan":
                    rank_data.append(row)
    
    return pd.DataFrame(rank_data)


def generate_comprehensive_report(all_stats, output_file=None):
    """
    Generate a comprehensive summary report with detailed statistics and recommendations

    Args:
        all_stats (dict): Dictionary containing statistics for all ranks and sources
        output_file (str, optional): Path to save the report
    """
    log_and_print(f"\n{'='*80}")
    log_and_print("📊 COMPREHENSIVE DATA QUALITY SUMMARY REPORT")
    log_and_print(f"{'='*80}")
    log_and_print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_and_print("Generating comprehensive data quality report", "info")

    report_lines = []
    report_lines.append("="*80)
    report_lines.append("📊 COMPREHENSIVE DATA QUALITY SUMMARY REPORT")
    report_lines.append("="*80)
    report_lines.append(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")

    # Overall statistics table with progress
    log_and_print(f"\n📈 OVERALL STATISTICS")
    log_and_print("-" * 90)
    header = f"{'Database':<12} {'Rank':<8} {'Total':<10} {'Complete':<10} {'Missing':<10} {'Missing %':<10} {'Status':<10}"
    log_and_print(header)
    log_and_print("-" * 90)

    report_lines.extend([
        "📈 OVERALL STATISTICS",
        "-" * 90,
        header,
        "-" * 90
    ])

    # Process statistics with progress bar
    total_entries = len(RANKS) * 2  # 2 sources per rank (GTDB and NCBI)
    with tqdm(total=total_entries, desc="Processing statistics", unit="entry") as pbar:
        for rank in RANKS:
            if rank in all_stats:
                for source in ["gtdb", "ncbi"]:
                    if source in all_stats[rank] and "error" not in all_stats[rank][source]:
                        stats = all_stats[rank][source]

                        total = stats['total_records']
                        complete = stats['complete_records']
                        missing = stats['missing_taxa']
                        missing_pct = stats['missing_percentage']

                        # Determine status
                        if missing_pct < 1.0:
                            status = "🟢 EXCELLENT"
                        elif missing_pct < 5.0:
                            status = "🟡 GOOD"
                        elif missing_pct < 10.0:
                            status = "🟠 FAIR"
                        else:
                            status = "🔴 POOR"

                        line = f"{source.upper():<12} {rank:<8} {total:<10,} {complete:<10,} {missing:<10,} {missing_pct:<10.1f} {status:<10}"
                        log_and_print(line)
                        report_lines.append(line)

                    elif source in all_stats[rank]:
                        line = f"{source.upper():<12} {rank:<8} {'ERROR':<10} {'N/A':<10} {'N/A':<10} {'N/A':<10} {'❌ ERROR':<10}"
                        log_and_print(line)
                        report_lines.append(line)

                    pbar.update(1)

    print("-" * 90)
    report_lines.append("-" * 90)

    # Data quality assessment
    print(f"\n🔍 DETAILED DATA QUALITY ASSESSMENT")
    report_lines.extend(["", "🔍 DETAILED DATA QUALITY ASSESSMENT"])

    quality_issues = []
    recommendations = []

    for rank in RANKS:
        if rank in all_stats:
            print(f"\n{rank.upper()} Level Analysis:")
            report_lines.append(f"\n{rank.upper()} Level Analysis:")

            for source in ["gtdb", "ncbi"]:
                if source in all_stats[rank] and "error" not in all_stats[rank][source]:
                    stats = all_stats[rank][source]

                    missing_pct = stats['missing_percentage']
                    total_records = stats['total_records']
                    complete_records = stats['complete_records']

                    # Quality assessment
                    if missing_pct < 1.0:
                        quality = "EXCELLENT"
                        emoji = "🟢"
                    elif missing_pct < 5.0:
                        quality = "GOOD"
                        emoji = "🟡"
                    elif missing_pct < 10.0:
                        quality = "FAIR"
                        emoji = "🟠"
                        quality_issues.append(f"{source.upper()} {rank}: {missing_pct:.1f}% missing data")
                    else:
                        quality = "POOR"
                        emoji = "🔴"
                        quality_issues.append(f"{source.upper()} {rank}: {missing_pct:.1f}% missing data (CRITICAL)")

                    # Coverage analysis
                    coverage = (complete_records / total_records) * 100

                    line = f"  {emoji} {source.upper()}: {quality} - {coverage:.1f}% coverage ({missing_pct:.1f}% missing)"
                    print(line)
                    report_lines.append(line)

                    # Generate recommendations
                    if missing_pct > 5.0:
                        recommendations.append(f"Review {source.upper()} {rank} data quality: {missing_pct:.1f}% missing taxa")

                elif source in all_stats[rank]:
                    error_msg = all_stats[rank][source]['error']
                    line = f"  ❌ {source.upper()}: ERROR - {error_msg}"
                    print(line)
                    report_lines.append(line)
                    recommendations.append(f"Fix {source.upper()} {rank} data loading issue: {error_msg}")

    # Summary and recommendations
    print(f"\n💡 RECOMMENDATIONS AND ACTION ITEMS")
    report_lines.extend(["", "💡 RECOMMENDATIONS AND ACTION ITEMS"])

    if recommendations:
        print(f"\n🔧 Immediate Actions Required:")
        report_lines.append("\n🔧 Immediate Actions Required:")
        for i, rec in enumerate(recommendations, 1):
            rec_line = f"  {i}. {rec}"
            print(rec_line)
            report_lines.append(rec_line)
    else:
        success_line = "✅ All databases show excellent data quality (< 5% missing data)"
        print(f"\n{success_line}")
        report_lines.append(f"\n{success_line}")

    if quality_issues:
        print(f"\n⚠️  Quality Issues Detected:")
        report_lines.append("\n⚠️  Quality Issues Detected:")
        for issue in quality_issues:
            issue_line = f"  • {issue}"
            print(issue_line)
            report_lines.append(issue_line)

    # File paths summary
    print(f"\n📁 DATA SOURCE PATHS")
    report_lines.extend(["", "📁 DATA SOURCE PATHS"])
    for rank in RANKS:
        if rank in all_stats:
            for source in ["gtdb", "ncbi"]:
                if source in all_stats[rank] and "file_path" in all_stats[rank][source]:
                    path_line = f"  {source.upper()} {rank}: {all_stats[rank][source]['file_path']}"
                    print(path_line)
                    report_lines.append(path_line)

    # Save report if requested
    if output_file:
        try:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'w') as f:
                f.write('\n'.join(report_lines))

            print(f"\n📄 Comprehensive report saved to: {output_path}")
            logging.info(f"Report saved to: {output_path}")

        except Exception as e:
            error_msg = f"Error saving report: {e}"
            print(f"\n❌ {error_msg}")
            logging.error(error_msg)


def check_data_consistency(sources, rank):
    """
    Check for data consistency issues between sources with progress tracking

    Args:
        sources (dict): Dictionary containing loaded dataframes
        rank (str): Taxonomic rank being checked
    """
    print(f"\n🔍 CROSS-SOURCE CONSISTENCY CHECK: {rank.upper()}")
    print("-" * 60)
    logging.info(f"Starting consistency check for {rank}")

    if len(sources) < 2:
        warning_msg = "Need at least 2 sources for consistency checking"
        print(f"⚠️  {warning_msg}")
        logging.warning(warning_msg)
        return

    # Check GTDB vs NCBI consistency
    if "gtdb" in sources and "ncbi" in sources:
        log_and_print("🔄 Analyzing GTDB vs NCBI consistency...")

        # Find common accessions with progress
        gtdb_accessions = set()
        ncbi_accessions = set()

        # Extract accessions with progress bars
        if 'accession_clean' in sources["gtdb"].columns:
            with tqdm(total=len(sources["gtdb"]), desc="Processing GTDB accessions", leave=False) as pbar:
                for acc in sources["gtdb"]['accession_clean'].dropna():
                    gtdb_accessions.add(acc.strip())  # Clean accessions
                    pbar.update(1)

        if 'accession_clean' in sources["ncbi"].columns:
            with tqdm(total=len(sources["ncbi"]), desc="Processing NCBI accessions", leave=False) as pbar:
                for acc in sources["ncbi"]['accession_clean'].dropna():
                    ncbi_accessions.add(acc.strip())  # Clean accessions
                    pbar.update(1)

        # Calculate overlap statistics
        common_accessions = gtdb_accessions.intersection(ncbi_accessions)
        gtdb_only = gtdb_accessions - ncbi_accessions
        ncbi_only = ncbi_accessions - gtdb_accessions

        # Log comprehensive accession statistics
        log_and_print(f"📊 ACCESSION OVERLAP ANALYSIS:")
        log_and_print(f"   🔵 GTDB total accessions: {len(gtdb_accessions):,}")
        log_and_print(f"   🟠 NCBI total accessions: {len(ncbi_accessions):,}")
        log_and_print(f"   🟢 Common accessions: {len(common_accessions):,}")
        log_and_print(f"   🔵 GTDB-only accessions: {len(gtdb_only):,}")
        log_and_print(f"   🟠 NCBI-only accessions: {len(ncbi_only):,}")

        if len(gtdb_accessions) > 0:
            gtdb_overlap_pct = (len(common_accessions) / len(gtdb_accessions)) * 100
            log_and_print(f"   📈 GTDB overlap percentage: {gtdb_overlap_pct:.1f}%")

        if len(ncbi_accessions) > 0:
            ncbi_overlap_pct = (len(common_accessions) / len(ncbi_accessions)) * 100
            log_and_print(f"   📈 NCBI overlap percentage: {ncbi_overlap_pct:.1f}%")

        if common_accessions:
            log_and_print(f"\n🔍 TAXONOMIC CONSISTENCY CHECK:")

            # Sample consistency check with progress
            sample_size = min(20, len(common_accessions))  # Increased sample size
            sample_accessions = list(common_accessions)[:sample_size]

            inconsistencies = []
            agreements = []

            with tqdm(sample_accessions, desc="Checking taxonomic consistency", unit="accession") as pbar:
                for acc in pbar:
                    # Clean accession for matching
                    acc_clean = acc.strip()
                    gtdb_row = sources["gtdb"][sources["gtdb"]['accession_clean'].str.strip() == acc_clean]
                    ncbi_row = sources["ncbi"][sources["ncbi"]['accession_clean'].str.strip() == acc_clean]

                    if not gtdb_row.empty and not ncbi_row.empty:
                        gtdb_taxon = gtdb_row[rank].iloc[0] if rank in gtdb_row.columns else "N/A"
                        ncbi_taxon = ncbi_row[rank].iloc[0] if rank in ncbi_row.columns else "N/A"

                        if pd.notna(gtdb_taxon) and pd.notna(ncbi_taxon):
                            if gtdb_taxon != ncbi_taxon:
                                inconsistency = {
                                    'accession': acc_clean,
                                    'gtdb_taxon': gtdb_taxon,
                                    'ncbi_taxon': ncbi_taxon,
                                    'rank': rank,
                                    'gtdb_domain': gtdb_row['domain'].iloc[0] if 'domain' in gtdb_row.columns else "N/A",
                                    'ncbi_domain': ncbi_row['domain'].iloc[0] if 'domain' in ncbi_row.columns else "N/A"
                                }
                                inconsistencies.append(inconsistency)
                                log_and_print(f"  ⚠️  {acc_clean}: GTDB='{gtdb_taxon}' vs NCBI='{ncbi_taxon}'")
                            else:
                                agreements.append(acc_clean)

            # Calculate consistency statistics
            total_compared = len(agreements) + len(inconsistencies)
            if total_compared > 0:
                agreement_pct = (len(agreements) / total_compared) * 100
                log_and_print(f"\n📈 TAXONOMIC CONSISTENCY STATISTICS:")
                log_and_print(f"   ✅ Agreements: {len(agreements):,} ({agreement_pct:.1f}%)")
                log_and_print(f"   ⚠️  Disagreements: {len(inconsistencies):,} ({100-agreement_pct:.1f}%)")
                log_and_print(f"   📊 Total compared: {total_compared:,} accessions")

            if inconsistencies:
                log_and_print(f"\n💾 Saving {len(inconsistencies)} taxonomic inconsistencies...")

                # Save inconsistencies to log file
                inconsistency_df = pd.DataFrame(inconsistencies)
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                inconsistency_file = script_dir / "logs" / f"taxonomic_inconsistencies_{rank}_{timestamp}.csv"
                inconsistency_df.to_csv(inconsistency_file, index=False)
                log_and_print(f"💾 Saved inconsistencies to: {inconsistency_file}")
            else:
                log_and_print("✅ No taxonomic inconsistencies found in sample")
        else:
            log_and_print("📊 No common accessions found between GTDB and NCBI")


def main():
    """Main function to run the enhanced data quality checker"""
    parser = argparse.ArgumentParser(
        description="Enhanced Data Quality Checker for Taxonomic Databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python enhanced_checker.py --rank all --detailed
  python enhanced_checker.py --rank phylum --output results.txt
  python enhanced_checker.py --detailed --log-level DEBUG
        """
    )

    parser.add_argument("--rank", choices=RANKS + ["all"], default="all",
                       help="Taxonomic rank to check (default: all)")
    parser.add_argument("--detailed", action="store_true",
                       help="Show detailed statistics and sample problematic entries")
    parser.add_argument("--output", type=str,
                       help="Save results to specified file")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       default="INFO", help="Logging level (default: INFO)")

    args = parser.parse_args()

    # Setup logging
    log_file = setup_logging(args.log_level)

    print("🔬 Enhanced Data Quality Checker for Taxonomic Databases")
    print("=" * 70)
    print(f"📝 Logging to: {log_file}")
    print(f"🎯 Target ranks: {args.rank}")
    print(f"📊 Detailed mode: {'ON' if args.detailed else 'OFF'}")

    logging.info("Starting enhanced data quality checker")
    logging.info(f"Arguments: rank={args.rank}, detailed={args.detailed}, output={args.output}")

    all_stats = {}
    ranks_to_check = RANKS if args.rank == "all" else [args.rank]

    # Process each rank with overall progress
    with tqdm(ranks_to_check, desc="Processing taxonomic ranks", unit="rank") as rank_pbar:
        for rank in rank_pbar:
            rank_pbar.set_description(f"Processing {rank}")

            sources, stats = load_source_files(rank, args.detailed, log_file)
            all_stats[rank] = stats

            if args.detailed:
                check_data_consistency(sources, rank)

            # Small delay for visual effect
            time.sleep(0.1)

    # Generate comprehensive report
    generate_comprehensive_report(all_stats, args.output)

    print(f"\n✅ Data quality check completed!")
    print(f"📝 Detailed logs saved to: {log_file}")
    logging.info("Enhanced data quality check completed successfully")


if __name__ == "__main__":
    main()
