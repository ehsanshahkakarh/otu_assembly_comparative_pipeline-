#!/usr/bin/env python3
"""
Check if any of the interesting unmapped entries at one taxonomic level
actually made it into the count files at other taxonomic levels.
This helps identify inconsistencies in the mapping process.
"""

import pandas as pd
from pathlib import Path
import logging
from datetime import datetime

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def find_latest_file(directory, pattern):
    """Find the most recent file matching the pattern"""
    files = list(directory.glob(pattern))
    if files:
        return max(files, key=lambda x: x.stat().st_mtime)
    return None

def load_count_files(base_dir):
    """Load the NCBI count files"""
    count_files = {}
    
    # Look for count files in the main ncbi_parse directory
    ncbi_parse_dir = base_dir.parent.parent
    
    # Common locations for count files
    possible_dirs = [
        ncbi_parse_dir,
        ncbi_parse_dir / "csv_ncbi",
        ncbi_parse_dir / "output",
        ncbi_parse_dir / "results"
    ]
    
    count_file_patterns = {
        'phylum': 'ncbi_phylum_with_accessions.csv',
        'family': 'ncbi_family_with_accessions.csv',
        'genus': 'ncbi_genus_with_accessions.csv'
    }
    
    for level, filename in count_file_patterns.items():
        found = False
        for dir_path in possible_dirs:
            if dir_path.exists():
                count_file = dir_path / filename
                if count_file.exists():
                    try:
                        df = pd.read_csv(count_file)
                        count_files[level] = df
                        logging.info(f"✅ Loaded {filename}: {len(df)} entries")
                        found = True
                        break
                    except Exception as e:
                        logging.warning(f"⚠️  Error loading {count_file}: {e}")
        
        if not found:
            logging.warning(f"⚠️  Could not find {filename}")
            count_files[level] = None
    
    return count_files

def check_taxids_in_counts(unmapped_taxids, count_files, unmapped_level):
    """Check if unmapped taxids appear in count files at other levels"""

    results = {}

    for count_level, count_df in count_files.items():
        if count_df is None:
            continue

        # Check if count file has taxid column
        if 'taxid' not in count_df.columns:
            logging.warning(f"⚠️  No 'taxid' column in {count_level} counts file")
            continue

        # Find matches
        count_taxids = set(count_df['taxid'].astype(str))
        unmapped_taxids_str = set(unmapped_taxids.astype(str))

        matches = unmapped_taxids_str.intersection(count_taxids)

        results[count_level] = {
            'matches': matches,
            'count': len(matches),
            'total_unmapped': len(unmapped_taxids),
            'total_in_counts': len(count_taxids)
        }

        if matches:
            logging.info(f"🚨 Found {len(matches)} unmapped {unmapped_level} taxids in {count_level} counts - INCONSISTENCY!")
            # Show a few examples
            example_matches = list(matches)[:3]
            logging.info(f"   Examples: {', '.join(example_matches)}")
        else:
            logging.info(f"✅ No unmapped {unmapped_level} taxids found in {count_level} counts")

    return results

def create_detailed_report(results, unmapped_level, output_dir):
    """Create detailed report of findings"""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = output_dir / f"unmapped_{unmapped_level}_in_counts_report_{timestamp}.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"UNMAPPED {unmapped_level.upper()} TAXIDS IN OTHER COUNT FILES ANALYSIS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        total_unmapped = 0
        total_found_elsewhere = 0
        
        for count_level, result in results.items():
            if result:
                total_unmapped = result['total_unmapped']
                matches_count = result['count']
                total_found_elsewhere += matches_count
                
                f.write(f"CHECKING {unmapped_level.upper()} UNMAPPED IN {count_level.upper()} COUNTS:\n")
                f.write("-" * 50 + "\n")
                f.write(f"Unmapped {unmapped_level} taxids: {total_unmapped:,}\n")
                f.write(f"Total {count_level} count entries: {result['total_in_counts']:,}\n")
                f.write(f"Matches found: {matches_count:,}\n")
                
                if matches_count > 0:
                    pct = (matches_count / total_unmapped) * 100
                    f.write(f"Percentage: {pct:.2f}% of unmapped {unmapped_level} entries\n")
                    f.write("⚠️  INCONSISTENCY DETECTED!\n")
                    f.write(f"These taxids failed {unmapped_level} mapping but succeeded {count_level} mapping\n")
                else:
                    f.write("✅ No inconsistencies found\n")
                f.write("\n")
        
        f.write("SUMMARY:\n")
        f.write("-" * 20 + "\n")
        if total_found_elsewhere > 0:
            f.write(f"🚨 TOTAL INCONSISTENCIES: {total_found_elsewhere:,} taxids\n")
            f.write("These entries failed mapping at one level but succeeded at another\n")
            f.write("This suggests potential issues in the mapping pipeline\n")
        else:
            f.write("✅ NO INCONSISTENCIES FOUND\n")
            f.write("Unmapped entries are consistently unmapped across all levels\n")
    
    logging.info(f"📄 Detailed report saved to: {report_file}")
    return report_file

def save_inconsistent_taxids(results, unmapped_df, unmapped_level, count_files, output_dir):
    """Save details of inconsistent taxids for further investigation"""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    for count_level, result in results.items():
        if result and result['matches']:
            # Get details of matching entries
            matching_taxids = [int(tid) for tid in result['matches']]
            
            # Get unmapped entries that have matches
            inconsistent_unmapped = unmapped_df[unmapped_df['taxid'].isin(matching_taxids)].copy()
            
            # Get corresponding count entries
            count_df = count_files[count_level]
            inconsistent_counts = count_df[count_df['taxid'].isin(matching_taxids)].copy()
            
            # Save both files
            unmapped_file = output_dir / f"inconsistent_{unmapped_level}_unmapped_but_in_{count_level}_{timestamp}.csv"
            counts_file = output_dir / f"inconsistent_{unmapped_level}_found_in_{count_level}_counts_{timestamp}.csv"
            
            inconsistent_unmapped.to_csv(unmapped_file, index=False)
            inconsistent_counts.to_csv(counts_file, index=False)
            
            logging.info(f"💾 Saved inconsistent entries:")
            logging.info(f"   Unmapped: {unmapped_file}")
            logging.info(f"   In counts: {counts_file}")

def main():
    setup_logging()
    logging.info("🚀 Checking if unmapped entries appear in count files at other levels")
    
    # File paths
    script_dir = Path(__file__).resolve().parent
    potential_loss_dir = script_dir / "potential_loss"
    
    # Load count files
    logging.info("📖 Loading NCBI count files...")
    count_files = load_count_files(script_dir)
    
    # Check phylum unmapped in family/genus counts
    logging.info("\n" + "="*60)
    logging.info("🔍 CHECKING PHYLUM UNMAPPED ENTRIES")
    logging.info("="*60)
    
    phylum_file = find_latest_file(potential_loss_dir, "interesting_unmapped_phylum_*.csv")
    if phylum_file:
        logging.info(f"📁 Loading: {phylum_file}")
        phylum_df = pd.read_csv(phylum_file)
        phylum_taxids = phylum_df['taxid']
        
        phylum_results = check_taxids_in_counts(phylum_taxids, count_files, 'phylum')
        phylum_report = create_detailed_report(phylum_results, 'phylum', script_dir)
        save_inconsistent_taxids(phylum_results, phylum_df, 'phylum', count_files, script_dir)
    else:
        logging.warning("⚠️  No phylum unmapped file found")
    
    # Check family unmapped in phylum/genus counts  
    logging.info("\n" + "="*60)
    logging.info("🔍 CHECKING FAMILY UNMAPPED ENTRIES")
    logging.info("="*60)
    
    family_file = find_latest_file(potential_loss_dir, "interesting_unmapped_family_*.csv")
    if family_file:
        logging.info(f"📁 Loading: {family_file}")
        family_df = pd.read_csv(family_file)
        family_taxids = family_df['taxid']
        
        family_results = check_taxids_in_counts(family_taxids, count_files, 'family')
        family_report = create_detailed_report(family_results, 'family', script_dir)
        save_inconsistent_taxids(family_results, family_df, 'family', count_files, script_dir)
    else:
        logging.warning("⚠️  No family unmapped file found")
    
    # Summary
    logging.info("\n" + "="*60)
    logging.info("📊 ANALYSIS COMPLETE")
    logging.info("="*60)
    logging.info("Check the generated reports for detailed findings")
    logging.info("Look for inconsistent_*.csv files if any inconsistencies were found")

if __name__ == "__main__":
    main()
