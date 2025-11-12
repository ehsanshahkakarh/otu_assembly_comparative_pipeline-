#!/usr/bin/env python3
"""
Analyze unmapped genus-level taxonomic entries and categorize them
to identify the interesting cases that represent potential data loss.
"""

import pandas as pd
from pathlib import Path
import logging
from datetime import datetime
from collections import defaultdict

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def categorize_unmapped_genus_entries(df):
    """Categorize unmapped genus entries into groups for filtering"""
    categories = {
        'viruses': [],
        'unclassified_bacteria': [],
        'unclassified_eukaryota': [],
        'unclassified_archaea': [],
        'environmental_samples': [],
        'incertae_sedis': [],
        'metagenomes': [],
        'no_genus_rank': [],
        'interesting_cases': []  # The ones you want to focus on
    }
    
    for _, row in df.iterrows():
        reason = str(row['reason']).lower()
        
        # Note: genus file only has taxid and reason columns, no lineage
        # Categorize based on reason content
        if 'virus' in reason or 'viral' in reason:
            categories['viruses'].append(row)
        elif 'metagenome' in reason:
            categories['metagenomes'].append(row)
        elif 'unclassified bacteria' in reason or 'bacteria' in reason and 'unclassified' in reason:
            categories['unclassified_bacteria'].append(row)
        elif 'unclassified eukaryota' in reason or 'eukaryota' in reason and 'unclassified' in reason:
            categories['unclassified_eukaryota'].append(row)
        elif 'unclassified archaea' in reason or 'archaea' in reason and 'unclassified' in reason:
            categories['unclassified_archaea'].append(row)
        elif ('environmental' in reason or 'uncultured' in reason or 
              'unidentified' in reason):
            categories['environmental_samples'].append(row)
        elif 'incertae sedis' in reason:
            categories['incertae_sedis'].append(row)
        elif 'no genus found' in reason or 'genus not found' in reason:
            # Since genus file lacks lineage info, treat these as interesting cases
            # They represent potential data loss that needs investigation
            categories['interesting_cases'].append(row)
        else:
            # Other specific issues
            categories['interesting_cases'].append(row)
    
    return categories

def analyze_genus_failure_reasons(interesting_cases):
    """Analyze reasons for genus mapping failures"""
    reason_patterns = defaultdict(int)
    
    for case in interesting_cases:
        reason = str(case['reason']).lower()
        
        # Categorize failure reasons
        if 'no genus found' in reason or 'genus not found' in reason:
            reason_patterns['Missing genus rank'] += 1
        elif 'multiple genus' in reason or 'multiple genera' in reason:
            reason_patterns['Multiple genus matches'] += 1
        elif 'genus not in lineage' in reason:
            reason_patterns['Genus not in lineage'] += 1
        elif 'invalid genus' in reason:
            reason_patterns['Invalid genus name'] += 1
        elif 'genus rank missing' in reason:
            reason_patterns['Genus rank missing from hierarchy'] += 1
        elif 'lineage too short' in reason:
            reason_patterns['Incomplete lineage'] += 1
        elif 'conflicting' in reason:
            reason_patterns['Conflicting taxonomic information'] += 1
        else:
            reason_patterns['Other mapping issues'] += 1
    
    return dict(reason_patterns)

def create_genus_summary_report(categories, reason_patterns, output_dir):
    """Create a comprehensive summary report for genus level"""
    
    total_unmapped = sum(len(entries) for entries in categories.values())
    interesting_count = len(categories['interesting_cases'])
    
    report_file = output_dir / f"genus_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NCBI GENUS-LEVEL TAXONOMIC MAPPING ERROR ANALYSIS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Total unmapped genus entries: {total_unmapped:,}\n\n")
        
        f.write("FILTERED CATEGORIES (Common/Expected Issues):\n")
        for category, entries in categories.items():
            if category != 'interesting_cases':
                pct = (len(entries) / total_unmapped) * 100
                f.write(f"  {category.replace('_', ' ').title()}: {len(entries):,} ({pct:.1f}%)\n")
        
        f.write(f"\n🎯 INTERESTING CASES TO EXAMINE: {interesting_count:,}\n")
        pct_interesting = (interesting_count / total_unmapped) * 100
        f.write(f"   ({pct_interesting:.1f}% of total unmapped entries)\n\n")
        
        f.write("Failure Reason Patterns in Interesting Cases:\n")
        for pattern, count in sorted(reason_patterns.items(), key=lambda x: x[1], reverse=True):
            f.write(f"  {pattern}: {count:,}\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("RECOMMENDATIONS:\n")
        f.write("=" * 80 + "\n")
        f.write("1. Focus analysis on the interesting cases - these represent\n")
        f.write("   potential data loss in genus-level classification\n\n")
        f.write("2. 'No genus rank' cases may need lineage reconstruction\n")
        f.write("   or mapping to family level instead\n\n")
        f.write("3. Viruses and metagenomes can likely be excluded from\n")
        f.write("   genus-level analysis as expected\n\n")
        f.write("4. Environmental samples represent uncultured organisms\n")
        f.write("   that may not have genus-level classification\n\n")
        f.write("5. Consider using taxonkit reformat with --fill-miss-rank\n")
        f.write("   to handle missing genus ranks\n\n")
    
    logging.info(f"📄 Genus analysis report saved to: {report_file}")
    return report_file

def save_interesting_genus_cases(categories, output_dir):
    """Save interesting genus cases to the potential_loss directory"""
    
    potential_loss_dir = output_dir / "potential_loss"
    potential_loss_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    interesting_file = potential_loss_dir / f"interesting_unmapped_genus_{timestamp}.csv"
    
    if categories['interesting_cases']:
        interesting_df = pd.DataFrame(categories['interesting_cases'])
        interesting_df.to_csv(interesting_file, index=False)
        logging.info(f"💾 Saved {len(categories['interesting_cases'])} interesting genus cases to: {interesting_file}")
    
    # Also save a summary of all categories
    summary_file = output_dir / f"genus_category_summary_{timestamp}.txt"
    
    with open(summary_file, 'w') as f:
        f.write("UNMAPPED GENUS CATEGORIES SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        total_entries = sum(len(entries) for entries in categories.values())
        
        for category, entries in categories.items():
            pct = (len(entries) / total_entries) * 100
            f.write(f"{category.replace('_', ' ').title()}: {len(entries):,} ({pct:.1f}%)\n")
            
            if category == 'interesting_cases' and entries:
                f.write("  Sample entries:\n")
                for i, entry in enumerate(entries[:5]):  # Show first 5
                    f.write(f"    {i+1}. TaxID {entry['taxid']}: {entry['reason']}\n")
                f.write("\n")
    
    logging.info(f"📄 Genus category summary saved to: {summary_file}")
    return interesting_file, summary_file

def main():
    setup_logging()
    logging.info("🚀 Starting genus-level unmapped taxonomic entries analysis")
    
    # File paths
    script_dir = Path(__file__).resolve().parent
    genus_file = script_dir / "unmapped_taxids_genus.csv"
    
    if not genus_file.exists():
        logging.error(f"Genus unmapped file not found: {genus_file}")
        return
    
    # Load and analyze genus data
    logging.info("🔍 Analyzing unmapped genus entries...")
    df = pd.read_csv(genus_file)
    logging.info(f"📖 Loaded {len(df)} unmapped genus entries")
    
    # Check columns available
    logging.info(f"📋 Available columns: {list(df.columns)}")
    
    # Categorize entries
    categories = categorize_unmapped_genus_entries(df)
    
    # Analyze failure reason patterns
    reason_patterns = analyze_genus_failure_reasons(categories['interesting_cases'])
    
    # Save interesting cases to potential_loss directory
    interesting_file, summary_file = save_interesting_genus_cases(categories, script_dir)
    
    # Create comprehensive report
    report_file = create_genus_summary_report(categories, reason_patterns, script_dir)
    
    # Print summary to console
    total_unmapped = len(df)
    interesting_count = len(categories['interesting_cases'])
    
    logging.info("=" * 60)
    logging.info("📊 GENUS-LEVEL UNMAPPED ENTRIES SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Total unmapped genus entries: {total_unmapped:,}")
    logging.info(f"🎯 Interesting genus cases: {interesting_count:,}")
    pct = (interesting_count / total_unmapped) * 100
    logging.info(f"   ({pct:.1f}% of unmapped - these are the ones to focus on)")
    
    logging.info("\n📋 FILTERED OUT (Common issues):")
    for category, entries in categories.items():
        if category != 'interesting_cases':
            pct_cat = (len(entries) / total_unmapped) * 100
            logging.info(f"   {category.replace('_', ' ').title()}: {len(entries):,} ({pct_cat:.1f}%)")
    
    logging.info("=" * 60)
    logging.info(f"📄 Detailed report: {report_file}")
    logging.info(f"💾 Interesting cases: {interesting_file}")

if __name__ == "__main__":
    main()
