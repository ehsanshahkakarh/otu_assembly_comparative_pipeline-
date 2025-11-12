#!/usr/bin/env python3
"""
Analyze unmapped family-level taxonomic entries and categorize them
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

def categorize_unmapped_family_entries(df):
    """Categorize unmapped family entries into groups for filtering"""
    categories = {
        'viruses': [],
        'unclassified_bacteria': [],
        'unclassified_eukaryota': [],
        'unclassified_archaea': [],
        'environmental_samples': [],
        'incertae_sedis': [],
        'metagenomes': [],
        'interesting_cases': []  # The ones you want to focus on
    }
    
    for _, row in df.iterrows():
        lineage = str(row['lineage']).lower()
        reason = str(row['reason']).lower()
        
        # Categorize based on lineage content
        if 'virus' in lineage or 'acellular root' in lineage:
            categories['viruses'].append(row)
        elif 'metagenome' in lineage or 'metagenomes' in lineage:
            categories['metagenomes'].append(row)
        elif 'unclassified bacteria' in lineage or ('bacteria' in lineage and 'unclassified' in lineage):
            categories['unclassified_bacteria'].append(row)
        elif 'unclassified eukaryota' in lineage or ('eukaryota' in lineage and 'unclassified' in lineage):
            categories['unclassified_eukaryota'].append(row)
        elif 'unclassified archaea' in lineage or ('archaea' in lineage and 'unclassified' in lineage):
            categories['unclassified_archaea'].append(row)
        elif ('environmental samples' in lineage or 'uncultured' in lineage or 
              'unidentified' in lineage):
            categories['environmental_samples'].append(row)
        elif 'incertae sedis' in lineage or 'incertae sedis' in reason:
            categories['incertae_sedis'].append(row)
        else:
            # These are the interesting cases you want to examine
            categories['interesting_cases'].append(row)
    
    return categories

def analyze_interesting_patterns(interesting_cases):
    """Analyze patterns in the interesting cases"""
    patterns = defaultdict(int)
    
    for case in interesting_cases:
        lineage = str(case['lineage']).lower()
        reason = str(case['reason']).lower()
        
        # Analyze reasons for failure
        if 'no family found' in reason:
            patterns['Missing family rank'] += 1
        elif 'multiple families' in reason:
            patterns['Multiple family matches'] += 1
        elif 'family not in lineage' in reason:
            patterns['Family not in lineage'] += 1
        
        # Analyze lineage patterns
        if 'proteobacteria' in lineage:
            patterns['Proteobacteria subgroups'] += 1
        elif 'firmicutes' in lineage or 'bacillota' in lineage:
            patterns['Firmicutes/Bacillota groups'] += 1
        elif 'bacteroidetes' in lineage or 'bacteroidota' in lineage:
            patterns['Bacteroidetes/Bacteroidota groups'] += 1
        elif 'actinobacteria' in lineage or 'actinobacteriota' in lineage:
            patterns['Actinobacteria groups'] += 1
        elif 'eukaryota' in lineage:
            if 'fungi' in lineage:
                patterns['Fungal groups'] += 1
            elif 'metazoa' in lineage:
                patterns['Animal groups'] += 1
            elif 'stramenopiles' in lineage:
                patterns['Stramenopile groups'] += 1
            else:
                patterns['Other eukaryotic groups'] += 1
        elif 'archaea' in lineage:
            patterns['Archaeal groups'] += 1
        else:
            patterns['Other complex cases'] += 1
    
    return dict(patterns)

def create_family_summary_report(categories, patterns, output_dir):
    """Create a comprehensive summary report for family level"""
    
    total_unmapped = sum(len(entries) for entries in categories.values())
    interesting_count = len(categories['interesting_cases'])
    
    report_file = output_dir / f"family_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NCBI FAMILY-LEVEL TAXONOMIC MAPPING ERROR ANALYSIS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Total unmapped family entries: {total_unmapped:,}\n\n")
        
        f.write("FILTERED CATEGORIES (Common/Expected Issues):\n")
        for category, entries in categories.items():
            if category != 'interesting_cases':
                pct = (len(entries) / total_unmapped) * 100
                f.write(f"  {category.replace('_', ' ').title()}: {len(entries):,} ({pct:.1f}%)\n")
        
        f.write(f"\n🎯 INTERESTING CASES TO EXAMINE: {interesting_count:,}\n")
        pct_interesting = (interesting_count / total_unmapped) * 100
        f.write(f"   ({pct_interesting:.1f}% of total unmapped entries)\n\n")
        
        f.write("Patterns in Interesting Cases:\n")
        for pattern, count in sorted(patterns.items(), key=lambda x: x[1], reverse=True):
            f.write(f"  {pattern}: {count:,}\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("RECOMMENDATIONS:\n")
        f.write("=" * 80 + "\n")
        f.write("1. Focus analysis on the interesting cases - these represent\n")
        f.write("   potential data loss in family-level classification\n\n")
        f.write("2. Viruses and metagenomes can likely be excluded from\n")
        f.write("   family-level analysis as expected\n\n")
        f.write("3. Incertae sedis cases may need special handling or\n")
        f.write("   mapping to higher taxonomic levels\n\n")
        f.write("4. Environmental samples represent uncultured organisms\n")
        f.write("   that may not have family-level classification\n\n")
    
    logging.info(f"📄 Family analysis report saved to: {report_file}")
    return report_file

def save_interesting_family_cases(categories, output_dir):
    """Save interesting family cases to the potential_loss directory"""
    
    potential_loss_dir = output_dir / "potential_loss"
    potential_loss_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    interesting_file = potential_loss_dir / f"interesting_unmapped_family_{timestamp}.csv"
    
    if categories['interesting_cases']:
        interesting_df = pd.DataFrame(categories['interesting_cases'])
        interesting_df.to_csv(interesting_file, index=False)
        logging.info(f"💾 Saved {len(categories['interesting_cases'])} interesting family cases to: {interesting_file}")
    
    # Also save a summary of all categories
    summary_file = output_dir / f"family_category_summary_{timestamp}.txt"
    
    with open(summary_file, 'w') as f:
        f.write("UNMAPPED FAMILY CATEGORIES SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        total_entries = sum(len(entries) for entries in categories.values())
        
        for category, entries in categories.items():
            pct = (len(entries) / total_entries) * 100
            f.write(f"{category.replace('_', ' ').title()}: {len(entries):,} ({pct:.1f}%)\n")
            
            if category == 'interesting_cases' and entries:
                f.write("  Sample entries:\n")
                for i, entry in enumerate(entries[:5]):  # Show first 5
                    f.write(f"    {i+1}. TaxID {entry['taxid']}: {entry['lineage'][:80]}...\n")
                f.write("\n")
    
    logging.info(f"📄 Family category summary saved to: {summary_file}")
    return interesting_file, summary_file

def main():
    setup_logging()
    logging.info("🚀 Starting family-level unmapped taxonomic entries analysis")
    
    # File paths
    script_dir = Path(__file__).resolve().parent
    family_file = script_dir / "unmapped_taxids_family.csv"
    
    if not family_file.exists():
        logging.error(f"Family unmapped file not found: {family_file}")
        return
    
    # Load and analyze family data
    logging.info("🔍 Analyzing unmapped family entries...")
    df = pd.read_csv(family_file)
    logging.info(f"📖 Loaded {len(df)} unmapped family entries")
    
    # Categorize entries
    categories = categorize_unmapped_family_entries(df)
    
    # Analyze interesting patterns
    patterns = analyze_interesting_patterns(categories['interesting_cases'])
    
    # Save interesting cases to potential_loss directory
    interesting_file, summary_file = save_interesting_family_cases(categories, script_dir)
    
    # Create comprehensive report
    report_file = create_family_summary_report(categories, patterns, script_dir)
    
    # Print summary to console
    total_unmapped = len(df)
    interesting_count = len(categories['interesting_cases'])
    
    logging.info("=" * 60)
    logging.info("📊 FAMILY-LEVEL UNMAPPED ENTRIES SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Total unmapped family entries: {total_unmapped:,}")
    logging.info(f"🎯 Interesting family cases: {interesting_count:,}")
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
