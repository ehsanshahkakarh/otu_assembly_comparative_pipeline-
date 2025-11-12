#!/usr/bin/env python3
"""
Analyze patterns in unmapped taxonomic entries to understand mapping failures
and suggest improvements to the mapping logic.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from collections import Counter, defaultdict
from datetime import datetime

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def categorize_unmapped_entries(df):
    """Categorize unmapped entries into groups for filtering"""

    # Initialize categories
    categories = {
        'viruses': [],
        'unclassified_bacteria': [],
        'unclassified_eukaryota': [],
        'unclassified_archaea': [],
        'environmental_samples': [],
        'interesting_cases': []  # The ones you want to focus on
    }

    for _, row in df.iterrows():
        lineage = str(row['lineage']).lower()

        # Categorize based on lineage content
        if 'virus' in lineage or 'acellular root' in lineage:
            categories['viruses'].append(row)
        elif 'unclassified bacteria' in lineage or ('bacteria' in lineage and 'unclassified' in lineage):
            categories['unclassified_bacteria'].append(row)
        elif 'unclassified eukaryota' in lineage or ('eukaryota' in lineage and 'unclassified' in lineage):
            categories['unclassified_eukaryota'].append(row)
        elif 'unclassified archaea' in lineage or ('archaea' in lineage and 'unclassified' in lineage):
            categories['unclassified_archaea'].append(row)
        elif ('environmental samples' in lineage or 'uncultured' in lineage or
              'unidentified' in lineage):
            categories['environmental_samples'].append(row)
        else:
            # These are the interesting cases you want to examine
            categories['interesting_cases'].append(row)

    return categories

def analyze_unmapped_phylum(file_path):
    """Analyze unmapped phylum entries with categorization"""
    logging.info("🔍 Analyzing unmapped phylum entries...")

    df = pd.read_csv(file_path)
    total_unmapped = len(df)

    # Categorize entries
    categories = categorize_unmapped_entries(df)

    # Count each category
    category_counts = {cat: len(entries) for cat, entries in categories.items()}

    # Analyze the interesting cases in more detail
    interesting_patterns = defaultdict(int)
    clade_instead_phylum = 0

    for entry in categories['interesting_cases']:
        lineage = str(entry['lineage']).lower()
        ranks = str(entry['ranks']).lower()

        # Check for clade instead of phylum
        if 'clade' in ranks and 'phylum' not in ranks:
            clade_instead_phylum += 1

        # Identify patterns in interesting cases
        if 'eukaryota' in lineage and 'sar' in lineage:
            interesting_patterns['SAR group eukaryotes'] += 1
        elif 'eukaryota' in lineage and 'opisthokonta' in lineage:
            interesting_patterns['Opisthokonta group'] += 1
        elif 'bacteria' in lineage and 'proteobacteria' in lineage:
            interesting_patterns['Proteobacteria subgroups'] += 1
        elif 'incertae sedis' in lineage:
            interesting_patterns['Incertae sedis (uncertain placement)'] += 1
        else:
            interesting_patterns['Other complex cases'] += 1

    # Generate summary
    results = {
        'total_unmapped': total_unmapped,
        'categories': category_counts,
        'interesting_patterns': dict(interesting_patterns),
        'clade_instead_phylum': clade_instead_phylum,
        'interesting_cases_count': len(categories['interesting_cases'])
    }
    
    return results, categories

def analyze_unmapped_family(file_path):
    """Analyze unmapped family entries"""
    logging.info("🔍 Analyzing unmapped family entries...")
    
    df = pd.read_csv(file_path)
    total_unmapped = len(df)
    
    # Analyze patterns
    incertae_sedis_count = 0
    no_rank_family = 0
    environmental_count = 0
    virus_count = 0
    
    domain_counts = defaultdict(int)
    
    for _, row in df.iterrows():
        lineage = str(row['lineage']).lower()
        ranks = str(row['ranks']).lower()
        
        # Count by domain
        if 'virus' in lineage or 'acellular root' in lineage:
            virus_count += 1
            domain_counts['Viruses'] += 1
        elif 'bacteria' in lineage:
            domain_counts['Bacteria'] += 1
        elif 'eukaryota' in lineage:
            domain_counts['Eukaryota'] += 1
        elif 'archaea' in lineage:
            domain_counts['Archaea'] += 1
        else:
            domain_counts['Other'] += 1
        
        # Check for specific issues
        if 'incertae sedis' in lineage:
            incertae_sedis_count += 1
        
        if 'environmental samples' in lineage or 'uncultured' in lineage:
            environmental_count += 1
        
        # Check if family rank is missing but should be there
        rank_list = ranks.split(';')
        if 'genus' in rank_list and 'family' not in rank_list:
            no_rank_family += 1
    
    results = {
        'total_unmapped': total_unmapped,
        'by_domain': dict(domain_counts),
        'issues': {
            'Incertae sedis (uncertain placement)': incertae_sedis_count,
            'Missing family rank': no_rank_family,
            'Environmental samples': environmental_count,
            'Viruses': virus_count
        }
    }
    
    return results

def analyze_unmapped_genus(file_path):
    """Analyze unmapped genus entries"""
    logging.info("🔍 Analyzing unmapped genus entries...")
    
    df = pd.read_csv(file_path)
    total_unmapped = len(df)
    
    # Note: genus file seems to have limited columns, just taxid and reason
    results = {
        'total_unmapped': total_unmapped,
        'note': 'Limited analysis - genus file contains only taxid and reason columns'
    }
    
    return results

def create_summary_report(phylum_results, family_results, genus_results, output_dir):
    """Create a comprehensive summary report focusing on interesting cases"""

    report_file = output_dir / f"unmapped_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NCBI TAXONOMIC MAPPING ERROR ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Phylum Analysis
        f.write("PHYLUM LEVEL MAPPING ISSUES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total unmapped entries: {phylum_results['total_unmapped']:,}\n\n")

        f.write("FILTERED CATEGORIES (Common/Expected Issues):\n")
        for category, count in phylum_results['categories'].items():
            if category != 'interesting_cases':
                pct = (count / phylum_results['total_unmapped']) * 100
                f.write(f"  {category.replace('_', ' ').title()}: {count:,} ({pct:.1f}%)\n")

        f.write(f"\n🎯 INTERESTING CASES TO EXAMINE: {phylum_results['interesting_cases_count']:,}\n")
        pct_interesting = (phylum_results['interesting_cases_count'] / phylum_results['total_unmapped']) * 100
        f.write(f"   ({pct_interesting:.1f}% of total unmapped entries)\n\n")

        f.write("Patterns in Interesting Cases:\n")
        for pattern, count in phylum_results['interesting_patterns'].items():
            f.write(f"  {pattern}: {count:,}\n")

        f.write(f"\nClade instead of phylum issues: {phylum_results['clade_instead_phylum']:,}\n\n")
        
        # Family Analysis
        f.write("FAMILY LEVEL MAPPING ISSUES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total unmapped entries: {family_results['total_unmapped']:,}\n\n")
        
        f.write("Distribution by Domain:\n")
        for domain, count in family_results['by_domain'].items():
            pct = (count / family_results['total_unmapped']) * 100
            f.write(f"  {domain}: {count:,} ({pct:.1f}%)\n")
        
        f.write("\nSpecific Issues:\n")
        for issue, count in family_results['issues'].items():
            pct = (count / family_results['total_unmapped']) * 100
            f.write(f"  {issue}: {count:,} ({pct:.1f}%)\n")
        
        f.write("\n")
        
        # Genus Analysis
        f.write("GENUS LEVEL MAPPING ISSUES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total unmapped entries: {genus_results['total_unmapped']:,}\n")
        f.write(f"Note: {genus_results['note']}\n\n")
        
        # Recommendations
        f.write("RECOMMENDATIONS FOR IMPROVEMENT\n")
        f.write("-" * 40 + "\n")
        f.write("1. VIRAL ENTRIES:\n")
        f.write("   - Consider creating virus-specific mapping logic\n")
        f.write("   - Map to higher taxonomic levels (realm, kingdom)\n")
        f.write("   - Create 'Viruses' as a pseudo-phylum category\n\n")
        
        f.write("2. EUKARYOTIC CLADE ISSUES:\n")
        f.write("   - Implement clade-to-phylum mapping rules\n")
        f.write("   - Use higher-level clades as phylum equivalents\n")
        f.write("   - Create mapping dictionary for common clades\n\n")
        
        f.write("3. ENVIRONMENTAL SAMPLES:\n")
        f.write("   - Consider excluding from analysis or\n")
        f.write("   - Create 'Environmental/Uncultured' categories\n")
        f.write("   - Use parent taxonomic information when available\n\n")
        
        f.write("4. INCERTAE SEDIS (UNCERTAIN PLACEMENT):\n")
        f.write("   - Use parent taxonomic level for classification\n")
        f.write("   - Create 'Unclassified [Parent]' categories\n\n")
        
        f.write("5. MISSING RANKS:\n")
        f.write("   - Implement rank interpolation logic\n")
        f.write("   - Use taxonkit reformat with --fill-miss-rank option\n")
        f.write("   - Create placeholder entries for missing ranks\n\n")
    
    logging.info(f"📄 Summary report saved to: {report_file}")
    return report_file

def save_interesting_cases(categories, output_dir, level="phylum"):
    """Save interesting cases to separate files for detailed examination"""

    interesting_file = output_dir / f"interesting_unmapped_{level}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

    if categories['interesting_cases']:
        interesting_df = pd.DataFrame(categories['interesting_cases'])
        interesting_df.to_csv(interesting_file, index=False)
        logging.info(f"💾 Saved {len(categories['interesting_cases'])} interesting cases to: {interesting_file}")

    # Also save a summary of all categories
    summary_file = output_dir / f"category_summary_{level}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    with open(summary_file, 'w') as f:
        f.write(f"UNMAPPED {level.upper()} CATEGORIES SUMMARY\n")
        f.write("=" * 50 + "\n\n")

        for category, entries in categories.items():
            f.write(f"{category.replace('_', ' ').title()}: {len(entries):,} entries\n")

            if category == 'interesting_cases' and entries:
                f.write("  Sample entries:\n")
                for i, entry in enumerate(entries[:5]):  # Show first 5
                    f.write(f"    {i+1}. TaxID {entry['taxid']}: {entry['lineage'][:100]}...\n")
                f.write("\n")

    logging.info(f"📄 Category summary saved to: {summary_file}")
    return interesting_file, summary_file

def main():
    setup_logging()
    logging.info("🚀 Starting unmapped taxonomic entries analysis")
    
    # Set up paths
    script_dir = Path(__file__).resolve().parent
    
    # File paths
    phylum_file = script_dir / "unmapped_taxids_phylum_improved.csv"
    family_file = script_dir / "unmapped_taxids_family.csv"
    genus_file = script_dir / "unmapped_taxids_genus.csv"
    
    # Check if files exist
    if not phylum_file.exists():
        phylum_file = script_dir / "unmapped_taxids_phylum.csv"
    
    results = {}
    
    # Analyze each level
    phylum_categories = None
    if phylum_file.exists():
        results['phylum'], phylum_categories = analyze_unmapped_phylum(phylum_file)
        # Save interesting cases for detailed examination
        save_interesting_cases(phylum_categories, script_dir, "phylum")
    else:
        logging.warning("Phylum unmapped file not found")
        results['phylum'] = {'total_unmapped': 0}

    if family_file.exists():
        results['family'] = analyze_unmapped_family(family_file)
    else:
        logging.warning("Family unmapped file not found")
        results['family'] = {'total_unmapped': 0}

    if genus_file.exists():
        results['genus'] = analyze_unmapped_genus(genus_file)
    else:
        logging.warning("Genus unmapped file not found")
        results['genus'] = {'total_unmapped': 0}
    
    # Create summary report
    report_file = create_summary_report(
        results['phylum'], 
        results['family'], 
        results['genus'], 
        script_dir
    )
    
    # Print summary to console
    logging.info("=" * 60)
    logging.info("📊 UNMAPPED ENTRIES SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Phylum level unmapped: {results['phylum']['total_unmapped']:,}")

    if 'interesting_cases_count' in results['phylum']:
        logging.info(f"🎯 Interesting phylum cases: {results['phylum']['interesting_cases_count']:,}")
        pct = (results['phylum']['interesting_cases_count'] / results['phylum']['total_unmapped']) * 100
        logging.info(f"   ({pct:.1f}% of unmapped - these are the ones to focus on)")

        logging.info("\n📋 FILTERED OUT (Common issues):")
        for category, count in results['phylum']['categories'].items():
            if category != 'interesting_cases':
                logging.info(f"   {category.replace('_', ' ').title()}: {count:,}")

    logging.info(f"\nFamily level unmapped: {results['family']['total_unmapped']:,}")
    logging.info(f"Genus level unmapped: {results['genus']['total_unmapped']:,}")
    logging.info("=" * 60)
    logging.info(f"📄 Detailed report: {report_file}")

if __name__ == "__main__":
    main()
