#!/usr/bin/env python3
"""
Fix the interesting unmapped cases by implementing improved mapping logic
for eukaryotic clades and other complex taxonomic hierarchies.
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

def create_clade_to_phylum_mapping():
    """Create mapping rules for eukaryotic clades to phylum equivalents"""
    
    clade_mappings = {
        # SAR supergroup mappings
        'stramenopiles': 'Stramenopiles',
        'alveolata': 'Alveolata', 
        'rhizaria': 'Rhizaria',
        'sar': 'SAR',
        
        # Archaeplastida mappings
        'viridiplantae': 'Chlorophyta',
        'rhodophyta': 'Rhodophyta',
        'glaucocystophyceae': 'Glaucophyta',
        
        # Opisthokonta mappings
        'fungi': 'Fungi',
        'metazoa': 'Metazoa',
        'choanoflagellata': 'Choanoflagellata',
        'ichthyosporea': 'Ichthyosporea',
        'filasterea': 'Filasterea',
        
        # Excavata mappings
        'discoba': 'Discoba',
        'metamonada': 'Metamonada',
        
        # Amoebozoa
        'amoebozoa': 'Amoebozoa',
        
        # CRuMs group
        'cryptista': 'Cryptista',
        'collodictyonidae': 'Collodictyonidae',
        
        # Specific problematic clades
        'ochrophyta': 'Ochrophyta',
        'bigyra': 'Bigyra',
        'ciliophora': 'Ciliophora',
        'dinoflagellata': 'Dinoflagellata',
        'apicomplexa': 'Apicomplexa'
    }
    
    return clade_mappings

def create_bacterial_mapping_rules():
    """Create mapping rules for bacterial groups with missing phylum"""
    
    bacterial_mappings = {
        # Proteobacteria subgroups
        'alphaproteobacteria': 'Pseudomonadota',
        'betaproteobacteria': 'Pseudomonadota', 
        'gammaproteobacteria': 'Pseudomonadota',
        'deltaproteobacteria': 'Pseudomonadota',
        'epsilonproteobacteria': 'Pseudomonadota',
        
        # Other bacterial groups
        'cyanobacteria': 'Cyanobacteriota',
        'actinobacteria': 'Actinobacteriota',
        'firmicutes': 'Bacillota',
        'bacteroidetes': 'Bacteroidota',
        'chloroflexi': 'Chloroflexota',
        'deinococcus-thermus': 'Deinococcota',
        'thermotogae': 'Thermotogota',
        'aquificae': 'Aquificota'
    }
    
    return bacterial_mappings

def fix_taxonomic_entry(row, clade_mappings, bacterial_mappings):
    """Fix a single taxonomic entry using mapping rules"""
    
    lineage = str(row['lineage']).lower()
    ranks = str(row['ranks']).lower()
    lineage_parts = lineage.split(';')
    rank_parts = ranks.split(';')
    
    # Try to find phylum using clade mappings
    for i, (lineage_part, rank_part) in enumerate(zip(lineage_parts, rank_parts)):
        lineage_part = lineage_part.strip()
        rank_part = rank_part.strip()
        
        # Check if this is a clade that can be mapped to phylum
        if rank_part == 'clade':
            for clade_key, phylum_name in clade_mappings.items():
                if clade_key in lineage_part:
                    return {
                        'taxid': row['taxid'],
                        'phylum': phylum_name,
                        'domain': determine_domain(lineage_parts),
                        'mapping_method': f'clade_mapping_{clade_key}',
                        'original_lineage': row['lineage']
                    }
        
        # Check bacterial mappings
        if 'bacteria' in lineage and rank_part in ['class', 'clade']:
            for bacterial_key, phylum_name in bacterial_mappings.items():
                if bacterial_key in lineage_part:
                    return {
                        'taxid': row['taxid'],
                        'phylum': phylum_name,
                        'domain': 'Bacteria',
                        'mapping_method': f'bacterial_mapping_{bacterial_key}',
                        'original_lineage': row['lineage']
                    }
    
    # Handle incertae sedis cases - use parent group
    if 'incertae sedis' in lineage:
        # Find the parent taxonomic group
        for i, (lineage_part, rank_part) in enumerate(zip(lineage_parts, rank_parts)):
            if rank_part.strip() in ['phylum', 'class', 'order'] and 'incertae sedis' not in lineage_part:
                parent_name = lineage_part.strip()
                return {
                    'taxid': row['taxid'],
                    'phylum': f'Unclassified_{parent_name}',
                    'domain': determine_domain(lineage_parts),
                    'mapping_method': 'incertae_sedis_parent',
                    'original_lineage': row['lineage']
                }
    
    # Handle cases where class can serve as phylum (for some eukaryotes)
    if 'eukaryota' in lineage:
        for i, (lineage_part, rank_part) in enumerate(zip(lineage_parts, rank_parts)):
            if rank_part.strip() == 'class' and i < 4:  # Early in hierarchy
                class_name = lineage_part.strip()
                # Some classes can serve as phylum equivalents
                if any(term in class_name.lower() for term in ['phyceae', 'mycetes', 'zoa']):
                    return {
                        'taxid': row['taxid'],
                        'phylum': class_name,
                        'domain': 'Eukaryota',
                        'mapping_method': 'class_as_phylum',
                        'original_lineage': row['lineage']
                    }
    
    return None

def determine_domain(lineage_parts):
    """Determine domain from lineage parts"""
    lineage_str = ';'.join(lineage_parts).lower()
    
    if 'eukaryota' in lineage_str:
        return 'Eukaryota'
    elif 'bacteria' in lineage_str:
        return 'Bacteria'
    elif 'archaea' in lineage_str:
        return 'Archaea'
    elif 'virus' in lineage_str:
        return 'Viruses'
    else:
        return 'Unknown'

def process_interesting_cases(interesting_file):
    """Process the interesting cases file and apply fixes"""
    
    logging.info(f"🔧 Processing interesting cases from: {interesting_file}")
    
    # Load the interesting cases
    df = pd.read_csv(interesting_file)
    logging.info(f"📖 Loaded {len(df)} interesting cases")
    
    # Create mapping rules
    clade_mappings = create_clade_to_phylum_mapping()
    bacterial_mappings = create_bacterial_mapping_rules()
    
    # Process each entry
    fixed_entries = []
    unfixed_entries = []
    
    for _, row in df.iterrows():
        fixed_entry = fix_taxonomic_entry(row, clade_mappings, bacterial_mappings)
        
        if fixed_entry:
            fixed_entries.append(fixed_entry)
        else:
            unfixed_entries.append(row.to_dict())
    
    logging.info(f"✅ Fixed {len(fixed_entries)} entries")
    logging.info(f"⚠️  Could not fix {len(unfixed_entries)} entries")
    
    return fixed_entries, unfixed_entries

def save_results(fixed_entries, unfixed_entries, output_dir):
    """Save the fixed and unfixed entries"""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save fixed entries
    if fixed_entries:
        fixed_file = output_dir / f"fixed_phylum_mappings_{timestamp}.csv"
        fixed_df = pd.DataFrame(fixed_entries)
        fixed_df.to_csv(fixed_file, index=False)
        logging.info(f"💾 Saved {len(fixed_entries)} fixed mappings to: {fixed_file}")
    
    # Save unfixed entries
    if unfixed_entries:
        unfixed_file = output_dir / f"still_unfixed_phylum_{timestamp}.csv"
        unfixed_df = pd.DataFrame(unfixed_entries)
        unfixed_df.to_csv(unfixed_file, index=False)
        logging.info(f"💾 Saved {len(unfixed_entries)} unfixed entries to: {unfixed_file}")
    
    # Create summary report
    summary_file = output_dir / f"fix_summary_{timestamp}.txt"
    with open(summary_file, 'w') as f:
        f.write("INTERESTING CASES FIX SUMMARY\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total interesting cases processed: {len(fixed_entries) + len(unfixed_entries)}\n")
        f.write(f"Successfully fixed: {len(fixed_entries)}\n")
        f.write(f"Still unfixed: {len(unfixed_entries)}\n")
        f.write(f"Success rate: {(len(fixed_entries) / (len(fixed_entries) + len(unfixed_entries))) * 100:.1f}%\n\n")
        
        if fixed_entries:
            f.write("MAPPING METHODS USED:\n")
            method_counts = {}
            for entry in fixed_entries:
                method = entry['mapping_method']
                method_counts[method] = method_counts.get(method, 0) + 1
            
            for method, count in sorted(method_counts.items()):
                f.write(f"  {method}: {count}\n")
    
    logging.info(f"📄 Summary report saved to: {summary_file}")
    
    return fixed_file if fixed_entries else None, unfixed_file if unfixed_entries else None

def main():
    setup_logging()
    logging.info("🚀 Starting fix for interesting unmapped cases")
    
    # Find the most recent interesting cases file
    script_dir = Path(__file__).resolve().parent
    interesting_files = list(script_dir.glob("interesting_unmapped_phylum_*.csv"))
    
    if not interesting_files:
        logging.error("No interesting cases file found!")
        return
    
    # Use the most recent file
    interesting_file = max(interesting_files, key=lambda x: x.stat().st_mtime)
    logging.info(f"📁 Using file: {interesting_file}")
    
    # Process the cases
    fixed_entries, unfixed_entries = process_interesting_cases(interesting_file)
    
    # Save results
    fixed_file, unfixed_file = save_results(fixed_entries, unfixed_entries, script_dir)
    
    # Final summary
    total_cases = len(fixed_entries) + len(unfixed_entries)
    success_rate = (len(fixed_entries) / total_cases) * 100 if total_cases > 0 else 0
    
    logging.info("=" * 60)
    logging.info("📊 FINAL RESULTS")
    logging.info("=" * 60)
    logging.info(f"Total interesting cases: {total_cases}")
    logging.info(f"Successfully fixed: {len(fixed_entries)}")
    logging.info(f"Still need attention: {len(unfixed_entries)}")
    logging.info(f"Success rate: {success_rate:.1f}%")
    logging.info("=" * 60)

if __name__ == "__main__":
    main()
