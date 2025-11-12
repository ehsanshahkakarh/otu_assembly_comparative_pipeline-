#!/usr/bin/env python3
"""
Improved family mapping script that addresses issues found in the original version.
This script uses taxonkit-like logic to ensure better lineage traversal and capture
the 3,085 recoverable family mappings identified in the sanity check.
"""

import pandas as pd
import subprocess
import sys
import os
from pathlib import Path
import logging
from tqdm import tqdm
import tempfile

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def extract_taxids_from_assembly(assembly_file):
    """Extract unique taxids from assembly file"""
    try:
        assembly_path = Path(assembly_file)
        if not assembly_path.exists():
            script_dir = Path(__file__).resolve().parent
            metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
            if metadata_assembly.exists():
                assembly_path = metadata_assembly
                logging.info(f"📁 Using assembly file from metadata folder: {assembly_path}")
            else:
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

        logging.info(f"📖 Reading assembly file: {assembly_path}")
        
        # Read with proper header handling
        df = pd.read_csv(assembly_path, sep="\t", low_memory=False, comment='#')
        
        # Handle different possible column names
        taxid_col = None
        for col in ['taxid', 'species_taxid', 'organism_taxid']:
            if col in df.columns:
                taxid_col = col
                break
        
        if taxid_col is None:
            # Try reading with skiprows if no taxid column found
            df = pd.read_csv(assembly_path, sep="\t", low_memory=False, skiprows=1)
            for col in ['taxid', 'species_taxid', 'organism_taxid']:
                if col in df.columns:
                    taxid_col = col
                    break
        
        if taxid_col is None:
            raise ValueError("Could not find taxid column in assembly file")
        
        # Extract unique taxids
        unique_taxids = df[taxid_col].dropna().astype(int).unique()
        logging.info(f"🔍 Extracted {len(unique_taxids):,} unique taxids from {len(df):,} assembly entries")
        
        return sorted(unique_taxids)
        
    except Exception as e:
        logging.error(f"Error reading assembly file: {e}")
        sys.exit(1)

def run_taxonkit_lineage_batch(taxids, output_dir, batch_size=10000):
    """Run taxonkit lineage on taxids in batches"""
    logging.info(f"🔍 Running taxonkit lineage on {len(taxids):,} taxids in batches of {batch_size:,}")
    
    # Create temporary files
    temp_dir = Path(output_dir)
    temp_dir.mkdir(exist_ok=True)
    
    all_results = []
    
    # Process in batches to avoid command line length limits
    for i in tqdm(range(0, len(taxids), batch_size), desc="Processing batches"):
        batch = taxids[i:i + batch_size]
        
        # Create temporary files for this batch
        temp_taxids_file = temp_dir / f"temp_taxids_batch_{i}.txt"
        temp_lineage_file = temp_dir / f"temp_lineages_batch_{i}.txt"
        
        try:
            # Write taxids to temporary file
            with open(temp_taxids_file, 'w') as f:
                for taxid in batch:
                    f.write(f"{taxid}\n")
            
            # Run taxonkit lineage command
            cmd = f"cat {temp_taxids_file} | taxonkit lineage -R > {temp_lineage_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                logging.warning(f"Taxonkit warning for batch {i}: {result.stderr}")
                # Continue processing even if some taxids fail
            
            # Parse results from this batch
            if temp_lineage_file.exists():
                with open(temp_lineage_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line:
                            all_results.append(line)
            
            # Clean up batch files
            temp_taxids_file.unlink(missing_ok=True)
            temp_lineage_file.unlink(missing_ok=True)
            
        except Exception as e:
            logging.warning(f"Error processing batch {i}: {e}")
            continue
    
    logging.info(f"✅ Processed {len(all_results):,} lineage results")
    return all_results

def parse_taxonkit_results(lineage_results):
    """Parse taxonkit lineage results and extract family information"""
    logging.info(f"📖 Parsing {len(lineage_results):,} lineage results")
    
    mapped_results = []
    unmapped_results = []
    
    for line in tqdm(lineage_results, desc="Parsing lineages"):
        parts = line.split('\t')
        if len(parts) >= 3:
            taxid = int(parts[0])
            lineage = parts[1] if len(parts) > 1 else ""
            ranks = parts[2] if len(parts) > 2 else ""
            
            # Extract family information
            family_name = ""
            domain = "Unknown"
            
            if lineage and ranks:
                lineage_parts = lineage.split(';')
                rank_parts = ranks.split(';')
                
                # Find family - first try 'family', then 'subfamily' as fallback
                for i, rank in enumerate(rank_parts):
                    if rank.lower().strip() == 'family' and i < len(lineage_parts):
                        family_name = lineage_parts[i].strip()
                        break

                # If no family found, try 'subfamily' as family-level classification
                if not family_name:
                    for i, rank in enumerate(rank_parts):
                        if rank.lower().strip() == 'subfamily' and i < len(lineage_parts):
                            family_name = lineage_parts[i].strip()
                            break
                
                # Determine domain using improved logic
                domain = determine_domain_from_lineage(lineage_parts, rank_parts)
            
            if family_name:
                mapped_results.append({
                    'taxid': taxid,
                    'family': family_name,
                    'domain': domain,
                    'lineage': lineage,
                    'ranks': ranks
                })
            else:
                unmapped_results.append({
                    'taxid': taxid,
                    'reason': 'No family found in lineage',
                    'lineage': lineage,
                    'ranks': ranks
                })
        else:
            # Handle malformed lines
            try:
                taxid = int(parts[0]) if parts else 0
                unmapped_results.append({
                    'taxid': taxid,
                    'reason': 'Malformed lineage result',
                    'lineage': '',
                    'ranks': ''
                })
            except:
                continue
    
    logging.info(f"✅ Found {len(mapped_results):,} entries with family information")
    logging.info(f"⚠️  Found {len(unmapped_results):,} entries without family information")
    
    return mapped_results, unmapped_results

def determine_domain_from_lineage(lineage_parts, rank_parts):
    """Determine domain from lineage with improved logic"""
    # First try to find explicit domain/superkingdom
    for i, rank in enumerate(rank_parts):
        if rank.lower().strip() in ['domain', 'superkingdom'] and i < len(lineage_parts):
            return lineage_parts[i].strip()
    
    # Check for viral indicators
    viral_indicators = [
        'viruses', 'virus', 'viral', 'viridae', 'viricota', 'viricetes',
        'orthornavirae', 'bamfordvirae', 'varidnaviria', 'riboviria',
        'duplodnaviria', 'monodnaviria', 'adnaviria', 'ribozyviria'
    ]
    
    for part in lineage_parts:
        if part and any(indicator in part.lower() for indicator in viral_indicators):
            return "Viruses"
    
    # Check for eukaryotic indicators
    eukaryotic_indicators = [
        'eukaryota', 'eukarya', 'fungi', 'metazoa', 'viridiplantae',
        'stramenopiles', 'alveolata', 'rhizaria', 'excavata', 'amoebozoa',
        'opisthokonta', 'archaeplastida'
    ]
    
    for part in lineage_parts:
        if part and any(indicator in part.lower() for indicator in eukaryotic_indicators):
            return "Eukaryota"
    
    # Check for archaea indicators
    archaea_indicators = ['archaea', 'archaeal']
    for part in lineage_parts:
        if part and any(indicator in part.lower() for indicator in archaea_indicators):
            return "Archaea"
    
    # Check for bacteria indicators or cellular organisms
    bacteria_indicators = ['bacteria', 'bacterial']
    for part in lineage_parts:
        if part and any(indicator in part.lower() for indicator in bacteria_indicators):
            return "Bacteria"
    
    # If we have cellular organisms, default to Bacteria
    if any('cellular organisms' in part.lower() for part in lineage_parts if part):
        return "Bacteria"
    
    return "Unknown"

def save_results(mapped_results, unmapped_results, output_dir):
    """Save the mapping results"""
    script_dir = Path(output_dir)
    error_dir = script_dir / "error_log"
    error_dir.mkdir(exist_ok=True)
    
    # Save mapped results
    output_file = script_dir / "taxid_to_family.csv"
    mapped_df = pd.DataFrame(mapped_results)
    
    if len(mapped_df) > 0:
        # Select only the columns we need for the final output
        final_df = mapped_df[['taxid', 'family', 'domain']].copy()
        final_df.to_csv(output_file, index=False)
        logging.info(f"💾 Saved {len(final_df):,} mapped entries to: {output_file}")
    else:
        logging.warning("No mapped results to save")
    
    # Save unmapped results with more detail
    unmapped_file = error_dir / "unmapped_taxids_family.csv"
    unmapped_df = pd.DataFrame(unmapped_results)
    
    if len(unmapped_df) > 0:
        unmapped_df.to_csv(unmapped_file, index=False)
        logging.info(f"💾 Saved {len(unmapped_df):,} unmapped entries to: {unmapped_file}")
    
    return output_file, unmapped_file

def compare_with_original(improved_file, original_backup_file):
    """Compare results with original mapping backup"""
    try:
        if not Path(original_backup_file).exists():
            logging.warning(f"Original backup file not found: {original_backup_file}")
            return

        # Load both files
        improved_df = pd.read_csv(improved_file)
        original_df = pd.read_csv(original_backup_file)

        # Compare coverage
        improved_taxids = set(improved_df['taxid'])
        original_taxids = set(original_df['taxid'])

        new_mappings = improved_taxids - original_taxids
        lost_mappings = original_taxids - improved_taxids

        logging.info("=" * 60)
        logging.info("📊 COMPARISON WITH ORIGINAL MAPPING")
        logging.info("=" * 60)
        logging.info(f"Original mappings: {len(original_taxids):,}")
        logging.info(f"Improved mappings: {len(improved_taxids):,}")
        logging.info(f"New mappings found: {len(new_mappings):,}")
        logging.info(f"Lost mappings: {len(lost_mappings):,}")
        logging.info(f"Net improvement: {len(new_mappings) - len(lost_mappings):,}")

        if new_mappings:
            logging.info(f"First 10 new taxids: {sorted(list(new_mappings))[:10]}")

    except Exception as e:
        logging.error(f"Error comparing with original: {e}")

def main():
    setup_logging()
    logging.info("🚀 Starting improved family mapping process using taxonkit")
    
    # Set up paths
    script_dir = Path(__file__).resolve().parent
    
    # Find assembly file
    metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
    if metadata_assembly.exists():
        assembly_file = metadata_assembly
    else:
        assembly_file = script_dir / "00assembly_summary_genbank.txt"
    
    # Extract taxids from assembly file
    taxids = extract_taxids_from_assembly(assembly_file)
    
    # Run taxonkit lineage in batches
    lineage_results = run_taxonkit_lineage_batch(taxids, script_dir)
    
    # Parse results
    mapped_results, unmapped_results = parse_taxonkit_results(lineage_results)
    
    # Save results
    output_file, unmapped_file = save_results(mapped_results, unmapped_results, script_dir)
    
    # Create backup of original file if it exists and compare
    original_backup_file = script_dir / "taxid_to_family_original_backup.csv"
    current_file = script_dir / "taxid_to_family.csv"

    # Create backup before overwriting (only if backup doesn't exist)
    if current_file.exists() and not original_backup_file.exists():
        import shutil
        shutil.copy2(current_file, original_backup_file)
        logging.info(f"📋 Created backup of original: {original_backup_file}")

    # Compare with backup if it exists
    if original_backup_file.exists():
        compare_with_original(output_file, original_backup_file)
    else:
        logging.info("No original backup file found for comparison")
    
    # Final statistics
    total_taxids = len(mapped_results) + len(unmapped_results)
    success_rate = (len(mapped_results) / total_taxids) * 100 if total_taxids > 0 else 0
    
    logging.info("=" * 60)
    logging.info("📊 FINAL STATISTICS")
    logging.info("=" * 60)
    logging.info(f"Total taxids processed: {total_taxids:,}")
    logging.info(f"Successfully mapped: {len(mapped_results):,}")
    logging.info(f"Unmapped: {len(unmapped_results):,}")
    logging.info(f"Success rate: {success_rate:.2f}%")
    logging.info(f"Results saved to: {output_file}")
    logging.info(f"Unmapped saved to: {unmapped_file}")
    logging.info("=" * 60)

if __name__ == "__main__":
    main()
