#!/usr/bin/env python3
"""
Test script to filter out bad phylum entries from merger results.

These are entries that exist only in GTDB due to lazy naming conventions
but have no verifiable phylum classification in NCBI taxonomy.

Process:
1. Identify suspicious entries (GTDB_Only with No_Match)
2. Get their accessions from the accession file
3. Look up original NCBI metadata for these accessions
4. Extract taxids and run taxonkit lineage -R
5. Check if they return verifiable phyla
6. Filter out entries with no verifiable phyla
"""

import pandas as pd
import subprocess
import tempfile
import os
import sys
from pathlib import Path
import logging
from typing import List, Dict, Set, Tuple

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_merger_results(results_file: Path) -> pd.DataFrame:
    """Load the merger results CSV file"""
    try:
        df = pd.read_csv(results_file)
        logger.info(f"Loaded {len(df)} entries from {results_file}")
        return df
    except Exception as e:
        logger.error(f"Error loading {results_file}: {e}")
        raise

def identify_suspicious_entries(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify suspicious entries that are GTDB_Only with No_Match.
    These are candidates for filtering.
    """
    suspicious = df[
        (df['Database_Presence'] == 'GTDB_Only') &
        (df['Match_Status'] == 'No_Match')
    ].copy()

    logger.info(f"Found {len(suspicious)} suspicious entries (GTDB_Only + No_Match)")
    return suspicious

def get_accessions_for_phyla(suspicious_df: pd.DataFrame, accession_file: Path) -> Dict[str, List[str]]:
    """
    Get accessions for suspicious phyla from the GTDB accession file.
    Returns dict: {phylum_name: [list_of_accessions]}
    """
    try:
        # Load GTDB accession file
        acc_df = pd.read_csv(accession_file)
        logger.info(f"Loaded {len(acc_df)} accessions from {accession_file}")
        
        phylum_to_accessions = {}
        
        for _, row in suspicious_df.iterrows():
            phylum_name = row['Phylum']
            
            # Find accessions for this phylum in GTDB data
            matching_accessions = acc_df[acc_df['phylum'] == phylum_name]['accession_clean'].tolist()
            
            if matching_accessions:
                phylum_to_accessions[phylum_name] = matching_accessions
                logger.info(f"Found {len(matching_accessions)} accessions for phylum '{phylum_name}'")
            else:
                logger.warning(f"No accessions found for phylum '{phylum_name}'")
        
        return phylum_to_accessions
        
    except Exception as e:
        logger.error(f"Error loading accession file {accession_file}: {e}")
        raise

def get_taxids_from_ncbi_metadata(accessions: List[str], ncbi_metadata_file: Path) -> Dict[str, List[str]]:
    """
    Get taxids for accessions from NCBI metadata file.
    Returns dict: {accession: [list_of_taxids]}
    """
    try:
        # Load NCBI metadata file (assembly summary)
        # Use low_memory=False to avoid dtype warnings
        ncbi_df = pd.read_csv(ncbi_metadata_file, sep='\t', skiprows=1, low_memory=False)
        
        # Rename the first column if it has the # prefix
        if ncbi_df.columns[0].startswith('#'):
            ncbi_df = ncbi_df.rename(columns={ncbi_df.columns[0]: ncbi_df.columns[0][1:]})
        
        logger.info(f"Loaded {len(ncbi_df)} entries from NCBI metadata")
        
        accession_to_taxids = {}
        
        for accession in accessions:
            # Find matching rows in NCBI metadata
            matching_rows = ncbi_df[ncbi_df['assembly_accession'] == accession]
            
            if not matching_rows.empty:
                # Get unique taxids for this accession
                taxids = matching_rows['taxid'].astype(str).unique().tolist()
                accession_to_taxids[accession] = taxids
                logger.debug(f"Accession {accession}: found taxids {taxids}")
            else:
                logger.warning(f"Accession {accession} not found in NCBI metadata")
        
        return accession_to_taxids
        
    except Exception as e:
        logger.error(f"Error processing NCBI metadata {ncbi_metadata_file}: {e}")
        raise

def run_taxonkit_lineage(taxids: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Run taxonkit lineage -R on taxids to get their lineages.
    Returns dict: {taxid: {rank: name}}
    """
    if not taxids:
        return {}
    
    try:
        # Create temporary file with taxids
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            for taxid in taxids:
                temp_file.write(f"{taxid}\n")
            temp_filename = temp_file.name
        
        # Run taxonkit lineage -R
        cmd = ['taxonkit', 'lineage', '-R', temp_filename]
        logger.info(f"Running: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse results
        taxid_lineages = {}
        
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 3:
                    taxid = parts[0]
                    lineage = parts[2] if len(parts) > 2 else ""
                    ranks = parts[1] if len(parts) > 1 else ""
                    
                    # Parse lineage into rank:name dict
                    lineage_dict = {}
                    if lineage and ranks:
                        lineage_names = lineage.split(';')
                        rank_names = ranks.split(';')
                        
                        for rank, name in zip(rank_names, lineage_names):
                            if rank and name:
                                lineage_dict[rank.strip()] = name.strip()
                    
                    taxid_lineages[taxid] = lineage_dict
                    logger.debug(f"Taxid {taxid}: {lineage_dict}")
        
        return taxid_lineages
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Taxonkit command failed: {e}")
        logger.error(f"Stderr: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error running taxonkit: {e}")
        raise
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filename)
        except:
            pass

def has_verifiable_phylum(lineage_dict: Dict[str, str]) -> bool:
    """
    Check if a lineage has a verifiable phylum.
    Returns True if phylum is present and not just a placeholder.
    """
    phylum = lineage_dict.get('phylum', '').strip()
    
    if not phylum:
        return False
    
    # Check for placeholder/invalid phylum names
    invalid_indicators = [
        'unclassified', 'unknown', 'environmental', 'uncultured',
        'candidate', 'candidatus', 'metagenome', 'synthetic'
    ]
    
    phylum_lower = phylum.lower()
    if any(indicator in phylum_lower for indicator in invalid_indicators):
        return False
    
    # If we get here, it seems like a real phylum
    return True

def filter_bad_entries(df: pd.DataFrame, bad_phyla: Set[str]) -> pd.DataFrame:
    """
    Filter out bad entries from the dataframe.
    """
    original_count = len(df)
    filtered_df = df[~df['Phylum'].isin(bad_phyla)].copy()
    filtered_count = len(filtered_df)

    logger.info(f"Filtered out {original_count - filtered_count} bad entries")
    logger.info(f"Remaining entries: {filtered_count}")

    return filtered_df

def main():
    """Main function to test the filtering process"""
    
    # File paths - adjust these as needed
    script_dir = Path(__file__).parent
    
    results_file = script_dir / "merge_outputs/phylum_triple_anchor_output/phylum_cleaned_results.csv"
    gtdb_accession_file = script_dir.parent / "gtdb_parse/csv_gtdb/gtdb_phylum_with_accessions.csv"
    ncbi_metadata_file = script_dir.parent / "ncbi_parse/metadata/00assembly_summary_genbank.txt"
    
    logger.info("=== Testing Bad Entry Filtering ===")
    logger.info(f"Results file: {results_file}")
    logger.info(f"GTDB accession file: {gtdb_accession_file}")
    logger.info(f"NCBI metadata file: {ncbi_metadata_file}")
    
    try:
        # Step 1: Load merger results
        df = load_merger_results(results_file)
        
        # Step 2: Identify suspicious entries
        suspicious_df = identify_suspicious_entries(df)
        
        if suspicious_df.empty:
            logger.info("No suspicious entries found!")
            return
        
        # Step 3: Get accessions for suspicious phyla
        logger.info("Getting accessions for suspicious phyla...")
        phylum_to_accessions = get_accessions_for_phyla(suspicious_df, gtdb_accession_file)
        
        # Step 4: Get taxids from NCBI metadata
        logger.info("Getting taxids from NCBI metadata...")
        all_accessions = []
        for accessions in phylum_to_accessions.values():
            all_accessions.extend(accessions)
        
        accession_to_taxids = get_taxids_from_ncbi_metadata(all_accessions, ncbi_metadata_file)
        
        # Step 5: Run taxonkit on all taxids
        logger.info("Running taxonkit lineage...")
        all_taxids = []
        for taxids in accession_to_taxids.values():
            all_taxids.extend(taxids)
        
        unique_taxids = list(set(all_taxids))
        logger.info(f"Running taxonkit on {len(unique_taxids)} unique taxids")
        
        taxid_lineages = run_taxonkit_lineage(unique_taxids)
        
        # Step 6: Check which phyla have no verifiable classification
        logger.info("Checking for verifiable phyla...")
        bad_phyla = set()
        
        for phylum_name, accessions in phylum_to_accessions.items():
            has_valid_phylum = False
            
            for accession in accessions:
                if accession in accession_to_taxids:
                    for taxid in accession_to_taxids[accession]:
                        if taxid in taxid_lineages:
                            if has_verifiable_phylum(taxid_lineages[taxid]):
                                has_valid_phylum = True
                                break
                    if has_valid_phylum:
                        break
            
            if not has_valid_phylum:
                bad_phyla.add(phylum_name)
                logger.info(f"❌ Phylum '{phylum_name}' has no verifiable classification - marking for removal")
            else:
                logger.info(f"✅ Phylum '{phylum_name}' has verifiable classification - keeping")
        
        # Step 7: Filter out bad entries
        if bad_phyla:
            logger.info(f"Filtering out {len(bad_phyla)} bad phyla: {bad_phyla}")
            filtered_df = filter_bad_entries(df, bad_phyla)
            
            # Save filtered results
            output_file = script_dir / "merge_outputs/phylum_triple_anchor_output/phylum_filtered_results.csv"
            filtered_df.to_csv(output_file, index=False)
            logger.info(f"Saved filtered results to: {output_file}")
        else:
            logger.info("No bad phyla found - no filtering needed")
        
        logger.info("=== Filtering test completed ===")
        
    except Exception as e:
        logger.error(f"Error in main process: {e}")
        raise

if __name__ == "__main__":
    main()
