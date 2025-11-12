#!/usr/bin/env python3
"""
Family-Level GTDB + NCBI Triple-Anchor Merger
=============================================

Combines GTDB and NCBI taxonomic data at the family level using a triple-anchor strategy:
1. Primary: Direct name matching between counts files (most comprehensive)
2. Secondary: Accession-based name mapping (bridges naming differences)
3. Tertiary: Lineage-based validation using taxid information

Key Approach:
- COUNTS FILES are the PRIMARY/REPRESENTATIVE datasets (comprehensive family data)
- ACCESSION FILES are MAPPING TOOLS (to bridge naming differences between databases)
- Goal: Use accession mapping to match MORE families from the counts files

Features:
- Counts files as authoritative data source for family names and counts
- Accession files used to discover naming mappings between databases
- Enhanced matching through accession-based name translation
- Lineage validation using taxid information
- Support for prokaryotic data (excludes eukaryotes and viruses)

Author: Enhanced Taxonomic Merger Team
Date: 2024
"""

import pandas as pd
import logging
import sys
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Set
# Removed unused tqdm import

# ============================================================================
# CANDIDAT FILTER CONFIGURATION - Simple output filter
# ============================================================================
# To EXCLUDE candidat taxa from output: uncomment the next line
EXCLUDE_CANDIDAT_TAXA = True
# To INCLUDE candidat taxa in output: comment out the line above
# ============================================================================

# Setup logging with proper directory structure
def setup_logging():
    """Setup logging with logs directory."""
    logs_dir = Path('logs')
    logs_dir.mkdir(exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(logs_dir / 'family_triple_anchor_merger.log')
        ]
    )
    return logging.getLogger(__name__)

logger = setup_logging()

def load_data_sources():
    """Load GTDB and NCBI data sources using cleaned counts files as authoritative source."""
    logger.info("Loading GTDB and NCBI data sources...")

    # Define file paths - use CLEANED counts files as authoritative source
    gtdb_counts_file = Path('../gtdb_parse/csv_gtdb/gtdb_family_counts.csv')
    gtdb_species_file = Path('../gtdb_parse/csv_gtdb/gtdb_family_species_counts.csv')
    ncbi_counts_file = Path('../ncbi_parse/csv_ncbi/ncbi_family_counts.csv')
    ncbi_species_file = Path('../ncbi_parse/csv_ncbi/ncbi_family_species_counts.csv')
    gtdb_accessions_file = Path('../gtdb_parse/csv_gtdb/gtdb_family_with_accessions.csv')
    ncbi_accessions_file = Path('../ncbi_parse/csv_ncbi/ncbi_family_with_accessions.csv')

    # Verify files exist
    for file_path in [gtdb_counts_file, gtdb_species_file, ncbi_counts_file, ncbi_species_file, gtdb_accessions_file, ncbi_accessions_file]:
        if not file_path.exists():
            raise FileNotFoundError(f"Required file not found: {file_path}")

    # Load cleaned counts files (authoritative family names)
    gtdb_counts_df = pd.read_csv(gtdb_counts_file)
    gtdb_species_df = pd.read_csv(gtdb_species_file)
    ncbi_counts_df = pd.read_csv(ncbi_counts_file)
    ncbi_species_df = pd.read_csv(ncbi_species_file)

    # Load accession files (for accession mapping)
    gtdb_accessions_df = pd.read_csv(gtdb_accessions_file)
    ncbi_accessions_df = pd.read_csv(ncbi_accessions_file)

    logger.info(f"Loaded GTDB cleaned families: {len(gtdb_counts_df):,}")
    logger.info(f"Loaded GTDB species counts: {len(gtdb_species_df):,}")
    logger.info(f"Loaded NCBI cleaned families: {len(ncbi_counts_df):,}")
    logger.info(f"Loaded NCBI species counts: {len(ncbi_species_df):,}")
    logger.info(f"Loaded GTDB accession entries: {len(gtdb_accessions_df):,}")
    logger.info(f"Loaded NCBI accession entries: {len(ncbi_accessions_df):,}")

    return gtdb_counts_df, gtdb_species_df, ncbi_counts_df, ncbi_species_df, gtdb_accessions_df, ncbi_accessions_df

def normalize_name(name):
    """Normalize taxonomic names for comparison."""
    if pd.isna(name):
        return ""

    # Convert to string and clean
    normalized = str(name).strip().lower()

    # Remove GTDB subdivision suffixes (e.g., Bacillota_A -> Bacillota)
    if '_' in normalized:
        base_name = normalized.split('_')[0]
        return base_name

    return normalized

def get_taxids_from_ncbi_metadata(accessions: List[str], ncbi_metadata_file: Path) -> Dict[str, List[str]]:
    """Get taxids for accessions from NCBI metadata file - VECTORIZED VERSION."""
    try:
        # Load NCBI metadata file (assembly summary)
        ncbi_df = pd.read_csv(ncbi_metadata_file, sep='\t', skiprows=1, low_memory=False)

        # Rename the first column if it has the # prefix
        if ncbi_df.columns[0].startswith('#'):
            ncbi_df = ncbi_df.rename(columns={ncbi_df.columns[0]: ncbi_df.columns[0][1:]})

        logger.info(f"Loaded {len(ncbi_df)} entries from NCBI metadata for taxid lookup")

        # VECTORIZED: Filter for matching accessions all at once
        accessions_set = set(accessions)
        matching_df = ncbi_df[ncbi_df['assembly_accession'].isin(accessions_set)].copy()

        if matching_df.empty:
            logger.warning("No matching accessions found in NCBI metadata")
            return {}

        # VECTORIZED: Group by accession and get unique taxids
        accession_to_taxids = (matching_df
                              .groupby('assembly_accession')['taxid']
                              .apply(lambda x: x.astype(str).unique().tolist())
                              .to_dict())

        logger.info(f"Found taxids for {len(accession_to_taxids)} out of {len(accessions)} accessions")
        return accession_to_taxids

    except Exception as e:
        logger.error(f"Error processing NCBI metadata {ncbi_metadata_file}: {e}")
        raise

def run_taxonkit_lineage(taxids: List[str]) -> Dict[str, Dict[str, str]]:
    """Run taxonkit lineage -R on taxids to get their lineages."""
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
        logger.info(f"Running taxonkit on {len(taxids)} taxids...")

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

        return taxid_lineages

    except subprocess.CalledProcessError as e:
        logger.error(f"Taxonkit command failed: {e}")
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

def has_verifiable_family(lineage_dict: Dict[str, str]) -> bool:
    """Check if a lineage has a verifiable family."""
    family = lineage_dict.get('family', '').strip()

    if not family:
        return False

    # Check for placeholder/invalid family names
    invalid_indicators = [
        'unclassified', 'unknown', 'environmental', 'uncultured',
        'candidate', 'candidatus', 'metagenome', 'synthetic'
    ]

    family_lower = family.lower()
    if any(indicator in family_lower for indicator in invalid_indicators):
        return False

    return True

def filter_non_verified_families(family_summary: Dict, gtdb_accessions_df: pd.DataFrame) -> Dict:
    """Filter out families that have no verifiable classification in NCBI taxonomy - VECTORIZED VERSION."""
    logger.info("Filtering out non-verified family entries (vectorized)...")

    # Identify suspicious entries (GTDB_Only with No_Match) using vectorized operations
    suspicious_families = [
        info['family'] for info in family_summary.values()
        if info.get('database_presence') == 'GTDB_Only' and info.get('match_status') == 'No_Match'
    ]

    if not suspicious_families:
        logger.info("No suspicious families found - no filtering needed")
        return family_summary

    logger.info(f"Found {len(suspicious_families)} suspicious families to verify")

    # VECTORIZED: Get accessions for all suspicious families at once
    suspicious_families_set = set(suspicious_families)
    gtdb_filtered = gtdb_accessions_df[
        gtdb_accessions_df['family'].isin(suspicious_families_set) &
        gtdb_accessions_df['accession_clean'].notna()
    ].copy()

    if gtdb_filtered.empty:
        logger.info("No accessions found for suspicious families")
        return family_summary

    # VECTORIZED: Group by family and limit accessions per family
    family_accessions_df = (gtdb_filtered
                           .groupby('family')['accession_clean']
                           .apply(lambda x: x.head(10).tolist())  # Limit to 10 per family
                           .reset_index())

    # Convert to dict for compatibility
    family_to_accessions = dict(zip(family_accessions_df['family'], family_accessions_df['accession_clean']))

    # Get all accessions for batch processing
    all_accessions = [acc for acc_list in family_to_accessions.values() for acc in acc_list]
    logger.info(f"Processing {len(all_accessions)} accessions for {len(family_to_accessions)} families")

    # Get taxids from NCBI metadata
    ncbi_metadata_file = Path('../ncbi_parse/metadata/00assembly_summary_genbank.txt')

    try:
        accession_to_taxids = get_taxids_from_ncbi_metadata(all_accessions, ncbi_metadata_file)

        # Get all unique taxids
        all_taxids = [taxid for taxids in accession_to_taxids.values() for taxid in taxids]
        unique_taxids = list(set(all_taxids))

        if not unique_taxids:
            logger.warning("No taxids found for suspicious families")
            return family_summary

        # Run taxonkit lineage
        taxid_lineages = run_taxonkit_lineage(unique_taxids)

        # VECTORIZED: Check which families have verifiable classification
        logger.info("Performing vectorized family verification...")

        # Create verification lookup for faster processing
        valid_taxids = {
            taxid for taxid, lineage in taxid_lineages.items()
            if has_verifiable_family(lineage)
        }

        bad_families = set()
        verified_count = 0

        for family_name, accessions in family_to_accessions.items():
            # Check if any accession for this family has a valid taxid
            has_valid_family = any(
                taxid in valid_taxids
                for accession in accessions
                if accession in accession_to_taxids
                for taxid in accession_to_taxids[accession]
            )

            if not has_valid_family:
                bad_families.add(family_name)
            else:
                verified_count += 1

        logger.info(f"Verification complete: {verified_count} verified, {len(bad_families)} non-verified")

        # Filter out bad families and save unverified entries for manual review
        if bad_families:
            # Create unverified entries list for manual review
            unverified_entries = []
            filtered_summary = {}

            for key, info in family_summary.items():
                if info['family'] in bad_families:
                    # Add to unverified list with additional context
                    unverified_entry = {
                        'family': info['family'],
                        'gtdb_domain': info.get('gtdb_domain', ''),
                        'ncbi_domain': info.get('ncbi_domain', ''),
                        'database_presence': info.get('database_presence', ''),
                        'match_status': info.get('match_status', ''),
                        'gtdb_genome_count': info.get('gtdb_genome_count', 0),
                        'gtdb_species_count': info.get('gtdb_species_count', 0),
                        'ncbi_genome_count': info.get('ncbi_genome_count', 0),
                        'ncbi_species_count': info.get('ncbi_species_count', 0),
                        'reason_filtered': 'No verifiable NCBI taxonomy classification',
                        'sample_accessions': ', '.join(family_to_accessions.get(info['family'], [])[:5])  # Show first 5 accessions
                    }
                    unverified_entries.append(unverified_entry)
                else:
                    filtered_summary[key] = info

            # Save unverified entries to CSV for manual review
            if unverified_entries:
                unverified_df = pd.DataFrame(unverified_entries)
                unverified_output_file = Path('merge_outputs/family_triple_anchor_output/family_unverified_entries.csv')
                unverified_output_file.parent.mkdir(parents=True, exist_ok=True)
                unverified_df.to_csv(unverified_output_file, index=False)
                logger.info(f"💾 Saved {len(unverified_entries)} unverified families to: {unverified_output_file}")
                logger.info(f"📋 Review these entries manually to check if any should be retained")

            logger.info(f"Filtered out {len(bad_families)} non-verified families")
            logger.info(f"Remaining families: {len(filtered_summary)} (was {len(family_summary)})")
            return filtered_summary
        else:
            logger.info("All suspicious families were verified - no filtering applied")
            return family_summary

    except Exception as e:
        logger.warning(f"Error during family verification: {e}")
        logger.warning("Continuing without filtering...")
        return family_summary

# Removed unused lineage compatibility function - archaic relic

def create_family_name_mappings(gtdb_counts_df, gtdb_species_df, ncbi_counts_df, ncbi_species_df):
    """Create mappings from raw family names to cleaned names using counts files with genome and species counts."""
    logger.info("Creating family name mappings from counts files with genome and species counts...")

    gtdb_name_map = {}
    ncbi_name_map = {}

    # Create GTDB species count lookup using vectorized operations
    gtdb_species_filtered = gtdb_species_df[gtdb_species_df['family'].notna()].copy()
    gtdb_species_filtered['family_norm'] = gtdb_species_filtered['family'].apply(normalize_name)
    gtdb_species_map = dict(zip(gtdb_species_filtered['family_norm'], gtdb_species_filtered['family_species_count'].fillna(0)))

    # Create GTDB name mapping using vectorized operations
    gtdb_filtered = gtdb_counts_df[
        (gtdb_counts_df['family'].notna()) &
        (~gtdb_counts_df['family'].str.lower().isin(['nan', 'none', ''])) &
        (~gtdb_counts_df['domain'].str.lower().isin(['eukaryota', 'eukaryotes', 'eukarya', 'viruses', 'viral']))
    ].copy()

    gtdb_filtered['family_norm'] = gtdb_filtered['family'].apply(normalize_name)

    for _, row in gtdb_filtered.iterrows():
        family_norm = row['family_norm']
        gtdb_name_map[family_norm] = {
            'cleaned_name': row['family'],
            'domain': row['domain'],
            'genome_count': row.get('family_counts', 0),
            'species_count': gtdb_species_map.get(family_norm, 0)
        }

    # Create NCBI species count lookup using vectorized operations
    ncbi_species_filtered = ncbi_species_df[ncbi_species_df['family'].notna()].copy()
    ncbi_species_filtered['family_norm'] = ncbi_species_filtered['family'].apply(normalize_name)
    ncbi_species_map = dict(zip(ncbi_species_filtered['family_norm'], ncbi_species_filtered['family_species_count'].fillna(0)))

    # Create NCBI name mapping using vectorized operations
    ncbi_filtered = ncbi_counts_df[
        (ncbi_counts_df['family'].notna()) &
        (~ncbi_counts_df['family'].str.lower().isin(['nan', 'none', ''])) &
        (~ncbi_counts_df['domain'].str.lower().isin(['eukaryota', 'eukaryotes', 'eukarya', 'viruses', 'viral'])) &
        (~ncbi_counts_df['family'].str.lower().str.contains('viric|virus|viral', na=False))
    ].copy()

    ncbi_filtered['family_norm'] = ncbi_filtered['family'].apply(normalize_name)

    for _, row in ncbi_filtered.iterrows():
        family_norm = row['family_norm']
        ncbi_name_map[family_norm] = {
            'cleaned_name': row['family'],
            'domain': row['domain'],
            'genome_count': row.get('family_counts', 0),
            'species_count': ncbi_species_map.get(family_norm, 0)
        }

    logger.info(f"GTDB name mappings created: {len(gtdb_name_map)}")
    logger.info(f"NCBI name mappings created: {len(ncbi_name_map)}")

    return gtdb_name_map, ncbi_name_map

# Removed unused create_accession_indices_with_cleaned_names function - archaic relic

# Removed unused perform_accession_based_matching function - archaic relic

# Removed unused add_name_based_matches function - archaic relic



# Removed unused add_unmatched_families function - archaic relic

def create_accession_based_name_mappings(gtdb_accessions_df, ncbi_accessions_df):
    """
    Create name mappings between GTDB and NCBI using accession files.
    This discovers how the same organisms are named differently in each database.
    """
    logger.info("Creating accession-based name mappings between databases...")

    # Filter and prepare GTDB data using vectorized operations
    gtdb_filtered = gtdb_accessions_df[
        (gtdb_accessions_df['accession_clean'].notna()) &
        (gtdb_accessions_df['family'].notna()) &
        (gtdb_accessions_df['domain'].str.lower().isin(['bacteria', 'archaea']))
    ].copy()

    # Create GTDB accession to family mapping
    gtdb_filtered['family_norm'] = gtdb_filtered['family'].apply(normalize_name)
    gtdb_acc_map = dict(zip(gtdb_filtered['accession_clean'], gtdb_filtered['family_norm']))

    logger.info(f"GTDB accession mappings created: {len(gtdb_acc_map):,}")

    # Filter and prepare NCBI data using vectorized operations
    ncbi_filtered = ncbi_accessions_df[
        (ncbi_accessions_df['accession_clean'].notna()) &
        (ncbi_accessions_df['family'].notna()) &
        (ncbi_accessions_df['domain'].str.lower().isin(['bacteria', 'archaea']))
    ].copy()

    # Create NCBI accession to family mapping
    ncbi_filtered['family_norm'] = ncbi_filtered['family'].apply(normalize_name)
    ncbi_acc_map = dict(zip(ncbi_filtered['accession_clean'], ncbi_filtered['family_norm']))

    logger.info(f"NCBI accession mappings created: {len(ncbi_acc_map):,}")

    # Find naming mappings through common accessions
    gtdb_to_ncbi_map = {}
    ncbi_to_gtdb_map = {}

    common_accessions = set(gtdb_acc_map.keys()) & set(ncbi_acc_map.keys())
    logger.info(f"Found {len(common_accessions):,} common accessions for name mapping")

    for acc in common_accessions:
        gtdb_name = gtdb_acc_map[acc]
        ncbi_name = ncbi_acc_map[acc]

        if gtdb_name != ncbi_name:  # Different names for same organism
            if gtdb_name not in gtdb_to_ncbi_map:
                gtdb_to_ncbi_map[gtdb_name] = set()
            if ncbi_name not in ncbi_to_gtdb_map:
                ncbi_to_gtdb_map[ncbi_name] = set()

            gtdb_to_ncbi_map[gtdb_name].add(ncbi_name)
            ncbi_to_gtdb_map[ncbi_name].add(gtdb_name)

    # Convert sets to most common mapping - OPTIMIZED VERSION
    gtdb_to_ncbi_final = {}
    ncbi_to_gtdb_final = {}

    # Pre-compute mapping counts for efficiency
    logger.info("Computing mapping frequencies (optimized)...")
    gtdb_ncbi_counts = {}
    ncbi_gtdb_counts = {}

    # Single pass through common accessions to count all mappings
    for acc in common_accessions:
        gtdb_name = gtdb_acc_map[acc]
        ncbi_name = ncbi_acc_map[acc]

        if gtdb_name != ncbi_name:
            # Count GTDB -> NCBI mappings
            if gtdb_name not in gtdb_ncbi_counts:
                gtdb_ncbi_counts[gtdb_name] = {}
            if ncbi_name not in gtdb_ncbi_counts[gtdb_name]:
                gtdb_ncbi_counts[gtdb_name][ncbi_name] = 0
            gtdb_ncbi_counts[gtdb_name][ncbi_name] += 1

            # Count NCBI -> GTDB mappings
            if ncbi_name not in ncbi_gtdb_counts:
                ncbi_gtdb_counts[ncbi_name] = {}
            if gtdb_name not in ncbi_gtdb_counts[ncbi_name]:
                ncbi_gtdb_counts[ncbi_name][gtdb_name] = 0
            ncbi_gtdb_counts[ncbi_name][gtdb_name] += 1

    # Find most common mappings using pre-computed counts
    for gtdb_name, ncbi_names in gtdb_to_ncbi_map.items():
        if gtdb_name in gtdb_ncbi_counts:
            # Use the NCBI name with highest count for this GTDB name
            gtdb_to_ncbi_final[gtdb_name] = max(ncbi_names, key=lambda x: gtdb_ncbi_counts[gtdb_name].get(x, 0))

    for ncbi_name, gtdb_names in ncbi_to_gtdb_map.items():
        if ncbi_name in ncbi_gtdb_counts:
            # Use the GTDB name with highest count for this NCBI name
            ncbi_to_gtdb_final[ncbi_name] = max(gtdb_names, key=lambda x: ncbi_gtdb_counts[ncbi_name].get(x, 0))

    logger.info(f"Created {len(gtdb_to_ncbi_final)} GTDB->NCBI name mappings")
    logger.info(f"Created {len(ncbi_to_gtdb_final)} NCBI->GTDB name mappings")

    return gtdb_to_ncbi_final, ncbi_to_gtdb_final

# Removed unused perform_accession_based_matching function - archaic relic

# Removed unused add_name_based_matches function - archaic relic

# Removed unused add_unmatched_families function - archaic relic

def generate_cleaned_outputs(family_summary, output_dir, gtdb_counts_df, ncbi_counts_df):
    """Generate final output files using cleaned data."""
    logger.info("Generating final output files from cleaned data...")

    # Convert to final clean format
    final_summary = []
    for _, info in family_summary.items():
        # Apply candidat filter if enabled
        if 'EXCLUDE_CANDIDAT_TAXA' in globals() and EXCLUDE_CANDIDAT_TAXA and info['family'].lower().startswith('candidat'):
            continue
        final_summary.append({
            'Family': info['family'],
            'GTDB_Domain': info['gtdb_domain'],
            'NCBI_Domain': info['ncbi_domain'],
            'Database_Presence': info['database_presence'],
            'Match_Status': info['match_status'],
            'Anchor_Type': info['anchor_type'],
            'GTDB_Genome_Count': info.get('gtdb_genome_count', 0),
            'GTDB_Species_Count': info.get('gtdb_species_count', 0),
            'NCBI_Genome_Count': info.get('ncbi_genome_count', 0),
            'NCBI_Species_Count': info.get('ncbi_species_count', 0)
        })

    # Save main results file
    if final_summary:
        family_df = pd.DataFrame(final_summary)
        # Sort by database presence and then by family name
        family_df = family_df.sort_values(['Database_Presence', 'Family'])
        family_df.to_csv(output_dir / 'family_cleaned_results.csv', index=False)
        logger.info(f"Main results saved with {len(family_df)} unique families")

    # Save separate GTDB and NCBI summaries for reference
    gtdb_summary = gtdb_counts_df[gtdb_counts_df['domain'].isin(['Bacteria', 'Archaea'])].copy()
    ncbi_summary = ncbi_counts_df[ncbi_counts_df['domain'].isin(['Bacteria', 'Archaea'])].copy()

    gtdb_summary.to_csv(output_dir / 'gtdb_family_reference.csv', index=False)
    ncbi_summary.to_csv(output_dir / 'ncbi_family_reference.csv', index=False)

# Removed unused log_cleaned_summary function - archaic relic

def log_corrected_summary(family_summary, gtdb_counts_df, ncbi_counts_df, direct_matches_count, mapped_matches_count):
    """Log comprehensive summary statistics for corrected approach."""
    logger.info(f"\n" + "="*70)
    logger.info("FAMILY TRIPLE-ANCHOR MERGER CORRECTED SUMMARY")
    logger.info("="*70)

    # Count prokaryotic entries in source data
    gtdb_prokaryotic = len(gtdb_counts_df[gtdb_counts_df['domain'].isin(['Bacteria', 'Archaea'])])
    ncbi_prokaryotic = len(ncbi_counts_df[ncbi_counts_df['domain'].isin(['Bacteria', 'Archaea'])])

    logger.info(f"Input Data (Prokaryotic Only - COUNTS FILES as PRIMARY DATA):")
    logger.info(f"  GTDB families: {gtdb_prokaryotic:,}")
    logger.info(f"  NCBI families: {ncbi_prokaryotic:,}")

    logger.info(f"\nMatching Strategy Results:")
    logger.info(f"  Direct name matches: {direct_matches_count}")
    logger.info(f"  Accession-mapped matches: {mapped_matches_count}")
    logger.info(f"  Total matches found: {direct_matches_count + mapped_matches_count}")

    logger.info(f"\nFinal Results:")
    logger.info(f"  Total unique families found: {len(family_summary)}")

    # Count by database presence
    both_count = sum(1 for info in family_summary.values() if info['database_presence'] == 'Both')
    gtdb_only_count = sum(1 for info in family_summary.values() if info['database_presence'] == 'GTDB_Only')
    ncbi_only_count = sum(1 for info in family_summary.values() if info['database_presence'] == 'NCBI_Only')

    logger.info(f"  Families in both databases: {both_count}")
    logger.info(f"  GTDB-only families: {gtdb_only_count}")
    logger.info(f"  NCBI-only families: {ncbi_only_count}")

    # Show top 10 families by presence with actual counts
    both_families = [(info['family'], info.get('gtdb_genome_count', 0), info.get('ncbi_genome_count', 0),
                   info.get('gtdb_species_count', 0), info.get('ncbi_species_count', 0), info.get('anchor_type', 'Unknown'))
                  for info in family_summary.values() if info['database_presence'] == 'Both']
    both_families.sort(key=lambda x: x[1] + x[2], reverse=True)

    logger.info(f"\nTop 10 Families Present in Both Databases (with genome and species counts):")
    for i, (family, gtdb_genome, ncbi_genome, gtdb_species, ncbi_species, anchor_type) in enumerate(both_families[:10], 1):
        # Check if this family has consolidation info
        family_info = None
        for info in family_summary.values():
            if info.get('family') == family:
                family_info = info
                break

        source_info = ""
        if family_info and 'gtdb_source_summary' in family_info:
            source_info = f" ({family_info['gtdb_source_summary']})"

        logger.info(f"  {i:2d}. {family}: GTDB={gtdb_genome:,}g/{gtdb_species:,}s, NCBI={ncbi_genome:,}g/{ncbi_species:,}s [{anchor_type}]{source_info}")

    logger.info(f"\nOutput Files:")
    logger.info(f"  Main results: family_cleaned_results.csv")
    logger.info(f"  GTDB reference: gtdb_family_reference.csv")
    logger.info(f"  NCBI reference: ncbi_family_reference.csv")
    logger.info(f"  All files saved to: family_triple_anchor_output/")
    logger.info("="*70)

def main():
    """Main execution function with corrected approach: counts files as primary data, accessions for mapping."""
    try:
        # Load data sources (counts files = primary data, accession files = mapping tools)
        gtdb_counts_df, gtdb_species_df, ncbi_counts_df, ncbi_species_df, gtdb_accessions_df, ncbi_accessions_df = load_data_sources()

        # Create name mappings from counts files (PRIMARY REPRESENTATIVE DATA)
        gtdb_name_map, ncbi_name_map = create_family_name_mappings(gtdb_counts_df, gtdb_species_df, ncbi_counts_df, ncbi_species_df)

        # Create accession-based name mappings (MAPPING TOOL to bridge naming differences)
        gtdb_to_ncbi_map, ncbi_to_gtdb_map = create_accession_based_name_mappings(gtdb_accessions_df, ncbi_accessions_df)

        # Perform enhanced matching using counts as primary data + accession mappings
        logger.info("Performing enhanced matching: counts files + accession-based name mappings...")

        # Strategy 1: Direct name matching (same normalized names) - Use NCBI names as preferred
        direct_matches = {}
        gtdb_families = set(gtdb_name_map.keys())
        ncbi_families = set(ncbi_name_map.keys())
        common_by_name = gtdb_families & ncbi_families

        for family_norm in common_by_name:
            gtdb_info = gtdb_name_map[family_norm]
            ncbi_info = ncbi_name_map[family_norm]

            direct_matches[family_norm] = {
                'family': ncbi_info['cleaned_name'],  # Use NCBI name as preferred
                'gtdb_domain': gtdb_info['domain'],
                'ncbi_domain': ncbi_info['domain'],
                'database_presence': 'Both',
                'match_status': 'Match',
                'anchor_type': 'Direct_Name',
                'gtdb_genome_count': gtdb_info['genome_count'],
                'gtdb_species_count': gtdb_info['species_count'],
                'ncbi_genome_count': ncbi_info['genome_count'],
                'ncbi_species_count': ncbi_info['species_count']
            }

        # Strategy 2: Accession-based name mapping (bidirectional)
        mapped_matches = {}
        already_matched = set(direct_matches.keys())

        # GTDB families that could map to NCBI via accession mappings - Use NCBI names as preferred
        for gtdb_norm, gtdb_info in gtdb_name_map.items():
            if gtdb_norm not in already_matched and gtdb_norm in gtdb_to_ncbi_map:
                ncbi_mapped_norm = gtdb_to_ncbi_map[gtdb_norm]
                if ncbi_mapped_norm in ncbi_name_map:
                    ncbi_info = ncbi_name_map[ncbi_mapped_norm]
                    mapped_matches[gtdb_norm] = {
                        'family': ncbi_info['cleaned_name'],  # Use NCBI name as preferred
                        'gtdb_domain': gtdb_info['domain'],
                        'ncbi_domain': ncbi_info['domain'],
                        'database_presence': 'Both',
                        'match_status': 'Match',
                        'anchor_type': 'Accession_Mapped_GTDB_to_NCBI',
                        'gtdb_genome_count': gtdb_info['genome_count'],
                        'gtdb_species_count': gtdb_info['species_count'],
                        'ncbi_genome_count': ncbi_info['genome_count'],
                        'ncbi_species_count': ncbi_info['species_count'],
                        'gtdb_original_name': gtdb_info['cleaned_name']  # Keep GTDB name for reference
                    }

        # NCBI families that could map to GTDB via accession mappings (reverse direction) - Use NCBI names as preferred
        for ncbi_norm, ncbi_info in ncbi_name_map.items():
            if ncbi_norm not in already_matched and ncbi_norm in ncbi_to_gtdb_map:
                gtdb_mapped_norm = ncbi_to_gtdb_map[ncbi_norm]
                if gtdb_mapped_norm in gtdb_name_map and gtdb_mapped_norm not in already_matched and gtdb_mapped_norm not in mapped_matches:
                    gtdb_info = gtdb_name_map[gtdb_mapped_norm]
                    mapped_matches[f"ncbi_mapped_{ncbi_norm}"] = {
                        'family': ncbi_info['cleaned_name'],  # Use NCBI name as preferred (consistent)
                        'gtdb_domain': gtdb_info['domain'],
                        'ncbi_domain': ncbi_info['domain'],
                        'database_presence': 'Both',
                        'match_status': 'Match',
                        'anchor_type': 'Accession_Mapped_NCBI_to_GTDB',
                        'gtdb_genome_count': gtdb_info['genome_count'],
                        'gtdb_species_count': gtdb_info['species_count'],
                        'ncbi_genome_count': ncbi_info['genome_count'],
                        'ncbi_species_count': ncbi_info['species_count'],
                        'gtdb_original_name': gtdb_info['cleaned_name']  # Keep GTDB name for reference
                    }

        # Strategy 3: Add unmatched families
        all_matched = set(direct_matches.keys()) | set(mapped_matches.keys())
        unmatched_families = {}

        # GTDB-only families
        for gtdb_norm, gtdb_info in gtdb_name_map.items():
            if gtdb_norm not in all_matched:
                unmatched_families[f"gtdb_{gtdb_norm}"] = {
                    'family': gtdb_info['cleaned_name'],
                    'gtdb_domain': gtdb_info['domain'],
                    'ncbi_domain': '',
                    'database_presence': 'GTDB_Only',
                    'match_status': 'No_Match',
                    'anchor_type': 'None',
                    'gtdb_genome_count': gtdb_info['genome_count'],
                    'gtdb_species_count': gtdb_info['species_count'],
                    'ncbi_genome_count': 0,
                    'ncbi_species_count': 0
                }

        # NCBI-only families - track which NCBI families were matched
        ncbi_matched_norms = set()
        for match in mapped_matches.values():
            # Since we're now using NCBI names as preferred, track by the family name in the match
            family_name = match['family']
            for ncbi_norm, ncbi_info in ncbi_name_map.items():
                if ncbi_info['cleaned_name'] == family_name:
                    ncbi_matched_norms.add(ncbi_norm)
                    break

        for ncbi_norm, ncbi_info in ncbi_name_map.items():
            if ncbi_norm not in common_by_name and ncbi_norm not in ncbi_matched_norms:
                unmatched_families[f"ncbi_{ncbi_norm}"] = {
                    'family': ncbi_info['cleaned_name'],
                    'gtdb_domain': '',
                    'ncbi_domain': ncbi_info['domain'],
                    'database_presence': 'NCBI_Only',
                    'match_status': 'No_Match',
                    'anchor_type': 'None',
                    'gtdb_genome_count': 0,
                    'gtdb_species_count': 0,
                    'ncbi_genome_count': ncbi_info['genome_count'],
                    'ncbi_species_count': ncbi_info['species_count']
                }

        # Consolidate results by NCBI family name to avoid duplicates
        logger.info("Consolidating results by NCBI family name to avoid duplicates...")

        consolidated_results = {}

        # Process all matches and consolidate by NCBI family name
        all_matches = {**direct_matches, **mapped_matches}

        for match_info in all_matches.values():
            ncbi_family_name = match_info['family']  # This is now the NCBI name

            if ncbi_family_name not in consolidated_results:
                # First occurrence of this NCBI family - use counts directly from source files
                consolidated_results[ncbi_family_name] = {
                    'family': ncbi_family_name,
                    'gtdb_domain': match_info['gtdb_domain'],
                    'ncbi_domain': match_info['ncbi_domain'],
                    'database_presence': 'Both',
                    'match_status': 'Match',
                    'anchor_types': [match_info['anchor_type']],
                    'gtdb_genome_count': match_info['gtdb_genome_count'],
                    'gtdb_species_count': match_info['gtdb_species_count'],
                    'ncbi_genome_count': match_info['ncbi_genome_count'],
                    'ncbi_species_count': match_info['ncbi_species_count'],
                    'gtdb_sources': [match_info.get('gtdb_original_name', ncbi_family_name)],
                    'mapping_count': 1
                }
            else:
                # Additional mapping to the same NCBI family - keep original counts, just track sources
                existing = consolidated_results[ncbi_family_name]
                # DO NOT modify counts - keep the original values from the first match
                existing['anchor_types'].append(match_info['anchor_type'])
                existing['gtdb_sources'].append(match_info.get('gtdb_original_name', ncbi_family_name))
                existing['mapping_count'] += 1

                # Update anchor type to show consolidation
                if 'Direct_Name' in existing['anchor_types']:
                    existing['anchor_type'] = 'Direct_Name+Accession_Mapped'
                else:
                    existing['anchor_type'] = 'Multiple_Accession_Mapped'

        # Clean up the anchor_types list and create summary
        for family_name, info in consolidated_results.items():
            unique_anchor_types = list(set(info['anchor_types']))
            if len(unique_anchor_types) == 1:
                info['anchor_type'] = unique_anchor_types[0]
            else:
                info['anchor_type'] = '+'.join(sorted(unique_anchor_types))

            # Create a summary of GTDB sources
            unique_sources = list(set(info['gtdb_sources']))
            info['gtdb_source_summary'] = f"{len(unique_sources)} GTDB sources" if len(unique_sources) > 1 else unique_sources[0]

            # Remove temporary fields
            del info['anchor_types']
            del info['gtdb_sources']

        # Add unmatched families
        family_summary = {**consolidated_results, **unmatched_families}

        # Filter out non-verified families (GTDB lazy naming convention entries)
        family_summary = filter_non_verified_families(family_summary, gtdb_accessions_df)

        # Log consolidation statistics
        original_matches = len(all_matches)
        consolidated_matches = len(consolidated_results)
        logger.info(f"Consolidated {original_matches} individual matches into {consolidated_matches} unique NCBI families")

        # Create output directory structure
        merge_outputs_dir = Path('merge_outputs')
        merge_outputs_dir.mkdir(exist_ok=True)
        output_dir = merge_outputs_dir / 'family_triple_anchor_output'
        output_dir.mkdir(exist_ok=True)

        # Generate final output files
        generate_cleaned_outputs(family_summary, output_dir, gtdb_counts_df, ncbi_counts_df)

        # Log comprehensive summary
        log_corrected_summary(family_summary, gtdb_counts_df, ncbi_counts_df, len(direct_matches), len(mapped_matches))

        logger.info("Corrected triple-anchor merger analysis complete!")

    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
