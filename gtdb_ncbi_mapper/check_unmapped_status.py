#!/usr/bin/env python3
"""
Check the status of unmapped GTDB accessions using NCBI datasets tool.

This script takes the unmapped GTDB accessions and queries NCBI to determine
why they're not in the assembly summary file (e.g., suppressed, replaced, etc.)
"""

import pandas as pd
import subprocess
import json
import logging
from pathlib import Path
from tqdm import tqdm
import time

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/unmapped_status_check.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def check_accession_status(accession: str, max_retries: int = 3) -> dict:
    """
    Check the status of a single accession using NCBI datasets.

    Args:
        accession: The accession to check
        max_retries: Maximum number of retry attempts

    Returns:
        Dictionary with status information
    """
    for attempt in range(max_retries):
        try:
            # Run datasets command using conda run to activate the ncbi_datasets environment
            cmd = ['conda', 'run', '-n', 'ncbi_datasets', 'datasets', 'summary', 'genome', 'accession', accession, '--as-json-lines']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0 and result.stdout.strip():
                # Parse JSON output
                data = json.loads(result.stdout)
                
                if 'reports' in data and len(data['reports']) > 0:
                    report = data['reports'][0]
                    assembly_info = report.get('assembly_info', {})
                    
                    return {
                        'accession': accession,
                        'status': 'found',
                        'assembly_status': assembly_info.get('assembly_status', 'unknown'),
                        'assembly_level': assembly_info.get('assembly_level', 'unknown'),
                        'bioproject': assembly_info.get('bioproject_accession', 'N/A'),
                        'refseq_accession': assembly_info.get('refseq_accession', 'N/A'),
                        'error': None
                    }
                else:
                    return {
                        'accession': accession,
                        'status': 'not_found',
                        'assembly_status': 'not_in_ncbi',
                        'assembly_level': 'N/A',
                        'bioproject': 'N/A',
                        'refseq_accession': 'N/A',
                        'error': 'No reports in response'
                    }
            else:
                # Check if it's a "not found" error
                if 'No data found' in result.stderr or 'not found' in result.stderr.lower():
                    return {
                        'accession': accession,
                        'status': 'not_found',
                        'assembly_status': 'not_in_ncbi',
                        'assembly_level': 'N/A',
                        'bioproject': 'N/A',
                        'refseq_accession': 'N/A',
                        'error': result.stderr.strip()
                    }
                else:
                    # Some other error, might retry
                    if attempt < max_retries - 1:
                        time.sleep(1)  # Wait before retry
                        continue
                    else:
                        return {
                            'accession': accession,
                            'status': 'error',
                            'assembly_status': 'query_failed',
                            'assembly_level': 'N/A',
                            'bioproject': 'N/A',
                            'refseq_accession': 'N/A',
                            'error': result.stderr.strip()
                        }
                        
        except subprocess.TimeoutExpired:
            if attempt < max_retries - 1:
                logger.warning(f"Timeout for {accession}, retrying...")
                continue
            else:
                return {
                    'accession': accession,
                    'status': 'error',
                    'assembly_status': 'timeout',
                    'assembly_level': 'N/A',
                    'bioproject': 'N/A',
                    'refseq_accession': 'N/A',
                    'error': 'Query timeout'
                }
        except Exception as e:
            return {
                'accession': accession,
                'status': 'error',
                'assembly_status': 'exception',
                'assembly_level': 'N/A',
                'bioproject': 'N/A',
                'refseq_accession': 'N/A',
                'error': str(e)
            }
    
    # Should not reach here
    return {
        'accession': accession,
        'status': 'error',
        'assembly_status': 'unknown',
        'assembly_level': 'N/A',
        'bioproject': 'N/A',
        'refseq_accession': 'N/A',
        'error': 'Max retries exceeded'
    }


def main():
    """Main execution."""
    import sys

    logger.info("="*70)
    logger.info("UNMAPPED ACCESSION STATUS CHECKER")
    logger.info("="*70)

    # Check command line arguments
    domain_filter = None
    if len(sys.argv) > 1:
        domain_filter = sys.argv[1].lower()
        if domain_filter not in ['archaea', 'bacteria']:
            logger.error(f"Invalid domain: {domain_filter}. Use 'archaea' or 'bacteria'")
            sys.exit(1)

    # Define paths
    output_dir = Path(__file__).parent / 'outputs'

    # Load unmapped files
    archaea_unmapped_file = output_dir / 'archaea_unmapped.csv'
    bacteria_unmapped_file = output_dir / 'bacteria_unmapped.csv'

    logger.info(f"Loading unmapped accessions...")
    archaea_unmapped = pd.read_csv(archaea_unmapped_file)
    bacteria_unmapped = pd.read_csv(bacteria_unmapped_file)

    logger.info(f"  Archaea unmapped: {len(archaea_unmapped):,}")
    logger.info(f"  Bacteria unmapped: {len(bacteria_unmapped):,}")

    # Filter by domain if specified
    if domain_filter == 'archaea':
        logger.info(f"\n🔬 Checking ARCHAEA only")
        all_unmapped = archaea_unmapped[['accession', 'gtdb_domain']].copy()
        output_suffix = '_archaea'
    elif domain_filter == 'bacteria':
        logger.info(f"\n🦠 Checking BACTERIA only")
        all_unmapped = bacteria_unmapped[['accession', 'gtdb_domain']].copy()
        output_suffix = '_bacteria'
    else:
        logger.info(f"\n🔬🦠 Checking ALL domains")
        all_unmapped = pd.concat([
            archaea_unmapped[['accession', 'gtdb_domain']],
            bacteria_unmapped[['accession', 'gtdb_domain']]
        ], ignore_index=True)
        output_suffix = ''

    logger.info(f"\nTotal unmapped accessions to check: {len(all_unmapped):,}")

    # Check status for each accession
    logger.info("\nChecking accession status with NCBI datasets...")
    results = []

    for idx, row in tqdm(all_unmapped.iterrows(), total=len(all_unmapped), desc="Checking accessions"):
        accession = row['accession']
        domain = row['gtdb_domain']

        status_info = check_accession_status(accession)
        status_info['gtdb_domain'] = domain
        results.append(status_info)

        # Rate limiting - be nice to NCBI
        time.sleep(0.5)

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Save results
    status_output = output_dir / f'unmapped_status_report{output_suffix}.csv'
    results_df.to_csv(status_output, index=False)
    logger.info(f"\n✅ Saved status report to: {status_output}")

    # Generate summary statistics
    logger.info("\n" + "="*70)
    logger.info("SUMMARY STATISTICS")
    logger.info("="*70)

    logger.info(f"\nTotal accessions checked: {len(results_df):,}")

    # Status breakdown
    logger.info("\nStatus breakdown:")
    status_counts = results_df['assembly_status'].value_counts()
    for status, count in status_counts.items():
        pct = count / len(results_df) * 100
        logger.info(f"  {status}: {count:,} ({pct:.1f}%)")

    # Domain breakdown
    logger.info("\nBy domain:")
    for domain in ['Archaea', 'Bacteria']:
        domain_df = results_df[results_df['gtdb_domain'] == domain]
        logger.info(f"\n  {domain} ({len(domain_df):,} total):")
        domain_status = domain_df['assembly_status'].value_counts()
        for status, count in domain_status.items():
            pct = count / len(domain_df) * 100
            logger.info(f"    {status}: {count:,} ({pct:.1f}%)")

    # Suppressed accessions
    suppressed = results_df[results_df['assembly_status'] == 'suppressed']
    if len(suppressed) > 0:
        logger.info(f"\n📋 Suppressed accessions: {len(suppressed):,}")
        suppressed_file = output_dir / f'unmapped_suppressed{output_suffix}.csv'
        suppressed.to_csv(suppressed_file, index=False)
        logger.info(f"   Saved to: {suppressed_file}")

    # Not found in NCBI
    not_found = results_df[results_df['assembly_status'] == 'not_in_ncbi']
    if len(not_found) > 0:
        logger.info(f"\n❌ Not found in NCBI: {len(not_found):,}")
        not_found_file = output_dir / f'unmapped_not_found{output_suffix}.csv'
        not_found.to_csv(not_found_file, index=False)
        logger.info(f"   Saved to: {not_found_file}")

    logger.info("\n" + "="*70)
    logger.info("STATUS CHECK COMPLETE!")
    logger.info("="*70)


if __name__ == '__main__':
    main()

