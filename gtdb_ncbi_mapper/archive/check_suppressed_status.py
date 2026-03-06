#!/usr/bin/env python3
"""
Check Suppressed Status of Unmapped Accessions
===============================================

Uses NCBI datasets CLI to check if unmapped GTDB accessions are suppressed in NCBI.
"""

import pandas as pd
import subprocess
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List
import time

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('logs/suppressed_check.log', mode='w')
    ]
)
logger = logging.getLogger(__name__)


def check_accession_status(accession: str) -> Dict:
    """
    Check the status of an accession using NCBI datasets.
    
    Returns:
        Dict with keys: accession, status, assembly_level, exists
    """
    try:
        cmd = f"datasets summary genome accession {accession} --as-json-lines"
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0 and result.stdout.strip():
            data = json.loads(result.stdout)
            
            if 'reports' in data and len(data['reports']) > 0:
                report = data['reports'][0]
                assembly_info = report.get('assembly_info', {})
                
                return {
                    'accession': accession,
                    'exists': True,
                    'status': assembly_info.get('assembly_status', 'unknown'),
                    'assembly_level': assembly_info.get('assembly_level', 'unknown'),
                    'assembly_name': assembly_info.get('assembly_name', 'unknown')
                }
        
        # If we get here, accession doesn't exist or error
        return {
            'accession': accession,
            'exists': False,
            'status': 'not_found',
            'assembly_level': 'N/A',
            'assembly_name': 'N/A'
        }
        
    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout checking {accession}")
        return {
            'accession': accession,
            'exists': None,
            'status': 'timeout',
            'assembly_level': 'N/A',
            'assembly_name': 'N/A'
        }
    except Exception as e:
        logger.warning(f"Error checking {accession}: {e}")
        return {
            'accession': accession,
            'exists': None,
            'status': 'error',
            'assembly_level': 'N/A',
            'assembly_name': 'N/A'
        }


def check_unmapped_accessions(unmapped_file: Path, output_file: Path, 
                              max_check: int = None, delay: float = 0.5):
    """
    Check status of unmapped accessions.
    
    Args:
        unmapped_file: Path to unmapped CSV file
        output_file: Path to save results
        max_check: Maximum number to check (None = all)
        delay: Delay between API calls in seconds
    """
    logger.info(f"Loading unmapped accessions from: {unmapped_file}")
    df = pd.read_csv(unmapped_file)
    
    accessions = df['accession'].tolist()
    total = len(accessions)
    
    if max_check:
        accessions = accessions[:max_check]
        logger.info(f"Checking first {max_check} of {total} accessions")
    else:
        logger.info(f"Checking all {total} accessions")
    
    results = []
    
    for i, accession in enumerate(accessions, 1):
        if i % 10 == 0:
            logger.info(f"Progress: {i}/{len(accessions)}")
        
        status_info = check_accession_status(accession)
        results.append(status_info)
        
        # Rate limiting
        time.sleep(delay)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Merge with original data
    final_df = pd.merge(df, results_df, on='accession', how='left')
    
    # Save results
    final_df.to_csv(output_file, index=False)
    logger.info(f"Results saved to: {output_file}")
    
    # Print summary
    logger.info("\n" + "="*70)
    logger.info("SUMMARY OF UNMAPPED ACCESSION STATUS")
    logger.info("="*70)
    
    if len(results_df) > 0:
        status_counts = results_df['status'].value_counts()
        logger.info("\nStatus breakdown:")
        for status, count in status_counts.items():
            pct = count / len(results_df) * 100
            logger.info(f"  {status}: {count} ({pct:.1f}%)")
    
    return results_df


def main():
    """Main execution."""
    logger.info("="*70)
    logger.info("CHECKING SUPPRESSED STATUS OF UNMAPPED ACCESSIONS")
    logger.info("="*70)
    
    output_dir = Path(__file__).parent / 'outputs'
    
    # Check archaea unmapped (small, check all)
    logger.info("\n### ARCHAEA UNMAPPED ###")
    archaea_results = check_unmapped_accessions(
        unmapped_file=output_dir / 'archaea_unmapped.csv',
        output_file=output_dir / 'archaea_unmapped_with_status.csv',
        max_check=None,  # Check all (only 24)
        delay=0.5
    )
    
    # Check bacteria unmapped (large, check sample)
    logger.info("\n### BACTERIA UNMAPPED (SAMPLE) ###")
    bacteria_results = check_unmapped_accessions(
        unmapped_file=output_dir / 'bacteria_unmapped.csv',
        output_file=output_dir / 'bacteria_unmapped_with_status_sample.csv',
        max_check=100,  # Check first 100
        delay=0.5
    )
    
    logger.info("\n" + "="*70)
    logger.info("STATUS CHECK COMPLETE!")
    logger.info("="*70)


if __name__ == '__main__':
    main()

