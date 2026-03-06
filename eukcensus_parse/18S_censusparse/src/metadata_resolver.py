#!/usr/bin/env python3
"""
Metadata-based taxonomic resolver for 18S census data.

This module resolves unmapped taxa by checking their division/family/genus
columns in the metadata and looking up those values in NCBI taxonomy.

Resolution hierarchy:
1. Check division column → lookup in NCBI
2. If division unmapped, check family column → lookup in NCBI
3. If family unmapped, check genus column → lookup in NCBI
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple
import pandas as pd

logger = logging.getLogger(__name__)


class MetadataResolver:
    """Resolve unmapped taxa using metadata division/family/genus columns"""
    
    def __init__(self, metadata_file: Path, taxonkit_db: str):
        """
        Initialize metadata resolver
        
        Args:
            metadata_file: Path to eukcensus metadata TSV
            taxonkit_db: Path to taxonkit database
        """
        self.metadata_file = metadata_file
        self.taxonkit_db = taxonkit_db
        self.metadata_df = None
        self._load_metadata()
    
    def _load_metadata(self):
        """Load metadata file"""
        logger.info(f"Loading metadata from {self.metadata_file}")
        self.metadata_df = pd.read_csv(
            self.metadata_file,
            sep='\t',
            usecols=['centroid', 'division', 'family', 'genus']
        )
        logger.info(f"Loaded {len(self.metadata_df)} metadata entries")
    
    def _lookup_taxid(self, taxon_name: str) -> Optional[Tuple[str, str, str]]:
        """
        Look up a taxon name in NCBI via taxonkit
        
        Args:
            taxon_name: Taxonomic name to look up
            
        Returns:
            Tuple of (taxid, name, rank) if found, None otherwise
        """
        if not taxon_name or pd.isna(taxon_name):
            return None
        
        # Skip obviously unmapped names
        if any(x in str(taxon_name) for x in ['.U.', '_X', '_XX', '_XXX', 'XXXX', '-lineage', '-Clade-']):
            return None
        
        try:
            cmd = f'echo "{taxon_name}" | taxonkit name2taxid --data-dir {self.taxonkit_db} | taxonkit lineage -r -n -L --data-dir {self.taxonkit_db}'
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 3 and parts[1]:  # Has TaxID
                    taxid = parts[1]
                    name = parts[2] if len(parts) > 2 else taxon_name
                    rank = parts[3] if len(parts) > 3 else "unknown"
                    logger.debug(f"Found: {taxon_name} → TaxID {taxid} ({rank})")
                    return (taxid, name, rank)
        except Exception as e:
            logger.debug(f"Lookup failed for {taxon_name}: {e}")
        
        return None
    
    def resolve_from_metadata(self, unmapped_taxon: str) -> Optional[Dict[str, str]]:
        """
        Resolve an unmapped taxon using its metadata columns
        
        Args:
            unmapped_taxon: The unmapped family/genus name
            
        Returns:
            Dict with resolution info if found, None otherwise
        """
        # Find rows where this taxon appears in family or genus column
        matches = self.metadata_df[
            (self.metadata_df['family'] == unmapped_taxon) |
            (self.metadata_df['genus'] == unmapped_taxon)
        ]
        
        if matches.empty:
            logger.debug(f"No metadata found for {unmapped_taxon}")
            return None
        
        # Get the first match's division/family/genus
        row = matches.iloc[0]
        division = row['division']
        family = row['family']
        genus = row['genus']
        
        logger.info(f"Checking metadata for {unmapped_taxon}: division={division}, family={family}, genus={genus}")
        
        # Try division first
        if division and division != unmapped_taxon:
            result = self._lookup_taxid(division)
            if result:
                taxid, name, rank = result
                return {
                    'parent_taxid': taxid,
                    'parent_name': name,
                    'parent_rank': rank,
                    'source': 'metadata_division',
                    'metadata_division': division
                }
        
        # Try family if division didn't work
        if family and family != unmapped_taxon:
            result = self._lookup_taxid(family)
            if result:
                taxid, name, rank = result
                return {
                    'parent_taxid': taxid,
                    'parent_name': name,
                    'parent_rank': rank,
                    'source': 'metadata_family',
                    'metadata_family': family
                }
        
        # Try genus as last resort
        if genus and genus != unmapped_taxon:
            result = self._lookup_taxid(genus)
            if result:
                taxid, name, rank = result
                return {
                    'parent_taxid': taxid,
                    'parent_name': name,
                    'parent_rank': rank,
                    'source': 'metadata_genus',
                    'metadata_genus': genus
                }
        
        logger.debug(f"No valid parent found in metadata for {unmapped_taxon}")
        return None

