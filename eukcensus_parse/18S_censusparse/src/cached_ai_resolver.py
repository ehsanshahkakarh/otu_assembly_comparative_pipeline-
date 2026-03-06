#!/usr/bin/env python3
"""
Cached AI-Assisted Taxonomic Resolver
======================================

Resolves unmapped taxonomic names using:
1. Cache-first approach (check ai_resolutions.json)
2. Manual database (known_parents.py)
3. AI web search (for new unmapped taxa)
4. Persistent caching for future runs
"""

from pathlib import Path
from typing import Dict, List, Optional
import logging
import subprocess
import sys
import os
import json

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from ai_resolution_cache import AIResolutionCache
from known_parents import KNOWN_PARENTS
from metadata_resolver import MetadataResolver

logger = logging.getLogger(__name__)


class CachedAIResolver:
    """AI resolver with persistent caching, metadata parsing, and manual database fallback"""

    def __init__(self, cache_file: Path, metadata_file: Path = None, taxonkit_db: str = None):
        self.cache = AIResolutionCache(cache_file)

        # Import manual database
        self.manual_db = KNOWN_PARENTS
        logger.info(f"Loaded manual database: {len(self.manual_db)} entries")

        # Initialize metadata resolver if metadata file provided
        self.metadata_resolver = None
        if metadata_file and taxonkit_db:
            self.metadata_resolver = MetadataResolver(metadata_file, taxonkit_db)
            logger.info(f"Initialized metadata resolver")
    
    def resolve_unmapped_taxa(
        self,
        unmapped_taxa: List[str],
        rank: str = "family"
    ) -> Dict[str, Dict]:
        """
        Resolve unmapped taxa using multi-tier approach:
        1. Check AI cache
        2. Check metadata (division/family/genus columns)
        3. Check manual database
        4. Perform AI research (placeholder for now)

        Args:
            unmapped_taxa: List of unmapped taxon names
            rank: Taxonomic rank (family or genus)

        Returns:
            Dictionary of resolutions
        """
        resolutions = {}
        needs_research = []

        logger.info(f"\n{'='*80}")
        logger.info(f"RESOLVING {len(unmapped_taxa)} UNMAPPED {rank.upper()} TAXA")
        logger.info(f"{'='*80}")

        # Step 1: Check AI cache
        logger.info(f"\n[Step 1/4] Checking AI cache...")
        for taxon in unmapped_taxa:
            cached = self.cache.get_resolution(taxon)
            if cached:
                logger.info(f"  ✅ Found in AI cache: {taxon}")
                resolutions[taxon] = cached

        # Step 2: Check metadata columns (division/family/genus)
        if self.metadata_resolver:
            logger.info(f"\n[Step 2/4] Checking metadata columns...")
            for taxon in unmapped_taxa:
                if taxon in resolutions:
                    continue  # Already found in cache

                metadata_result = self.metadata_resolver.resolve_from_metadata(taxon)
                if metadata_result:
                    logger.info(f"  ✅ Found in metadata: {taxon} → {metadata_result['parent_name']} (TaxID: {metadata_result['parent_taxid']}) via {metadata_result['source']}")

                    # Get full lineage from taxonkit
                    lineage_info = self._get_lineage_from_taxonkit(
                        metadata_result['parent_taxid'],
                        taxon,
                        rank
                    )

                    if lineage_info:
                        resolutions[taxon] = {
                            "parent_taxid": metadata_result['parent_taxid'],
                            "parent_name": metadata_result['parent_name'],
                            "lineage": lineage_info["lineage"],
                            "lineage_ranks": lineage_info["lineage_ranks"],
                            "lineage_taxids": lineage_info["lineage_taxids"],
                            "rank": rank,
                            "source": metadata_result['source'],
                            "confidence": 0.9,
                            "research_notes": f"Resolved from metadata {metadata_result['source']}"
                        }
        else:
            logger.info(f"\n[Step 2/4] Metadata resolver not initialized (skipping)")

        # Step 3: Check manual database
        logger.info(f"\n[Step 3/4] Checking manual database...")
        logger.info(f"  Manual DB has {len(self.manual_db)} entries")
        logger.info(f"  Sample manual DB keys: {list(self.manual_db.keys())[:5]}")

        for taxon in unmapped_taxa:
            if taxon in resolutions:
                continue  # Already found in cache

            logger.debug(f"  Checking: {taxon} in manual DB...")

            if taxon in self.manual_db:
                parent_taxid, parent_name, notes, db_rank = self.manual_db[taxon]
                logger.info(f"  ✅ Found in manual DB: {taxon} → {parent_name} (TaxID: {parent_taxid})")

                # Get lineage from taxonkit
                lineage_info = self._get_lineage_from_taxonkit(parent_taxid, taxon, db_rank)

                if lineage_info:
                    resolutions[taxon] = {
                        "parent_taxid": parent_taxid,
                        "parent_name": parent_name,
                        "lineage": lineage_info["lineage"],
                        "lineage_ranks": lineage_info["lineage_ranks"],
                        "lineage_taxids": lineage_info["lineage_taxids"],
                        "rank": db_rank,
                        "source": "manual-database",
                        "confidence": 1.0,
                        "research_notes": notes,
                        "validated": True
                    }
                    logger.info(f"     Lineage: {lineage_info['lineage']}")
                else:
                    logger.warning(f"     Failed to get lineage from taxonkit for {taxon}")
        
        # Step 4: Identify taxa needing research
        logger.info(f"\n[Step 4/4] Identifying taxa needing AI research...")
        for taxon in unmapped_taxa:
            if taxon not in resolutions:
                needs_research.append(taxon)
                logger.info(f"  🔍 Needs research: {taxon}")

        # Summary
        logger.info(f"\n{'='*80}")
        logger.info(f"RESOLUTION SUMMARY")
        logger.info(f"{'='*80}")
        logger.info(f"Total unmapped: {len(unmapped_taxa)}")
        logger.info(f"  From AI cache: {sum(1 for r in resolutions.values() if r['source'] == 'AI-web-search')}")
        logger.info(f"  From metadata: {sum(1 for r in resolutions.values() if 'metadata' in r['source'])}")
        logger.info(f"  From manual DB: {sum(1 for r in resolutions.values() if r['source'] == 'manual-database')}")
        logger.info(f"  Need research: {len(needs_research)}")

        if needs_research:
            logger.info(f"\n⚠️  {len(needs_research)} taxa need AI web research")
            logger.info(f"   Generating research prompts...")

            # Generate AI research prompts for unmapped taxa
            self.research_taxa_with_ai(needs_research, rank=rank, interactive=False)

            logger.info(f"\n📝 NEXT STEPS:")
            logger.info(f"   1. Use the AI prompts above to research each taxon")
            logger.info(f"   2. Add findings to the cache using add_manual_resolution()")
            logger.info(f"   3. Re-run this script to apply cached resolutions")

        return resolutions
    
    def _get_lineage_from_taxonkit(
        self,
        parent_taxid: str,
        taxon_name: str,
        rank: str
    ) -> Optional[Dict]:
        """Get lineage from taxonkit and append taxon name"""
        try:
            # Use taxonkit with built-in taxdump (no external database needed)
            env = os.environ.copy()

            result = subprocess.run(
                ['taxonkit', 'lineage', '-R', '-t'],
                input=parent_taxid,
                capture_output=True,
                text=True,
                env=env
            )

            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 3:
                    lineage = parts[1]
                    lineage_taxids = parts[2]
                    lineage_ranks = parts[3] if len(parts) > 3 else ""

                    # Append taxon name
                    full_lineage = f"{lineage};{taxon_name}"
                    full_taxids = f"{lineage_taxids};NA"
                    full_ranks = f"{lineage_ranks};{rank}"

                    return {
                        "lineage": full_lineage,
                        "lineage_taxids": full_taxids,
                        "lineage_ranks": full_ranks
                    }
        except Exception as e:
            logger.error(f"Error getting lineage for {taxon_name}: {e}")

        return None

    def _ai_search_for_parent(
        self,
        taxon_name: str,
        rank: str = "family"
    ) -> Optional[Dict]:
        """
        Use AI to search for parent taxon information.

        This creates a prompt file for the AI to research the taxon and find:
        1. The parent taxon name
        2. The parent taxon NCBI TaxID
        3. Confidence level
        4. Research notes

        Args:
            taxon_name: Name of the unmapped taxon
            rank: Taxonomic rank (family or genus)

        Returns:
            Dictionary with parent info or None if not found
        """
        logger.info(f"\n🤖 AI SEARCH NEEDED: {taxon_name}")
        logger.info(f"   Please research this {rank} and provide parent taxon information")
        logger.info(f"   The AI will be prompted to search for taxonomic information")

        # Create a prompt for the AI
        prompt = f"""
TAXONOMIC RESEARCH REQUEST
==========================

Taxon Name: {taxon_name}
Rank: {rank}

TASK:
Please search for information about this taxonomic name and identify:
1. What is the parent taxon (the taxonomic group this belongs to)?
2. What is the NCBI TaxID of the parent taxon?
3. What is the full taxonomic lineage?

SEARCH STRATEGY:
- Search for "{taxon_name}" + "taxonomy" + "18S rRNA"
- Look for scientific publications, taxonomic databases
- Check if this is an environmental clade, phylogenetic group, or formal taxon
- For environmental clades (e.g., "MAST-12", "Dino-Group-II"), find the broader taxonomic group

REQUIRED OUTPUT FORMAT (JSON):
{{
    "taxon_name": "{taxon_name}",
    "parent_taxid": "NCBI_TAXID_HERE",
    "parent_name": "PARENT_TAXON_NAME",
    "confidence": 0.0-1.0,
    "research_notes": "Brief explanation of findings and sources",
    "lineage": "Full lineage from parent (will be fetched from taxonkit)",
    "rank": "{rank}"
}}

If you cannot find reliable information, set confidence to 0.0 and explain in research_notes.
"""

        logger.info(f"\n{'='*80}")
        logger.info("AI PROMPT:")
        logger.info(f"{'='*80}")
        logger.info(prompt)
        logger.info(f"{'='*80}\n")

        # For now, return None to indicate manual intervention needed
        # In a full implementation, this would call an AI API
        return None

    def research_taxa_with_ai(
        self,
        taxa_list: List[str],
        rank: str = "family",
        interactive: bool = True
    ) -> Dict[str, Dict]:
        """
        Research multiple taxa using AI assistance.

        Args:
            taxa_list: List of taxon names to research
            rank: Taxonomic rank
            interactive: If True, prompt user for each taxon

        Returns:
            Dictionary of resolutions
        """
        resolutions = {}

        logger.info(f"\n{'='*80}")
        logger.info(f"AI RESEARCH MODE: {len(taxa_list)} taxa")
        logger.info(f"{'='*80}")

        for taxon in taxa_list:
            logger.info(f"\n--- Researching: {taxon} ---")

            if interactive:
                # Prompt user to provide information
                logger.info(f"\nPlease research '{taxon}' and provide the following:")
                logger.info(f"  1. Parent taxon name")
                logger.info(f"  2. Parent NCBI TaxID")
                logger.info(f"  3. Confidence (0.0-1.0)")
                logger.info(f"  4. Research notes")
                logger.info(f"\nEnter 'skip' to skip this taxon")
                logger.info(f"Enter 'quit' to stop research session")

                # This would be replaced with actual AI API call or user input
                # For now, just log the prompt
                self._ai_search_for_parent(taxon, rank)
            else:
                # Non-interactive mode: just generate prompts
                self._ai_search_for_parent(taxon, rank)

        return resolutions

    def add_manual_resolution(
        self,
        taxon_name: str,
        parent_taxid: str,
        parent_name: str,
        confidence: float,
        research_notes: str,
        rank: str = "family",
        validated: bool = False
    ) -> bool:
        """
        Manually add a resolution to the cache after research.

        Args:
            taxon_name: Name of the unmapped taxon
            parent_taxid: NCBI TaxID of the parent taxon
            parent_name: Name of the parent taxon
            confidence: Confidence level (0.0-1.0)
            research_notes: Notes about the research/sources
            rank: Taxonomic rank
            validated: Whether this has been validated

        Returns:
            True if successful, False otherwise
        """
        try:
            # Get lineage from taxonkit
            lineage_info = self._get_lineage_from_taxonkit(parent_taxid, taxon_name, rank)

            if not lineage_info:
                logger.error(f"Failed to get lineage for parent TaxID: {parent_taxid}")
                return False

            # Add to cache
            self.cache.add_resolution(
                taxon_name=taxon_name,
                parent_taxid=parent_taxid,
                parent_name=parent_name,
                lineage=lineage_info["lineage"],
                lineage_ranks=lineage_info["lineage_ranks"],
                lineage_taxids=lineage_info["lineage_taxids"],
                confidence=confidence,
                research_notes=research_notes,
                rank=rank,
                validated=validated
            )

            logger.info(f"✅ Added resolution to cache: {taxon_name} → {parent_name}")
            return True

        except Exception as e:
            logger.error(f"Error adding manual resolution: {e}")
            return False

