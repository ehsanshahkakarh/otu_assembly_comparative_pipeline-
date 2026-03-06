#!/usr/bin/env python3
"""
Cached AI-Assisted Taxonomic Resolver for 16S Census Parser
============================================================

Resolves unmapped prokaryotic taxonomic names using:
1. Cache-first approach (check ai_resolutions_16S.json)
2. Manual database (known_parents.py)
3. AI web search (for new unmapped taxa)
4. Persistent caching for future runs

Focus: Bacteria and Archaea (16S rRNA)
"""

from pathlib import Path
from typing import Dict, List, Optional
import logging
import subprocess
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from ai_resolution_cache import AIResolutionCache, get_default_cache_path
from known_parents import KNOWN_PARENTS

logger = logging.getLogger(__name__)


class CachedAIResolver:
    """AI resolver with persistent caching and manual database fallback for prokaryotes"""

    def __init__(self, cache_file: Path = None):
        # Use default cache path if not provided
        if cache_file is None:
            cache_file = get_default_cache_path()
        
        self.cache = AIResolutionCache(cache_file)

        # Import manual database
        self.manual_db = KNOWN_PARENTS
        logger.info(f"Loaded manual database: {len(self.manual_db)} entries")
    
    def resolve_unmapped_taxa(
        self,
        unmapped_taxa: List[str],
        rank: str = "family"
    ) -> Dict[str, Dict]:
        """
        Resolve unmapped taxa using multi-tier approach:
        1. Check AI cache
        2. Check manual database
        3. Perform AI research (placeholder for now)

        Args:
            unmapped_taxa: List of unmapped taxon names
            rank: Taxonomic rank (phylum, family, or genus)

        Returns:
            Dictionary of resolutions
        """
        resolutions = {}
        needs_research = []

        logger.info(f"\n{'='*80}")
        logger.info(f"RESOLVING {len(unmapped_taxa)} UNMAPPED {rank.upper()} TAXA (16S - Prokaryotes)")
        logger.info(f"{'='*80}")

        # Step 1: Check AI cache
        logger.info(f"\n[Step 1/3] Checking AI cache...")
        for taxon in unmapped_taxa:
            cached = self.cache.get_resolution(taxon)
            if cached:
                logger.info(f"  ✅ Found in AI cache: {taxon}")
                resolutions[taxon] = cached

        # Step 2: Check manual database
        logger.info(f"\n[Step 2/3] Checking manual database...")
        logger.info(f"  Manual DB has {len(self.manual_db)} entries")

        for taxon in unmapped_taxa:
            if taxon in resolutions:
                continue  # Already found in cache

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
        
        # Step 3: Identify taxa needing research
        logger.info(f"\n[Step 3/3] Identifying taxa needing AI research...")
        for taxon in unmapped_taxa:
            if taxon not in resolutions:
                needs_research.append(taxon)
                logger.info(f"  🔍 Needs research: {taxon}")

        # Summary
        logger.info(f"\n{'='*80}")
        logger.info(f"RESOLUTION SUMMARY")
        logger.info(f"{'='*80}")
        logger.info(f"Total unmapped: {len(unmapped_taxa)}")
        logger.info(f"  From AI cache: {sum(1 for r in resolutions.values() if r.get('source') == 'AI-web-search')}")
        logger.info(f"  From manual DB: {sum(1 for r in resolutions.values() if r.get('source') == 'manual-database')}")
        logger.info(f"  Need research: {len(needs_research)}")

        if needs_research:
            logger.info(f"\n⚠️  {len(needs_research)} taxa need AI web research")
            logger.info(f"   Generating research prompts...")

            # Generate AI research prompts for unmapped taxa
            self._generate_research_prompts(needs_research, rank=rank)

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
            # Use taxonkit with built-in taxdump
            result = subprocess.run(
                ['taxonkit', 'lineage', '-R', '-t'],
                input=parent_taxid,
                capture_output=True,
                text=True
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

    def _generate_research_prompts(self, taxa_list: List[str], rank: str = "family"):
        """Generate AI research prompts for unmapped taxa"""
        logger.info(f"\n{'='*80}")
        logger.info(f"AI RESEARCH PROMPTS FOR {len(taxa_list)} UNMAPPED TAXA")
        logger.info(f"{'='*80}")

        for taxon in taxa_list:
            logger.info(f"\n🔬 RESEARCH NEEDED: {taxon} ({rank})")
            logger.info(f"   Context: Prokaryotic taxon (Bacteria or Archaea) from 16S rRNA survey")
            logger.info(f"   Task: Find parent taxon and NCBI TaxID")
            logger.info(f"   Possible issues:")
            logger.info(f"     - GTDB name not in NCBI (e.g., Lokiarchaeia → Asgardarchaeota)")
            logger.info(f"     - CPR group (Candidate Phyla Radiation)")
            logger.info(f"     - Misspelling or alternative name")
            logger.info(f"     - Environmental clade without formal taxonomy")
            logger.info(f"   Search terms: '{taxon}', '{taxon} taxonomy', '{taxon} NCBI'")
            logger.info(f"   ---")

    def add_manual_resolution(
        self,
        taxon_name: str,
        parent_taxid: str,
        parent_name: str,
        confidence: float,
        research_notes: str,
        rank: str = "family",
        validated: bool = False
    ):
        """
        Manually add a resolution to the cache after AI research.

        Args:
            taxon_name: Name of the unmapped taxon
            parent_taxid: NCBI TaxID of the parent taxon
            parent_name: Name of the parent taxon
            confidence: Confidence level (0.0-1.0)
            research_notes: Notes about the research/resolution
            rank: Taxonomic rank
            validated: Whether this has been validated
        """
        # Get lineage from taxonkit
        lineage_info = self._get_lineage_from_taxonkit(parent_taxid, taxon_name, rank)

        if lineage_info:
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
            logger.info(f"✅ Added resolution: {taxon_name} → {parent_name}")
        else:
            logger.error(f"❌ Failed to get lineage for {taxon_name}")

    def get_cache_statistics(self) -> Dict:
        """Get statistics about the cache"""
        return self.cache.get_statistics()


if __name__ == "__main__":
    # Test the resolver
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Create resolver
    resolver = CachedAIResolver()

    # Test with some unmapped taxa from the log
    test_taxa = [
        "Lokiarchaeia",
        "Microgenomatia",
        "ABY",
        "WWE",
        "Candidatus Gracilibacteria"
    ]

    logger.info("Testing 16S AI Resolver")
    logger.info("=" * 80)

    # Resolve taxa
    resolutions = resolver.resolve_unmapped_taxa(test_taxa, rank="family")

    # Print results
    logger.info("\n" + "=" * 80)
    logger.info("RESOLUTION RESULTS")
    logger.info("=" * 80)
    for taxon, resolution in resolutions.items():
        logger.info(f"\n{taxon}:")
        logger.info(f"  Parent: {resolution['parent_name']} (TaxID: {resolution['parent_taxid']})")
        logger.info(f"  Source: {resolution['source']}")
        logger.info(f"  Confidence: {resolution['confidence']}")
        logger.info(f"  Lineage: {resolution['lineage']}")

    # Print cache statistics
    stats = resolver.get_cache_statistics()
    logger.info("\n" + "=" * 80)
    logger.info("CACHE STATISTICS")
    logger.info("=" * 80)
    for key, value in stats.items():
        logger.info(f"  {key}: {value}")

