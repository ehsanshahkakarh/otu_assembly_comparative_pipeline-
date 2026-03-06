#!/usr/bin/env python3
"""
AI Resolution Cache Manager for 16S Census Parser
==================================================

Manages persistent caching of AI-researched taxonomic resolutions for prokaryotes.
Implements cache-first approach: check cache before doing AI web search.

This module is part of the 16S census parser pipeline and provides
persistent storage for AI-researched taxonomic resolutions for Bacteria and Archaea.
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, List
import logging

logger = logging.getLogger(__name__)


def get_default_cache_path() -> Path:
    """Get the default cache file path for 16S"""
    return Path(__file__).parent.parent.parent / "cache" / "ai_resolutions_16S.json"


class AIResolutionCache:
    """Manages persistent caching of AI-researched taxonomic resolutions for prokaryotes"""
    
    def __init__(self, cache_file: Path):
        self.cache_file = cache_file
        self.cache = self._load_cache()
    
    def _load_cache(self) -> Dict:
        """Load cache from JSON file"""
        if self.cache_file.exists():
            logger.info(f"Loading existing cache from: {self.cache_file}")
            with open(self.cache_file, 'r') as f:
                return json.load(f)
        else:
            logger.info(f"Creating new cache at: {self.cache_file}")
            return {
                "metadata": {
                    "created": datetime.now().isoformat(),
                    "last_updated": datetime.now().isoformat(),
                    "total_resolutions": 0,
                    "version": "1.0",
                    "description": "AI-researched taxonomic resolutions for 16S census data (Bacteria and Archaea)"
                },
                "resolutions": {}
            }
    
    def _save_cache(self):
        """Save cache to JSON file"""
        self.cache["metadata"]["last_updated"] = datetime.now().isoformat()
        self.cache["metadata"]["total_resolutions"] = len(self.cache["resolutions"])
        
        self.cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.cache_file, 'w') as f:
            json.dump(self.cache, f, indent=2)
        
        logger.info(f"Cache saved: {len(self.cache['resolutions'])} resolutions")
    
    def get_resolution(self, taxon_name: str) -> Optional[Dict]:
        """Get cached resolution for a taxon"""
        return self.cache["resolutions"].get(taxon_name)
    
    def has_resolution(self, taxon_name: str) -> bool:
        """Check if taxon has cached resolution"""
        return taxon_name in self.cache["resolutions"]
    
    def add_resolution(
        self,
        taxon_name: str,
        parent_taxid: str,
        parent_name: str,
        lineage: str,
        lineage_ranks: str,
        lineage_taxids: str,
        confidence: float,
        research_notes: str,
        rank: str = "family",
        validated: bool = False
    ):
        """Add a new resolution to cache"""
        self.cache["resolutions"][taxon_name] = {
            "parent_taxid": parent_taxid,
            "parent_name": parent_name,
            "lineage": lineage,
            "lineage_ranks": lineage_ranks,
            "lineage_taxids": lineage_taxids,
            "rank": rank,
            "source": "AI-web-search",
            "confidence": confidence,
            "research_date": datetime.now().isoformat()[:10],
            "research_notes": research_notes,
            "validated": validated,
            "validator": None,
            "validation_date": None
        }
        self._save_cache()
        logger.info(f"Added to cache: {taxon_name} → {parent_name} (confidence: {confidence:.2f})")
    
    def validate_resolution(self, taxon_name: str, validator: str = "human"):
        """Mark a resolution as validated"""
        if taxon_name in self.cache["resolutions"]:
            self.cache["resolutions"][taxon_name]["validated"] = True
            self.cache["resolutions"][taxon_name]["validator"] = validator
            self.cache["resolutions"][taxon_name]["validation_date"] = datetime.now().isoformat()[:10]
            self._save_cache()
            logger.info(f"Validated: {taxon_name}")
    
    def get_unvalidated(self) -> List[str]:
        """Get list of unvalidated resolutions"""
        return [
            name for name, data in self.cache["resolutions"].items()
            if not data.get("validated", False)
        ]
    
    def get_validated(self) -> List[str]:
        """Get list of validated resolutions"""
        return [
            name for name, data in self.cache["resolutions"].items()
            if data.get("validated", False)
        ]
    
    def get_statistics(self) -> Dict:
        """Get cache statistics"""
        total = len(self.cache["resolutions"])
        validated = sum(1 for r in self.cache["resolutions"].values() if r.get("validated", False))
        
        # Confidence distribution
        confidences = [r.get("confidence", 0) for r in self.cache["resolutions"].values()]
        avg_confidence = sum(confidences) / len(confidences) if confidences else 0
        
        return {
            "total_resolutions": total,
            "validated": validated,
            "unvalidated": total - validated,
            "average_confidence": avg_confidence,
            "last_updated": self.cache["metadata"]["last_updated"]
        }
    
    def export_to_known_parents_format(self) -> Dict[str, tuple]:
        """
        Export validated resolutions to known_parents.py format
        
        Returns:
            Dictionary in format: {taxon_name: (parent_taxid, parent_name, notes, rank)}
        """
        known_parents = {}
        
        for taxon_name, data in self.cache["resolutions"].items():
            if data.get("validated", False):
                known_parents[taxon_name] = (
                    data["parent_taxid"],
                    data["parent_name"],
                    data["research_notes"],
                    data.get("rank", "family")
                )
        
        return known_parents

