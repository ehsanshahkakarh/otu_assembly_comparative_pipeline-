#!/usr/bin/env python3
"""
Test script to verify enhanced taxonkit processing with system version.

This script tests the problematic taxa mentioned by the user to verify
that the enhanced parser with system taxonkit will resolve them correctly.
"""

import subprocess
import os
from collections import defaultdict

def test_single_taxid_with_fallbacks(taxon_name, env):
    """
    Test the enhanced fallback logic for a single taxon name.
    """
    print(f"\n🔍 Testing: {taxon_name}")
    
    fallback_strategies = []
    
    # Strategy 1: Original name (especially important for Candidatus)
    fallback_strategies.append(("original", taxon_name))
    
    # Strategy 2: For Candidatus names, try without the prefix
    if taxon_name.startswith("Candidatus "):
        candidatus_stripped = taxon_name[len("Candidatus "):]
        fallback_strategies.append(("candidatus_stripped", candidatus_stripped))
    
    # Strategy 3: Extract meaningful taxonomic part for complex uncultured names
    if 'uncultured' in taxon_name.lower():
        parts = taxon_name.replace('_', ' ').split()
        generic_terms = ['uncultured', 'bacterium', 'organism', 'eukaryote', 'sp', 'species']
        meaningful_parts = []
        
        for part in parts:
            if part.lower() not in generic_terms and len(part) > 2:
                if part[0].isupper() and part[1:].islower():
                    meaningful_parts.append(part)
                elif part.endswith('bacteria') or part.endswith('proteobacteria'):
                    meaningful_parts.append(part)
        
        if meaningful_parts:
            meaningful_part = meaningful_parts[-1] if len(meaningful_parts) == 1 else ' '.join(meaningful_parts[:2])
            fallback_strategies.append(("meaningful_part", meaningful_part))
    
    # Try each strategy
    for method, name_to_try in fallback_strategies:
        if not name_to_try or name_to_try.lower() in ['uncultured', 'unknown', 'environmental']:
            print(f"   ❌ {method}: '{name_to_try}' - skipped (generic term)")
            continue
            
        try:
            result = subprocess.run(
                ["taxonkit", "name2taxid"],
                input=name_to_try,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
                cwd=".",
                timeout=5
            )
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 2 and parts[1] != "0" and parts[1].strip():
                    print(f"   ✅ {method}: '{name_to_try}' → taxid {parts[1]}")
                    return parts[1], method
                else:
                    print(f"   ❌ {method}: '{name_to_try}' → no taxid")
                    
        except Exception as e:
            print(f"   ❌ {method}: '{name_to_try}' → error: {e}")
            continue
    
    print(f"   ❌ All strategies failed for '{taxon_name}'")
    return "NA", "failed_all_strategies"

def main():
    """Test problematic taxa from the user's log."""
    
    print("🧪 Testing Enhanced Taxonkit Processing")
    print("=" * 50)
    
    # Use default environment
    env = os.environ.copy()
    
    # Test cases from the user's unmapped log
    test_cases = [
        # Procabacter issue
        "Procabacter",
        
        # Candidatus issues
        "Candidatus Sumerlaea",
        "Candidatus Cardinium", 
        "Candidatus Blochmannia",
        "Candidatus Ishikawaella",
        "Candidatus Atelocyanobacterium",
        "Candidatus Lokiarchaeum",
        "Candidatus Phytoplasma",
        "Candidatus Baumannia",
        
        # Complex uncultured names with heavy underscores
        "uncultured_Alphaproteobacteria_bacterium",
        "uncultured_Rickettsia_sp",
        "uncultured_Rickettsiales_bacterium",
        "uncultured_deep-sea_bacterium",
        "uncultured_proteobacterium",
        
        # Other problematic names
        "Raphid-pennate_X_sp.",
        "Polar-centric-Mediophyceae_X_sp.",
        "Embryophyceae_XXX_sp.",
    ]
    
    results = defaultdict(list)
    
    for taxon_name in test_cases:
        taxid, method = test_single_taxid_with_fallbacks(taxon_name, env)
        results[method].append((taxon_name, taxid))
    
    # Summary
    print("\n📊 SUMMARY")
    print("=" * 50)
    
    total_tested = len(test_cases)
    successful = sum(len(names) for method, names in results.items() if method != "failed_all_strategies")
    failed = len(results.get("failed_all_strategies", []))
    
    print(f"Total tested: {total_tested}")
    print(f"Successfully resolved: {successful} ({successful/total_tested*100:.1f}%)")
    print(f"Failed to resolve: {failed} ({failed/total_tested*100:.1f}%)")
    
    print("\nSuccess by method:")
    for method, names in results.items():
        if method != "failed_all_strategies":
            print(f"  {method}: {len(names)} taxa")
    
    if results.get("failed_all_strategies"):
        print(f"\nStill failing:")
        for name, _ in results["failed_all_strategies"]:
            print(f"  - {name}")

if __name__ == "__main__":
    main()
