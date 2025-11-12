#!/usr/bin/env python3
"""
Demonstration of enhanced organelle handling for EukCensus data.

This script shows how the enhanced parser handles real organelle entries
from the EukCensus dataset, demonstrating the improvements over the original parser.
"""

import sys
sys.path.append('.')

from enhanced_eukcensus_parser import (
    detect_organelle_type,
    clean_organelle_taxon_name,
    extract_genus_from_species,
    extract_appropriate_rank_name
)

def demo_real_organelle_entries():
    """Demonstrate handling of real organelle entries from the dataset."""
    
    print("🧬 Enhanced Organelle Handling Demonstration")
    print("=" * 50)
    print("Using real entries from EukCensus 16S dataset\n")
    
    # Real organelle entries from the dataset
    real_entries = [
        "Annulohypoxylon_stygium.Mitochondria",
        "Pyronema_omphalodes.Mitochondria", 
        "Leptographium_lundbergii.Mitochondria",
        "Aspergillus_nidulans_FGSC_A4.Mitochondria",
        "uncultured_bacterium.Mitochondria"
    ]
    
    print("Processing real organelle entries from the dataset:")
    print("-" * 55)
    
    for entry in real_entries:
        print(f"\n📋 Original entry: {entry}")
        
        # Step 1: Detect organelle
        is_organelle, org_type, base_name = detect_organelle_type(entry)
        print(f"   🔍 Organelle detected: {is_organelle} ({org_type})")
        print(f"   🧹 Base organism: {base_name}")
        
        # Step 2: Clean the name
        cleaned = clean_organelle_taxon_name(entry)
        print(f"   ✨ Cleaned name: {cleaned}")
        
        # Step 3: Extract genus for genus-level parsing
        genus = extract_genus_from_species(cleaned)
        print(f"   🏷️  Genus extracted: {genus}")
        
        # Step 4: Show rank-appropriate processing
        for rank in ['genus', 'family', 'phylum']:
            appropriate = extract_appropriate_rank_name(entry, rank)
            status = "✅ KEEP" if appropriate else "❌ FILTER"
            result = appropriate if appropriate else "FILTERED OUT"
            print(f"   📊 {rank.capitalize()} level: {result} {status}")
        
        print()

def compare_with_original():
    """Compare enhanced vs original handling."""
    
    print("\n" + "=" * 60)
    print("🔄 Comparison: Enhanced vs Original Processing")
    print("=" * 60)
    
    test_entry = "Vitis_vinifera:plas.Chloroplast"
    
    print(f"Test entry: {test_entry}")
    print()
    
    # Enhanced processing
    print("🚀 Enhanced Parser:")
    is_organelle, org_type, base_name = detect_organelle_type(test_entry)
    cleaned = clean_organelle_taxon_name(test_entry)
    genus = extract_genus_from_species(cleaned)
    
    print(f"   • Organelle detection: {is_organelle} ({org_type})")
    print(f"   • Cleaned name: {cleaned}")
    print(f"   • Genus extraction: {genus}")
    print(f"   • Result for genus parsing: {genus} ✅")
    
    print()
    print("📜 Original Parser (simulated):")
    print(f"   • Basic cleaning: {test_entry.replace('_', ' ')}")
    print(f"   • No organelle detection")
    print(f"   • No genus extraction")
    print(f"   • Result: Vitis vinifera:plas.Chloroplast ❌")

def show_taxonomic_filtering():
    """Demonstrate taxonomic rank filtering."""
    
    print("\n" + "=" * 60)
    print("🎯 Taxonomic Rank Filtering Demonstration")
    print("=" * 60)
    
    test_cases = [
        ("Annulohypoxylon_stygium.Mitochondria", "Species with organelle"),
        ("Bacillus_subtilis", "Species without organelle"),
        ("Proteobacteria", "Phylum-level entry"),
        ("Enterobacteriaceae", "Family-level entry")
    ]
    
    for entry, description in test_cases:
        print(f"\n📋 {description}: {entry}")
        
        for target_rank in ['genus', 'family', 'phylum']:
            result = extract_appropriate_rank_name(entry, target_rank)
            if result:
                print(f"   ✅ {target_rank.capitalize()}: {result}")
            else:
                print(f"   ❌ {target_rank.capitalize()}: FILTERED (inappropriate rank)")

def main():
    """Run the demonstration."""
    
    demo_real_organelle_entries()
    compare_with_original()
    show_taxonomic_filtering()
    
    print("\n" + "=" * 60)
    print("✅ Demonstration Complete!")
    print("=" * 60)
    print("\nKey Benefits of Enhanced Parser:")
    print("• 🧬 Proper organelle detection and cleaning")
    print("• 🏷️  Accurate genus extraction from species names")
    print("• 🎯 Intelligent taxonomic rank filtering")
    print("• 📊 Better data quality for downstream analysis")
    print("• 🔍 Detailed logging and traceability")
    
    print(f"\nTo run the enhanced parser on the full dataset:")
    print(f"python enhanced_eukcensus_parser.py eukcensus_16S.clusters.97.tsv")

if __name__ == "__main__":
    main()