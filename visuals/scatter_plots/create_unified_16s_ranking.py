#!/usr/bin/env python3
"""
Create Unified 16S Prokaryotic Ranking from TaxID-Based Results
==============================================================
Parse the TaxID-based 16S results and create a unified ranking of all prokaryotic groups.
"""

import re
import os
from datetime import datetime

def parse_taxonomic_groups_from_file(filepath):
    """Parse taxonomic group data from the TaxID-based results file."""
    
    if not os.path.exists(filepath):
        print(f"❌ File not found: {filepath}")
        return {}
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find all taxonomic group sections
    group_pattern = r'🏷️\s+TAXONOMIC GROUP:\s+([^(]+)\s+\(TaxID:\s+(\d+)\)'
    groups = {}
    
    # Split content by taxonomic groups
    sections = re.split(r'🏷️\s+TAXONOMIC GROUP:', content)
    
    for section in sections[1:]:  # Skip first empty section
        lines = section.strip().split('\n')
        if not lines:
            continue
            
        # Extract group name and TaxID from first line
        first_line = lines[0]
        match = re.search(r'([^(]+)\s+\(TaxID:\s+(\d+)\)', first_line)
        if not match:
            continue
            
        group_name = match.group(1).strip()
        taxid = int(match.group(2))
        
        # Extract genome count and diversity info
        genomes = 0
        genera_count = 0
        species_count = 0
        top_genus = ""
        concentration = 0.0
        concentration_level = ""
        
        for line in lines:
            # Total genomes
            if "Total genomes:" in line:
                match = re.search(r'Total genomes:\s+(\d+)', line)
                if match:
                    genomes = int(match.group(1))
            
            # Diversity info
            if "Genera:" in line and "Species:" in line:
                match = re.search(r'Genera:\s+(\d+).*Species:\s+(\d+)', line)
                if match:
                    genera_count = int(match.group(1))
                    species_count = int(match.group(2))
            
            # Top genus info
            if "TOP GENERA (by genome count):" in line:
                # Find the next line with genus info
                line_idx = lines.index(line)
                if line_idx + 1 < len(lines):
                    next_line = lines[line_idx + 1]
                    genus_match = re.search(r'([^:]+):\s+\d+\s+genomes\s+\(([0-9.]+)%\)', next_line)
                    if genus_match:
                        top_genus = genus_match.group(1).strip()
                        concentration = float(genus_match.group(2))
            
            # Concentration level
            if "🔥 HIGHLY CONCENTRATED" in line:
                concentration_level = "🔥 HIGHLY CONCENTRATED"
            elif "📊 MODERATELY CONCENTRATED" in line:
                concentration_level = "📊 MODERATELY CONCENTRATED"
            elif "🌳 DIVERSE" in line:
                concentration_level = "🌳 DIVERSE"
        
        groups[group_name] = {
            'taxid': taxid,
            'genomes': genomes,
            'genera': genera_count,
            'species': species_count,
            'top_genus': top_genus,
            'concentration': concentration,
            'concentration_level': concentration_level
        }
    
    return groups

def create_unified_ranking(filepath):
    """Create unified ranking from TaxID-based 16S results."""
    
    print(f"🔍 Parsing TaxID-based results from: {filepath}")
    
    # Parse the file
    groups = parse_taxonomic_groups_from_file(filepath)
    
    if not groups:
        print("❌ No taxonomic groups found in file")
        return
    
    print(f"📋 Found {len(groups)} taxonomic groups")
    
    # Sort by genome count
    sorted_groups = sorted(groups.items(), key=lambda x: x[1]['genomes'], reverse=True)
    
    # Calculate totals
    total_genomes = sum(data['genomes'] for _, data in sorted_groups)
    total_genera = sum(data['genera'] for _, data in sorted_groups)
    total_species = sum(data['species'] for _, data in sorted_groups)
    
    # Create unified ranking content
    content = []
    content.append("# UNIFIED 16S PROKARYOTIC OVERREPRESENTED TAXONOMIC GROUPS")
    content.append("# TaxID-Based Extraction Results")
    content.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    content.append("# " + "="*70)
    content.append("")
    content.append("🎯 UNIFIED RANKING OF ALL PROKARYOTIC OVERREPRESENTED TAXONOMIC GROUPS:")
    content.append("="*90)
    content.append(f"Total taxonomic groups: {len(groups)}")
    content.append(f"Total genomes: {total_genomes:,}")
    content.append(f"Total genera: {total_genera:,}")
    content.append(f"Total species: {total_species:,}")
    content.append(f"Average genomes per group: {total_genomes/len(groups):,.1f}")
    content.append("")
    content.append("📊 COMPLETE RANKING (Bacteria + Archaea combined, by genome count):")
    content.append("")
    
    for i, (group_name, data) in enumerate(sorted_groups, 1):
        content.append(f"   [{i:2d}] {group_name}: {data['genomes']:,} genomes (TaxID: {data['taxid']})")
        content.append(f"        Diversity: {data['genera']} genera, {data['species']} species")
        if data['top_genus']:
            content.append(f"        Top genus: {data['top_genus']} ({data['concentration']:.1f}% of group genomes)")
        content.append(f"        {data['concentration_level']}")
        content.append("")
    
    # Add concentration analysis
    content.append("🔍 CONCENTRATION ANALYSIS:")
    content.append("="*50)
    
    highly_concentrated = [name for name, data in sorted_groups if "🔥 HIGHLY CONCENTRATED" in data['concentration_level']]
    moderately_concentrated = [name for name, data in sorted_groups if "📊 MODERATELY CONCENTRATED" in data['concentration_level']]
    diverse = [name for name, data in sorted_groups if "🌳 DIVERSE" in data['concentration_level']]
    
    content.append(f"🔥 Highly concentrated groups ({len(highly_concentrated)}): {', '.join(highly_concentrated)}")
    content.append(f"📊 Moderately concentrated groups ({len(moderately_concentrated)}): {', '.join(moderately_concentrated)}")
    content.append(f"🌳 Diverse groups ({len(diverse)}): {', '.join(diverse)}")
    content.append("")
    
    # Top 5 analysis
    content.append("🏆 TOP 5 OVERREPRESENTED PROKARYOTIC TAXONOMIC GROUPS:")
    content.append("="*60)
    
    for i, (group_name, data) in enumerate(sorted_groups[:5], 1):
        pct_of_total = (data['genomes'] / total_genomes) * 100
        content.append(f"{i}. **{group_name}**: {data['genomes']:,} genomes ({pct_of_total:.1f}% of all overrepresented prokaryotic genomes)")
        content.append(f"   - Diversity: {data['genera']} genera across {data['species']} species")
        content.append(f"   - Dominated by: {data['top_genus']} ({data['concentration']:.1f}% of group)")
        content.append(f"   - Classification: {data['concentration_level']}")
        content.append("")
    
    # Write to file
    output_file = filepath.replace('.txt', '_unified_ranking.txt')
    with open(output_file, 'w') as f:
        f.write('\n'.join(content))
    
    print(f"✅ Created unified ranking: {output_file}")
    
    # Also append to the original file
    with open(filepath, 'a') as f:
        f.write("\n\n" + "="*90 + "\n")
        f.write("UNIFIED PROKARYOTIC RANKING (BACTERIA + ARCHAEA COMBINED)\n")
        f.write("="*90 + "\n\n")
        f.write('\n'.join(content[6:]))  # Skip header lines
    
    print(f"✅ Appended unified ranking to: {filepath}")
    
    return sorted_groups

def main():
    """Create unified ranking for 16S prokaryotic groups."""
    
    filepath = "family_diversity_results/16S_overrepresented_family_diversity.txt"
    
    if not os.path.exists(filepath):
        print(f"❌ File not found: {filepath}")
        return
    
    print("🚀 CREATING UNIFIED 16S PROKARYOTIC RANKING")
    print("="*60)
    
    results = create_unified_ranking(filepath)
    
    if results:
        print(f"\n🎯 SUCCESS: Created unified ranking for {len(results)} prokaryotic taxonomic groups")
        print(f"📊 Total genomes: {sum(data['genomes'] for _, data in results):,}")
        
        # Show top 3
        print(f"\n🏆 TOP 3 OVERREPRESENTED PROKARYOTIC GROUPS:")
        for i, (name, data) in enumerate(results[:3], 1):
            print(f"   {i}. {name}: {data['genomes']:,} genomes")

if __name__ == "__main__":
    main()
