#!/usr/bin/env python3
"""
Create Simple Unified 16S Ranking from TaxID Results
===================================================
Extract the ranking sections and combine them into one unified ranking.
"""

import re
import os
from datetime import datetime

def extract_rankings_from_file(filepath):
    """Extract ranking data from both bacteria and archaea sections."""
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find both ranking sections
    bacteria_pattern = r'📊 TAXONOMIC GROUP RANKING \(by genome count\):(.*?)(?=📋 Loaded \d+ overrepresented archaea|$)'
    archaea_pattern = r'📊 TAXONOMIC GROUP RANKING \(by genome count\):(.*?)(?=\n\n|$)'
    
    bacteria_match = re.search(bacteria_pattern, content, re.DOTALL)
    archaea_matches = list(re.finditer(archaea_pattern, content, re.DOTALL))
    
    all_groups = []
    
    # Extract bacteria groups
    if bacteria_match:
        bacteria_section = bacteria_match.group(1)
        bacteria_groups = extract_groups_from_section(bacteria_section, "Bacteria")
        all_groups.extend(bacteria_groups)
        print(f"📊 Found {len(bacteria_groups)} bacteria groups")
    
    # Extract archaea groups (take the last match which should be the archaea ranking)
    if len(archaea_matches) >= 2:
        archaea_section = archaea_matches[-1].group(1)
        archaea_groups = extract_groups_from_section(archaea_section, "Archaea")
        all_groups.extend(archaea_groups)
        print(f"📊 Found {len(archaea_groups)} archaea groups")
    
    return all_groups

def extract_groups_from_section(section, domain):
    """Extract individual group data from a ranking section."""
    
    groups = []
    
    # Split by lines and find group entries
    lines = section.strip().split('\n')
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Look for ranking entries like "[ 1] Bartonellaceae: 2,191 genomes"
        match = re.match(r'\[\s*\d+\]\s+([^:]+):\s+([\d,]+)\s+genomes', line)
        if match:
            group_name = match.group(1).strip()
            genomes_str = match.group(2).replace(',', '')
            genomes = int(genomes_str)
            
            # Extract additional info from next lines
            diversity_info = ""
            top_genus_info = ""
            concentration_level = ""
            
            # Look ahead for diversity and top genus info
            for j in range(i+1, min(i+4, len(lines))):
                if j < len(lines):
                    next_line = lines[j].strip()
                    
                    if "Diversity:" in next_line:
                        diversity_info = next_line.replace("Diversity:", "").strip()
                    elif "Top genus:" in next_line:
                        top_genus_info = next_line.replace("Top genus:", "").strip()
                    elif any(marker in next_line for marker in ["🔥 HIGHLY CONCENTRATED", "📊 MODERATELY CONCENTRATED", "🌳 DIVERSE"]):
                        concentration_level = next_line
            
            groups.append({
                'name': group_name,
                'genomes': genomes,
                'domain': domain,
                'diversity': diversity_info,
                'top_genus': top_genus_info,
                'concentration': concentration_level
            })
        
        i += 1
    
    return groups

def create_unified_ranking_simple(filepath):
    """Create a simple unified ranking."""
    
    print(f"🔍 Extracting rankings from: {filepath}")
    
    # Extract all groups
    all_groups = extract_rankings_from_file(filepath)
    
    if not all_groups:
        print("❌ No groups found")
        return
    
    print(f"📋 Total groups found: {len(all_groups)}")
    
    # Sort by genome count
    sorted_groups = sorted(all_groups, key=lambda x: x['genomes'], reverse=True)
    
    # Calculate totals
    total_genomes = sum(g['genomes'] for g in sorted_groups)
    bacteria_groups = [g for g in sorted_groups if g['domain'] == 'Bacteria']
    archaea_groups = [g for g in sorted_groups if g['domain'] == 'Archaea']
    
    bacteria_genomes = sum(g['genomes'] for g in bacteria_groups)
    archaea_genomes = sum(g['genomes'] for g in archaea_groups)
    
    # Create content
    content = []
    content.append("# UNIFIED 16S PROKARYOTIC OVERREPRESENTED TAXONOMIC GROUPS")
    content.append("# TaxID-Based Extraction Results - Complete Ranking")
    content.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    content.append("# " + "="*70)
    content.append("")
    content.append("🎯 COMPLETE UNIFIED RANKING (ALL PROKARYOTIC OVERREPRESENTED GROUPS):")
    content.append("="*80)
    content.append(f"📊 SUMMARY STATISTICS:")
    content.append(f"   Total taxonomic groups: {len(sorted_groups)}")
    content.append(f"   Total genomes: {total_genomes:,}")
    content.append(f"   Average genomes per group: {total_genomes/len(sorted_groups):,.1f}")
    content.append("")
    content.append(f"   Bacteria: {len(bacteria_groups)} groups, {bacteria_genomes:,} genomes ({bacteria_genomes/total_genomes*100:.1f}%)")
    content.append(f"   Archaea: {len(archaea_groups)} groups, {archaea_genomes:,} genomes ({archaea_genomes/total_genomes*100:.1f}%)")
    content.append("")
    content.append("🏆 UNIFIED RANKING (by genome count):")
    content.append("")
    
    for i, group in enumerate(sorted_groups, 1):
        pct_of_total = (group['genomes'] / total_genomes) * 100
        domain_icon = "🦠" if group['domain'] == "Bacteria" else "🔥"
        
        content.append(f"   [{i:2d}] {domain_icon} {group['name']}: {group['genomes']:,} genomes ({pct_of_total:.1f}%)")
        if group['diversity']:
            content.append(f"        {group['diversity']}")
        if group['top_genus']:
            content.append(f"        {group['top_genus']}")
        if group['concentration']:
            content.append(f"        {group['concentration']}")
        content.append("")
    
    # Top 5 analysis
    content.append("🎯 TOP 5 ANALYSIS:")
    content.append("="*40)
    
    top5_genomes = sum(g['genomes'] for g in sorted_groups[:5])
    top5_pct = (top5_genomes / total_genomes) * 100
    
    content.append(f"Top 5 groups represent {top5_genomes:,} genomes ({top5_pct:.1f}% of all overrepresented prokaryotic genomes)")
    content.append("")
    
    for i, group in enumerate(sorted_groups[:5], 1):
        pct = (group['genomes'] / total_genomes) * 100
        content.append(f"{i}. {group['name']} ({group['domain']}): {group['genomes']:,} genomes ({pct:.1f}%)")
    
    # Write to file
    output_file = filepath.replace('.txt', '_UNIFIED.txt')
    with open(output_file, 'w') as f:
        f.write('\n'.join(content))
    
    print(f"✅ Created unified ranking: {output_file}")
    
    return sorted_groups

def main():
    """Create simple unified ranking."""
    
    filepath = "family_diversity_results/16S_overrepresented_family_diversity.txt"
    
    print("🚀 CREATING SIMPLE UNIFIED 16S RANKING")
    print("="*50)
    
    results = create_unified_ranking_simple(filepath)
    
    if results:
        total_genomes = sum(g['genomes'] for g in results)
        print(f"\n🎯 SUCCESS!")
        print(f"📊 {len(results)} groups, {total_genomes:,} total genomes")
        
        print(f"\n🏆 TOP 5:")
        for i, group in enumerate(results[:5], 1):
            pct = (group['genomes'] / total_genomes) * 100
            print(f"   {i}. {group['name']}: {group['genomes']:,} genomes ({pct:.1f}%)")

if __name__ == "__main__":
    main()
