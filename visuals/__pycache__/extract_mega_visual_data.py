#!/usr/bin/env python3
"""
Extract and organize all source data from mega stacked visuals
Created: 2025-08-18
Purpose: Export all data used in mega visuals to organized .txt file for external research
"""

import pandas as pd
import os
from datetime import datetime

def load_16s_data():
    """Load 16S data from merged files"""
    data_dir = "../Eukcensus_merge/merged_output/16s_merged/results"
    
    datasets = {}
    levels = ["phylum", "family"]
    domains = ["bacteria", "archaea"]
    
    for domain in domains:
        datasets[domain] = {}
        for level in levels:
            filename = f"16s_ncbi_merged_clean_{domain}_{level}.csv"
            filepath = os.path.join(data_dir, filename)
            
            if os.path.exists(filepath):
                df = pd.read_csv(filepath)
                # Calculate ratios
                df['novelty_ratio'] = df['census_otu_count'] / df['ncbi_species_count']
                df['coverage_factor'] = df['ncbi_species_count'] / df['census_otu_count']
                df['weighted_novelty'] = df['novelty_ratio'] * (df['ncbi_species_count'] + 1).apply(lambda x: x**0.5)
                df['weighted_coverage'] = df['coverage_factor'] * (df['census_otu_count'] + 1).apply(lambda x: x**0.5)
                
                # Filter for meaningful data (factors > 1.0)
                meaningful = df[(df['novelty_ratio'] > 1.0) | (df['coverage_factor'] > 1.0)]
                
                datasets[domain][level] = {
                    'all_data': df,
                    'meaningful': meaningful,
                    'high_novelty': df[df['novelty_ratio'] > 1.0].sort_values('novelty_ratio', ascending=False),
                    'high_coverage': df[df['coverage_factor'] > 1.0].sort_values('coverage_factor', ascending=False)
                }
                print(f"Loaded {domain} {level}: {len(df)} total, {len(meaningful)} meaningful")
            else:
                print(f"Warning: {filepath} not found")
    
    return datasets

def load_18s_data():
    """Load 18S data from merged files"""
    data_dir = "../Eukcensus_merge/merged_output/18s_merged/results"
    
    datasets = {}
    levels = ["phylum", "family"]
    
    for level in levels:
        filename = f"18s_ncbi_merged_clean_{level}.csv"
        filepath = os.path.join(data_dir, filename)
        
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            # Calculate ratios
            df['novelty_ratio'] = df['census_otu_count'] / df['ncbi_species_count']
            df['coverage_factor'] = df['ncbi_species_count'] / df['census_otu_count']
            df['weighted_novelty'] = df['novelty_ratio'] * (df['ncbi_species_count'] + 1).apply(lambda x: x**0.5)
            df['weighted_coverage'] = df['coverage_factor'] * (df['census_otu_count'] / 1).apply(lambda x: x**0.5)
            
            # Filter for meaningful data (factors > 1.0)
            meaningful = df[(df['novelty_ratio'] > 1.0) | (df['coverage_factor'] > 1.0)]
            
            datasets[level] = {
                'all_data': df,
                'meaningful': meaningful,
                'high_novelty': df[df['novelty_ratio'] > 1.0].sort_values('novelty_ratio', ascending=False),
                'high_coverage': df[df['coverage_factor'] > 1.0].sort_values('coverage_factor', ascending=False)
            }
            print(f"Loaded 18S {level}: {len(df)} total, {len(meaningful)} meaningful")
        else:
            print(f"Warning: {filepath} not found")
    
    return datasets

def format_data_section(data, title, max_entries=50):
    """Format a data section for the output file"""
    if data.empty:
        return f"\n{title}\n{'='*len(title)}\nNo data available\n"
    
    output = f"\n{title}\n{'='*len(title)}\n"
    output += f"Total entries: {len(data)}\n\n"
    
    # Show top entries (limited to max_entries)
    display_data = data.head(max_entries)
    
    for idx, row in display_data.iterrows():
        taxon = row.iloc[0] if hasattr(row.iloc[0], 'strip') else str(row.iloc[0])
        census_count = int(row['census_otu_count']) if 'census_otu_count' in row else 0
        ncbi_count = int(row['ncbi_species_count']) if 'ncbi_species_count' in row else 0
        novelty = round(row['novelty_ratio'], 2) if 'novelty_ratio' in row else 0
        coverage = round(row['coverage_factor'], 2) if 'coverage_factor' in row else 0
        
        output += f"{taxon:<30} | Census: {census_count:>6} | NCBI: {ncbi_count:>6} | Nov: {novelty:>6.2f}x | Cov: {coverage:>6.2f}x\n"
    
    if len(data) > max_entries:
        output += f"\n... and {len(data) - max_entries} more entries\n"
    
    return output

def main():
    print("Extracting mega visual source data...")
    
    # Load all data
    data_16s = load_16s_data()
    data_18s = load_18s_data()
    
    # Create output file
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_file = "mega_visual_source_data.txt"
    
    with open(output_file, 'w') as f:
        f.write(f"MEGA VISUAL SOURCE DATA EXPORT\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"{'='*80}\n")
        f.write(f"This file contains all source data used in the mega stacked visualizations\n")
        f.write(f"organized by domain (16S: Bacteria/Archaea, 18S: Eukaryotes) and factors.\n")
        f.write(f"{'='*80}\n")
        
        # 16S Data Section
        f.write(f"\n\n16S rRNA DATA (BACTERIA & ARCHAEA)\n")
        f.write(f"{'='*80}\n")
        
        for domain in ['bacteria', 'archaea']:
            if domain in data_16s:
                f.write(f"\n\n{domain.upper()} DATA\n")
                f.write(f"{'-'*50}\n")
                
                for level in ['phylum', 'family']:
                    if level in data_16s[domain]:
                        data = data_16s[domain][level]
                        
                        f.write(format_data_section(
                            data['high_novelty'], 
                            f"{domain.title()} {level.title()} - HIGH NOVELTY (>1.0x)"
                        ))
                        
                        f.write(format_data_section(
                            data['high_coverage'], 
                            f"{domain.title()} {level.title()} - HIGH COVERAGE (>1.0x)"
                        ))
        
        # 18S Data Section
        f.write(f"\n\n\n18S rRNA DATA (EUKARYOTES)\n")
        f.write(f"{'='*80}\n")
        
        for level in ['phylum', 'family']:
            if level in data_18s:
                data = data_18s[level]
                
                level_name = "Divisions" if level == "phylum" else "Family"
                
                f.write(format_data_section(
                    data['high_novelty'], 
                    f"Eukaryote {level_name} - HIGH NOVELTY (>1.0x)"
                ))
                
                f.write(format_data_section(
                    data['high_coverage'], 
                    f"Eukaryote {level_name} - HIGH COVERAGE (>1.0x)"
                ))
        
        # Summary Statistics
        f.write(f"\n\n\nSUMMARY STATISTICS\n")
        f.write(f"{'='*80}\n")
        
        # 16S Summary
        total_16s_taxa = 0
        total_16s_meaningful = 0
        for domain in data_16s:
            for level in data_16s[domain]:
                total_16s_taxa += len(data_16s[domain][level]['all_data'])
                total_16s_meaningful += len(data_16s[domain][level]['meaningful'])
        
        f.write(f"16S rRNA (Bacteria + Archaea):\n")
        f.write(f"  Total taxa: {total_16s_taxa}\n")
        f.write(f"  Meaningful taxa (>1.0x factor): {total_16s_meaningful}\n")
        
        # 18S Summary
        total_18s_taxa = 0
        total_18s_meaningful = 0
        for level in data_18s:
            total_18s_taxa += len(data_18s[level]['all_data'])
            total_18s_meaningful += len(data_18s[level]['meaningful'])
        
        f.write(f"\n18S rRNA (Eukaryotes):\n")
        f.write(f"  Total taxa: {total_18s_taxa}\n")
        f.write(f"  Meaningful taxa (>1.0x factor): {total_18s_meaningful}\n")
        
        f.write(f"\nGRAND TOTAL:\n")
        f.write(f"  All taxa: {total_16s_taxa + total_18s_taxa}\n")
        f.write(f"  All meaningful: {total_16s_meaningful + total_18s_meaningful}\n")
        
        f.write(f"\n\nNOTES:\n")
        f.write(f"- Novelty Ratio: Census OTU Count / NCBI Species Count\n")
        f.write(f"- Coverage Factor: NCBI Species Count / Census OTU Count\n")
        f.write(f"- Values >1.0x indicate higher representation in that database\n")
        f.write(f"- Data sorted by factor value (highest first)\n")
        f.write(f"- Limited to top 50 entries per category for readability\n")
    
    print(f"✅ Data exported to: {output_file}")
    print(f"   Total file size: {os.path.getsize(output_file)} bytes")

if __name__ == "__main__":
    main()
