#!/usr/bin/env python3
"""
Extract all taxa data from individual scatter plot CSV outputs
Created: 2025-08-18
Purpose: Export all taxa names, ranks, classifications, and factors for external research
"""

import pandas as pd
import os
import glob
from datetime import datetime

def find_csv_files():
    """Find all CSV output files from individual scatter plots"""
    csv_files = {}
    
    # Look in final_visualizations directory
    viz_dir = "final_visualizations"
    if os.path.exists(viz_dir):
        # 16S files
        csv_files['16s_bacteria_phylum'] = glob.glob(os.path.join(viz_dir, "16s_bacteria_phylum_top_*.csv"))
        csv_files['16s_bacteria_family'] = glob.glob(os.path.join(viz_dir, "16s_bacteria_family_top_*.csv"))
        csv_files['16s_bacteria_genus'] = glob.glob(os.path.join(viz_dir, "16s_bacteria_genus_top_*.csv"))
        csv_files['16s_archaea_phylum'] = glob.glob(os.path.join(viz_dir, "16s_archaea_phylum_top_*.csv"))
        csv_files['16s_archaea_family'] = glob.glob(os.path.join(viz_dir, "16s_archaea_family_top_*.csv"))
        
        # 18S files
        csv_files['18s_eukaryote_division'] = glob.glob(os.path.join(viz_dir, "18s_eukaryote_division_top_*.csv"))
        csv_files['18s_eukaryote_family'] = glob.glob(os.path.join(viz_dir, "18s_eukaryote_family_top_*.csv"))
        csv_files['18s_eukaryote_genus'] = glob.glob(os.path.join(viz_dir, "18s_eukaryote_genus_top_*.csv"))
        csv_files['18s_eukaryote_phylum'] = glob.glob(os.path.join(viz_dir, "18s_eukaryote_phylum_top_*.csv"))
    
    return csv_files

def load_and_process_csv(filepath, dataset_name):
    """Load and process a single CSV file"""
    try:
        df = pd.read_csv(filepath)
        
        # Add dataset information
        df['Dataset'] = dataset_name
        df['Source_File'] = os.path.basename(filepath)
        
        # Determine factor type from filename
        if 'novelty' in filepath.lower():
            df['Factor_Type'] = 'Novelty'
        elif 'coverage' in filepath.lower():
            df['Factor_Type'] = 'Coverage'
        else:
            df['Factor_Type'] = 'Unknown'
        
        # Determine domain and taxonomic level
        if '16s_bacteria' in dataset_name:
            df['Domain'] = 'Bacteria'
            df['rRNA'] = '16S'
        elif '16s_archaea' in dataset_name:
            df['Domain'] = 'Archaea'
            df['rRNA'] = '16S'
        elif '18s_eukaryote' in dataset_name:
            df['Domain'] = 'Eukaryota'
            df['rRNA'] = '18S'
        
        # Determine taxonomic level
        if 'phylum' in dataset_name or 'division' in dataset_name:
            df['Taxonomic_Level'] = 'Phylum/Division'
        elif 'family' in dataset_name:
            df['Taxonomic_Level'] = 'Family'
        elif 'genus' in dataset_name:
            df['Taxonomic_Level'] = 'Genus'
        
        return df
        
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def format_taxa_section(data, title):
    """Format a section of taxa data"""
    if data.empty:
        return f"\n{title}\n{'='*len(title)}\nNo data available\n"
    
    output = f"\n{title}\n{'='*len(title)}\n"
    output += f"Total entries: {len(data)}\n\n"
    
    # Column headers
    output += f"{'Taxon':<35} | {'Level':<15} | {'Domain':<12} | {'Nov.Ratio':<10} | {'Cov.Factor':<10} | {'Census':<8} | {'NCBI':<8}\n"
    output += f"{'-'*35} | {'-'*15} | {'-'*12} | {'-'*10} | {'-'*10} | {'-'*8} | {'-'*8}\n"
    
    for idx, row in data.iterrows():
        taxon = str(row.get('Taxon', row.iloc[0]))[:34]
        level = str(row.get('Taxonomic_Level', 'Unknown'))[:14]
        domain = str(row.get('Domain', 'Unknown'))[:11]
        
        # Handle different column names for ratios
        novelty = 0
        coverage = 0
        census = 0
        ncbi = 0
        
        # Try different possible column names
        for col in row.index:
            if 'novelty' in col.lower() and 'ratio' in col.lower():
                novelty = round(float(row[col]), 2) if pd.notna(row[col]) else 0
            elif 'coverage' in col.lower() and 'factor' in col.lower():
                coverage = round(float(row[col]), 2) if pd.notna(row[col]) else 0
            elif 'census' in col.lower() and 'count' in col.lower():
                census = int(row[col]) if pd.notna(row[col]) else 0
            elif 'ncbi' in col.lower() and ('species' in col.lower() or 'count' in col.lower()):
                ncbi = int(row[col]) if pd.notna(row[col]) else 0
        
        output += f"{taxon:<35} | {level:<15} | {domain:<12} | {novelty:<10.2f} | {coverage:<10.2f} | {census:<8} | {ncbi:<8}\n"
    
    return output

def main():
    print("Extracting all taxa data from scatter plot outputs...")
    
    # Find all CSV files
    csv_files = find_csv_files()
    
    all_data = []
    file_count = 0
    
    # Process each category
    for dataset_name, file_list in csv_files.items():
        for filepath in file_list:
            if os.path.exists(filepath):
                df = load_and_process_csv(filepath, dataset_name)
                if df is not None:
                    all_data.append(df)
                    file_count += 1
                    print(f"Loaded: {os.path.basename(filepath)} ({len(df)} entries)")
    
    if not all_data:
        print("No CSV files found! Make sure scatter plots have been generated.")
        return
    
    # Combine all data
    combined_data = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal entries loaded: {len(combined_data)} from {file_count} files")
    
    # Create output file
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_file = "all_taxa_research_data.txt"
    
    with open(output_file, 'w') as f:
        f.write(f"ALL TAXA RESEARCH DATA EXPORT\n")
        f.write(f"Generated: {timestamp}\n")
        f.write(f"{'='*100}\n")
        f.write(f"This file contains all taxa from individual scatter plot analyses\n")
        f.write(f"with their taxonomic ranks, classifications, and novelty/coverage factors.\n")
        f.write(f"Source: Individual scatter plot CSV output files\n")
        f.write(f"{'='*100}\n")
        
        # Organize by domain and level
        domains = ['Bacteria', 'Archaea', 'Eukaryota']
        levels = ['Phylum/Division', 'Family', 'Genus']
        
        for domain in domains:
            domain_data = combined_data[combined_data['Domain'] == domain]
            if not domain_data.empty:
                f.write(f"\n\n{domain.upper()} TAXA\n")
                f.write(f"{'='*80}\n")
                
                for level in levels:
                    level_data = domain_data[domain_data['Taxonomic_Level'] == level]
                    if not level_data.empty:
                        # Separate by factor type
                        novelty_data = level_data[level_data['Factor_Type'] == 'Novelty'].sort_values(
                            level_data.columns[level_data.columns.str.contains('novelty', case=False)].tolist()[0] 
                            if any(level_data.columns.str.contains('novelty', case=False)) else level_data.columns[1], 
                            ascending=False
                        )
                        coverage_data = level_data[level_data['Factor_Type'] == 'Coverage'].sort_values(
                            level_data.columns[level_data.columns.str.contains('coverage', case=False)].tolist()[0] 
                            if any(level_data.columns.str.contains('coverage', case=False)) else level_data.columns[2], 
                            ascending=False
                        )
                        
                        if not novelty_data.empty:
                            f.write(format_taxa_section(novelty_data, f"{domain} {level} - HIGH NOVELTY"))
                        
                        if not coverage_data.empty:
                            f.write(format_taxa_section(coverage_data, f"{domain} {level} - HIGH COVERAGE"))
        
        # Summary by domain
        f.write(f"\n\n\nSUMMARY BY DOMAIN\n")
        f.write(f"{'='*80}\n")
        
        for domain in domains:
            domain_data = combined_data[combined_data['Domain'] == domain]
            if not domain_data.empty:
                f.write(f"\n{domain}:\n")
                f.write(f"  Total unique taxa: {len(domain_data['Taxon'].unique()) if 'Taxon' in domain_data.columns else len(domain_data.iloc[:, 0].unique())}\n")
                f.write(f"  Phylum/Division level: {len(domain_data[domain_data['Taxonomic_Level'] == 'Phylum/Division'])}\n")
                f.write(f"  Family level: {len(domain_data[domain_data['Taxonomic_Level'] == 'Family'])}\n")
                f.write(f"  Genus level: {len(domain_data[domain_data['Taxonomic_Level'] == 'Genus'])}\n")
        
        f.write(f"\n\nCOLUMN DEFINITIONS:\n")
        f.write(f"- Taxon: Scientific name of the taxonomic group\n")
        f.write(f"- Level: Taxonomic level (Phylum/Division, Family, Genus)\n")
        f.write(f"- Domain: Bacteria, Archaea, or Eukaryota\n")
        f.write(f"- Nov.Ratio: Novelty Ratio (Census OTUs / NCBI Species)\n")
        f.write(f"- Cov.Factor: Coverage Factor (NCBI Species / Census OTUs)\n")
        f.write(f"- Census: Number of OTUs in environmental census data\n")
        f.write(f"- NCBI: Number of species in NCBI GenBank\n")
        f.write(f"\nValues >1.0 indicate higher representation in that database\n")
    
    print(f"✅ All taxa data exported to: {output_file}")
    print(f"   File contains {len(combined_data)} total entries")
    print(f"   Unique taxa: {len(combined_data.iloc[:, 0].unique())}")
    
    # Also create a simple CSV for easy import
    csv_output = "all_taxa_research_data.csv"
    combined_data.to_csv(csv_output, index=False)
    print(f"✅ CSV version saved to: {csv_output}")

if __name__ == "__main__":
    main()
