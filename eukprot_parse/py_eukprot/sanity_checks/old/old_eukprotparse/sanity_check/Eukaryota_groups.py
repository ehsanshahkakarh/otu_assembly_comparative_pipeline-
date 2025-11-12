#!/usr/bin/env python3
"""
Analyze taxonomic groups in the EukProt database.
This script analyzes the distribution of taxonomic groups in the EukProt database
and generates statistics on taxonomic coverage.
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
import logging
from datetime import datetime

# Set up logging
log_dir = Path(__file__).parent / "logs"
log_dir.mkdir(exist_ok=True)
log_file = log_dir / f"eukaryota_groups_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

def load_eukprot_data(csv_dir):
    """
    Load EukProt taxonomic data from CSV files.
    
    Args:
        csv_dir (Path): Directory containing EukProt CSV files
        
    Returns:
        dict: Dictionary containing DataFrames for phyla, families, and genera
    """
    data = {}
    
    # Define files to load
    files = {
        "phylum": csv_dir / "eukprot_phylum_counts.csv",
        "family": csv_dir / "eukprot_family_counts.csv",
        "genus": csv_dir / "eukprot_genus_counts.csv"
    }
    
    # Load each file
    for rank, file_path in files.items():
        if file_path.exists():
            try:
                df = pd.read_csv(file_path)
                data[rank] = df
                logging.info(f"Loaded {len(df)} {rank} records from {file_path.name}")
            except Exception as e:
                logging.error(f"Error loading {file_path}: {e}")
        else:
            logging.warning(f"File not found: {file_path}")
    
    return data

def analyze_taxonomic_groups(data):
    """
    Analyze taxonomic groups in the EukProt database.
    
    Args:
        data (dict): Dictionary containing DataFrames for phyla, families, and genera
        
    Returns:
        dict: Dictionary containing analysis results
    """
    results = {}
    
    for rank, df in data.items():
        if df is None or df.empty:
            logging.warning(f"No data available for {rank}")
            continue
        
        # Get total number of genomes
        if "eukprot_genome_count" in df.columns:
            total_genomes = df["eukprot_genome_count"].sum()
        else:
            total_genomes = 0
            logging.warning(f"Column 'eukprot_genome_count' not found in {rank} data")
        
        # Get total number of taxa
        total_taxa = len(df)
        
        # Get top 10 taxa by genome count
        if "eukprot_genome_count" in df.columns:
            top_taxa = df.sort_values("eukprot_genome_count", ascending=False).head(10)
        else:
            top_taxa = df.head(10)
            logging.warning(f"Column 'eukprot_genome_count' not found in {rank} data")
        
        # Calculate coverage of top 10 taxa
        if "eukprot_genome_count" in df.columns:
            top_10_coverage = top_taxa["eukprot_genome_count"].sum() / total_genomes * 100 if total_genomes > 0 else 0
        else:
            top_10_coverage = 0
        
        # Store results
        results[rank] = {
            "total_genomes": total_genomes,
            "total_taxa": total_taxa,
            "top_taxa": top_taxa,
            "top_10_coverage": top_10_coverage
        }
        
        logging.info(f"Analysis for {rank}:")
        logging.info(f"  Total genomes: {total_genomes}")
        logging.info(f"  Total taxa: {total_taxa}")
        logging.info(f"  Top 10 taxa coverage: {top_10_coverage:.2f}%")
    
    return results

def plot_taxonomic_distribution(results, output_dir):
    """
    Plot taxonomic distribution for each rank.
    
    Args:
        results (dict): Dictionary containing analysis results
        output_dir (Path): Directory to save plots
    """
    output_dir.mkdir(exist_ok=True)
    
    for rank, result in results.items():
        if "top_taxa" not in result or result["top_taxa"] is None or result["top_taxa"].empty:
            logging.warning(f"No data available for plotting {rank} distribution")
            continue
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Get data for plotting
        if "eukprot_genome_count" in result["top_taxa"].columns and rank in result["top_taxa"].columns:
            taxa = result["top_taxa"][rank]
            counts = result["top_taxa"]["eukprot_genome_count"]
        else:
            logging.warning(f"Required columns not found in {rank} data for plotting")
            continue
        
        # Create bar plot
        plt.bar(range(len(taxa)), counts, color='skyblue')
        plt.xticks(range(len(taxa)), taxa, rotation=45, ha='right')
        plt.xlabel(f'{rank.capitalize()} Name')
        plt.ylabel('Number of Genomes')
        plt.title(f'Top 10 {rank.capitalize()} in EukProt Database')
        plt.tight_layout()
        
        # Save plot
        output_file = output_dir / f"eukprot_{rank}_distribution.png"
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logging.info(f"Saved plot to {output_file}")

def generate_report(results, output_dir):
    """
    Generate a report of the analysis results.
    
    Args:
        results (dict): Dictionary containing analysis results
        output_dir (Path): Directory to save report
    """
    output_dir.mkdir(exist_ok=True)
    report_file = output_dir / "eukprot_taxonomic_analysis.md"
    
    with open(report_file, 'w') as f:
        f.write("# EukProt Taxonomic Analysis\n\n")
        f.write(f"Report generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for rank, result in results.items():
            f.write(f"## {rank.capitalize()} Analysis\n\n")
            f.write(f"- Total genomes: {result['total_genomes']}\n")
            f.write(f"- Total {rank} taxa: {result['total_taxa']}\n")
            f.write(f"- Top 10 {rank} taxa coverage: {result['top_10_coverage']:.2f}%\n\n")
            
            f.write(f"### Top 10 {rank.capitalize()} Taxa\n\n")
            
            if "top_taxa" in result and result["top_taxa"] is not None and not result["top_taxa"].empty:
                if "eukprot_genome_count" in result["top_taxa"].columns and rank in result["top_taxa"].columns:
                    f.write("| Rank | Taxon Name | Genome Count | Percentage |\n")
                    f.write("|------|-----------|--------------|------------|\n")
                    
                    for i, (_, row) in enumerate(result["top_taxa"].iterrows(), 1):
                        taxon = row[rank]
                        count = row["eukprot_genome_count"]
                        percentage = count / result["total_genomes"] * 100 if result["total_genomes"] > 0 else 0
                        f.write(f"| {i} | {taxon} | {count} | {percentage:.2f}% |\n")
                else:
                    f.write("Required columns not found in data for reporting\n")
            else:
                f.write("No data available for reporting\n")
            
            f.write("\n")
        
        f.write("## Summary\n\n")
        f.write("This analysis provides an overview of the taxonomic distribution in the EukProt database.\n")
        f.write("The results show the diversity of eukaryotic taxa represented in the database and the coverage of the top taxa.\n")
    
    logging.info(f"Generated report: {report_file}")
    return report_file

def main():
    """Main function to analyze EukProt taxonomic groups."""
    # Set up paths
    script_dir = Path(__file__).parent.resolve()
    eukprot_dir = script_dir.parent
    output_dir = script_dir / "output"
    output_dir.mkdir(exist_ok=True)
    
    logging.info("Starting EukProt taxonomic group analysis")
    
    # Load EukProt data
    data = load_eukprot_data(eukprot_dir)
    
    # Analyze taxonomic groups
    results = analyze_taxonomic_groups(data)
    
    # Plot taxonomic distribution
    plot_taxonomic_distribution(results, output_dir)
    
    # Generate report
    report_file = generate_report(results, output_dir)
    
    logging.info(f"Analysis complete. Report saved to {report_file}")

if __name__ == "__main__":
    main()
