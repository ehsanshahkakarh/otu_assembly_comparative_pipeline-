import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import matplotlib.colors as mcolors

# Set style
plt.style.use('default')

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths relative to the script location
csv_output_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "csv_output"))
domain_mismatches_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "domain_taxonomic_tables", "mismatches"))
visualization_dir = script_dir  # Current directory

# Define input files for mismatches
family_mismatches_file = os.path.join(csv_output_dir, "family_mismatches.csv")
bacteria_family_mismatches_file = os.path.join(domain_mismatches_dir, "bacteria_family_mismatches.csv")
archaea_family_mismatches_file = os.path.join(domain_mismatches_dir, "archaea_family_mismatches.csv")

# Check if input files exist and use fallback paths if needed
def check_file_exists(file_path, description):
    if not os.path.exists(file_path):
        # Try absolute path as fallback
        if "bacteria" in file_path:
            fallback_path = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/domain_taxonomic_tables/mismatches/bacteria_family_mismatches.csv"
        elif "archaea" in file_path:
            fallback_path = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/domain_taxonomic_tables/mismatches/archaea_family_mismatches.csv"
        else:
            fallback_path = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/csv_output/family_mismatches.csv"
        
        if not os.path.exists(fallback_path):
            print(f"Error: {description} file not found: {file_path}")
            print(f"Fallback path also not found: {fallback_path}")
            return None
        return fallback_path
    return file_path

# Check all files
family_mismatches_file = check_file_exists(family_mismatches_file, "Family mismatches")
bacteria_family_mismatches_file = check_file_exists(bacteria_family_mismatches_file, "Bacteria family mismatches")
archaea_family_mismatches_file = check_file_exists(archaea_family_mismatches_file, "Archaea family mismatches")

# Create a figure with 2 subplots (one for each domain)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16))

# Function to analyze and visualize family mismatches
def analyze_family_mismatches(file_path, ax, domain_name):
    if file_path is None:
        ax.text(0.5, 0.5, f"No {domain_name} family mismatch data available", 
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        ax.axis('off')
        return None
    
    print(f"Loading {domain_name} family mismatches from: {file_path}")
    try:
        mismatches_data = pd.read_csv(file_path)
        
        # Check if the file is empty or has no data
        if mismatches_data.empty:
            ax.text(0.5, 0.5, f"No {domain_name} family mismatch data available (empty file)", 
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.axis('off')
            return None
            
        # Print the first few rows to understand the structure
        print(f"\n{domain_name} family mismatches - first few rows:")
        print(mismatches_data.head())
        
        # Print column names
        print(f"\n{domain_name} family mismatches - columns:")
        print(mismatches_data.columns.tolist())
        
        # Check if the necessary columns exist
        required_cols = ['family_gtdb', 'family_ncbi', 'gtdb_genome_count', 'ncbi_genome_count']
        missing_cols = [col for col in required_cols if col not in mismatches_data.columns]
        
        if missing_cols:
            print(f"Warning: Missing required columns in {domain_name} family mismatches: {missing_cols}")
            ax.text(0.5, 0.5, f"Cannot visualize {domain_name} family mismatches\nMissing columns: {', '.join(missing_cols)}", 
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.axis('off')
            return None
        
        # Get top families by genome count
        top_gtdb_families = mismatches_data.groupby('family_gtdb')['gtdb_genome_count'].sum().nlargest(10).index.tolist()
        top_ncbi_families = mismatches_data.groupby('family_ncbi')['ncbi_genome_count'].sum().nlargest(10).index.tolist()
        
        # Create a matrix of mismatches between top families
        mismatch_matrix = pd.DataFrame(0, index=top_gtdb_families, columns=top_ncbi_families)
        
        for _, row in mismatches_data.iterrows():
            if row['family_gtdb'] in top_gtdb_families and row['family_ncbi'] in top_ncbi_families:
                mismatch_matrix.loc[row['family_gtdb'], row['family_ncbi']] += row['gtdb_genome_count']
        
        # Normalize by row sums for better visualization
        row_sums = mismatch_matrix.sum(axis=1)
        normalized_matrix = mismatch_matrix.div(row_sums, axis=0).fillna(0)
        
        # Create heatmap
        im = ax.imshow(normalized_matrix, cmap='YlOrRd')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Proportion of GTDB genomes', rotation=270, labelpad=15)
        
        # Add labels
        ax.set_xticks(np.arange(len(normalized_matrix.columns)))
        ax.set_yticks(np.arange(len(normalized_matrix.index)))
        ax.set_xticklabels(normalized_matrix.columns, rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(normalized_matrix.index, fontsize=9)
        
        # Add title
        ax.set_title(f"{domain_name} Family Reclassification: GTDB vs NCBI", fontsize=14)
        
        # Add axis labels
        ax.set_xlabel('NCBI Classification', fontsize=12)
        ax.set_ylabel('GTDB Classification', fontsize=12)
        
        # Add grid lines
        ax.set_xticks(np.arange(-.5, len(normalized_matrix.columns), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(normalized_matrix.index), 1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=1)
        
        # Calculate mismatch statistics
        total_genomes = mismatches_data['gtdb_genome_count'].sum()
        mismatched_genomes = mismatches_data[mismatches_data['family_gtdb'] != mismatches_data['family_ncbi']]['gtdb_genome_count'].sum()
        mismatch_percent = (mismatched_genomes / total_genomes * 100) if total_genomes > 0 else 0
        
        # Add text with statistics
        ax.text(0.5, -0.15, 
                f"Total genomes: {total_genomes:,}\n"
                f"Mismatched genomes: {mismatched_genomes:,} ({mismatch_percent:.1f}%)",
                ha='center', va='center', transform=ax.transAxes, fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        return {
            'total_genomes': total_genomes,
            'mismatched_genomes': mismatched_genomes,
            'mismatch_percent': mismatch_percent,
            'top_gtdb_families': top_gtdb_families,
            'top_ncbi_families': top_ncbi_families
        }
    
    except Exception as e:
        print(f"Error processing {domain_name} family mismatches: {str(e)}")
        ax.text(0.5, 0.5, f"Error processing {domain_name} family mismatches:\n{str(e)}", 
                ha='center', va='center', fontsize=12, transform=ax.transAxes)
        ax.axis('off')
        return None

# Analyze and visualize family mismatches for each domain
bacteria_stats = analyze_family_mismatches(bacteria_family_mismatches_file, ax1, "Bacteria")
archaea_stats = analyze_family_mismatches(archaea_family_mismatches_file, ax2, "Archaea")

# Add overall title
plt.suptitle('Prokaryotic Family Reclassification: GTDB vs NCBI', fontsize=16, y=0.92)

# Add explanation text
fig.text(0.5, 0.01, 
        "Heatmaps show the proportion of genomes classified differently between GTDB and NCBI at the family level.\n"
        "Each row represents a GTDB family, and each column represents an NCBI family.\n"
        "Darker colors indicate a higher proportion of genomes with that particular reclassification.",
        ha='center', fontsize=12,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.tight_layout(rect=[0, 0.03, 1, 0.90])

# Save figure
output_path = os.path.join(visualization_dir, "prokaryote_family_mismatches.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nProkaryotic Family Reclassification Summary:")
if bacteria_stats:
    print(f"Bacteria: {bacteria_stats['mismatched_genomes']:,} of {bacteria_stats['total_genomes']:,} genomes ({bacteria_stats['mismatch_percent']:.1f}%) have different family classifications")
if archaea_stats:
    print(f"Archaea: {archaea_stats['mismatched_genomes']:,} of {archaea_stats['total_genomes']:,} genomes ({archaea_stats['mismatch_percent']:.1f}%) have different family classifications")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
