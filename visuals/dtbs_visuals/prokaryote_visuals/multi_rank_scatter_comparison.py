import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Set style
plt.style.use('default')

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths relative to the script location
merged_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables"))
visualization_dir = script_dir  # Current directory

# Define input files for prokaryotes (bacteria and archaea)
prokaryote_phylum_file = os.path.join(merged_tables_dir, "merged_taxonomic_tables", "prokaryote_phylum.csv")
prokaryote_family_file = os.path.join(merged_tables_dir, "merged_taxonomic_tables", "prokaryote_family.csv")
prokaryote_genus_file = os.path.join(merged_tables_dir, "merged_taxonomic_tables", "prokaryote_genus.csv")

# Check if input files exist and use fallback paths if needed
def check_file_exists(file_path, rank):
    if not os.path.exists(file_path):
        # Try absolute path as fallback
        fallback_path = f"/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/merged_taxonomic_tables/prokaryote_{rank}.csv"
        if not os.path.exists(fallback_path):
            print(f"Error: Prokaryote {rank} file not found: {file_path}")
            sys.exit(1)
        return fallback_path
    return file_path

# Check all files
prokaryote_phylum_file = check_file_exists(prokaryote_phylum_file, "phylum")
prokaryote_family_file = check_file_exists(prokaryote_family_file, "family")
prokaryote_genus_file = check_file_exists(prokaryote_genus_file, "genus")

# Function to process data for a specific taxonomic rank
def process_prokaryote_data(file_path, rank_name):
    print(f"Loading prokaryote {rank_name} data from: {file_path}")

    # Load the data
    df = pd.read_csv(file_path)

    # Filter for prokaryotes only (Bacteria and Archaea)
    df = df[df['domain'].isin(['Bacteria', 'Archaea'])]

    # Extract GTDB and NCBI counts directly from the file
    # The prokaryote files already have gtdb_genome_count and ncbi_genome_count columns
    data = df[['taxon_name', 'gtdb_genome_count', 'ncbi_genome_count']].copy()

    # Rename columns to match expected format
    data.columns = [rank_name, 'GTDB_Count', 'NCBI_Count']

    # Fill NaN values with 0
    data = data.fillna(0)

    return data

# Process data for each taxonomic rank
phylum_data = process_prokaryote_data(prokaryote_phylum_file, "phylum")
family_data = process_prokaryote_data(prokaryote_family_file, "family")
genus_data = process_prokaryote_data(prokaryote_genus_file, "genus")

# Create a figure with 3 subplots (one for each taxonomic rank)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))

# Function to create a scatter plot for a specific taxonomic rank
def create_scatter_plot(ax, data, rank_name):
    # Calculate fold change for coloring
    fold_change = np.log2((data["NCBI_Count"] + 1) / (data["GTDB_Count"] + 1))

    # Create scatter plot
    scatter = ax.scatter(data["GTDB_Count"], data["NCBI_Count"],
                         alpha=0.8,
                         s=80,
                         c=fold_change,
                         cmap='RdBu_r',
                         vmin=-4,
                         vmax=4,
                         label=f'{rank_name}s')

    # Add labels for points that deviate significantly from the diagonal
    for idx, row in data.iterrows():
        ratio = (row["NCBI_Count"] + 1) / (row["GTDB_Count"] + 1)
        # Adjust thresholds to catch extreme cases while avoiding medium values
        if ((ratio > 3 or ratio < 1/3) and max(row["NCBI_Count"], row["GTDB_Count"]) > 10000) or \
           (ratio > 10 or ratio < 1/10):  # Catch extreme ratios regardless of count
            label = f"{row[rank_name]}\n({ratio:.1f}x)"
            ax.annotate(label,
                        (row["GTDB_Count"], row["NCBI_Count"]),
                        xytext=(10, 10),
                        textcoords='offset points',
                        fontsize=8,
                        bbox=dict(facecolor='white', edgecolor='black', alpha=0.7),
                        arrowprops=dict(arrowstyle='->', color='gray'))

    # Log scale for both axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Diagonal line
    max_val = max(data["GTDB_Count"].max(), data["NCBI_Count"].max())
    min_val = min(data["GTDB_Count"].min(), data["NCBI_Count"].min())
    min_val = max(0.1, min_val)  # Ensure min_val is positive for log scale
    ax.plot([min_val, max_val], [min_val, max_val],
             linestyle="--", color="gray", alpha=0.8, label="Equal representation")

    # Labels and title
    ax.set_xlabel("GTDB Genome Count (log scale)")
    ax.set_ylabel("NCBI Genome Count (log scale)")
    ax.set_title(f"GTDB vs. NCBI Genome Counts per {rank_name.capitalize()}", pad=20)

    # Add grid with higher contrast
    ax.grid(True, which="major", ls="-", alpha=0.3, color='gray')
    ax.grid(True, which="minor", ls=":", alpha=0.2, color='gray')

    # Add legend
    ax.legend(loc='upper left')

    # Calculate statistics
    total_taxa = len(data)
    gtdb_only = len(data[(data['GTDB_Count'] > 0) & (data['NCBI_Count'] == 0)])
    ncbi_only = len(data[(data['NCBI_Count'] > 0) & (data['GTDB_Count'] == 0)])
    both = len(data[(data['GTDB_Count'] > 0) & (data['NCBI_Count'] > 0)])

    # Add text with statistics
    ax.text(0.05, 0.95,
            f"Total {rank_name}s: {total_taxa}\n"
            f"GTDB only: {gtdb_only} ({gtdb_only/total_taxa*100:.1f}%)\n"
            f"NCBI only: {ncbi_only} ({ncbi_only/total_taxa*100:.1f}%)\n"
            f"Both: {both} ({both/total_taxa*100:.1f}%)",
            transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
            verticalalignment='top')

    return scatter

# Create scatter plots for each taxonomic rank
scatter1 = create_scatter_plot(ax1, phylum_data, "phylum")
scatter2 = create_scatter_plot(ax2, family_data, "family")
scatter3 = create_scatter_plot(ax3, genus_data, "genus")

# Add a single colorbar for all subplots
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(scatter1, cax=cbar_ax)
cbar.set_label('Log2 Fold Change (NCBI/GTDB)\n← GTDB enriched | NCBI enriched →',
               rotation=270, labelpad=25)

# Add overall title
plt.suptitle('Multi-Rank Comparison: GTDB vs NCBI Genome Counts', fontsize=16, y=0.98)

# Add explanation text
fig.text(0.5, 0.01,
        "Scatter plots show the comparison of genome counts between GTDB and NCBI databases across taxonomic ranks.\n"
        "Points above the diagonal line have more genomes in NCBI, while points below have more genomes in GTDB.\n"
        "Color indicates the log2 fold change between databases (blue = GTDB enriched, red = NCBI enriched).",
        ha='center', fontsize=12,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.subplots_adjust(left=0.05, right=0.9, bottom=0.1, top=0.9, wspace=0.3)

# Save figure
output_path = os.path.join(visualization_dir, "multi_rank_scatter_comparison.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nMulti-Rank Comparison Summary:")
print(f"Phylum level: {len(phylum_data)} total phyla")
print(f"  - GTDB only: {len(phylum_data[(phylum_data['GTDB_Count'] > 0) & (phylum_data['NCBI_Count'] == 0)])}")
print(f"  - NCBI only: {len(phylum_data[(phylum_data['NCBI_Count'] > 0) & (phylum_data['GTDB_Count'] == 0)])}")
print(f"  - Both: {len(phylum_data[(phylum_data['GTDB_Count'] > 0) & (phylum_data['NCBI_Count'] > 0)])}")

print(f"\nFamily level: {len(family_data)} total families")
print(f"  - GTDB only: {len(family_data[(family_data['GTDB_Count'] > 0) & (family_data['NCBI_Count'] == 0)])}")
print(f"  - NCBI only: {len(family_data[(family_data['NCBI_Count'] > 0) & (family_data['GTDB_Count'] == 0)])}")
print(f"  - Both: {len(family_data[(family_data['GTDB_Count'] > 0) & (family_data['NCBI_Count'] > 0)])}")

print(f"\nGenus level: {len(genus_data)} total genera")
print(f"  - GTDB only: {len(genus_data[(genus_data['GTDB_Count'] > 0) & (genus_data['NCBI_Count'] == 0)])}")
print(f"  - NCBI only: {len(genus_data[(genus_data['NCBI_Count'] > 0) & (genus_data['GTDB_Count'] == 0)])}")
print(f"  - Both: {len(genus_data[(genus_data['GTDB_Count'] > 0) & (genus_data['NCBI_Count'] > 0)])}")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
