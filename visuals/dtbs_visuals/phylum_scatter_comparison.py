import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths relative to the script location
merged_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables"))
visualization_dir = script_dir  # Current directory

# Define input file
merged_phylum_file = os.path.join(merged_tables_dir, "merged_taxonomic_tables", "merged_phylum.csv")

# Check if input file exists
if not os.path.exists(merged_phylum_file):
    # Try absolute path as fallback
    merged_phylum_file = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/merged_taxonomic_tables/merged_phylum.csv"
    if not os.path.exists(merged_phylum_file):
        print(f"Error: Merged phylum file not found: {merged_phylum_file}")
        sys.exit(1)

print(f"Loading data from: {merged_phylum_file}")

# Load the data
df = pd.read_csv(merged_phylum_file)

# Process the data to get counts by database
# First, filter for prokaryotes only (Bacteria and Archaea)
df = df[df['domain'].isin(['Bacteria', 'Archaea'])]

# Create separate dataframes for GTDB and NCBI
gtdb_data = df[df['source'].str.contains('GTDB')].copy()
ncbi_data = df[df['source'].str.contains('NCBI')].copy()

# Group by phylum and sum genome counts
gtdb_counts = gtdb_data.groupby('taxon_name')['genome_count'].sum().reset_index()
ncbi_counts = ncbi_data.groupby('taxon_name')['genome_count'].sum().reset_index()

# Merge the counts
merged_counts = pd.merge(gtdb_counts, ncbi_counts, on='taxon_name', how='outer', suffixes=('_GTDB', '_NCBI'))
merged_counts = merged_counts.fillna(0)
merged_counts.columns = ['phylum', 'GTDB_Count', 'NCBI_Count']

# Use the merged counts for the scatter plot
df = merged_counts

# Create figure with controlled size
plt.figure(figsize=(10, 8))

# Calculate fold change for coloring
fold_change = np.log2((df["NCBI_Count"] + 1) / (df["GTDB_Count"] + 1))

# Create a more saturated colormap
scatter = plt.scatter(df["GTDB_Count"], df["NCBI_Count"],
                     alpha=1.0,
                     s=100,
                     c=fold_change,
                     cmap='RdBu_r',
                     vmin=-4,
                     vmax=4,
                     label='Phyla')

# Add labels for points that deviate significantly from the diagonal
for idx, row in df.iterrows():
    ratio = (row["NCBI_Count"] + 1) / (row["GTDB_Count"] + 1)
    # Adjust thresholds to catch extreme cases while avoiding medium values
    if ((ratio > 3 or ratio < 1/3) and max(row["NCBI_Count"], row["GTDB_Count"]) > 10000) or \
       (ratio > 10 or ratio < 1/10):  # Catch extreme ratios regardless of count
        label = f"{row['phylum']}\n({ratio:.1f}x)"
        plt.annotate(label,
                    (row["GTDB_Count"], row["NCBI_Count"]),
                    xytext=(10, 10),
                    textcoords='offset points',
                    fontsize=8,
                    bbox=dict(facecolor='white', edgecolor='black', alpha=0.7),
                    arrowprops=dict(arrowstyle='->', color='gray'))

# Add colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('Log2 Fold Change (NCBI/GTDB)\n← GTDB enriched | NCBI enriched →',
               rotation=270, labelpad=25)

# Log scale for both axes
plt.xscale("log")
plt.yscale("log")

# Diagonal line
max_val = max(df["GTDB_Count"].max(), df["NCBI_Count"].max())
min_val = min(df["GTDB_Count"].min(), df["NCBI_Count"].min())
plt.plot([min_val, max_val], [min_val, max_val],
         linestyle="--", color="gray", alpha=0.8, label="Equal representation")

# Labels and title
plt.xlabel("GTDB Genome Count (log scale)")
plt.ylabel("NCBI Genome Count (log scale)")
plt.title("GTDB vs. NCBI Genome Counts per Phylum", pad=20)

# Add grid with higher contrast
plt.grid(True, which="major", ls="-", alpha=0.3, color='gray')
plt.grid(True, which="minor", ls=":", alpha=0.2, color='gray')

# Add legend
plt.legend(loc='upper left')

# Adjust layout
plt.tight_layout()

# Save with reasonable DPI
plt.savefig("gtdb_vs_ncbi_scatter.png", dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
