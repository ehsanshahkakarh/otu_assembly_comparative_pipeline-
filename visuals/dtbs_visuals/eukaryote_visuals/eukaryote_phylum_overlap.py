import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# Set style
plt.style.use('default')

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths relative to the script location
eukaryota_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "eukaryota_taxonomic_tables"))
visualization_dir = script_dir  # Current directory

# Define input file
phylum_file = os.path.join(eukaryota_tables_dir, "eukaryota_phylum.csv")

# Check if input file exists
if not os.path.exists(phylum_file):
    # Try absolute path as fallback
    phylum_file = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/eukaryota_taxonomic_tables/eukaryota_phylum.csv"
    if not os.path.exists(phylum_file):
        print(f"Error: File not found: {phylum_file}")
        sys.exit(1)

# Load data
print(f"Loading data from: {phylum_file}")
phylum_data = pd.read_csv(phylum_file)

# Create sets of taxa in each database
ncbi_phyla = set(phylum_data[phylum_data['ncbi_genome_count'] > 0]['taxon_name'])
eukprot_phyla = set(phylum_data[phylum_data['eukprot_genome_count'] > 0]['taxon_name'])

# Calculate sizes for overlap analysis
only_ncbi = len(ncbi_phyla - eukprot_phyla)
only_eukprot = len(eukprot_phyla - ncbi_phyla)
both = len(ncbi_phyla.intersection(eukprot_phyla))

# Calculate totals
total_phyla = len(phylum_data['taxon_name'].unique())
total_ncbi = len(ncbi_phyla)
total_eukprot = len(eukprot_phyla)

# Calculate percentages
ncbi_percent = (total_ncbi / total_phyla * 100) if total_phyla > 0 else 0
eukprot_percent = (total_eukprot / total_phyla * 100) if total_phyla > 0 else 0

# Create a dictionary for the bar plot
counts = {
    "NCBI only": only_ncbi,
    "Both": both,
    "EukProt only": only_eukprot
}

# Define colors for the bars
colors = {
    "NCBI only": "#1f77b4",       # muted blue
    "Both": "#9467bd",            # purple
    "EukProt only": "#2ca02c"     # green
}

# Create the plot
plt.figure(figsize=(9, 6))
bars = plt.bar(counts.keys(), counts.values(), color=[colors[k] for k in counts])

# Add titles and labels
plt.title("Phylum-Level Overlap: NCBI vs EukProt", fontsize=16)
plt.ylabel("Number of Unique Phyla", fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Annotate bars with values
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height + 1, f"{height}", 
             ha='center', va='bottom', fontsize=12)

# Add summary box
plt.figtext(0.5, 0.01,
        f'Total unique phyla (combined): {total_phyla} | '
        f'NCBI: {total_ncbi} ({ncbi_percent:.1f}%) | '
        f'EukProt: {total_eukprot} ({eukprot_percent:.1f}%)',
        ha='center', fontsize=11,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='gray', alpha=0.9))

# Save and show the plot
output_path = os.path.join(visualization_dir, "eukaryote_phylum_overlap_barplot.png")
plt.tight_layout()
plt.savefig(output_path, dpi=300)
plt.show()

# Print summary statistics
print("\nEukaryote Phylum Coverage Summary:")
print(f"Total unique phyla: {total_phyla}")
print(f"NCBI-only phyla: {only_ncbi}")
print(f"EukProt-only phyla: {only_eukprot}")
print(f"Phyla in both databases: {both}")
print(f"NCBI total phyla: {total_ncbi} ({ncbi_percent:.1f}%)")
print(f"EukProt total phyla: {total_eukprot} ({eukprot_percent:.1f}%)")
print(f"Saved to: {output_path}")
