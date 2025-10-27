import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

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

# Create figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Function to create a pie chart for a specific database
def create_pie_chart(ax, data, database, title):
    # Filter data for the specific database
    if database == 'NCBI':
        db_data = data[data['ncbi_genome_count'] > 0].copy()
        db_data['count'] = db_data['ncbi_genome_count']
    else:  # EukProt
        db_data = data[data['eukprot_genome_count'] > 0].copy()
        db_data['count'] = db_data['eukprot_genome_count']
    
    # Sort by count
    db_data = db_data.sort_values('count', ascending=False)
    
    # Get top 10 phyla and group the rest as "Other"
    top_phyla = db_data.head(10).copy()
    other_phyla = db_data.iloc[10:].copy() if len(db_data) > 10 else None
    
    # Prepare data for pie chart
    labels = list(top_phyla['taxon_name'])
    sizes = list(top_phyla['count'])
    
    # Add "Other" category if needed
    if other_phyla is not None and len(other_phyla) > 0:
        labels.append('Other')
        sizes.append(other_phyla['count'].sum())
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie(
        sizes, 
        labels=None,  # We'll add a legend instead
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops=dict(width=0.5, edgecolor='w'),
        textprops={'fontsize': 12}
    )
    
    # Make the percentage labels more readable
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    # Add a legend
    legend_labels = [f"{label} ({size})" for label, size in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    
    # Add title
    ax.set_title(f"{title} (n={sum(sizes)})", fontsize=16, pad=20)
    
    return db_data

# Create pie charts for each database
ncbi_data = create_pie_chart(ax1, phylum_data, 'NCBI', 'NCBI Eukaryote Phyla Distribution')
eukprot_data = create_pie_chart(ax2, phylum_data, 'EukProt', 'EukProt Eukaryote Phyla Distribution')

# Add overall title
plt.suptitle('Eukaryote Phyla Distribution Comparison', fontsize=18, y=0.98)

# Add summary statistics
ncbi_total = ncbi_data['count'].sum()
eukprot_total = eukprot_data['count'].sum()
ncbi_unique = len(ncbi_data)
eukprot_unique = len(eukprot_data)

summary_text = (
    f"NCBI: {ncbi_total} genomes across {ncbi_unique} phyla\n"
    f"EukProt: {eukprot_total} genomes across {eukprot_unique} phyla\n"
    f"Shared phyla: {len(set(ncbi_data['taxon_name']).intersection(set(eukprot_data['taxon_name'])))}"
)

fig.text(0.5, 0.01, summary_text, ha='center', fontsize=12, 
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Save figure
output_path = os.path.join(visualization_dir, "eukaryote_phylum_distribution.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nEukaryote Phylum Distribution Summary:")
print(f"NCBI: {ncbi_total} genomes across {ncbi_unique} phyla")
print(f"EukProt: {eukprot_total} genomes across {eukprot_unique} phyla")

# Print top 5 phyla for each database
print("\nTop 5 phyla in NCBI:")
for i, (_, row) in enumerate(ncbi_data.head(5).iterrows(), 1):
    percent = (row['count'] / ncbi_total * 100)
    print(f"{i}. {row['taxon_name']}: {int(row['count'])} genomes ({percent:.1f}%)")

print("\nTop 5 phyla in EukProt:")
for i, (_, row) in enumerate(eukprot_data.head(5).iterrows(), 1):
    percent = (row['count'] / eukprot_total * 100)
    print(f"{i}. {row['taxon_name']}: {int(row['count'])} genomes ({percent:.1f}%)")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
