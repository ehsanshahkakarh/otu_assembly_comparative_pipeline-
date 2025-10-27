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

# Define input files
phylum_file = os.path.join(eukaryota_tables_dir, "eukaryota_phylum.csv")
family_file = os.path.join(eukaryota_tables_dir, "eukaryota_family.csv")
genus_file = os.path.join(eukaryota_tables_dir, "eukaryota_genus.csv")

# Check if input files exist
for file_path in [phylum_file, family_file, genus_file]:
    if not os.path.exists(file_path):
        # Try absolute path as fallback
        base_name = os.path.basename(file_path)
        fallback_path = f"/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/eukaryota_taxonomic_tables/{base_name}"
        if not os.path.exists(fallback_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
        else:
            if file_path == phylum_file:
                phylum_file = fallback_path
            elif file_path == family_file:
                family_file = fallback_path
            else:
                genus_file = fallback_path

# Load data
print(f"Loading phylum data from: {phylum_file}")
phylum_data = pd.read_csv(phylum_file)

print(f"Loading family data from: {family_file}")
family_data = pd.read_csv(family_file)

print(f"Loading genus data from: {genus_file}")
genus_data = pd.read_csv(genus_file)

# Calculate diversity metrics for each taxonomic level
def calculate_diversity(data):
    # Count unique taxa in each database
    ncbi_taxa = len(data[data['ncbi_genome_count'] > 0])
    eukprot_taxa = len(data[data['eukprot_genome_count'] > 0])
    shared_taxa = len(set(data[data['ncbi_genome_count'] > 0]['taxon_name']).intersection(
                      set(data[data['eukprot_genome_count'] > 0]['taxon_name'])))
    
    # Calculate total genomes in each database
    ncbi_genomes = data['ncbi_genome_count'].sum()
    eukprot_genomes = data['eukprot_genome_count'].sum()
    
    # Calculate average genomes per taxon
    ncbi_avg = ncbi_genomes / ncbi_taxa if ncbi_taxa > 0 else 0
    eukprot_avg = eukprot_genomes / eukprot_taxa if eukprot_taxa > 0 else 0
    
    return {
        'ncbi_taxa': ncbi_taxa,
        'eukprot_taxa': eukprot_taxa,
        'shared_taxa': shared_taxa,
        'ncbi_genomes': ncbi_genomes,
        'eukprot_genomes': eukprot_genomes,
        'ncbi_avg': ncbi_avg,
        'eukprot_avg': eukprot_avg
    }

# Calculate diversity metrics for each taxonomic level
phylum_diversity = calculate_diversity(phylum_data)
family_diversity = calculate_diversity(family_data)
genus_diversity = calculate_diversity(genus_data)

# Create a figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# Function to create a bar chart comparing taxonomic diversity
def create_diversity_chart(ax, ncbi_count, eukprot_count, shared_count, title):
    # Data for plotting
    categories = ['NCBI', 'EukProt', 'Shared']
    counts = [ncbi_count, eukprot_count, shared_count]
    
    # Create bar chart
    bars = ax.bar(categories, counts, color=['#1f77b4', '#2ca02c', '#9467bd'])
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom')
    
    # Add title and labels
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Number of Taxa', fontsize=12)
    
    # Add percentage of shared taxa
    if ncbi_count > 0 and eukprot_count > 0:
        ncbi_percent = (shared_count / ncbi_count * 100)
        eukprot_percent = (shared_count / eukprot_count * 100)
        ax.text(0.5, 0.05, 
                f"Shared: {ncbi_percent:.1f}% of NCBI\n{eukprot_percent:.1f}% of EukProt",
                transform=ax.transAxes, ha='center', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

# Create diversity charts for each taxonomic level
create_diversity_chart(ax1, 
                      phylum_diversity['ncbi_taxa'], 
                      phylum_diversity['eukprot_taxa'], 
                      phylum_diversity['shared_taxa'], 
                      'Phylum-Level Diversity')

create_diversity_chart(ax2, 
                      family_diversity['ncbi_taxa'], 
                      family_diversity['eukprot_taxa'], 
                      family_diversity['shared_taxa'], 
                      'Family-Level Diversity')

create_diversity_chart(ax3, 
                      genus_diversity['ncbi_taxa'], 
                      genus_diversity['eukprot_taxa'], 
                      genus_diversity['shared_taxa'], 
                      'Genus-Level Diversity')

# Add overall title
plt.suptitle('Eukaryote Taxonomic Diversity Comparison', fontsize=16, y=0.98)

# Add summary statistics
summary_text = (
    f"NCBI: {int(phylum_diversity['ncbi_genomes'])} genomes across {phylum_diversity['ncbi_taxa']} phyla, "
    f"{family_diversity['ncbi_taxa']} families, and {genus_diversity['ncbi_taxa']} genera\n"
    f"EukProt: {int(phylum_diversity['eukprot_genomes'])} genomes across {phylum_diversity['eukprot_taxa']} phyla, "
    f"{family_diversity['eukprot_taxa']} families, and {genus_diversity['eukprot_taxa']} genera"
)

fig.text(0.5, 0.01, summary_text, ha='center', fontsize=12, 
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Save figure
output_path = os.path.join(visualization_dir, "eukaryote_taxonomic_diversity.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nEukaryote Taxonomic Diversity Summary:")
print(f"NCBI: {int(phylum_diversity['ncbi_genomes'])} genomes")
print(f"  - {phylum_diversity['ncbi_taxa']} phyla (avg. {phylum_diversity['ncbi_avg']:.1f} genomes/phylum)")
print(f"  - {family_diversity['ncbi_taxa']} families (avg. {family_diversity['ncbi_avg']:.1f} genomes/family)")
print(f"  - {genus_diversity['ncbi_taxa']} genera (avg. {genus_diversity['ncbi_avg']:.1f} genomes/genus)")
print(f"\nEukProt: {int(phylum_diversity['eukprot_genomes'])} genomes")
print(f"  - {phylum_diversity['eukprot_taxa']} phyla (avg. {phylum_diversity['eukprot_avg']:.1f} genomes/phylum)")
print(f"  - {family_diversity['eukprot_taxa']} families (avg. {family_diversity['eukprot_avg']:.1f} genomes/family)")
print(f"  - {genus_diversity['eukprot_taxa']} genera (avg. {genus_diversity['eukprot_avg']:.1f} genomes/genus)")
print(f"\nShared taxa:")
print(f"  - {phylum_diversity['shared_taxa']} phyla")
print(f"  - {family_diversity['shared_taxa']} families")
print(f"  - {genus_diversity['shared_taxa']} genera")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
