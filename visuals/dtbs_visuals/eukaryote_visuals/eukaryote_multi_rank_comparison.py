import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from matplotlib.gridspec import GridSpec

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
for file_path, rank in zip([phylum_file, family_file, genus_file], ["phylum", "family", "genus"]):
    if not os.path.exists(file_path):
        # Try absolute path as fallback
        fallback_path = f"/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/eukaryota_taxonomic_tables/eukaryota_{rank}.csv"
        if not os.path.exists(fallback_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
        else:
            if rank == "phylum":
                phylum_file = fallback_path
            elif rank == "family":
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

# Create a figure with a complex layout using GridSpec
fig = plt.figure(figsize=(18, 12))
gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1])

# Function to create donut charts for overlap visualization
def create_overlap_donut(ax, data, rank_name):
    # Create sets of taxa in each database
    ncbi_taxa = set(data[data['ncbi_genome_count'] > 0]['taxon_name'])
    eukprot_taxa = set(data[data['eukprot_genome_count'] > 0]['taxon_name'])
    
    # Calculate sizes for overlap analysis
    only_ncbi = len(ncbi_taxa - eukprot_taxa)
    only_eukprot = len(eukprot_taxa - ncbi_taxa)
    both = len(ncbi_taxa.intersection(eukprot_taxa))
    
    # Create data for the donut chart
    sizes = [only_ncbi, both, only_eukprot]
    labels = ['NCBI only', 'Both', 'EukProt only']
    colors = ['#1f77b4', '#9467bd', '#2ca02c']  # blue, purple, green
    
    # Create donut chart
    wedges, texts, autotexts = ax.pie(
        sizes, 
        labels=None,  # We'll add a legend instead
        autopct='%1.1f%%',
        startangle=90,
        colors=colors,
        wedgeprops=dict(width=0.5, edgecolor='w'),
        textprops={'fontsize': 10}
    )
    
    # Make the percentage labels more readable
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)
    
    # Add a legend
    ax.legend(wedges, [f"{label} ({size})" for label, size in zip(labels, sizes)], 
              loc="center", bbox_to_anchor=(0.5, 0), fontsize=9)
    
    # Add title
    ax.set_title(f"{rank_name.capitalize()} Overlap", fontsize=12)
    
    # Calculate total counts
    total_taxa = len(data['taxon_name'].unique())
    total_ncbi = len(ncbi_taxa)
    total_eukprot = len(eukprot_taxa)
    
    # Add text with counts
    ax.text(0.5, -0.15, f"Total: {total_taxa}", 
            ha='center', va='center', transform=ax.transAxes, fontsize=9)
    
    return {
        'only_ncbi': only_ncbi,
        'only_eukprot': only_eukprot,
        'both': both,
        'total': total_taxa,
        'total_ncbi': total_ncbi,
        'total_eukprot': total_eukprot
    }

# Function to create bar charts for count comparison
def create_count_bars(ax, data, rank_name):
    # Calculate total genomes in each database
    ncbi_genomes = data['ncbi_genome_count'].sum()
    eukprot_genomes = data['eukprot_genome_count'].sum()
    
    # Create bar chart
    bars = ax.bar(['NCBI', 'EukProt'], [ncbi_genomes, eukprot_genomes], 
                  color=['#1f77b4', '#2ca02c'])
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom', fontsize=9)
    
    # Add title and labels
    ax.set_title(f"{rank_name.capitalize()} Genome Counts", fontsize=12)
    ax.set_ylabel('Number of Genomes', fontsize=10)
    
    # Calculate average genomes per taxon
    ncbi_taxa = len(data[data['ncbi_genome_count'] > 0])
    eukprot_taxa = len(data[data['eukprot_genome_count'] > 0])
    
    ncbi_avg = ncbi_genomes / ncbi_taxa if ncbi_taxa > 0 else 0
    eukprot_avg = eukprot_genomes / eukprot_taxa if eukprot_taxa > 0 else 0
    
    # Add text with averages
    ax.text(0.5, -0.15, 
            f"Avg. genomes per {rank_name}:\nNCBI: {ncbi_avg:.1f}, EukProt: {eukprot_avg:.1f}", 
            ha='center', va='center', transform=ax.transAxes, fontsize=9)
    
    return {
        'ncbi_genomes': ncbi_genomes,
        'eukprot_genomes': eukprot_genomes,
        'ncbi_avg': ncbi_avg,
        'eukprot_avg': eukprot_avg
    }

# Function to create top taxa bar charts
def create_top_taxa_bars(ax, data, database, rank_name, top_n=5):
    # Filter and sort data
    if database == 'NCBI':
        filtered_data = data[data['ncbi_genome_count'] > 0].copy()
        filtered_data = filtered_data.sort_values('ncbi_genome_count', ascending=False).head(top_n)
        count_col = 'ncbi_genome_count'
        color = '#1f77b4'  # blue
    else:  # EukProt
        filtered_data = data[data['eukprot_genome_count'] > 0].copy()
        filtered_data = filtered_data.sort_values('eukprot_genome_count', ascending=False).head(top_n)
        count_col = 'eukprot_genome_count'
        color = '#2ca02c'  # green
    
    # Create horizontal bar chart
    bars = ax.barh(filtered_data['taxon_name'], filtered_data[count_col], color=color, alpha=0.7)
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
                f'{int(width)}', ha='left', va='center', fontsize=9)
    
    # Add title and labels
    ax.set_title(f"Top {top_n} {rank_name.capitalize()}s in {database}", fontsize=12)
    ax.set_xlabel('Number of Genomes', fontsize=10)
    
    # Adjust y-axis labels
    ax.tick_params(axis='y', labelsize=9)
    
    return filtered_data

# Create overlap donut charts (first row)
phylum_overlap = create_overlap_donut(fig.add_subplot(gs[0, 0]), phylum_data, 'phylum')
family_overlap = create_overlap_donut(fig.add_subplot(gs[0, 1]), family_data, 'family')
genus_overlap = create_overlap_donut(fig.add_subplot(gs[0, 2]), genus_data, 'genus')

# Create count comparison bar charts (second row)
phylum_counts = create_count_bars(fig.add_subplot(gs[1, 0]), phylum_data, 'phylum')
family_counts = create_count_bars(fig.add_subplot(gs[1, 1]), family_data, 'family')
genus_counts = create_count_bars(fig.add_subplot(gs[1, 2]), genus_data, 'genus')

# Create top taxa bar charts (third row)
# For phylum, show top taxa for both NCBI and EukProt side by side
ax_phylum_top = fig.add_subplot(gs[2, 0])
ax_phylum_top.set_title("Top Phyla Comparison", fontsize=12)
ax_phylum_top.axis('off')

# Create a nested gridspec for the phylum comparison
gs_phylum = GridSpec(2, 1, figure=fig, height_ratios=[1, 1])
gs_phylum.update(left=0.05, right=0.32, bottom=0.05, top=0.28)

top_phyla_ncbi = create_top_taxa_bars(fig.add_subplot(gs_phylum[0, 0]), phylum_data, 'NCBI', 'phylum')
top_phyla_eukprot = create_top_taxa_bars(fig.add_subplot(gs_phylum[1, 0]), phylum_data, 'EukProt', 'phylum')

# For family, show top NCBI families
top_families_ncbi = create_top_taxa_bars(fig.add_subplot(gs[2, 1]), family_data, 'NCBI', 'family')

# For genus, show top EukProt genera
top_genera_eukprot = create_top_taxa_bars(fig.add_subplot(gs[2, 2]), genus_data, 'EukProt', 'genus')

# Add overall title
plt.suptitle('Eukaryote Taxonomic Comparison Across Ranks: NCBI vs EukProt', fontsize=16, y=0.98)

# Add summary statistics as text
summary_text = (
    f"NCBI: {int(phylum_counts['ncbi_genomes'])} genomes across {phylum_overlap['total_ncbi']} phyla, "
    f"{family_overlap['total_ncbi']} families, and {genus_overlap['total_ncbi']} genera\n"
    f"EukProt: {int(phylum_counts['eukprot_genomes'])} genomes across {phylum_overlap['total_eukprot']} phyla, "
    f"{family_overlap['total_eukprot']} families, and {genus_overlap['total_eukprot']} genera\n"
    f"Shared taxa: {phylum_overlap['both']} phyla, {family_overlap['both']} families, and {genus_overlap['both']} genera"
)

fig.text(0.5, 0.01, summary_text, ha='center', fontsize=12, 
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Save figure
output_path = os.path.join(visualization_dir, "eukaryote_multi_rank_comparison.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nEukaryote Multi-Rank Comparison Summary:")
print(f"NCBI: {int(phylum_counts['ncbi_genomes'])} genomes")
print(f"  - {phylum_overlap['total_ncbi']} phyla (avg. {phylum_counts['ncbi_avg']:.1f} genomes/phylum)")
print(f"  - {family_overlap['total_ncbi']} families (avg. {family_counts['ncbi_avg']:.1f} genomes/family)")
print(f"  - {genus_overlap['total_ncbi']} genera (avg. {genus_counts['ncbi_avg']:.1f} genomes/genus)")
print(f"\nEukProt: {int(phylum_counts['eukprot_genomes'])} genomes")
print(f"  - {phylum_overlap['total_eukprot']} phyla (avg. {phylum_counts['eukprot_avg']:.1f} genomes/phylum)")
print(f"  - {family_overlap['total_eukprot']} families (avg. {family_counts['eukprot_avg']:.1f} genomes/family)")
print(f"  - {genus_overlap['total_eukprot']} genera (avg. {genus_counts['eukprot_avg']:.1f} genomes/genus)")
print(f"\nShared taxa:")
print(f"  - {phylum_overlap['both']} phyla ({(phylum_overlap['both']/phylum_overlap['total']*100):.1f}% of total)")
print(f"  - {family_overlap['both']} families ({(family_overlap['both']/family_overlap['total']*100):.1f}% of total)")
print(f"  - {genus_overlap['both']} genera ({(genus_overlap['both']/genus_overlap['total']*100):.1f}% of total)")

print(f"\nTop phyla in NCBI:")
for i, row in enumerate(top_phyla_ncbi.itertuples(), 1):
    print(f"{i}. {row.taxon_name}: {int(row.ncbi_genome_count)} genomes")

print(f"\nTop phyla in EukProt:")
for i, row in enumerate(top_phyla_eukprot.itertuples(), 1):
    print(f"{i}. {row.taxon_name}: {int(row.eukprot_genome_count)} genomes")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
