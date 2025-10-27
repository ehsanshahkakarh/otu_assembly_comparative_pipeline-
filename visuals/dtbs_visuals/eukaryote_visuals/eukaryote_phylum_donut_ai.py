import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

# Set style
plt.style.use('default')

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths
eukaryota_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "eukaryota_taxonomic_tables"))
visualization_dir = script_dir

# Input file
phylum_file = os.path.join(eukaryota_tables_dir, "eukaryota_phylum.csv")

# Check file existence
if not os.path.exists(phylum_file):
    phylum_file = "/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/eukaryota_taxonomic_tables/eukaryota_phylum.csv"
    if not os.path.exists(phylum_file):
        print(f"Error: File not found: {phylum_file}")
        sys.exit(1)

print(f"Loading data from: {phylum_file}")
phylum_data = pd.read_csv(phylum_file)

# Define sets for overlap
ncbi_phyla = set(phylum_data[phylum_data['ncbi_genome_count'] > 0]['taxon_name'])
eukprot_phyla = set(phylum_data[phylum_data['eukprot_genome_count'] > 0]['taxon_name'])

# Overlap stats
only_ncbi = len(ncbi_phyla - eukprot_phyla)
only_eukprot = len(eukprot_phyla - ncbi_phyla)
both = len(ncbi_phyla.intersection(eukprot_phyla))

# Totals and percents
total_phyla = len(phylum_data['taxon_name'].unique())
total_ncbi = len(ncbi_phyla)
total_eukprot = len(eukprot_phyla)
ncbi_percent = (total_ncbi / total_phyla * 100) if total_phyla else 0
eukprot_percent = (total_eukprot / total_phyla * 100) if total_phyla else 0

# Overlap counts
counts = {
    "NCBI only": only_ncbi,
    "Both": both,
    "EukProt only": only_eukprot
}

colors = {
    "NCBI only": "#1f77b4",       # Blue
    "Both": "#a05195",            # Purple
    "EukProt only": "#2ca02c"     # Green
}

# Setup figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True, gridspec_kw={'width_ratios': [1, 1]})

# --- Donut Chart ---
labels = list(counts.keys())
sizes = list(counts.values())
donut_colors = [colors[label] for label in labels]

wedges, texts, autotexts = ax1.pie(
    sizes,
    labels=labels,
    autopct='%1.1f%%',
    startangle=90,
    pctdistance=0.85,
    colors=donut_colors,
    wedgeprops=dict(width=0.4, edgecolor='white')
)

for autotext in autotexts:
    autotext.set_fontsize(12)
    autotext.set_weight('bold')
    autotext.set_color('white')

ax1.set_title("Eukaryote Phylum Coverage\nNCBI vs EukProt", fontsize=16, pad=10)
centre_circle = plt.Circle((0, 0), 0.55, fc='white')
ax1.add_artist(centre_circle)
ax1.axis('equal')

# --- Overlap Table ---
overlap_phyla = ncbi_phyla.intersection(eukprot_phyla)
overlap_data = phylum_data[phylum_data['taxon_name'].isin(overlap_phyla)].copy()
overlap_data['total_genomes'] = overlap_data['ncbi_genome_count'] + overlap_data['eukprot_genome_count']
overlap_data = overlap_data.sort_values('total_genomes', ascending=False)
total_ncbi_genomes = phylum_data['ncbi_genome_count'].sum()
total_eukprot_genomes = phylum_data['eukprot_genome_count'].sum()
overlap_data['ncbi_percent'] = (overlap_data['ncbi_genome_count'] / total_ncbi_genomes * 100).round(2)
overlap_data['eukprot_percent'] = (overlap_data['eukprot_genome_count'] / total_eukprot_genomes * 100).round(2)

if len(overlap_data) > 0:
    overlap_data_trimmed = overlap_data.head(10)
    table_data = []
    for _, row in overlap_data_trimmed.iterrows():
        table_data.append([
            row['taxon_name'],
            f"{int(row['ncbi_genome_count'])} ({row['ncbi_percent']}%)",
            f"{int(row['eukprot_genome_count'])} ({row['eukprot_percent']}%)"
        ])
    
    table = ax2.table(
        cellText=table_data,
        colLabels=["Phylum", "NCBI Count (%)", "EukProt Count (%)"],
        cellLoc='center',
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.5)

    for (i, j), cell in table.get_celld().items():
        if i == 0:
            cell.set_fontsize(12)
            cell.set_text_props(weight='bold')
            cell.set_facecolor("#f0f0f0")

    ax2.set_title("Top Shared Phyla (by Total Genomes)", fontsize=16, pad=10)
    ax2.axis('off')
else:
    ax2.text(0.5, 0.5, "No shared phyla", ha='center', va='center', fontsize=14)
    ax2.axis('off')

# Summary text
summary_text = (
    f"Total unique phyla: {total_phyla} | "
    f"NCBI: {total_ncbi} ({ncbi_percent:.1f}%) | "
    f"EukProt: {total_eukprot} ({eukprot_percent:.1f}%) | "
    f"Overlap: {both} ({(both / total_phyla * 100):.1f}%)"
)

fig.text(0.5, 0.02, summary_text, ha='center', fontsize=12,
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='gray', alpha=0.9))

# Save figure
output_path = os.path.join(visualization_dir, "eukaryote_phylum_donut_overlap.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

# Console summary
print("\nEukaryote Phylum Coverage Summary:")
print(f"Total unique phyla: {total_phyla}")
print(f"NCBI-only phyla: {only_ncbi}")
print(f"EukProt-only phyla: {only_eukprot}")
print(f"Phyla in both databases: {both}")
print(f"NCBI total phyla: {total_ncbi} ({ncbi_percent:.1f}%)")
print(f"EukProt total phyla: {total_eukprot} ({eukprot_percent:.1f}%)")
if len(overlap_phyla) > 0:
    print("\nPhyla in both databases:")
    for i, phylum in enumerate(sorted(overlap_phyla), 1):
        row = phylum_data[phylum_data['taxon_name'] == phylum].iloc[0]
        print(f"{i}. {phylum} - NCBI: {row['ncbi_genome_count']}, EukProt: {row['eukprot_genome_count']}")
print(f"\nSaved to: {output_path}")
