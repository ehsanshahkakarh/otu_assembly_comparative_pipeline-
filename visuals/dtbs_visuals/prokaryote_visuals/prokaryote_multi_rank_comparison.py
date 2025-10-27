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
merged_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "domain_taxonomic_tables"))
visualization_dir = script_dir  # Current directory

# Define input files for bacteria
bacteria_phylum_file = os.path.join(merged_tables_dir, "bacteria_phylum.csv")
bacteria_family_file = os.path.join(merged_tables_dir, "bacteria_family.csv")
bacteria_genus_file = os.path.join(merged_tables_dir, "bacteria_genus.csv")

# Define input files for archaea
archaea_phylum_file = os.path.join(merged_tables_dir, "archaea_phylum.csv")
archaea_family_file = os.path.join(merged_tables_dir, "archaea_family.csv")
archaea_genus_file = os.path.join(merged_tables_dir, "archaea_genus.csv")

# Check if input files exist and use fallback paths if needed
def check_file_exists(file_path, domain, rank):
    if not os.path.exists(file_path):
        # Try absolute path as fallback
        fallback_path = f"/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/merged_tables/domain_taxonomic_tables/{domain}_{rank}.csv"
        if not os.path.exists(fallback_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
        return fallback_path
    return file_path

# Check all files
bacteria_phylum_file = check_file_exists(bacteria_phylum_file, "bacteria", "phylum")
bacteria_family_file = check_file_exists(bacteria_family_file, "bacteria", "family")
bacteria_genus_file = check_file_exists(bacteria_genus_file, "bacteria", "genus")
archaea_phylum_file = check_file_exists(archaea_phylum_file, "archaea", "phylum")
archaea_family_file = check_file_exists(archaea_family_file, "archaea", "family")
archaea_genus_file = check_file_exists(archaea_genus_file, "archaea", "genus")

# Load data for bacteria
print(f"Loading bacteria phylum data from: {bacteria_phylum_file}")
bacteria_phylum_data = pd.read_csv(bacteria_phylum_file)

print(f"Loading bacteria family data from: {bacteria_family_file}")
bacteria_family_data = pd.read_csv(bacteria_family_file)

print(f"Loading bacteria genus data from: {bacteria_genus_file}")
bacteria_genus_data = pd.read_csv(bacteria_genus_file)

# Load data for archaea
print(f"Loading archaea phylum data from: {archaea_phylum_file}")
archaea_phylum_data = pd.read_csv(archaea_phylum_file)

print(f"Loading archaea family data from: {archaea_family_file}")
archaea_family_data = pd.read_csv(archaea_family_file)

print(f"Loading archaea genus data from: {archaea_genus_file}")
archaea_genus_data = pd.read_csv(archaea_genus_file)

# Function to create overlap donut charts
def create_overlap_donut(ax, data, rank_name):
    # Create sets of taxa in each database
    gtdb_taxa = set(data[data['source'].str.contains('GTDB')]['taxon_name'])
    ncbi_taxa = set(data[data['source'].str.contains('NCBI')]['taxon_name'])
    
    # Calculate sizes for overlap analysis
    only_gtdb = len(gtdb_taxa - ncbi_taxa)
    only_ncbi = len(ncbi_taxa - gtdb_taxa)
    both = len(gtdb_taxa.intersection(ncbi_taxa))
    
    # Create data for the donut chart
    sizes = [only_gtdb, both, only_ncbi]
    labels = ['GTDB only', 'Both', 'NCBI only']
    colors = ['#ff7f0e', '#9467bd', '#1f77b4']  # orange, purple, blue
    
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
    total_gtdb = len(gtdb_taxa)
    total_ncbi = len(ncbi_taxa)
    
    # Add text with counts
    ax.text(0.5, -0.15, f"Total: {total_taxa}", 
            ha='center', va='center', transform=ax.transAxes, fontsize=9)
    
    return {
        'only_gtdb': only_gtdb,
        'only_ncbi': only_ncbi,
        'both': both,
        'total': total_taxa,
        'total_gtdb': total_gtdb,
        'total_ncbi': total_ncbi
    }

# Function to create genome count bar charts
def create_count_bars(ax, data, rank_name):
    # Calculate total genomes in each database
    gtdb_data = data[data['source'].str.contains('GTDB')]
    ncbi_data = data[data['source'].str.contains('NCBI')]
    
    gtdb_genomes = gtdb_data['genome_count'].sum() if not gtdb_data.empty else 0
    ncbi_genomes = ncbi_data['genome_count'].sum() if not ncbi_data.empty else 0
    
    # Create bar chart
    bars = ax.bar(['GTDB', 'NCBI'], [gtdb_genomes, ncbi_genomes], 
                  color=['#ff7f0e', '#1f77b4'])  # orange, blue
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom', fontsize=9)
    
    # Add title and labels
    ax.set_title(f"{rank_name.capitalize()} Genome Counts", fontsize=12)
    ax.set_ylabel('Number of Genomes', fontsize=10)
    
    # Calculate average genomes per taxon
    gtdb_taxa = len(gtdb_data['taxon_name'].unique()) if not gtdb_data.empty else 0
    ncbi_taxa = len(ncbi_data['taxon_name'].unique()) if not ncbi_data.empty else 0
    
    gtdb_avg = gtdb_genomes / gtdb_taxa if gtdb_taxa > 0 else 0
    ncbi_avg = ncbi_genomes / ncbi_taxa if ncbi_taxa > 0 else 0
    
    # Add text with averages
    ax.text(0.5, -0.15, 
            f"Avg. genomes per {rank_name}:\nGTDB: {gtdb_avg:.1f}, NCBI: {ncbi_avg:.1f}", 
            ha='center', va='center', transform=ax.transAxes, fontsize=9)
    
    return {
        'gtdb_genomes': gtdb_genomes,
        'ncbi_genomes': ncbi_genomes,
        'gtdb_avg': gtdb_avg,
        'ncbi_avg': ncbi_avg
    }

# Function to create top taxa bar charts
def create_top_taxa_bars(ax, data, database, rank_name, top_n=5):
    # Filter and sort data
    if database == 'GTDB':
        filtered_data = data[data['source'].str.contains('GTDB')].copy()
        color = '#ff7f0e'  # orange
    else:  # NCBI
        filtered_data = data[data['source'].str.contains('NCBI')].copy()
        color = '#1f77b4'  # blue
    
    # Group by taxon name and sum genome counts
    grouped_data = filtered_data.groupby('taxon_name')['genome_count'].sum().reset_index()
    
    # Sort and get top N
    top_data = grouped_data.sort_values('genome_count', ascending=False).head(top_n)
    
    # Create horizontal bar chart
    bars = ax.barh(top_data['taxon_name'], top_data['genome_count'], color=color, alpha=0.7)
    
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
    
    return top_data

# Create figures for bacteria and archaea
def create_domain_dashboard(domain, phylum_data, family_data, genus_data):
    # Create figure with GridSpec
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1])
    
    # Create overlap donut charts (first row)
    phylum_overlap = create_overlap_donut(fig.add_subplot(gs[0, 0]), phylum_data, 'phylum')
    family_overlap = create_overlap_donut(fig.add_subplot(gs[0, 1]), family_data, 'family')
    genus_overlap = create_overlap_donut(fig.add_subplot(gs[0, 2]), genus_data, 'genus')
    
    # Create count comparison bar charts (second row)
    phylum_counts = create_count_bars(fig.add_subplot(gs[1, 0]), phylum_data, 'phylum')
    family_counts = create_count_bars(fig.add_subplot(gs[1, 1]), family_data, 'family')
    genus_counts = create_count_bars(fig.add_subplot(gs[1, 2]), genus_data, 'genus')
    
    # Create top taxa bar charts (third row)
    # For phylum, show top taxa for both GTDB and NCBI side by side
    ax_phylum_top = fig.add_subplot(gs[2, 0])
    ax_phylum_top.set_title("Top Phyla Comparison", fontsize=12)
    ax_phylum_top.axis('off')
    
    # Create a nested gridspec for the phylum comparison
    gs_phylum = GridSpec(2, 1, figure=fig, height_ratios=[1, 1])
    gs_phylum.update(left=0.05, right=0.32, bottom=0.05, top=0.28)
    
    top_phyla_gtdb = create_top_taxa_bars(fig.add_subplot(gs_phylum[0, 0]), phylum_data, 'GTDB', 'phylum')
    top_phyla_ncbi = create_top_taxa_bars(fig.add_subplot(gs_phylum[1, 0]), phylum_data, 'NCBI', 'phylum')
    
    # For family, show top GTDB families
    top_families_gtdb = create_top_taxa_bars(fig.add_subplot(gs[2, 1]), family_data, 'GTDB', 'family')
    
    # For genus, show top NCBI genera
    top_genera_ncbi = create_top_taxa_bars(fig.add_subplot(gs[2, 2]), genus_data, 'NCBI', 'genus')
    
    # Add overall title
    plt.suptitle(f'{domain.capitalize()} Taxonomic Comparison Across Ranks: GTDB vs NCBI', fontsize=16, y=0.98)
    
    # Add summary statistics as text
    summary_text = (
        f"GTDB: {int(phylum_counts['gtdb_genomes'])} genomes across {phylum_overlap['total_gtdb']} phyla, "
        f"{family_overlap['total_gtdb']} families, and {genus_overlap['total_gtdb']} genera\n"
        f"NCBI: {int(phylum_counts['ncbi_genomes'])} genomes across {phylum_overlap['total_ncbi']} phyla, "
        f"{family_overlap['total_ncbi']} families, and {genus_overlap['total_ncbi']} genera\n"
        f"Shared taxa: {phylum_overlap['both']} phyla, {family_overlap['both']} families, and {genus_overlap['both']} genera"
    )
    
    fig.text(0.5, 0.01, summary_text, ha='center', fontsize=12, 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    
    # Save figure
    output_path = os.path.join(visualization_dir, f"{domain}_multi_rank_comparison.png")
    print(f"Saving figure to: {output_path}")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Print summary statistics
    print(f"\n{domain.capitalize()} Multi-Rank Comparison Summary:")
    print(f"GTDB: {int(phylum_counts['gtdb_genomes'])} genomes")
    print(f"  - {phylum_overlap['total_gtdb']} phyla (avg. {phylum_counts['gtdb_avg']:.1f} genomes/phylum)")
    print(f"  - {family_overlap['total_gtdb']} families (avg. {family_counts['gtdb_avg']:.1f} genomes/family)")
    print(f"  - {genus_overlap['total_gtdb']} genera (avg. {genus_counts['gtdb_avg']:.1f} genomes/genus)")
    print(f"\nNCBI: {int(phylum_counts['ncbi_genomes'])} genomes")
    print(f"  - {phylum_overlap['total_ncbi']} phyla (avg. {phylum_counts['ncbi_avg']:.1f} genomes/phylum)")
    print(f"  - {family_overlap['total_ncbi']} families (avg. {family_counts['ncbi_avg']:.1f} genomes/family)")
    print(f"  - {genus_overlap['total_ncbi']} genera (avg. {genus_counts['ncbi_avg']:.1f} genomes/genus)")
    print(f"\nShared taxa:")
    print(f"  - {phylum_overlap['both']} phyla ({(phylum_overlap['both']/phylum_overlap['total']*100):.1f}% of total)")
    print(f"  - {family_overlap['both']} families ({(family_overlap['both']/family_overlap['total']*100):.1f}% of total)")
    print(f"  - {genus_overlap['both']} genera ({(genus_overlap['both']/genus_overlap['total']*100):.1f}% of total)")
    
    print(f"\nTop phyla in GTDB:")
    for i, row in enumerate(top_phyla_gtdb.itertuples(), 1):
        print(f"{i}. {row.taxon_name}: {int(row.genome_count)} genomes")
    
    print(f"\nTop phyla in NCBI:")
    for i, row in enumerate(top_phyla_ncbi.itertuples(), 1):
        print(f"{i}. {row.taxon_name}: {int(row.genome_count)} genomes")
    
    print(f"\nSaved to: {output_path}")
    
    return fig

# Create dashboards for bacteria and archaea
bacteria_fig = create_domain_dashboard('bacteria', bacteria_phylum_data, bacteria_family_data, bacteria_genus_data)
archaea_fig = create_domain_dashboard('archaea', archaea_phylum_data, archaea_family_data, archaea_genus_data)

# Show the plots
plt.show()
