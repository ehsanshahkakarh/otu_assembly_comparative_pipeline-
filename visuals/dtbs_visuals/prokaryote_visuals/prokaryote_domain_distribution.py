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
merged_tables_dir = os.path.normpath(os.path.join(script_dir, "..", "merged_tables", "domain_taxonomic_tables"))
visualization_dir = script_dir  # Current directory

# Define input files for bacteria and archaea
bacteria_phylum_file = os.path.join(merged_tables_dir, "bacteria_phylum.csv")
archaea_phylum_file = os.path.join(merged_tables_dir, "archaea_phylum.csv")

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

# Check files
bacteria_phylum_file = check_file_exists(bacteria_phylum_file, "bacteria", "phylum")
archaea_phylum_file = check_file_exists(archaea_phylum_file, "archaea", "phylum")

# Load data
print(f"Loading bacteria phylum data from: {bacteria_phylum_file}")
bacteria_phylum_data = pd.read_csv(bacteria_phylum_file)

print(f"Loading archaea phylum data from: {archaea_phylum_file}")
archaea_phylum_data = pd.read_csv(archaea_phylum_file)

# Create a figure with 2x2 subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 14))

# 1. Domain Distribution Comparison (Pie Charts)
def create_domain_distribution_pies(ax1, ax2, bacteria_data, archaea_data):
    # Calculate total genomes for each database and domain
    bacteria_gtdb = bacteria_data[bacteria_data['source'].str.contains('GTDB')]['genome_count'].sum()
    bacteria_ncbi = bacteria_data[bacteria_data['source'].str.contains('NCBI')]['genome_count'].sum()
    archaea_gtdb = archaea_data[archaea_data['source'].str.contains('GTDB')]['genome_count'].sum()
    archaea_ncbi = archaea_data[archaea_data['source'].str.contains('NCBI')]['genome_count'].sum()
    
    # Create data for pie charts
    gtdb_sizes = [bacteria_gtdb, archaea_gtdb]
    ncbi_sizes = [bacteria_ncbi, archaea_ncbi]
    labels = ['Bacteria', 'Archaea']
    colors = ['#ff9999', '#66b3ff']
    
    # Create pie chart for GTDB
    ax1.pie(gtdb_sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90,
            wedgeprops=dict(width=0.5, edgecolor='w'))
    ax1.set_title('GTDB Domain Distribution', fontsize=14)
    
    # Create pie chart for NCBI
    ax2.pie(ncbi_sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90,
            wedgeprops=dict(width=0.5, edgecolor='w'))
    ax2.set_title('NCBI Domain Distribution', fontsize=14)
    
    # Add text with total counts
    ax1.text(0, -1.2, f"Total: {bacteria_gtdb + archaea_gtdb:,} genomes\nBacteria: {bacteria_gtdb:,}\nArchaea: {archaea_gtdb:,}",
             ha='center', fontsize=12)
    ax2.text(0, -1.2, f"Total: {bacteria_ncbi + archaea_ncbi:,} genomes\nBacteria: {bacteria_ncbi:,}\nArchaea: {archaea_ncbi:,}",
             ha='center', fontsize=12)
    
    return {
        'bacteria_gtdb': bacteria_gtdb,
        'bacteria_ncbi': bacteria_ncbi,
        'archaea_gtdb': archaea_gtdb,
        'archaea_ncbi': archaea_ncbi
    }

# 2. Taxonomic Diversity Comparison (Bar Charts)
def create_taxonomic_diversity_bars(ax, bacteria_data, archaea_data, database):
    # Filter data for the specified database
    if database == 'GTDB':
        bacteria_filtered = bacteria_data[bacteria_data['source'].str.contains('GTDB')]
        archaea_filtered = archaea_data[archaea_data['source'].str.contains('GTDB')]
        color = '#ff7f0e'  # orange
    else:  # NCBI
        bacteria_filtered = bacteria_data[bacteria_data['source'].str.contains('NCBI')]
        archaea_filtered = archaea_data[archaea_data['source'].str.contains('NCBI')]
        color = '#1f77b4'  # blue
    
    # Count unique phyla
    bacteria_phyla = len(bacteria_filtered['taxon_name'].unique())
    archaea_phyla = len(archaea_filtered['taxon_name'].unique())
    
    # Create bar chart
    domains = ['Bacteria', 'Archaea']
    phyla_counts = [bacteria_phyla, archaea_phyla]
    
    bars = ax.bar(domains, phyla_counts, color=color, alpha=0.7)
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom', fontsize=12)
    
    # Add title and labels
    ax.set_title(f"{database} Phylum Diversity by Domain", fontsize=14)
    ax.set_ylabel('Number of Unique Phyla', fontsize=12)
    
    # Calculate average genomes per phylum
    bacteria_genomes = bacteria_filtered['genome_count'].sum()
    archaea_genomes = archaea_filtered['genome_count'].sum()
    
    bacteria_avg = bacteria_genomes / bacteria_phyla if bacteria_phyla > 0 else 0
    archaea_avg = archaea_genomes / archaea_phyla if archaea_phyla > 0 else 0
    
    # Add text with averages
    ax.text(0.5, 0.05, 
            f"Avg. genomes per phylum:\nBacteria: {bacteria_avg:.1f}\nArchaea: {archaea_avg:.1f}", 
            ha='center', va='center', transform=ax.transAxes, fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    return {
        'bacteria_phyla': bacteria_phyla,
        'archaea_phyla': archaea_phyla,
        'bacteria_genomes': bacteria_genomes,
        'archaea_genomes': archaea_genomes,
        'bacteria_avg': bacteria_avg,
        'archaea_avg': archaea_avg
    }

# 3. Top Phyla Comparison (Horizontal Bar Charts)
def create_top_phyla_comparison(ax, data, domain, top_n=10):
    # Filter data for the domain
    domain_data = data.copy()
    
    # Create separate dataframes for GTDB and NCBI
    gtdb_data = domain_data[domain_data['source'].str.contains('GTDB')].copy()
    ncbi_data = domain_data[domain_data['source'].str.contains('NCBI')].copy()
    
    # Group by taxon name and sum genome counts
    gtdb_grouped = gtdb_data.groupby('taxon_name')['genome_count'].sum().reset_index()
    ncbi_grouped = ncbi_data.groupby('taxon_name')['genome_count'].sum().reset_index()
    
    # Sort and get top N
    gtdb_top = gtdb_grouped.sort_values('genome_count', ascending=False).head(top_n)
    ncbi_top = ncbi_grouped.sort_values('genome_count', ascending=False).head(top_n)
    
    # Get union of top phyla from both databases
    top_phyla = pd.concat([gtdb_top['taxon_name'], ncbi_top['taxon_name']]).unique()
    
    # Limit to top_n phyla
    if len(top_phyla) > top_n:
        top_phyla = top_phyla[:top_n]
    
    # Create dataframe for plotting
    plot_data = pd.DataFrame({'taxon_name': top_phyla})
    plot_data['GTDB'] = plot_data['taxon_name'].map(gtdb_grouped.set_index('taxon_name')['genome_count']).fillna(0)
    plot_data['NCBI'] = plot_data['taxon_name'].map(ncbi_grouped.set_index('taxon_name')['genome_count']).fillna(0)
    
    # Sort by total genome count
    plot_data['total'] = plot_data['GTDB'] + plot_data['NCBI']
    plot_data = plot_data.sort_values('total', ascending=True)
    
    # Create horizontal bar chart
    y_pos = np.arange(len(plot_data))
    width = 0.35
    
    ax.barh(y_pos - width/2, plot_data['GTDB'], width, label='GTDB', color='#ff7f0e', alpha=0.7)
    ax.barh(y_pos + width/2, plot_data['NCBI'], width, label='NCBI', color='#1f77b4', alpha=0.7)
    
    # Add labels and title
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_data['taxon_name'])
    ax.set_xlabel('Number of Genomes', fontsize=12)
    ax.set_title(f'Top {domain} Phyla Comparison', fontsize=14)
    ax.legend()
    
    # Add grid lines
    ax.grid(axis='x', linestyle='--', alpha=0.7)
    
    return plot_data

# Create the visualizations
domain_stats = create_domain_distribution_pies(ax1, ax2, bacteria_phylum_data, archaea_phylum_data)
gtdb_diversity = create_taxonomic_diversity_bars(ax3, bacteria_phylum_data, archaea_phylum_data, 'GTDB')
ncbi_diversity = create_taxonomic_diversity_bars(ax4, bacteria_phylum_data, archaea_phylum_data, 'NCBI')

# Add overall title
plt.suptitle('Prokaryotic Domain Distribution and Diversity: GTDB vs NCBI', fontsize=16, y=0.98)

# Add summary statistics as text
summary_text = (
    f"GTDB: {domain_stats['bacteria_gtdb'] + domain_stats['archaea_gtdb']:,} genomes "
    f"({domain_stats['bacteria_gtdb']:,} Bacteria, {domain_stats['archaea_gtdb']:,} Archaea) "
    f"across {gtdb_diversity['bacteria_phyla'] + gtdb_diversity['archaea_phyla']} phyla\n"
    f"NCBI: {domain_stats['bacteria_ncbi'] + domain_stats['archaea_ncbi']:,} genomes "
    f"({domain_stats['bacteria_ncbi']:,} Bacteria, {domain_stats['archaea_ncbi']:,} Archaea) "
    f"across {ncbi_diversity['bacteria_phyla'] + ncbi_diversity['archaea_phyla']} phyla"
)

fig.text(0.5, 0.01, summary_text, ha='center', fontsize=12, 
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.9))

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Save figure
output_path = os.path.join(visualization_dir, "prokaryote_domain_distribution.png")
print(f"Saving figure to: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Print summary statistics
print("\nProkaryotic Domain Distribution and Diversity Summary:")
print(f"GTDB: {domain_stats['bacteria_gtdb'] + domain_stats['archaea_gtdb']:,} genomes")
print(f"  - Bacteria: {domain_stats['bacteria_gtdb']:,} genomes ({domain_stats['bacteria_gtdb']/(domain_stats['bacteria_gtdb'] + domain_stats['archaea_gtdb'])*100:.1f}%)")
print(f"  - Archaea: {domain_stats['archaea_gtdb']:,} genomes ({domain_stats['archaea_gtdb']/(domain_stats['bacteria_gtdb'] + domain_stats['archaea_gtdb'])*100:.1f}%)")
print(f"  - {gtdb_diversity['bacteria_phyla']} bacterial phyla (avg. {gtdb_diversity['bacteria_avg']:.1f} genomes/phylum)")
print(f"  - {gtdb_diversity['archaea_phyla']} archaeal phyla (avg. {gtdb_diversity['archaea_avg']:.1f} genomes/phylum)")

print(f"\nNCBI: {domain_stats['bacteria_ncbi'] + domain_stats['archaea_ncbi']:,} genomes")
print(f"  - Bacteria: {domain_stats['bacteria_ncbi']:,} genomes ({domain_stats['bacteria_ncbi']/(domain_stats['bacteria_ncbi'] + domain_stats['archaea_ncbi'])*100:.1f}%)")
print(f"  - Archaea: {domain_stats['archaea_ncbi']:,} genomes ({domain_stats['archaea_ncbi']/(domain_stats['bacteria_ncbi'] + domain_stats['archaea_ncbi'])*100:.1f}%)")
print(f"  - {ncbi_diversity['bacteria_phyla']} bacterial phyla (avg. {ncbi_diversity['bacteria_avg']:.1f} genomes/phylum)")
print(f"  - {ncbi_diversity['archaea_phyla']} archaeal phyla (avg. {ncbi_diversity['archaea_avg']:.1f} genomes/phylum)")

print(f"\nSaved to: {output_path}")

# Show the plot
plt.show()
