#!/usr/bin/env python3
"""
16S Prokaryotic Sankey Diagram with COLORED FLOWS (Python version)
==================================================================

Creates a Sankey diagram showing the flow from:
NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species

For prokaryotic (16S) data - Bacteria and Archaea
Uses shared taxonomic color mapping for FLOWS (not just nodes)
Shows top phyla by total representation + Other

Author: Enhanced Sankey Team
Date: 2026-01-28
"""

import pandas as pd
import plotly.graph_objects as go
import yaml
from pathlib import Path

print("=== 16S Prokaryotic Sankey Diagram with Colored Flows (Python) ===\n")

# Configuration: Set which domain to process
PROCESS_DOMAIN = "Bacteria"  # Change to "Archaea" for archaea plot
print(f"Processing domain: {PROCESS_DOMAIN}\n")

# Load color configuration
color_config_path = Path("../../shared_config/taxonomic_color_mapping.yaml")
if not color_config_path.exists():
    raise FileNotFoundError(f"Cannot find color config at: {color_config_path}")

with open(color_config_path, 'r') as f:
    color_config = yaml.safe_load(f)

# Load alluvial filtering configuration
alluvial_config_path = Path("../config/alluvial_filtering_config.yaml")
if not alluvial_config_path.exists():
    raise FileNotFoundError(f"Cannot find alluvial config at: {alluvial_config_path}")

with open(alluvial_config_path, 'r') as f:
    alluvial_config = yaml.safe_load(f)

# Get domain-specific configuration
if PROCESS_DOMAIN == "Bacteria":
    domain_config = alluvial_config['bacteria_16s']
    top_n = domain_config['filtering']['hybrid']['max_total_taxa']
    color_key = 'bacteria_colors'
elif PROCESS_DOMAIN == "Archaea":
    domain_config = alluvial_config['archaea_16s']
    top_n = domain_config['filtering']['top_n']
    color_key = 'archaea_colors'
else:
    raise ValueError("Invalid domain. Must be 'Bacteria' or 'Archaea'")

print(f"Top N phyla to show: {top_n}\n")

# Load the 16S data
merged_data_path = Path("../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv")
if not merged_data_path.exists():
    raise FileNotFoundError(f"Cannot find 16S merged data at: {merged_data_path}")

print("Loading 16S data...")
data_16s = pd.read_csv(merged_data_path)
print(f"16S data loaded: {len(data_16s)} rows\n")

# Filter for matched phyla for specific domain
matched_data = data_16s[
    (data_16s['match_status'] == 'matched') &
    (data_16s['census_size_count'].notna()) &
    (data_16s['census_otu_count'].notna()) &
    (data_16s['ncbi_genome_count'].notna()) &
    (data_16s['ncbi_species_count'].notna()) &
    (data_16s['phylum'].notna()) &
    (data_16s['phylum'] != 'N/A') &
    (data_16s['phylum'] != '') &
    (data_16s['domain'] == PROCESS_DOMAIN)
].copy()

print(f"Matched {PROCESS_DOMAIN.lower()} phyla found: {len(matched_data)}\n")

# Add total representation across all 4 nodes and get top N
matched_data['total_representation'] = (
    matched_data['ncbi_genome_count'] +
    matched_data['census_size_count'] +
    matched_data['census_otu_count'] +
    matched_data['ncbi_species_count']
)

top_phyla = matched_data.nlargest(top_n, 'total_representation')
print(f"Top {top_n} phyla selected: {len(top_phyla)}\n")

# Calculate total counts
total_genome_count = matched_data['ncbi_genome_count'].sum()
total_size_count = matched_data['census_size_count'].sum()
total_otu_count = matched_data['census_otu_count'].sum()
total_species_count = matched_data['ncbi_species_count'].sum()

# Calculate "Other" counts
other_genome_count = total_genome_count - top_phyla['ncbi_genome_count'].sum()
other_size_count = total_size_count - top_phyla['census_size_count'].sum()
other_otu_count = total_otu_count - top_phyla['census_otu_count'].sum()
other_species_count = total_species_count - top_phyla['ncbi_species_count'].sum()

print("Other counts (remaining phyla):")
print(f"  Other Genomes: {other_genome_count:,}")
print(f"  Other Sequences: {other_size_count:,}")
print(f"  Other OTUs: {other_otu_count:,}")
print(f"  Other Species: {other_species_count:,}\n")

# Get phyla names and assign colors
phyla_names = list(top_phyla['phylum']) + ['Other']
domain_colors = color_config[color_key]
other_color = color_config['special_colors']['Other']

# Assign colors
phyla_colors = []
for phylum in phyla_names:
    if phylum == 'Other':
        phyla_colors.append(other_color)
    elif phylum in domain_colors:
        phyla_colors.append(domain_colors[phylum])
    else:
        # Fallback color
        fallback_colors = color_config['fallback_colors'][PROCESS_DOMAIN.lower()]
        idx = len(phyla_colors) % len(fallback_colors)
        phyla_colors.append(fallback_colors[idx])

print("Assigning colors from shared taxonomic color mapping...")
print("\nCOLOR MAPPING SUMMARY:")
print("=" * 60)
for i, (phylum, color) in enumerate(zip(phyla_names, phyla_colors), 1):
    print(f"{i:2d}. {phylum:30s} {color}")
print("=" * 60)
print()

# Create nodes
# Node structure: Each phylum appears in each stage
stages = ["Genomes", "Sequences", "OTUs", "Species"]
node_labels = []
node_colors = []
node_x = []
node_y = []

# Position nodes
x_positions = [0.01, 0.33, 0.66, 0.99]

# Create node mapping
node_map = {}
node_id = 0

for stage_idx, stage in enumerate(stages):
    for phylum_idx, phylum in enumerate(phyla_names):
        node_labels.append(f"{phylum}")
        node_colors.append(phyla_colors[phylum_idx])
        node_x.append(x_positions[stage_idx])
        node_y.append(phylum_idx / len(phyla_names))
        node_map[(phylum, stage)] = node_id
        node_id += 1

print(f"Created {len(node_labels)} nodes\n")

# Create links with colors
sources = []
targets = []
values = []
link_colors = []

# Helper function to add link
def add_link(phylum, source_stage, target_stage, value, color):
    if value > 0:
        sources.append(node_map[(phylum, source_stage)])
        targets.append(node_map[(phylum, target_stage)])
        values.append(value)
        # Make link color semi-transparent
        if color.startswith('#'):
            # Convert hex to rgba with 0.4 opacity
            link_colors.append(f"rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.4)")
        else:
            link_colors.append(color)

# Create links for top phyla
for idx, row in top_phyla.iterrows():
    phylum = row['phylum']
    color = phyla_colors[phyla_names.index(phylum)]

    # Genomes -> Sequences
    add_link(phylum, "Genomes", "Sequences", row['census_size_count'], color)

    # Sequences -> OTUs
    add_link(phylum, "Sequences", "OTUs", row['census_otu_count'], color)

    # OTUs -> Species
    add_link(phylum, "OTUs", "Species", row['ncbi_species_count'], color)

# Add "Other" links
add_link("Other", "Genomes", "Sequences", other_size_count, other_color)
add_link("Other", "Sequences", "OTUs", other_otu_count, other_color)
add_link("Other", "OTUs", "Species", other_species_count, other_color)

print(f"Created {len(sources)} links with individual colors\n")

# Create the Sankey diagram
print("Creating Sankey diagram with colored flows...\n")

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="white", width=0.5),
        label=node_labels,
        color=node_colors,
        x=node_x,
        y=node_y
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors  # This colors the flows!
    )
)])

# Update layout
fig.update_layout(
    title=dict(
        text=f"{PROCESS_DOMAIN} Sankey Diagram with Colored Flows<br>" +
             f"<sub>Genomes → Sequences → OTUs → Species | Top {top_n} Phyla + Other</sub>",
        font=dict(size=20)
    ),
    font=dict(size=12, family="Arial"),
    height=800,
    width=1400,
    plot_bgcolor='white',
    paper_bgcolor='white'
)

# Save the Sankey diagram
output_file = f"sankey_16s_{PROCESS_DOMAIN.lower()}_colored_flows.html"
fig.write_html(output_file)
print(f"Sankey diagram saved: {output_file}\n")

# Also save as static image if kaleido is available
try:
    png_file = f"sankey_16s_{PROCESS_DOMAIN.lower()}_colored_flows.png"
    fig.write_image(png_file, width=1400, height=800)
    print(f"Static PNG version saved: {png_file}\n")
except Exception as e:
    print(f"Note: Could not save PNG (kaleido not available): {e}\n")

# Display summary statistics
print("=" * 70)
print(f"🧬 {PROCESS_DOMAIN.upper()} SANKEY DIAGRAM WITH COLORED FLOWS - SUMMARY")
print("=" * 70)
print()
print("📊 DATASET OVERVIEW:")
print(f"   • Domain: {PROCESS_DOMAIN}")
print(f"   • Total phyla analyzed: {len(top_phyla)} + Other")
print("   • Data flow: Genomes → Sequences → OTUs → Species")
print("   • Interactive HTML output with colored flows")
print("   • Colors from shared taxonomic color mapping")
print()
print("📈 TOTAL COUNTS:")
print(f"   • Total genomes: {total_genome_count:,}")
print(f"   • Total sequences: {total_size_count:,}")
print(f"   • Total OTUs: {total_otu_count:,}")
print(f"   • Total species: {total_species_count:,}")
print()
print("🎨 VISUALIZATION FEATURES:")
print("   • Interactive Sankey diagram with hover details")
print("   • COLORED FLOWS matching phylum colors")
print(f"   • Color-coded phyla with {len(phyla_colors)} distinct colors")
print("   • Flow widths represent absolute counts")
print("   • HTML output for web viewing and sharing")
print()
print("📋 FLOW SUMMARY:")
print("-" * 80)
print(f"{'Phylum':<30s} {'Genomes':>12s} {'Sequences':>12s} {'OTUs':>12s} {'Species':>12s} {'Color':>10s}")
print("-" * 80)

for idx, row in top_phyla.iterrows():
    phylum = row['phylum']
    color = phyla_colors[phyla_names.index(phylum)]
    print(f"{phylum[:30]:<30s} {row['ncbi_genome_count']:>12,} {row['census_size_count']:>12,} "
          f"{row['census_otu_count']:>12,} {row['ncbi_species_count']:>12,} {color:>10s}")

print(f"{'Other':<30s} {other_genome_count:>12,} {other_size_count:>12,} "
      f"{other_otu_count:>12,} {other_species_count:>12,} {other_color:>10s}")
print("-" * 80)
print()
print(f"✅ {PROCESS_DOMAIN.upper()} Sankey diagram with colored flows completed!")
print(f"   📁 HTML file: {output_file}")
print("   💡 Open in web browser to view interactive diagram")
print("   🎨 Flows are colored to match their phylum colors")
print("=" * 70)

