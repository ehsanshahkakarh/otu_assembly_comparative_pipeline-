#!/usr/bin/env python3
"""
EukCensus-EukProt Overlap Visualization Script

This script creates visualizations to represent the taxonomic overlap
between EukCensus divisions and EukProt species data.

Usage:
    python visualize_overlap.py

Input:
    - merged_lineages.csv
    - species_matches_summary.csv

Output:
    - Multiple visualization files (PNG format)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def setup_plot_style():
    """Set up consistent plot styling."""
    plt.style.use('default')
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12

def create_species_count_barplot(df):
    """Create a horizontal bar plot showing EukProt species count per division."""
    print("📊 Creating species count bar plot...")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Sort by species count for better visualization
    df_sorted = df.sort_values('eukprot_count', ascending=True)

    # Create horizontal bar plot
    bars = ax.barh(df_sorted['eukcensus_taxon_name'], df_sorted['eukprot_count'],
                   color=plt.cm.viridis(np.linspace(0, 1, len(df_sorted))))

    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, df_sorted['eukprot_count'])):
        ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2,
                f'{count}', ha='left', va='center', fontweight='bold')

    ax.set_xlabel('Number of EukProt Species')
    ax.set_ylabel('EukCensus Division')
    ax.set_title('EukProt Species Coverage by EukCensus Division')
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig('species_count_by_division.png', dpi=300, bbox_inches='tight')
    plt.close()


def create_taxonomic_hierarchy_plot(df):
    """Create a hierarchical visualization showing taxonomic relationships."""
    print("📊 Creating taxonomic hierarchy plot...")

    fig, ax = plt.subplots(figsize=(14, 10))

    # Parse lineages to get supergroups
    supergroups = {}
    for idx, row in df.iterrows():
        lineage = row['eukcensus_lineage'].split(';')
        if len(lineage) >= 3:
            supergroup = lineage[2] if lineage[2] != 'Eukaryota' else (lineage[3] if len(lineage) > 3 else 'Other')
        else:
            supergroup = 'Other'

        if supergroup not in supergroups:
            supergroups[supergroup] = []
        supergroups[supergroup].append({
            'division': row['eukcensus_taxon_name'],
            'species_count': row['eukprot_count'],
            'member_size': row['eukcensus_member_size']
        })

    # Create nested bar plot
    y_pos = 0
    colors = plt.cm.Set3(np.linspace(0, 1, len(supergroups)))

    for i, (supergroup, divisions) in enumerate(supergroups.items()):
        # Plot supergroup header
        ax.barh(y_pos, max([d['species_count'] for d in divisions]),
                height=0.3, color=colors[i], alpha=0.3,
                label=f'{supergroup} (supergroup)')
        ax.text(5, y_pos, supergroup, fontweight='bold', fontsize=12, va='center')
        y_pos += 0.5

        # Plot divisions within supergroup
        for division in sorted(divisions, key=lambda x: x['species_count'], reverse=True):
            ax.barh(y_pos, division['species_count'], height=0.6,
                   color=colors[i], alpha=0.8)
            ax.text(division['species_count'] + 5, y_pos,
                   f"{division['division']} ({division['species_count']} spp.)",
                   va='center', fontsize=10)
            y_pos += 0.8

        y_pos += 0.5  # Space between supergroups

    ax.set_xlabel('Number of EukProt Species')
    ax.set_title('Taxonomic Hierarchy: EukProt Species Distribution')
    ax.set_ylim(-0.5, y_pos)
    ax.set_yticks([])
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig('taxonomic_hierarchy.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_coverage_donut_chart(df):
    """Create a nested donut chart showing EukCensus-EukProt overlap with atlas for small groups."""
    print("📊 Creating coverage donut chart...")

    fig = plt.figure(figsize=(20, 12))

    # Create custom layout: main donut chart and legend area
    ax_donut = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    ax_legend = plt.subplot2grid((2, 3), (0, 2), colspan=1, rowspan=2)

    # Calculate coverage statistics
    total_divisions = len(df)
    matched_divisions = len(df[df['eukprot_count'] > 0])
    unmatched_divisions = total_divisions - matched_divisions

    # Calculate percentages
    pct_matched = (matched_divisions / total_divisions) * 100
    pct_unmatched = (unmatched_divisions / total_divisions) * 100

    # Get species distribution by major taxonomic groups
    major_groups = {}
    for idx, row in df.iterrows():
        lineage = row['eukcensus_lineage'].split(';')
        if len(lineage) >= 3:
            group = lineage[2] if lineage[2] != 'Eukaryota' else (lineage[3] if len(lineage) > 3 else 'Other')
        else:
            group = 'Other'

        if group not in major_groups:
            major_groups[group] = 0
        major_groups[group] += row['eukprot_count']

    # Sort groups by species count
    sorted_groups = sorted(major_groups.items(), key=lambda x: x[1], reverse=True)
    total_species = sum(major_groups.values())

    # Separate large and small groups (threshold: 5% of total for cleaner visualization)
    threshold_percent = 5.0
    large_groups = []
    small_groups = []

    for group, count in sorted_groups:
        percentage = (count / total_species) * 100
        if percentage >= threshold_percent:
            large_groups.append((group, count, percentage))
        else:
            small_groups.append((group, count, percentage))

    # Combine small groups into "Others" category
    if small_groups:
        others_count = sum([count for _, count, _ in small_groups])
        others_percentage = (others_count / total_species) * 100
        large_groups.append(("Others", others_count, others_percentage))

    # Outer ring: Division coverage (EukCensus perspective)
    division_sizes = [matched_divisions, unmatched_divisions]
    division_labels = [f'Matched Divisions\n{matched_divisions} divisions\n({pct_matched:.1f}%)',
                      f'Unmatched Divisions\n{unmatched_divisions} divisions\n({pct_unmatched:.1f}%)']
    division_colors = ['#2ecc71', '#e74c3c']

    ax_donut.pie(division_sizes, labels=division_labels, colors=division_colors,
                autopct='', startangle=90, wedgeprops=dict(width=0.3, edgecolor='w'))

    # Inner ring: Species distribution by taxonomic groups
    species_sizes = [count for _, count, _ in large_groups]
    species_labels = []

    # Create cleaner labels for inner ring
    for group, count, percentage in large_groups:
        if group == "Others":
            species_labels.append('')  # No label for Others in inner ring
        else:
            species_labels.append(f'{group}\n({percentage:.1f}%)')

    # Use distinct colors for taxonomic groups
    species_colors = plt.cm.Set3(np.linspace(0, 1, len(large_groups)))

    ax_donut.pie(species_sizes, labels=species_labels, colors=species_colors,
                autopct='', startangle=90, wedgeprops=dict(width=0.3, edgecolor='w'), radius=0.7)

    # Add center circle for donut effect
    centre_circle = plt.Circle((0, 0), 0.4, fc='white')
    ax_donut.add_patch(centre_circle)

    # Add title
    ax_donut.set_title('EukCensus-EukProt Overlap Analysis\nOuter: Division Coverage | Inner: Species Distribution',
                      fontsize=14, fontweight='bold')

    # Create legend/atlas for taxonomic groups
    ax_legend.axis('off')
    ax_legend.set_title('Taxonomic Groups Atlas', fontsize=12, fontweight='bold', pad=20)

    # Display major groups information
    legend_text = []
    legend_text.append("Major Taxonomic Groups:\n")

    for group, count, percentage in large_groups:
        if group != "Others":
            legend_text.append(f"• {group}")
            legend_text.append(f"  {count:,} species ({percentage:.1f}%)")
            legend_text.append("")  # Add spacing

    # Add small groups information if any
    if small_groups:
        legend_text.append("Small Groups (in 'Others'):")
        legend_text.append("")

        for group, count, percentage in small_groups:
            legend_text.append(f"• {group}: {count} species ({percentage:.1f}%)")

        others_count = sum([count for _, count, _ in small_groups])
        others_percentage = (others_count / total_species) * 100
        legend_text.append("")
        legend_text.append(f"Total 'Others': {others_count:,} species ({others_percentage:.1f}%)")

    # Display legend text
    legend_str = '\n'.join(legend_text)
    ax_legend.text(0.05, 0.95, legend_str, transform=ax_legend.transAxes,
                  fontsize=10, verticalalignment='top',
                  bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))

    # Add overall statistics
    stats_text = (f"📊 Summary Statistics:\n\n"
                 f"Total Divisions: {total_divisions}\n"
                 f"Matched Divisions: {matched_divisions} ({pct_matched:.1f}%)\n"
                 f"Total Species: {total_species:,}\n"
                 f"Taxonomic Groups: {len(sorted_groups)}")

    ax_legend.text(0.05, 0.25, stats_text, transform=ax_legend.transAxes,
                  fontsize=11, fontweight='bold',
                  bbox=dict(boxstyle="round,pad=0.4", facecolor="lightblue", alpha=0.7))

    plt.tight_layout()
    plt.savefig('coverage_overview.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_stats_table(df):
    """Create a summary statistics visualization."""
    print("📊 Creating summary statistics table...")

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('tight')
    ax.axis('off')

    # Calculate summary statistics
    stats = {
        'Total EukCensus Divisions': len(df),
        'Divisions with EukProt Matches': len(df[df['eukprot_count'] > 0]),
        'Total EukProt Species': df['eukprot_count'].sum(),
        'Average Species per Division': f"{df['eukprot_count'].mean():.1f}",
        'Max Species in Single Division': df['eukprot_count'].max(),
        'Min Species in Single Division': df['eukprot_count'].min(),
        'Average Match Score': f"{df['match_score'].mean():.1f}",
        'Perfect Matches (Score ≥90)': len(df[df['match_score'] >= 90]),
        'Good Matches (Score 70-89)': len(df[(df['match_score'] >= 70) & (df['match_score'] < 90)]),
        'Moderate Matches (Score <70)': len(df[df['match_score'] < 70])
    }

    # Create table
    table_data = [[k, v] for k, v in stats.items()]
    table = ax.table(cellText=table_data,
                    colLabels=['Metric', 'Value'],
                    cellLoc='left',
                    loc='center',
                    colWidths=[0.7, 0.3])

    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 2)

    # Style the table
    for i in range(len(table_data) + 1):
        for j in range(2):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#3498db')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#ecf0f1' if i % 2 == 0 else 'white')

    ax.set_title('EukCensus-EukProt Overlap Summary Statistics',
                fontsize=16, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig('summary_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("🎨 Starting EukCensus-EukProt overlap visualization...")

    # Setup
    setup_plot_style()

    # Check input files
    merged_file = "merged_lineages.csv"
    if not os.path.exists(merged_file):
        print(f"❌ Error: {merged_file} not found")
        return

    # Read data
    print(f"📖 Reading {merged_file}...")
    df = pd.read_csv(merged_file)
    print(f"   Loaded {len(df)} division entries")

    # Create visualizations
    create_species_count_barplot(df)
    create_taxonomic_hierarchy_plot(df)
    create_coverage_donut_chart(df)
    create_summary_stats_table(df)

    print("\n📁 Visualization files created:")
    print("   species_count_by_division.png")
    print("   taxonomic_hierarchy.png")
    print("   coverage_overview.png")
    print("   summary_statistics.png")
    print("✨ Visualization completed successfully!")

if __name__ == "__main__":
    main()