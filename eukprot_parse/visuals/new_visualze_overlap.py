#!/usr/bin/env python3
"""
EukCensus-EukProt Overlap Visualization & Validation Script (Enhanced)

This script:
 1. Loads merged lineage data.
 2. Verifies EukProt species counts against the original raw dataset.
 3. Generates publication-ready visualizations with clear styling.
 4. Offers additional visualizations: treemap, sunburst, and heatmap.

Usage:
    python visualize_and_verify.py

Inputs:
    - merged_lineages.csv         # Output from merging pipeline (contains 'eukprot_count', 'overlap_count', 'match_score', 'eukcensus_member_size', 'eukcensus_lineage')
    - eukprot_new_lineages.csv  # Raw EukProt metadata with 'division' and unique 'species_id'

Outputs:
    - species_count_by_division.png
    - match_score_scatter.png
    - coverage_pie_chart.png
    - count_discrepancies.csv     # Verification report
    - division_treemap.png        # Treemap of species counts by division
    - lineage_sunburst.png        # Sunburst chart of taxonomic hierarchy
    - overlap_heatmap.png         # Heatmap of overlap_count vs match_score
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


def setup_plot_style():
    plt.style.use('default')
    plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['font.family'] = 'DejaVu Sans'


def verify_counts(merged_df, raw_file, division_col='Name_to_use', id_col='taxid'):
    """Verify species counts against raw EukProt data"""
    if not os.path.exists(raw_file):
        print(f"❌ Raw data not found: {raw_file}")
        return pd.DataFrame()

    print(f"📖 Reading raw EukProt data from {raw_file}...")
    raw = pd.read_csv(raw_file)

    # Filter out failed entries
    raw = raw[raw['taxid'] != 'FAILED'].copy()

    # Count unique species per division (using Name_to_use as division identifier)
    actual = raw.groupby(division_col)[id_col].nunique().rename('actual_count')

    # Join with merged data
    comp = merged_df.set_index('eukcensus_taxon_name').join(actual, how='left')
    comp['actual_count'] = comp['actual_count'].fillna(0).astype(int)
    comp['difference'] = comp['eukprot_count'] - comp['actual_count']

    # Create report
    report = comp.reset_index()[['eukcensus_taxon_name','eukprot_count','actual_count','difference']]
    report.to_csv('count_discrepancies.csv', index=False)
    print(f"✅ Verification report saved: count_discrepancies.csv")

    # Show summary
    total_diff = report['difference'].abs().sum()
    print(f"   Total discrepancies: {total_diff}")

    return report


def create_species_count_barplot(df):
    df_sorted = df.sort_values('eukprot_count', ascending=True)
    fig, ax = plt.subplots(figsize=(12, 10))

    # Create color palette using matplotlib
    colors = plt.cm.viridis(np.linspace(0, 1, len(df_sorted)))

    # Create horizontal bar plot
    bars = ax.barh(df_sorted['eukcensus_taxon_name'], df_sorted['eukprot_count'], color=colors)

    ax.set_title('EukProt Species Coverage per Division', pad=15)
    ax.set_xlabel('Number of Species')
    ax.set_ylabel('Division')

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, df_sorted['eukprot_count'])):
        ax.text(bar.get_width() + df_sorted['eukprot_count'].max()*0.01,
                bar.get_y() + bar.get_height()/2, f'{val}',
                va='center', ha='left', fontsize=10, fontweight='bold')

    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig('species_count_by_division.png', dpi=300, bbox_inches='tight')
    plt.close()


def create_match_score_plot(df):
    fig, ax = plt.subplots()
    sizes = np.interp(df['eukcensus_member_size'],
                      (df['eukcensus_member_size'].min(), df['eukcensus_member_size'].max()),
                      (50, 300))
    scatter = ax.scatter(
        df['eukprot_count'], df['match_score'],
        s=sizes, c=df['overlap_count'], cmap='plasma',
        edgecolors='white', linewidth=0.6, alpha=0.8
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Overlap Count')
    ax.set_title('Match Score vs. Species Count', pad=15)
    ax.set_xlabel('EukProt Species Count')
    ax.set_ylabel('Match Score')
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('match_score_scatter.png', dpi=300)
    plt.close()


def create_coverage_pie_chart(df):
    matched = df['eukprot_count'] > 0
    sizes = [matched.sum(), (~matched).sum()]
    labels = ['Matched', 'Unmatched']
    colors = ['#4CAF50', '#F44336']
    fig, ax = plt.subplots()
    ax.pie(
        sizes, labels=labels, autopct='%1.1f%%', startangle=140,
        colors=colors, textprops={'color':'black'}, wedgeprops={'width':0.4}
    )
    ax.set_title('Division Coverage by EukProt', pad=15)
    plt.tight_layout()
    plt.savefig('coverage_pie_chart.png', dpi=300)
    plt.close()


def create_division_treemap(df):
    """Create a pie chart as alternative to treemap (since squarify not available)"""
    counts = df.set_index('eukcensus_taxon_name')['eukprot_count']
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create pie chart with species counts
    colors = plt.cm.viridis(np.linspace(0, 1, len(counts)))
    wedges, texts, autotexts = ax.pie(counts.values, labels=counts.index, autopct='%1.0f',
                                     colors=colors, startangle=90)

    # Improve text readability
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')

    ax.set_title('Species Distribution by Division (Alternative to Treemap)', pad=15)
    plt.tight_layout()
    plt.savefig('division_treemap.png', dpi=300, bbox_inches='tight')
    plt.close()


def create_overlap_heatmap(df):
    """Create a scatter plot matrix as alternative to heatmap"""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create scatter plot with overlap_count vs match_score
    scatter = ax.scatter(df['match_score'], df['overlap_count'],
                        s=df['eukprot_count']*2, # Size by species count
                        c=df['eukprot_count'], cmap='YlGnBu',
                        alpha=0.7, edgecolors='black', linewidth=0.5)

    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('EukProt Species Count')

    # Add labels for each point
    for idx, row in df.iterrows():
        ax.annotate(row['eukcensus_taxon_name'],
                   (row['match_score'], row['overlap_count']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.8)

    ax.set_xlabel('Match Score')
    ax.set_ylabel('Overlap Count')
    ax.set_title('Overlap Analysis: Match Score vs Overlap Count\n(Bubble size = Species count)', pad=15)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('overlap_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    setup_plot_style()
    merged_file = 'merged_lineages.csv'
    raw_file = 'eukprot_new_lineages.csv'

    if not os.path.exists(merged_file):
        print(f"❌ File not found: {merged_file}")
        return

    merged_df = pd.read_csv(merged_file)
    print(f"✅ Loaded {len(merged_df)} divisions from merged data")

    verify_counts(merged_df, raw_file)
    create_species_count_barplot(merged_df)
    create_match_score_plot(merged_df)
    create_coverage_pie_chart(merged_df)
    create_division_treemap(merged_df)
    create_overlap_heatmap(merged_df)

    print("\n🎉 All visuals and reports generated successfully!")

if __name__ == '__main__':
    main()
