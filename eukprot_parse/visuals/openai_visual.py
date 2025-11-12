#!/usr/bin/env python3
"""
EukCensus-EukProt Overlap Visualization Script (Refined)

This updated script enhances the visual presentation of the taxonomic overlap
between EukCensus divisions and EukProt species data.

Usage:
    python visualize_overlap.py

Input:
    - merged_lineages.csv

Output:
    - High-resolution and publication-ready visualization files (PNG format)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


def setup_plot_style():
    sns.set_context("notebook")
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['font.family'] = 'DejaVu Sans'


def create_species_barplot(df):
    df_sorted = df.sort_values('eukprot_count', ascending=True)

    fig, ax = plt.subplots()
    palette = sns.color_palette("viridis", len(df_sorted))
    
    sns.barplot(
        y='eukcensus_taxon_name', x='eukprot_count', data=df_sorted,
        palette=palette, ax=ax
    )

    ax.set_title('EukProt Species Coverage per EukCensus Division')
    ax.set_xlabel('Number of EukProt Species')
    ax.set_ylabel('EukCensus Division')

    for i, val in enumerate(df_sorted['eukprot_count']):
        ax.text(val + 1, i, f'{val}', va='center')

    plt.tight_layout()
    plt.savefig('species_count_by_division.png', dpi=300)
    plt.close()


def create_match_score_plot(df):
    fig, ax = plt.subplots()
    scatter = ax.scatter(
        df['eukprot_count'], df['match_score'],
        s=np.sqrt(df['eukcensus_member_size'])*2,
        c=df['overlap_count'], cmap='plasma', edgecolor='k', alpha=0.75
    )

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Overlap Count')

    ax.set_title('Match Score vs EukProt Species Count')
    ax.set_xlabel('EukProt Species Count')
    ax.set_ylabel('Match Score')

    plt.tight_layout()
    plt.savefig('match_score_scatter.png', dpi=300)
    plt.close()


def create_coverage_pie_chart(df):
    matched = df[df['eukprot_count'] > 0]
    unmatched = df[df['eukprot_count'] == 0]
    
    labels = ['Matched Divisions', 'Unmatched Divisions']
    sizes = [len(matched), len(unmatched)]
    colors = ['#4CAF50', '#F44336']

    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, autopct='%1.1f%%', startangle=140,
        colors=colors, textprops=dict(color="w"), wedgeprops=dict(width=0.5)
    )

    for text in texts:
        text.set_color('black')

    ax.set_title('EukCensus Division Coverage by EukProt')
    plt.tight_layout()
    plt.savefig('coverage_pie_chart.png', dpi=300)
    plt.close()


def main():
    setup_plot_style()
    filepath = "merged_lineages.csv"

    if not os.path.exists(filepath):
        print(f"❌ File not found: {filepath}")
        return

    df = pd.read_csv(filepath)
    print(f"✅ Loaded {len(df)} entries from {filepath}")

    create_species_barplot(df)
    create_match_score_plot(df)
    create_coverage_pie_chart(df)

    print("\n🎉 Visualizations created:")
    print("- species_count_by_division.png")
    print("- match_score_scatter.png")
    print("- coverage_pie_chart.png")


if __name__ == "__main__":
    main()