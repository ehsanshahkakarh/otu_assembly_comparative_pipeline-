#!/usr/bin/env python3
"""
16S Census + GTDB Visualization Script
======================================

Creates three key visualizations comparing 16S census data with GTDB genomic data:
1. Environment vs Genomic Comparison - Side-by-side bar chart
2. Coverage Bar Chart - Horizontal bars showing genomic coverage of OTUs
3. OTU Coverage Chart - Detailed coverage breakdown with color coding

Usage:
    python create_16s_gtdb_visualizations.py

Input:
    - merged_output/16s_gtdb_merged_phylum.csv
    - merged_output/16s_gtdb_merged_family.csv
    - merged_output/16s_gtdb_merged_genus.csv

Output:
    - visualizations/16s_gtdb_environment_vs_genomic_comparison_{level}.png
    - visualizations/16s_gtdb_coverage_bar_chart_{level}.png
    - visualizations/16s_gtdb_otu_coverage_chart_{level}.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality figures
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# Configuration
FIGURE_SIZE = (12, 8)
DPI = 300
FONT_SIZE = 12
TITLE_SIZE = 16
LABEL_SIZE = 14

# Color schemes
COVERAGE_COLORS = {
    'no_coverage': '#1b5e20',      # Dark Green - No coverage (0%)
    'low_coverage': '#ff9800',     # Orange - Low coverage (<5%)
    'moderate_coverage': '#ffc107', # Amber - Moderate coverage (5-20%)
    'good_coverage': '#4caf50',     # Green - Good coverage (20-50%)
    'high_coverage': '#2196f3'     # Blue - High coverage (>50%)
}

def setup_output_directory():
    """Create output directory for visualizations."""
    output_dir = Path('visualizations')
    output_dir.mkdir(exist_ok=True)
    return output_dir

def load_merged_data():
    """Load merged 16S + GTDB data."""
    print("Loading merged 16S + GTDB data...")
    
    data = {}
    levels = ['phylum', 'family', 'genus']
    
    for level in levels:
        file_path = Path(f'merged_output/16s_gtdb_merged_{level}.csv')
        if file_path.exists():
            df = pd.read_csv(file_path)
            # Filter to only include entries present in both datasets
            df_both = df[df['in_both'] == 1].copy()
            data[level] = df_both
            print(f"  {level.capitalize()}: {len(df_both)} entries with both datasets")
        else:
            print(f"  Warning: {file_path} not found")
    
    return data

def get_coverage_color(coverage):
    """Get color based on coverage percentage."""
    if coverage == 0:
        return COVERAGE_COLORS['no_coverage']
    elif coverage < 5:
        return COVERAGE_COLORS['low_coverage']
    elif coverage < 20:
        return COVERAGE_COLORS['moderate_coverage']
    elif coverage < 50:
        return COVERAGE_COLORS['good_coverage']
    elif coverage < 100:
        return COVERAGE_COLORS['high_coverage']
    else:
        return COVERAGE_COLORS['high_coverage']  # Very high coverage (≥100%)

def create_environment_vs_genomic_comparison(df, level, output_dir):
    """Create overlapping bar chart comparing environmental occurrence vs genomic representation."""
    print(f"Creating environment vs genomic comparison chart for {level} level...")

    fig, ax = plt.subplots(figsize=(16, 10))

    # Sort by census occurrence count for better visualization
    df_sorted = df.sort_values('census_occurrence_count', ascending=False)

    # Take top 15 for readability
    df_top = df_sorted.head(15)

    # Create positions for bars
    x_pos = range(len(df_top))
    width = 0.35

    # Create bars
    bars1 = ax.bar([x - width/2 for x in x_pos], df_top['census_occurrence_count'],
                   width, label='Environmental Occurrence (16S Census OTUs)',
                   color='lightgreen', alpha=0.8, edgecolor='darkgreen', linewidth=0.5)

    bars2 = ax.bar([x + width/2 for x in x_pos], df_top['gtdb_species_count'],
                   width, label='Genomic Representation (GTDB Species)',
                   color='darkgreen', alpha=0.8, edgecolor='black', linewidth=0.5)

    # Customize chart
    ax.set_xlabel(f'Prokaryotic {level.capitalize()}s (Top 15 by Environmental OTU Count)',
                  fontsize=LABEL_SIZE, weight='bold')
    ax.set_ylabel('Count', fontsize=LABEL_SIZE, weight='bold')
    ax.set_title(f'Environmental vs Genomic Representation: {level.capitalize()} Level\n' +
                f'Comparing 16S Census OTU Count with GTDB Species Count',
                fontsize=TITLE_SIZE, weight='bold')

    # Set x-axis labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df_top['taxon_name'], rotation=45, ha='right', fontsize=14)

    # Add legend
    ax.legend(loc='upper right', fontsize=12)

    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='y')

    # Calculate y-axis limit to reduce white space
    all_values = list(df_top['census_occurrence_count']) + list(df_top['gtdb_species_count'])

    # Set specific y-axis limit for phylum level, use dynamic scaling for others
    if level == 'phylum':
        y_limit = 40000 # Fixed limit for phylum level
    else:
        q75 = np.percentile(all_values, 75)
        y_limit = max(q75 * 1.5, max(all_values) * 0.4)

    # Set y-axis limit
    ax.set_ylim(0, y_limit)

    # Add value labels on bars with special handling for tall bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:  # Only label non-zero bars
                bar_center = bar.get_x() + bar.get_width()/2.
                if height <= y_limit:
                    # Normal label for bars within limit
                    ax.text(bar_center, height + y_limit * 0.01,
                           f'{int(height)}',
                           ha='center', va='bottom', fontsize=9, weight='bold')
                else:
                    # Yellow highlighted label for tall bars that exceed limit
                    ax.text(bar_center, y_limit * 0.95,
                           f'{int(height)}',
                           ha='center', va='center', fontsize=9, weight='bold',
                           bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.8))

    # Add note about clipped bars if any exist
    clipped_bars = [height for height in all_values if height > y_limit]
    if clipped_bars:
        ax.text(0.02, 0.98, f'Note: {len(clipped_bars)} bar(s) exceed y-axis limit',
                transform=ax.transAxes, fontsize=10, va='top', ha='left',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_dir / f'16s_gtdb_environment_vs_genomic_comparison_{level}.png', dpi=DPI, bbox_inches='tight')
    plt.close()

def create_coverage_bar_chart(df, level, output_dir):
    """Create horizontal bar chart showing genomic coverage of OTUs with optimized y-axis scaling and outlier handling."""
    print(f"Creating coverage bar chart for {level} level...")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Calculate coverage percentage
    df['coverage_pct'] = (df['gtdb_species_count'] / df['census_occurrence_count'] * 100).fillna(0)

    # Sort by coverage for better visualization
    df_sorted = df.sort_values('coverage_pct', ascending=True)

    # Take top 20 for readability
    df_top = df_sorted.tail(20)

    # Set fixed axis limit for consistency with NCBI charts
    axis_limit = 2000

    # Determine which entries are outliers
    df_top['clipped'] = df_top['coverage_pct'] > axis_limit
    display_values = np.minimum(df_top['coverage_pct'], axis_limit)

    # Assign colors
    colors = [get_coverage_color(cov) for cov in df_top['coverage_pct']]

    # Plot horizontal bars
    bars = ax.barh(df_top['taxon_name'], display_values, color=colors, alpha=0.8, edgecolor='black')

    # Annotate bars
    for i, (bar, actual_value, taxon) in enumerate(zip(bars, df_top['coverage_pct'], df_top['taxon_name'])):
        bar_center = bar.get_y() + bar.get_height() / 2
        if actual_value <= axis_limit:
            ax.text(actual_value + axis_limit * 0.01, bar_center,
                    f'{actual_value:.1f}%', va='center', fontsize=10, weight='bold')
        else:
            ax.text(axis_limit * 0.95, bar_center,
                    f'{actual_value:.1f}%', va='center', fontsize=10, weight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.7))

    # Label axes and title
    ax.set_xlim(0, axis_limit)
    ax.set_xlabel('Genomic Coverage of OTUs (%)', fontsize=LABEL_SIZE, weight='bold')
    ax.set_ylabel(f'Prokaryotic {level.capitalize()}', fontsize=LABEL_SIZE, weight='bold')
    ax.set_title(f'GTDB Genomic Coverage of 16S Census OTUs: {level.capitalize()} Level\n'
                 f'Percentage Coverage = (GTDB Species / 16S OTUs) × 100',
                 fontsize=TITLE_SIZE, weight='bold')

    # Create coverage legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['no_coverage'], label='No Coverage (0%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['low_coverage'], label='Low Coverage (<5%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['moderate_coverage'], label='Moderate Coverage (5-20%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['good_coverage'], label='Good Coverage (20-50%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['high_coverage'], label='High Coverage (≥50%)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

    # Add outlier notice
    outliers = df_top[df_top['clipped']]
    if not outliers.empty:
        ax.text(0.02, 0.98, f'Note: {len(outliers)} outlier(s) clipped for resolution',
                transform=ax.transAxes, fontsize=10, va='top', ha='left',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_dir / f'16s_gtdb_coverage_bar_chart_{level}.png', dpi=DPI, bbox_inches='tight')
    plt.close()

def create_otu_coverage_chart(df, level, output_dir):
    """Create vertical bar chart showing OTU coverage with optimized scaling and clipped outliers."""
    print(f"Creating OTU coverage chart for {level} level...")

    fig, ax = plt.subplots(figsize=(16, 10))

    df_sorted = df.sort_values('census_occurrence_count', ascending=False).head(20)
    df_sorted['coverage_pct'] = (df_sorted['gtdb_species_count'] / df_sorted['census_occurrence_count'] * 100).fillna(0)

    # More aggressive scaling to minimize white space - focus on bulk of data
    coverage_values = df_sorted['coverage_pct'].values
    q80 = np.percentile(coverage_values, 80)
    # Use 80th percentile + 20% padding, minimum 25% for context
    axis_limit = max(min(q80 * 1.2, max(coverage_values) * 0.6), 25)

    df_sorted['clipped'] = df_sorted['coverage_pct'] > axis_limit
    display_values = np.minimum(df_sorted['coverage_pct'], axis_limit)

    x_pos = range(len(df_sorted))
    colors = [get_coverage_color(c) for c in df_sorted['coverage_pct']]
    bars = ax.bar(x_pos, display_values, color=colors, edgecolor='black', alpha=0.8)

    for i, (bar, actual_value) in enumerate(zip(bars, df_sorted['coverage_pct'])):
        bar_center = bar.get_x() + bar.get_width() / 2
        if actual_value <= axis_limit:
            ax.text(bar_center, actual_value + axis_limit * 0.01,
                    f"{actual_value:.1f}%", ha='center', va='bottom', fontsize=10, weight='bold')
        else:
            clipped_height = axis_limit * 0.95
            ax.text(bar_center, clipped_height + axis_limit * 0.01,
                    f"{actual_value:.1f}%", ha='center', va='bottom', fontsize=10, weight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.8))

    ax.set_ylim(0, axis_limit)
    ax.set_xlabel(f'Prokaryotic {level.capitalize()}s (Top 20 by Environmental OTU Count)', fontsize=LABEL_SIZE, weight='bold')
    ax.set_ylabel('OTU Coverage (%)', fontsize=LABEL_SIZE, weight='bold')
    ax.set_title(f'Coverage = (GTDB Species / 16S OTUs) × 100\n'
                 f'GTDB Genomic Coverage of 16S Census OTUs: {level.capitalize()} Level\n'
                 f'Percentage of 16S Census OTUs with GTDB Genomic Representation', fontsize=TITLE_SIZE, weight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(df_sorted['taxon_name'], rotation=45, ha='right', fontsize=10)
    ax.grid(True, axis='y', linestyle='--', alpha=0.3)

    outliers = df_sorted[df_sorted['clipped']]
    if not outliers.empty:
        outlier_names = ', '.join(outliers['taxon_name'].head(3))
        if len(outliers) > 3:
            outlier_names += f' (+{len(outliers) - 3} more)'
        ax.text(0.02, 0.98, f'Note: {len(outliers)} outlier(s) clipped for resolution\nOutliers: {outlier_names}',
                transform=ax.transAxes, fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle="round,pad=0.4", facecolor='lightyellow', alpha=0.9))

    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['no_coverage'], label='No Coverage (0%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['low_coverage'], label='Low Coverage (<5%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['moderate_coverage'], label='Moderate Coverage (5-20%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['good_coverage'], label='Good Coverage (20-50%)'),
        plt.Rectangle((0,0),1,1, facecolor=COVERAGE_COLORS['high_coverage'], label='High Coverage (≥50%)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / f'16s_gtdb_otu_coverage_chart_{level}.png', dpi=DPI, bbox_inches='tight')
    plt.close()

def main():
    """Main function to generate all 16S + GTDB visualizations."""
    print("="*60)
    print("16S Census vs GTDB Visualization Generator")
    print("="*60)

    try:
        # Setup
        output_dir = setup_output_directory()
        print(f"Output directory: {output_dir}")

        # Load data
        data = load_merged_data()

        if not data:
            print("No merged data found. Please run 16s_gtdb_merger.py first.")
            return

        # Generate visualizations for each level
        for level, df in data.items():
            print(f"\nGenerating visualizations for {level} level...")
            print(f"  Data points: {len(df)}")

            create_environment_vs_genomic_comparison(df, level, output_dir)
            create_coverage_bar_chart(df, level, output_dir)
            create_otu_coverage_chart(df, level, output_dir)

        print("\n" + "="*60)
        print("VISUALIZATION GENERATION COMPLETE!")
        print("="*60)
        print(f"\nGenerated files in {output_dir}:")

        # List all generated files
        for png_file in sorted(output_dir.glob("16s_gtdb_*.png")):
            print(f"  - {png_file.name}")

        print(f"\nTotal 16S + GTDB files generated: {len(list(output_dir.glob('16s_gtdb_*.png')))}")

        print("\nVisualization types generated for each level:")
        print("  1. Environment vs Genomic comparison - 16S OTU count vs GTDB species count")
        print("  2. Coverage bar chart - Horizontal bars showing genomic coverage percentages")
        print("  3. OTU Coverage chart - Detailed coverage breakdown with color coding")

    except Exception as e:
        print(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()
