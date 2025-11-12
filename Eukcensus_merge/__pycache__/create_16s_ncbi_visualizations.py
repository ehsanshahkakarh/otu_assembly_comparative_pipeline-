#!/usr/bin/env python3
"""
16S Census + NCBI Visualization Generator
========================================

Creates comprehensive visualizations comparing 16S census data with NCBI genomic data.
Generates multiple chart types for different taxonomic levels.

Author: Enhanced EukCensus Parser Team
Date: 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style and parameters
plt.style.use('default')

# Font sizes
TITLE_SIZE = 16
LABEL_SIZE = 14
LEGEND_SIZE = 12
TICK_SIZE = 10

class SixteenSNCBIVisualizer:
    """
    Creates visualizations for 16S census + NCBI merged data.
    """
    
    def __init__(self, data_dir: str = "merged_output", output_dir: str = "visualizations"):
        """
        Initialize the visualizer.
        
        Args:
            data_dir: Directory containing merged CSV files
            output_dir: Directory for output visualizations
        """
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Taxonomic levels to process
        self.taxonomic_levels = ['phylum', 'family', 'genus']
        
    def load_merged_data(self, level: str) -> pd.DataFrame:
        """
        Load merged data for a specific taxonomic level.
        
        Args:
            level: Taxonomic level ('phylum', 'family', 'genus')
            
        Returns:
            Merged DataFrame
        """
        file_path = self.data_dir / f"16s_ncbi_merged_{level}.csv"
        if not file_path.exists():
            logger.warning(f"File not found: {file_path}")
            return pd.DataFrame()
        
        df = pd.read_csv(file_path)
        logger.info(f"Loaded {level} data: {len(df)} entries")
        return df
    
    def create_environment_vs_genomic_comparison(self, df: pd.DataFrame, level: str) -> None:
        """
        Create side-by-side bar chart comparing environmental vs genomic representation.
        
        Args:
            df: Merged DataFrame
            level: Taxonomic level
        """
        # Filter to entries present in both datasets for meaningful comparison
        df_both = df[df['in_both'] == 1].copy()
        
        if df_both.empty:
            logger.warning(f"No overlapping data for {level} level")
            return
        
        # Sort by NCBI species count to show top genomically represented taxa
        df_sorted = df.sort_values('ncbi_species_count', ascending=False)

        # Take top 10 for readability (focusing on NCBI representation)
        df_top = df_sorted.head(10)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Set up bar positions
        x_pos = np.arange(len(df_top))
        width = 0.35
        
        # Create bars for census occurrence count and NCBI species count
        bars1 = ax.bar([x - width/2 for x in x_pos], df_top['census_occurrence_count'], 
                       width, label='Environmental Representation (16S Census OTUs)',
                       color='darkblue', alpha=0.8, edgecolor='black', linewidth=0.5)
        bars2 = ax.bar([x + width/2 for x in x_pos], df_top['ncbi_species_count'],
                       width, label='Genomic Representation (NCBI Species)',
                       color='darkgreen', alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Customize chart
        ax.set_xlabel(f'Prokaryotic {level.capitalize()}s (Top 10 by NCBI Species Count)',
                      fontsize=LABEL_SIZE, weight='bold')
        ax.set_title(f'Environmental vs Genomic Representation: {level.capitalize()} Level\n' +
                    f'Comparing 16S Census OTU Count with NCBI Species Count',
                    fontsize=TITLE_SIZE, weight='bold')
        
        # Set x-axis labels
        ax.set_xticks(x_pos)
        ax.set_xticklabels(df_top['taxon_name'], rotation=45, ha='right', fontsize=14)
        
        # Calculate y-axis limit to reduce white space
        all_values = list(df_top['census_occurrence_count']) + list(df_top['ncbi_species_count'])
        
        # Set specific y-axis limit for phylum level, use dynamic scaling for others
        if level == 'phylum':
            y_limit = 40000  # Fixed limit for phylum level
        else:
            q75 = np.percentile(all_values, 75)
            y_limit = max(q75 * 1.5, max(all_values) * 0.4)
        
        # Set y-axis limit
        ax.set_ylim(0, y_limit)
        
        # Add value labels on bars with special handling for tall bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > y_limit * 0.9:  # If bar is very tall
                    # Place label inside the bar
                    ax.text(bar.get_x() + bar.get_width()/2., height * 0.85,
                           f'{int(height):,}', ha='center', va='center', 
                           fontsize=9, weight='bold', color='white')
                    # Add clipping indicator
                    ax.text(bar.get_x() + bar.get_width()/2., y_limit * 0.95,
                           '↑', ha='center', va='center', 
                           fontsize=12, weight='bold', color='red')
                else:
                    # Place label above the bar
                    ax.text(bar.get_x() + bar.get_width()/2., height + y_limit * 0.01,
                           f'{int(height):,}', ha='center', va='bottom', fontsize=9)
        
        # Customize y-axis
        ax.set_ylabel('Count', fontsize=LABEL_SIZE, weight='bold')
        ax.tick_params(axis='y', labelsize=TICK_SIZE)
        ax.tick_params(axis='x', labelsize=TICK_SIZE)
        
        # Add legend
        ax.legend(fontsize=LEGEND_SIZE, loc='upper right')
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3, axis='y')
        
        # Tight layout
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / f"16s_ncbi_environment_vs_genomic_comparison_{level}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved environment vs genomic comparison chart: {output_file}")
    
    def create_coverage_bar_chart(self, df: pd.DataFrame, level: str) -> None:
        """
        Create horizontal bar chart showing genomic coverage percentages.

        Args:
            df: Merged DataFrame
            level: Taxonomic level
        """
        # Filter to entries present in census (environmental data)
        df_census = df[df['in_census'] == 1].copy()

        if df_census.empty:
            logger.warning(f"No census data for {level} level")
            return

        # Calculate coverage percentage on census data only
        df_census['coverage_pct'] = (df_census['ncbi_species_count'] / df_census['census_occurrence_count'] * 100).fillna(0)

        # Filter out entries with zero coverage for better visualization
        df_with_coverage = df_census[df_census['coverage_pct'] > 0].copy()

        if df_with_coverage.empty:
            logger.warning(f"No entries with coverage > 0 for {level} level")
            # Create a chart showing the most abundant taxa even if they have 0 coverage
            df_sorted = df_census.sort_values('census_occurrence_count', ascending=False)
            df_top = df_sorted.head(15)  # Top 15 by environmental abundance
        else:
            # Sort by coverage percentage in descending order for better visualization
            df_sorted = df_with_coverage.sort_values('coverage_pct', ascending=False)
            df_top = df_sorted.head(15)  # Top 15 by coverage percentage

        # Create figure with larger size for better readability
        fig, ax = plt.subplots(figsize=(14, 10))

        # Create color coding based on coverage levels (matching GTDB thresholds)
        colors = []
        for pct in df_top['coverage_pct']:
            if pct > 50:
                colors.append('#2196f3')    # Blue - High coverage (>50%)
            elif pct >= 20:
                colors.append('#4caf50')    # Green - Good coverage (20-50%)
            elif pct >= 5:
                colors.append('#ffc107')    # Amber - Moderate coverage (5-20%)
            elif pct > 0:
                colors.append('#ff9800')    # Orange - Low coverage (<5%)
            else:
                colors.append('#1b5e20')    # Dark Green - No coverage (0%)

        # Create horizontal bars
        y_positions = range(len(df_top))
        bars = ax.barh(y_positions, df_top['coverage_pct'],
                      color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

        # Customize chart
        ax.set_yticks(y_positions)
        ax.set_yticklabels(df_top['taxon_name'], fontsize=TICK_SIZE)
        ax.set_xlabel('Coverage Percentage (%)', fontsize=LABEL_SIZE, weight='bold')
        ax.set_title(f'NCBI Genomic Coverage of 16S Census {level.capitalize()}s\n'
                     f'Coverage = (NCBI Species / 16S OTUs) × 100 | Sorted by Coverage %',
                     fontsize=TITLE_SIZE, weight='bold')

        # Set fixed x-axis limit for consistency with GTDB charts
        ax.set_xlim(0, 2000)

        # Add value labels on bars
        for i, (bar, coverage) in enumerate(zip(bars, df_top['coverage_pct'])):
            width = bar.get_width()
            if width > 0:
                # Position label based on bar width
                if width < max_coverage * 0.1:  # For very small bars
                    ax.text(width + max_coverage * 0.02, bar.get_y() + bar.get_height()/2,
                           f'{width:.1f}%', ha='left', va='center', fontsize=9)
                else:
                    ax.text(width + max_coverage * 0.01, bar.get_y() + bar.get_height()/2,
                           f'{width:.1f}%', ha='left', va='center', fontsize=9)
            else:
                # For zero coverage, show "0%" inside the chart area
                ax.text(max_coverage * 0.02, bar.get_y() + bar.get_height()/2,
                       '0.0%', ha='left', va='center', fontsize=9, color='gray')

        # Add legend for color coding (matching GTDB thresholds)
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#2196f3', label='High Coverage (>50%)'),
            plt.Rectangle((0,0),1,1, facecolor='#4caf50', label='Good Coverage (20-50%)'),
            plt.Rectangle((0,0),1,1, facecolor='#ffc107', label='Moderate Coverage (5-20%)'),
            plt.Rectangle((0,0),1,1, facecolor='#ff9800', label='Low Coverage (<5%)'),
            plt.Rectangle((0,0),1,1, facecolor='#1b5e20', label='No Coverage (0%)')
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=LEGEND_SIZE-1)

        # Add grid for better readability
        ax.grid(True, alpha=0.3, axis='x')

        # Invert y-axis to show highest coverage at top
        ax.invert_yaxis()

        # Tight layout
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / f"16s_ncbi_coverage_bar_chart_{level}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Saved coverage bar chart: {output_file}")
        logger.info(f"Chart shows {len(df_top)} {level}s with coverage range: {df_top['coverage_pct'].min():.1f}% - {df_top['coverage_pct'].max():.1f}%")
    
    def create_otu_coverage_chart(self, df: pd.DataFrame, level: str) -> None:
        """
        Create detailed OTU coverage chart with color coding.
        
        Args:
            df: Merged DataFrame
            level: Taxonomic level
        """
        # Filter to entries present in census
        df_census = df[df['in_census'] == 1].copy()
        
        if df_census.empty:
            logger.warning(f"No census data for {level} level")
            return
        
        # Sort by census occurrence count (environmental abundance)
        df_sorted = df_census.sort_values('census_occurrence_count', ascending=True)
        df_sorted['coverage_pct'] = (df_sorted['ncbi_species_count'] / df_sorted['census_occurrence_count'] * 100).fillna(0)
        
        # Take a reasonable number for visualization
        df_plot = df_sorted.tail(30)  # Top 30 by environmental abundance
        
        # Create figure
        fig, ax = plt.subplots(figsize=(14, 12))
        
        # Color coding based on coverage levels (matching GTDB thresholds)
        colors = []
        for pct in df_plot['coverage_pct']:
            if pct > 50:
                colors.append('#2196f3')    # Blue - High coverage (>50%)
            elif pct >= 20:
                colors.append('#4caf50')    # Green - Good coverage (20-50%)
            elif pct >= 5:
                colors.append('#ffc107')    # Amber - Moderate coverage (5-20%)
            elif pct > 0:
                colors.append('#ff9800')    # Orange - Low coverage (<5%)
            else:
                colors.append('#1b5e20')    # Dark Green - No coverage (0%)
        
        # Create horizontal bars
        bars = ax.barh(range(len(df_plot)), df_plot['coverage_pct'], 
                      color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Customize chart
        ax.set_yticks(range(len(df_plot)))
        ax.set_yticklabels(df_plot['taxon_name'], fontsize=TICK_SIZE)
        ax.set_xlabel('Coverage Percentage (%)', fontsize=LABEL_SIZE, weight='bold')
        ax.set_title(f'Coverage = (NCBI Species / 16S OTUs) × 100\n'
                     f'NCBI Genomic Coverage of 16S Census OTUs: {level.capitalize()} Level\n'
                     f'Percentage of 16S Census OTUs with NCBI Genomic Representation', fontsize=TITLE_SIZE, weight='bold')
        
        # Add value labels on bars
        for i, bar in enumerate(bars):
            width = bar.get_width()
            if width > 0:
                ax.text(width + 2, bar.get_y() + bar.get_height()/2,
                       f'{width:.1f}%', ha='left', va='center', fontsize=9)
        
        # Add legend for color coding (matching GTDB thresholds)
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#2196f3', label='High Coverage (>50%)'),
            plt.Rectangle((0,0),1,1, facecolor='#4caf50', label='Good Coverage (20-50%)'),
            plt.Rectangle((0,0),1,1, facecolor='#ffc107', label='Moderate Coverage (5-20%)'),
            plt.Rectangle((0,0),1,1, facecolor='#ff9800', label='Low Coverage (<5%)'),
            plt.Rectangle((0,0),1,1, facecolor='#1b5e20', label='No Coverage (0%)')
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=LEGEND_SIZE)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='x')
        
        # Tight layout
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / f"16s_ncbi_otu_coverage_chart_{level}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved OTU coverage chart: {output_file}")

    def generate_all_visualizations(self) -> None:
        """
        Generate all visualization types for all taxonomic levels.
        """
        logger.info("Starting 16S + NCBI visualization generation...")

        for level in self.taxonomic_levels:
            logger.info(f"Processing {level} level...")

            # Load data
            df = self.load_merged_data(level)
            if df.empty:
                logger.warning(f"No data available for {level} level")
                continue

            # Generate visualizations
            try:
                self.create_environment_vs_genomic_comparison(df, level)
                self.create_coverage_bar_chart(df, level)
                self.create_otu_coverage_chart(df, level)
                logger.info(f"Completed visualizations for {level} level")
            except Exception as e:
                logger.error(f"Error creating visualizations for {level}: {e}")

        logger.info("16S + NCBI visualization generation completed!")

        # Print summary
        print("\n16S + NCBI VISUALIZATION SUMMARY")
        print("=" * 35)
        print("Generated visualizations for each taxonomic level:")
        print("  1. Environment vs Genomic comparison - Top 10 by NCBI species count")
        print("  2. Coverage bar chart - Horizontal bars showing genomic coverage percentages")
        print("  3. OTU Coverage chart - Detailed coverage breakdown with color coding")
        print(f"\nAll visualizations saved to: {self.output_dir}")


def main():
    """
    Main function to generate 16S + NCBI visualizations.
    """
    import sys

    # Default paths
    data_dir = "merged_output"
    output_dir = "visualizations"

    # Allow command line arguments
    if len(sys.argv) >= 2:
        data_dir = sys.argv[1]
    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]

    # Create visualizer and generate plots
    visualizer = SixteenSNCBIVisualizer(data_dir, output_dir)
    visualizer.generate_all_visualizations()


if __name__ == "__main__":
    main()
