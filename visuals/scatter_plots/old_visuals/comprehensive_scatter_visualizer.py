 #!/usr/bin/env python3
"""
Comprehensive Scatter Plot Visualizer for Merged Taxonomic Data

This script creates high-quality scatter plots for all merged datasets:
- 16S and 18S data from NCBI and GTDB (prokaryote_eukcensus)
- EukProt 18S merger data
- Organized by taxonomic levels (phylum, family, genus)
- Separated by domains (Bacteria, Archaea, Eukaryota)

INPUT FILES:
===========
Prokaryote Census Data (16S/18S):
- prokaryote_eukcensus/merged_output/16s_ncbi_merged_phylum.csv
- prokaryote_eukcensus/merged_output/16s_ncbi_merged_family.csv
- prokaryote_eukcensus/merged_output/16s_ncbi_merged_genus.csv
- prokaryote_eukcensus/merged_output/16s_gtdb_merged_phylum.csv
- prokaryote_eukcensus/merged_output/16s_gtdb_merged_family.csv
- prokaryote_eukcensus/merged_output/16s_gtdb_merged_genus.csv

Eukaryote Data (18S):
- eukprot_18S_merger/results/division_coverage_summary.csv (phylum level)
- eukprot_18S_merger/results/division_analysis_summary.csv (backup)
- eukprot_18S_merger/results/family_coverage_summary.csv
- eukprot_18S_merger/results/genus_coverage_summary.csv 
- prokaryote_eukcensus/merged_output/18s_ncbi_merged_phylum.csv
- prokaryote_eukcensus/merged_output/18s_ncbi_merged_family.csv
- prokaryote_eukcensus/merged_output/18s_ncbi_merged_genus.csv

Expected Columns in Input Files:
- Census_OTU_Count: Environmental census OTU counts
- NCBI_Species_Count / GTDB_Species_Count: Genomic database species counts
- EukProt_Species_Count: EukProt species counts (for 18S data)
- Coverage_Percentage: Coverage percentage calculation
- taxon_name / taxon_name_clean: Taxonomic names
- in_both: Flag indicating presence in both datasets

OUTPUT FILES:
============
Directory Structure: new_visualizations/
├── 16S/
│   ├── ncbi_16s_bacteria_phylum_comparison.png
│   ├── ncbi_16s_bacteria_family_comparison.png
│   ├── ncbi_16s_bacteria_genus_comparison.png
│   ├── ncbi_16s_archaea_phylum_comparison.png
│   ├── ncbi_16s_archaea_family_comparison.png
│   ├── ncbi_16s_archaea_genus_comparison.png
│   ├── gtdb_16s_bacteria_phylum_comparison.png
│   ├── gtdb_16s_bacteria_family_comparison.png
│   ├── gtdb_16s_bacteria_genus_comparison.png
│   ├── gtdb_16s_archaea_phylum_comparison.png
│   ├── gtdb_16s_archaea_family_comparison.png
│   └── gtdb_16s_archaea_genus_comparison.png
├── 18S/
│   ├── ncbi_18s_eukaryota_phylum_comparison.png
│   ├── ncbi_18s_eukaryota_family_comparison.png
│   ├── ncbi_18s_eukaryota_genus_comparison.png
│   ├── eukprot_18s_eukaryota_phylum_comparison.png
│   ├── eukprot_18s_eukaryota_family_comparison.png
│   └── eukprot_18s_eukaryota_genus_comparison.png
├── by_domain/
│   ├── bacteria_subset_plots.png
│   ├── archaea_subset_plots.png
│   └── eukaryota_subset_plots.png
└── by_database/
    ├── ncbi_subset_plots.png
    ├── gtdb_subset_plots.png
    └── eukprot_subset_plots.png

Log Files:
- comprehensive_visualization.log

FEATURES:
========
- Professional matplotlib styling with configurable parameters
- Log-scale scatter plots with species richness color coding
- Novelty factor circle sizing (OTU count / species count ratio)
  * Larger circles indicate higher novelty (more environmental diversity than genomic representation)
  * Smaller circles indicate lower novelty (well-represented taxa in genomic databases)
- Top 15 taxa labeling with intelligent positioning
- Domain-specific color schemes (Bacteria: red, Archaea: purple, Eukaryota: green)
- Database-specific color schemes (NCBI: blue, GTDB: orange, EukProt: green)
- High-resolution output (300 DPI) suitable for publication
- Robust path handling for different working directories
- Comprehensive logging and error handling

Author: Generated for comprehensive taxonomic analysis
Date: 2025-01-23
Updated: 2025-06-30 - Implemented novelty factor circle sizing (Frederik Schulz suggestion)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from pathlib import Path
import warnings
import logging
import argparse
import json

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class VisualizationConfig:
    """Configuration class for consistent, professional-quality plotting."""

    def __init__(self):
        # Canvas settings
        self.figure_width = 12    # Adjusted to reduce excess width
        self.figure_height = 8    # Slightly more compact
        self.dpi = 300            # High resolution for publication

        # Font settings
        self.font_size = 16       # Base font for axis ticks, grid labels
        self.label_size = 20      # Axis labels
        self.title_size = 24      # Main title
        self.legend_size = 14     # Legend text
        self.annotation_size = 12 # Inline text/annotations

        # Scale options
        self.x_scale = 'symlog'      # Options: 'linear', 'log', 'symlog', 'logit'
        self.y_scale = 'symlog'

        # Plot element styling
        self.alpha = 0.8                  # Transparency for points/lines
        self.point_size = 100             # Suitable for scatter plots
        self.edge_width = 0.8             # Circle border width
        self.line_width = 2.5             # Line thickness for plots
        self.grid_alpha = 0.3             # Grid transparency

        # Color themes for databases
        self.database_colors = {
            'NCBI': '#1f77b4',    # Blue
            'GTDB': '#ff7f0e',    # Orange
            'EukProt': '#2ca02c'  # Green
        }

        # Color themes for domains
        self.domain_colors = {
            'Bacteria': '#e74c3c',     # Red
            'Archaea': '#9b59b6',      # Purple
            'Eukaryota': '#27ae60'     # Green
        }

        # EukProt representation levels
        self.representation_colors = {
            'No_Representation': '#d32f2f',       # Strong red
            'Low_Representation': '#ff9800',      # Orange
            'Moderate_Representation': '#ffc107', # Amber
            'Good_Representation': '#4caf50',     # Green
            'High_Representation': '#2196f3',     # Blue
            'No_Genomic_Data': '#9e9e9e',         # Gray
            'Under_Represented': '#f44336',       # Red
            'Proportional': '#4caf50',            # Green
            'Over_Represented': '#2196f3'         # Blue
        }

    def apply_matplotlib_style(self, plt):
        """Apply professional formatting to matplotlib plots."""
        import matplotlib.font_manager as fm

        # Professional font options (in order of preference)
        professional_fonts = [
            'DejaVu Sans',     # Clean, professional, widely available on Linux
            'Cantarell',       # Modern, clean (GNOME default)
            'Droid Sans',      # Clean, readable
            'Liberation Sans', # Open source alternative to Arial
            'Arial',           # Clean, professional (if available)
            'Helvetica',       # Classic, clean (if available)
            'sans-serif'       # System default fallback
        ]

        # Find the first available font
        available_fonts = [f.name for f in fm.fontManager.ttflist]
        selected_font = 'sans-serif'  # Default fallback

        for font in professional_fonts:
            if font in available_fonts or font == 'sans-serif':
                selected_font = font
                break

        plt.rcParams.update({
            'font.family': [selected_font],
            'font.weight': 'normal',
            'font.size': self.font_size,
            'figure.dpi': self.dpi,
            'axes.titlesize': self.title_size,
            'axes.labelsize': self.label_size,
            'axes.linewidth': 1.2,
            'xtick.labelsize': self.font_size,
            'ytick.labelsize': self.font_size,
            'legend.fontsize': self.legend_size,
            'legend.frameon': False,
            'grid.alpha': self.grid_alpha,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'text.color': '#2c3e50',      # Professional dark blue-gray
            'axes.labelcolor': '#2c3e50',
            'xtick.color': '#34495e',     # Slightly lighter for ticks
            'ytick.color': '#34495e',
            'axes.edgecolor': '#7f8c8d',  # Light gray for axes
        })

        logging.info(f"Using font: {selected_font}")

    def load_from_file(self, config_file):
        """Load configuration from JSON file."""
        try:
            with open(config_file, 'r') as f:
                config_data = json.load(f)

            # Update attributes from config file
            for key, value in config_data.items():
                if hasattr(self, key):
                    setattr(self, key, value)

            logging.info(f"Configuration loaded from {config_file}")
        except FileNotFoundError:
            logging.warning(f"Config file {config_file} not found, using defaults")
        except Exception as e:
            logging.error(f"Error loading config file: {e}")

    def save_to_file(self, config_file):
        """Save current configuration to JSON file."""
        config_data = {
            'figure_width': self.figure_width,
            'figure_height': self.figure_height,
            'dpi': self.dpi,
            'font_size': self.font_size,
            'label_size': self.label_size,
            'title_size': self.title_size,
            'legend_size': self.legend_size,
            'annotation_size': self.annotation_size,
            'x_scale': self.x_scale,
            'y_scale': self.y_scale,
            'alpha': self.alpha,
            'point_size': self.point_size,
            'edge_width': self.edge_width,
            'grid_alpha': self.grid_alpha,
            'line_width': self.line_width,
            'database_colors': self.database_colors,
            'domain_colors': self.domain_colors,
            'representation_colors': self.representation_colors
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        logging.info(f"Configuration saved to {config_file}")

    @property
    def figure_size(self):
        """Return figure size as tuple."""
        return (self.figure_width, self.figure_height)

# Global configuration instance
config = VisualizationConfig()

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('comprehensive_visualization.log'),
            logging.StreamHandler()
        ]
    )

def create_output_directory():
    """Create output directory for visualizations."""
    output_dir = Path('new_visualizations')
    output_dir.mkdir(exist_ok=True)

    # Create subdirectories for 16S and 18S
    for subdir in ['16S', '18S']:
        (output_dir / subdir).mkdir(exist_ok=True)

    return output_dir

def load_prokaryote_data():
    """Load prokaryote census merged data (16S and 18S, NCBI and GTDB)."""
    # Try multiple possible base paths
    possible_base_paths = [
        Path('../Eukcensus_merge/merged_output'),  # From visuals/ directory
        Path('Eukcensus_merge/merged_output'),     # From OTU_97_eukcensus_merger/ directory
        Path('../../Eukcensus_merge/merged_output'),
        Path('OTU_97_eukcensus_merger/Eukcensus_merge/merged_output'),
        Path('../OTU_97_eukcensus_merger/Eukcensus_merge/merged_output')
    ]

    base_path = None
    for path in possible_base_paths:
        if path.exists():
            base_path = path
            logging.info(f"Found prokaryote data directory: {base_path}")
            break

    if base_path is None:
        logging.warning("Could not find prokaryote data directory")
        return {}

    data = {}

    # 16S data
    for db in ['ncbi', 'gtdb']:
        for level in ['phylum', 'family', 'genus']:
            file_path = base_path / f'16s_{db}_merged_{level}.csv'
            if file_path.exists():
                df = pd.read_csv(file_path)
                data[f'16S_{db}_{level}'] = df
                logging.info(f"Loaded {file_path}: {len(df)} entries")
            else:
                logging.debug(f"File not found: {file_path}")

    # 18S data (only NCBI available)
    for level in ['phylum', 'family', 'genus']:
        file_path = base_path / f'18s_ncbi_merged_{level}.csv'
        if file_path.exists():
            df = pd.read_csv(file_path)
            data[f'18S_ncbi_{level}'] = df
            logging.info(f"Loaded {file_path}: {len(df)} entries")
        else:
            logging.debug(f"File not found: {file_path}")

    return data

def load_eukprot_data():
    """Load EukProt 18S merger data."""
    # Try multiple possible base paths - EukProt data is in the same merged_output directory
    possible_base_paths = [
        Path('../Eukcensus_merge/merged_output'),  # From visuals/ directory
        Path('Eukcensus_merge/merged_output'),     # From OTU_97_eukcensus_merger/ directory
        Path('../../Eukcensus_merge/merged_output'),
        Path('OTU_97_eukcensus_merger/Eukcensus_merge/merged_output'),
        Path('../OTU_97_eukcensus_merger/Eukcensus_merge/merged_output')
    ]

    data = {}

    # Find the base path
    base_path = None
    for path in possible_base_paths:
        if path.exists():
            base_path = path
            logging.info(f"Found EukProt data directory: {base_path}")
            break

    if base_path is None:
        logging.warning("Could not find EukProt data directory")
        return {}

    # Division (phylum) level - EukProt files are directly in merged_output
    division_file = base_path / '18s_eukprot_merged_division.csv'
    if division_file.exists():
        df = pd.read_csv(division_file)
        data['EukProt_18S_phylum'] = df
        logging.info(f"Loaded EukProt phylum/division data from {division_file.name}: {len(df)} entries")
    else:
        logging.warning("No EukProt phylum/division data found")

    # Family level
    family_file = base_path / '18s_eukprot_merged_family.csv'
    if family_file.exists():
        df = pd.read_csv(family_file)
        data['EukProt_18S_family'] = df
        logging.info(f"Loaded EukProt family data from {family_file.name}: {len(df)} entries")
    else:
        logging.debug(f"EukProt family file not found: {family_file}")

    # Genus level
    genus_file = base_path / '18s_eukprot_merged_genus.csv'
    if genus_file.exists():
        df = pd.read_csv(genus_file)
        data['EukProt_18S_genus'] = df
        logging.info(f"Loaded EukProt genus data from {genus_file.name}: {len(df)} entries")
    else:
        logging.debug(f"EukProt genus file not found: {genus_file}")

    return data

def create_subset_scatter_plots(df, dataset_name, level, output_dir):
    """Create simple scatter plots for each dataset comparison."""
    logging.info(f"Creating scatter plot for {dataset_name}")

    # Filter data with both census and database representation
    df_both = df[df['in_both'] == 1].copy()

    # Apply different Census_OTU_Count cutoffs based on taxonomic level
    if level == 'phylum':
        cutoff = 100  # Most restrictive for phylum
    elif level == 'family':
        cutoff = 10   # Moderate for family
    else:  # genus
        cutoff = 2    # Minimal cutoff for genus (filter out very low counts)

    # Apply the cutoff (all levels now have cutoffs)
    initial_count = len(df_both)
    df_both = df_both[df_both['Census_OTU_Count'] >= cutoff].copy()
    filtered_count = initial_count - len(df_both)
    if filtered_count > 0:
        logging.info(f"Filtered out {filtered_count} taxa with Census_OTU_Count < {cutoff} from {dataset_name} ({level} level)")
    else:
        logging.info(f"No taxa filtered from {dataset_name} ({level} level) - all had Census_OTU_Count >= {cutoff}")

    # Apply 100% coverage cutoff to filter out unrealistic coverage values
    if 'Coverage_Percentage' in df_both.columns:
        initial_coverage_count = len(df_both)
        df_both = df_both[df_both['Coverage_Percentage'] <= 100.0].copy()
        filtered_coverage_count = initial_coverage_count - len(df_both)
        if filtered_coverage_count > 0:
            logging.info(f"Filtered out {filtered_coverage_count} taxa with >100% coverage from {dataset_name}")

    if len(df_both) == 0:
        logging.warning(f"No meaningful data found for {dataset_name} after filtering")
        return

    # Determine database and gene from dataset name
    parts = dataset_name.split('_')
    gene = parts[0]  # 16S or 18S
    database = parts[1].upper()  # NCBI or GTDB

    # Determine the correct species count column name based on database
    if database == 'GTDB':
        species_count_col = 'GTDB_Species_Count'
    else:  # NCBI
        species_count_col = 'NCBI_Species_Count'

    # Check if the column exists
    if species_count_col not in df_both.columns:
        logging.error(f"Column '{species_count_col}' not found in {dataset_name}. Available columns: {list(df_both.columns)}")
        return

    # Create gene-specific subdirectory and filename
    gene_dir = output_dir / gene
    filename = f'{gene.lower()}_{database.lower()}_{level}.png'

    # Create single scatter plot
    create_single_scatter_plot(
        df_both, gene, database, level, species_count_col,
        gene_dir / filename,
        f'{gene} EukCensus vs {database}: {level.capitalize()} Level'
    )

def create_single_scatter_plot(df, gene, database, level, species_count_col, output_path, title):
    """Create a single scatter plot with specified parameters."""

    # Create figure with minimal white space
    _, ax = plt.subplots(figsize=(config.figure_width + 1, config.figure_height))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)  # Minimize white space

    # Use coverage percentage for summary statistics
    if 'Coverage_Percentage' in df.columns:
        coverage_values = df['Coverage_Percentage']
    else:
        # Fallback to species count if coverage not available
        coverage_values = df[species_count_col]

    # Calculate novelty factor (ratio of OTU count to species count) for point sizes
    # Higher novelty = larger circles (more OTUs relative to known species)
    novelty_factors = []
    for _, row in df.iterrows():
        otu_count = row['Census_OTU_Count']
        species_count = row[species_count_col]

        if species_count == 0:
            # No genomic data available - assign moderate novelty
            novelty_factors.append(1.0)
        else:
            # Calculate novelty as OTU/Species ratio
            novelty = otu_count / species_count
            novelty_factors.append(novelty)

    # Add novelty factors to dataframe for sorting
    df_with_novelty = df.copy()
    df_with_novelty['Novelty_Factor'] = novelty_factors

    # Get top 10 taxa by novelty factor (instead of coverage)
    df_sorted = df_with_novelty.nlargest(10, 'Novelty_Factor')

    novelty_factors = pd.Series(novelty_factors)
    min_novelty = novelty_factors.min()
    max_novelty = novelty_factors.max()

    # Create point sizes based on novelty factor
    # Higher novelty = larger circles (indicating more environmental diversity than genomic representation)
    point_sizes = []
    for novelty in novelty_factors:
        if novelty == 0:
            point_sizes.append(50)  # Very small for zero novelty
        else:
            # Scale novelty from 100 to 800 (reasonable range for scatter plots)
            if max_novelty > min_novelty:
                normalized = (novelty - min_novelty) / (max_novelty - min_novelty)
                point_sizes.append(100 + normalized * 700)
            else:
                point_sizes.append(300)  # Default size if all novelty values are the same

    # Create professional color palette for top 10 taxa
    # Beautiful yellow → green → blue → violet gradient spectrum
    top_10_colors = [
        '#FFEB3B',  # Bright yellow
        '#CDDC39',  # Yellow-green
        '#8BC34A',  # Light green
        '#4CAF50',  # Green
        '#009688',  # Teal
        '#00BCD4',  # Cyan
        '#2196F3',  # Blue
        '#3F51B5',  # Indigo
        '#673AB7',  # Deep purple
        '#9C27B0'   # Violet
    ]

    # Create color mapping for all points
    point_colors = []
    top_10_names = df_sorted.get('taxon_name_clean', df_sorted.get('taxon_name', df_sorted.index)).tolist()

    for _, row in df_with_novelty.iterrows():
        taxon_name = row.get('taxon_name_clean', row.get('taxon_name', 'Unknown'))
        if taxon_name in top_10_names:
            # Assign specific color to top 10 taxa
            idx = top_10_names.index(taxon_name)
            point_colors.append(top_10_colors[idx])
        else:
            # Gray for all other taxa
            point_colors.append('#cccccc')

    # Create scatter plot with individual colors for top 10 novelty taxa
    scatter = ax.scatter(
        df['Census_OTU_Count'],
        df[species_count_col],
        c=point_colors,  # Individual colors for top 10, gray for others
        s=point_sizes,   # Variable sizes based on novelty factor
        alpha=config.alpha + 0.1,  # Slightly more opaque for vibrant colors
        edgecolors='black',
        linewidth=config.edge_width
    )

    # No graph annotations - taxa will be identified by colors in the legend

    # Add diagonal line for proportional representation
    if len(df) > 0:
        max_val = max(df['Census_OTU_Count'].max(), df[species_count_col].max())
        min_val = min(df['Census_OTU_Count'].min(), df[species_count_col].min())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.7, linewidth=config.line_width,
                label='Proportional Representation')

    # Customize plot with much larger labels
    ax.set_xlabel(f'{gene} EukCensus OTU Count', fontsize=config.label_size + 8, weight='bold', labelpad=20)
    ax.set_ylabel(f'{database} Species Count', fontsize=config.label_size + 8, weight='bold', labelpad=20)
    ax.set_title(title, fontsize=config.title_size + 6, weight='bold', pad=25)

    # Make axis tick labels much larger
    ax.tick_params(axis='both', which='major', labelsize=config.font_size + 2)
    ax.tick_params(axis='both', which='minor', labelsize=config.font_size)

    # Use configurable scales
    ax.set_xscale(config.x_scale)
    ax.set_yscale(config.y_scale)

    ax.grid(True, alpha=config.grid_alpha)

    # No colorbar - using individual colors for top 10 taxa

    # Legend shows top 10 highest novelty taxa
    legend_count = 10
    legend_title = "Top 10 Highest Novelty Taxa:"

    # Create colored legend with markers for top 10 novelty taxa
    legend_elements = []
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        if i < legend_count:
            taxon_name = row.get('taxon_name_clean', row.get('taxon_name', f'Taxon_{i+1}'))
            novelty = row['Novelty_Factor']
            color = top_10_colors[i]

            # Create legend entry with colored marker
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                        markerfacecolor=color, markersize=10,
                                        label=f'{taxon_name}: {novelty:.1f}x'))

    # Add legend in top left corner
    legend = ax.legend(handles=legend_elements, loc='upper left',
                      fontsize=config.legend_size - 2,
                      frameon=True, fancybox=True, shadow=True,
                      title=legend_title, title_fontsize=config.legend_size)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.95)

    # Add properly scaled circle size legend in bottom right
    size_legend_y = 0.30
    ax.text(0.98, size_legend_y + 0.12, 'Circle Size = Novelty Factor', transform=ax.transAxes,
            fontsize=config.annotation_size + 1, weight='bold',
            verticalalignment='bottom', horizontalalignment='right')

    # Use actual novelty values from the data for truly representative legend
    # Get the actual novelty values from the top 10 taxa shown in the colored legend
    top_novelties = df_sorted['Novelty_Factor'].head(10).sort_values(ascending=False)

    # Select 4 representative values from the actual top taxa
    if len(top_novelties) >= 4:
        legend_novelties = [
            top_novelties.iloc[0],  # Highest novelty from top 10
            top_novelties.iloc[2],  # 3rd highest
            top_novelties.iloc[5],  # 6th highest
            top_novelties.iloc[-1]  # Lowest of top 10
        ]
    else:
        # Fallback for small datasets
        legend_novelties = top_novelties.tolist()
        while len(legend_novelties) < 4:
            legend_novelties.append(legend_novelties[-1])

    # Calculate sizes using the EXACT same formula as the main plot
    legend_sizes = []
    for novelty in legend_novelties:
        if novelty == 0:
            legend_sizes.append(50)
        else:
            # Use the exact same scaling formula as the main scatter plot
            if max_novelty > min_novelty:
                normalized = (novelty - min_novelty) / (max_novelty - min_novelty)
                legend_sizes.append(100 + normalized * 700)
            else:
                legend_sizes.append(300)

    # Add example circles with proper scaling for legend
    for i, (size, novelty) in enumerate(zip(legend_sizes, legend_novelties)):
        y_pos = size_legend_y - i*0.08  # Spacing for circles

        # Scale down by a smaller factor to make legend circles more visible
        # Use a scale factor that makes the circles clearly visible and proportional
        legend_display_size = size / 1.3  # Just a smidgen larger for perfect visibility

        ax.scatter(0.94, y_pos, s=legend_display_size, c='#666666', alpha=0.8,
                  transform=ax.transAxes, edgecolors='black', linewidth=0.8)
        ax.text(0.91, y_pos, f'{novelty:.1f}x', transform=ax.transAxes,
                fontsize=config.annotation_size, va='center', ha='right', weight='bold')

    # Add summary statistics in bottom left
    total_taxa = len(df)
    correlation = df['Census_OTU_Count'].corr(df[species_count_col])
    avg_coverage = coverage_values.mean()
    avg_novelty = novelty_factors.mean()

    stats_text = f'Total Taxa: {total_taxa}\nAvg Coverage: {avg_coverage:.1f}%\nAvg Novelty: {avg_novelty:.2f}x\nCorrelation: {correlation:.3f}'
    ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=config.annotation_size + 2,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            weight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()

    logging.info(f"Saved {output_path.name}")

def create_eukprot_subset_plots(df, dataset_name, level, output_dir):
    """Create simple scatter plot for EukProt data."""
    logging.info(f"Creating EukProt scatter plot for {dataset_name}")

    # Filter data with EukProt matches
    df_with_data = df[df['EukProt_Species_Count'] > 0].copy()

    # Apply different Census_OTU_Count cutoffs based on taxonomic level
    if level == 'phylum':
        cutoff = 100  # Most restrictive for phylum
    elif level == 'family':
        cutoff = 10   # Moderate for family
    else:  # genus
        cutoff = 2    # Minimal cutoff for genus (filter out very low counts)

    # Apply the cutoff (all levels now have cutoffs)
    initial_count_low = len(df_with_data)
    df_with_data = df_with_data[df_with_data['Census_OTU_Count'] >= cutoff].copy()
    filtered_count_low = initial_count_low - len(df_with_data)
    if filtered_count_low > 0:
        logging.info(f"Filtered out {filtered_count_low} taxa with Census_OTU_Count < {cutoff} from {dataset_name} ({level} level)")
    else:
        logging.info(f"No taxa filtered from {dataset_name} ({level} level) - all had Census_OTU_Count >= {cutoff}")

    # Apply 100% coverage cutoff to filter out taxa with unrealistically high coverage
    # (e.g., Rhodophyta with few OTU entries but high EukProt representation)
    coverage_col = 'Coverage_Numeric' if 'Coverage_Numeric' in df_with_data.columns else 'Coverage_Percentage'
    if coverage_col in df_with_data.columns:
        initial_count = len(df_with_data)
        df_with_data = df_with_data[df_with_data[coverage_col] <= 100.0].copy()
        filtered_count = initial_count - len(df_with_data)
        if filtered_count > 0:
            logging.info(f"Filtered out {filtered_count} taxa with >100% coverage from {dataset_name}")

    if len(df_with_data) == 0:
        logging.warning(f"No meaningful EukProt data found for {dataset_name} after filtering")
        return

    # Create 18S subdirectory and filename: 18s_eukprot_level.png
    gene_dir = output_dir / '18S'
    filename = f'18s_eukprot_{level}.png'

    # Create single scatter plot
    create_eukprot_single_plot(
        df_with_data, level,
        gene_dir / filename,
        f'18S EukCensus vs EukProt: {level.capitalize()} Level'
    )

def create_eukprot_single_plot(df, level, output_path, title):
    """Create a single EukProt scatter plot with specified parameters."""

    # Create figure with minimal white space
    _, ax = plt.subplots(figsize=(config.figure_width + 1, config.figure_height))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)  # Minimize white space

    # Determine the correct coverage column (phylum data has Coverage_Numeric)
    coverage_col = 'Coverage_Numeric' if 'Coverage_Numeric' in df.columns else 'Coverage_Percentage'

    # Get coverage values for coloring
    if coverage_col == 'Coverage_Numeric':
        coverage_values = df[coverage_col]
    else:
        # For family/genus data, Coverage_Percentage should be numeric
        coverage_values = df[coverage_col]

    # Calculate novelty factor (ratio of OTU count to species count) for point sizes
    # Higher novelty = larger circles (more OTUs relative to known species)
    novelty_factors = []
    for _, row in df.iterrows():
        otu_count = row['Census_OTU_Count']
        species_count = row['EukProt_Species_Count']

        if species_count == 0:
            # No genomic data available - assign moderate novelty
            novelty_factors.append(1.0)
        else:
            # Calculate novelty as OTU/Species ratio
            novelty = otu_count / species_count
            novelty_factors.append(novelty)

    # Add novelty factors to dataframe for sorting
    df_with_novelty = df.copy()
    df_with_novelty['Novelty_Factor'] = novelty_factors

    # Get top 10 taxa by novelty factor (instead of coverage)
    df_sorted = df_with_novelty.nlargest(10, 'Novelty_Factor')

    novelty_factors = pd.Series(novelty_factors)
    min_novelty = novelty_factors.min()
    max_novelty = novelty_factors.max()

    # Create point sizes based on novelty factor
    # Higher novelty = larger circles (indicating more environmental diversity than genomic representation)
    point_sizes = []
    for novelty in novelty_factors:
        if novelty == 0:
            point_sizes.append(50)  # Very small for zero novelty
        else:
            # Scale novelty from 100 to 800 (reasonable range for scatter plots)
            if max_novelty > min_novelty:
                normalized = (novelty - min_novelty) / (max_novelty - min_novelty)
                point_sizes.append(100 + normalized * 700)
            else:
                point_sizes.append(300)  # Default size if all novelty values are the same

    # Create professional color palette for top 10 taxa
    # Beautiful yellow → green → blue → violet gradient spectrum
    top_10_colors = [
        '#FFEB3B',  # Bright yellow
        '#CDDC39',  # Yellow-green
        '#8BC34A',  # Light green
        '#4CAF50',  # Green
        '#009688',  # Teal
        '#00BCD4',  # Cyan
        '#2196F3',  # Blue
        '#3F51B5',  # Indigo
        '#673AB7',  # Deep purple
        '#9C27B0'   # Violet
    ]

    # Create color mapping for all points
    point_colors = []
    top_10_names = df_sorted['Division'].tolist() if 'Division' in df_sorted.columns else df_sorted.index.tolist()

    for _, row in df_with_novelty.iterrows():
        taxon_name = row['Division'] if 'Division' in row else 'Unknown'
        if taxon_name in top_10_names:
            # Assign specific color to top 10 taxa
            idx = top_10_names.index(taxon_name)
            point_colors.append(top_10_colors[idx])
        else:
            # Gray for all other taxa
            point_colors.append('#cccccc')

    # Create scatter plot with individual colors for top 10 novelty taxa
    ax.scatter(
        df['Census_OTU_Count'],
        df['EukProt_Species_Count'],
        s=point_sizes,   # Variable sizes based on novelty factor
        c=point_colors,  # Individual colors for top 10, gray for others
        alpha=config.alpha + 0.1,  # Slightly more opaque for vibrant colors
        edgecolors='black',
        linewidth=config.edge_width
    )

    # No graph annotations - taxa will be identified by colors in the legend

    # Add diagonal line for proportional representation
    if len(df) > 0:
        max_val = max(df['Census_OTU_Count'].max(), df['EukProt_Species_Count'].max())
        min_val = min(df['Census_OTU_Count'].min(), df['EukProt_Species_Count'].min())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.7, linewidth=config.line_width,
                label='Proportional Representation')

    # Customize main plot with much larger labels
    ax.set_xlabel('18S EukCensus OTU Count', fontsize=config.label_size + 8, weight='bold', labelpad=20)
    ax.set_ylabel('EukProt Species Count', fontsize=config.label_size + 8, weight='bold', labelpad=20)
    ax.set_title(title, fontsize=config.title_size + 6, weight='bold', pad=25)

    # Make axis tick labels much larger
    ax.tick_params(axis='both', which='major', labelsize=config.font_size + 2)
    ax.tick_params(axis='both', which='minor', labelsize=config.font_size)

    # Use configurable scales
    ax.set_xscale(config.x_scale)
    ax.set_yscale(config.y_scale)

    ax.grid(True, alpha=config.grid_alpha)

    # No colorbar - using individual colors for top 10 taxa

    # Legend shows top 10 highest novelty taxa
    legend_count = 10
    legend_title = "Top 10 Highest Novelty Taxa:"

    # Create colored legend with markers for top 10 novelty taxa
    legend_elements = []
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        if i < legend_count:
            taxon_name = row['Division'] if 'Division' in row else f'Taxon_{i+1}'
            novelty = row['Novelty_Factor']
            color = top_10_colors[i]

            # Create legend entry with colored marker
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                        markerfacecolor=color, markersize=10,
                                        label=f'{taxon_name}: {novelty:.1f}x'))

    # Add legend in top left corner
    legend = ax.legend(handles=legend_elements, loc='upper left',
                      fontsize=config.legend_size - 2,
                      frameon=True, fancybox=True, shadow=True,
                      title=legend_title, title_fontsize=config.legend_size)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.95)

    # Add properly scaled circle size legend in bottom right
    size_legend_y = 0.30
    ax.text(0.98, size_legend_y + 0.12, 'Circle Size = Novelty Factor', transform=ax.transAxes,
            fontsize=config.annotation_size + 1, weight='bold',
            verticalalignment='bottom', horizontalalignment='right')

    # Use actual novelty values from the data for truly representative legend
    # Get the actual novelty values from the top 10 taxa shown in the colored legend
    top_novelties = df_sorted['Novelty_Factor'].head(10).sort_values(ascending=False)

    # Select 4 representative values from the actual top taxa
    if len(top_novelties) >= 4:
        legend_novelties = [
            top_novelties.iloc[0],  # Highest novelty from top 10
            top_novelties.iloc[2],  # 3rd highest
            top_novelties.iloc[5],  # 6th highest
            top_novelties.iloc[-1]  # Lowest of top 10
        ]
    else:
        # Fallback for small datasets
        legend_novelties = top_novelties.tolist()
        while len(legend_novelties) < 4:
            legend_novelties.append(legend_novelties[-1])

    # Calculate sizes using the EXACT same formula as the main plot
    legend_sizes = []
    for novelty in legend_novelties:
        if novelty == 0:
            legend_sizes.append(50)
        else:
            # Use the exact same scaling formula as the main scatter plot
            if max_novelty > min_novelty:
                normalized = (novelty - min_novelty) / (max_novelty - min_novelty)
                legend_sizes.append(100 + normalized * 700)
            else:
                legend_sizes.append(300)

    # Add example circles with proper scaling for legend
    for i, (size, novelty) in enumerate(zip(legend_sizes, legend_novelties)):
        y_pos = size_legend_y - i*0.08  # Spacing for circles

        # Scale down by a smaller factor to make legend circles more visible
        # Use a scale factor that makes the circles clearly visible and proportional
        legend_display_size = size / 1.3  # Just a smidgen larger for perfect visibility

        ax.scatter(0.94, y_pos, s=legend_display_size, c='#666666', alpha=0.8,
                  transform=ax.transAxes, edgecolors='black', linewidth=0.8)
        ax.text(0.91, y_pos, f'{novelty:.1f}x', transform=ax.transAxes,
                fontsize=config.annotation_size, va='center', ha='right', weight='bold')

    # Add summary statistics in bottom left
    total_taxa = len(df)
    correlation = df['Census_OTU_Count'].corr(df['EukProt_Species_Count'])
    avg_coverage = coverage_values.mean()
    avg_novelty = novelty_factors.mean()

    stats_text = f'Total Taxa: {total_taxa}\nAvg Coverage: {avg_coverage:.1f}%\nAvg Novelty: {avg_novelty:.2f}x\nCorrelation: {correlation:.3f}'
    ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=config.annotation_size + 2,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            weight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()

    logging.info(f"Saved {output_path.name}")



def create_default_config():
    """Create a default configuration file."""
    config_file = 'visualization_config.json'
    if not Path(config_file).exists():
        config.save_to_file(config_file)
        print(f"Created default configuration file: {config_file}")
        print("You can edit this file to customize visualization parameters.")
    else:
        print(f"Configuration file {config_file} already exists.")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate comprehensive scatter plot visualizations')
    parser.add_argument('--config', '-c', type=str, default='visualization_config.json',
                       help='Path to configuration file (default: visualization_config.json)')
    parser.add_argument('--create-config', action='store_true',
                       help='Create a default configuration file and exit')
    parser.add_argument('--font-size', type=int, help='Override font size')
    parser.add_argument('--label-size', type=int, help='Override label font size')
    parser.add_argument('--title-size', type=int, help='Override title font size')
    parser.add_argument('--x-scale', choices=['linear', 'log', 'symlog', 'logit'],
                       help='Override X-axis scale')
    parser.add_argument('--y-scale', choices=['linear', 'log', 'symlog', 'logit'],
                       help='Override Y-axis scale')
    parser.add_argument('--figure-width', type=float, help='Override figure width')
    parser.add_argument('--figure-height', type=float, help='Override figure height')
    parser.add_argument('--dpi', type=int, help='Override figure DPI')

    return parser.parse_args()

def main():
    """Main function to generate all scatter plots."""
    print("="*80)
    print("COMPREHENSIVE SCATTER PLOT VISUALIZER")
    print("="*80)

    # Parse command line arguments
    args = parse_arguments()

    # Handle config creation
    if args.create_config:
        create_default_config()
        return

    setup_logging()

    # Load configuration from file if it exists
    if Path(args.config).exists():
        config.load_from_file(args.config)
        logging.info(f"Loaded configuration from {args.config}")
    else:
        logging.info("Using default configuration")

    # Override config with command line arguments
    if args.font_size:
        config.font_size = args.font_size
    if args.label_size:
        config.label_size = args.label_size
    if args.title_size:
        config.title_size = args.title_size
    if args.x_scale:
        config.x_scale = args.x_scale
    if args.y_scale:
        config.y_scale = args.y_scale
    if args.figure_width:
        config.figure_width = args.figure_width
    if args.figure_height:
        config.figure_height = args.figure_height
    if args.dpi:
        config.dpi = args.dpi

    # Print current configuration
    print(f"\nCurrent Configuration:")
    print(f"  Figure size: {config.figure_width} x {config.figure_height}")
    print(f"  DPI: {config.dpi}")
    print(f"  Font sizes - Main: {config.font_size}, Label: {config.label_size}, Title: {config.title_size}")
    print(f"  Scales - X: {config.x_scale}, Y: {config.y_scale}")
    print(f"  Point size: {config.point_size}, Alpha: {config.alpha}")

    try:
        # Apply matplotlib style configuration
        config.apply_matplotlib_style(plt)
        logging.info("Applied professional matplotlib styling")

        # Setup output directory
        output_dir = create_output_directory()
        logging.info(f"Output directory: {output_dir}")

        # Load all datasets
        logging.info("Loading prokaryote census data...")
        prokaryote_data = load_prokaryote_data()

        logging.info("Loading EukProt data...")
        eukprot_data = load_eukprot_data()

        # Generate subset scatter plots
        logging.info("Generating subset scatter plots...")
        for dataset_name, df in prokaryote_data.items():
            parts = dataset_name.split('_')
            level = parts[-1]  # phylum, family, or genus
            create_subset_scatter_plots(df, dataset_name, level, output_dir)

        for dataset_name, df in eukprot_data.items():
            parts = dataset_name.split('_')
            level = parts[-1]  # phylum, family, or genus
            create_eukprot_subset_plots(df, dataset_name, level, output_dir)

        logging.info("Subset scatter plot generation completed!")

        print("\n" + "="*80)
        print("SUBSET VISUALIZATION GENERATION COMPLETE!")
        print("="*80)

        # List generated files
        total_files = 0
        for subdir in ['16S', '18S']:
            subdir_path = output_dir / subdir
            files = list(subdir_path.glob("*.png"))
            print(f"\n{subdir} Visualizations ({len(files)} files):")
            for file in sorted(files):
                print(f"  - {file.name}")
            total_files += len(files)

        print(f"\nTotal files generated: {total_files}")

        # Print summary statistics
        print("\nDATASET SUMMARY:")
        print(f"Prokaryote datasets loaded: {len(prokaryote_data)}")
        print(f"EukProt datasets loaded: {len(eukprot_data)}")

        print("\nPlot organization:")
        print("  - 16S/: 16S EukCensus vs NCBI and GTDB (phylum, family, genus)")
        print("  - 18S/: 18S EukCensus vs NCBI and EukProt (phylum, family, genus)")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()
