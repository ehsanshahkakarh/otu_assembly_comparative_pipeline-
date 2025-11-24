#!/usr/bin/env python3
"""
Color Palette Visualizer for Comprehensive Visual Config
========================================================
Reads the comprehensive_visual_config.yaml file and creates a visual representation
of all color schemes organized by domain.

Created: 2025-11-14
"""

import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path

def load_config(config_path):
    """Load the YAML configuration file."""
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def apply_color_exclusions(config):
    """Apply color exclusions without modifying the original YAML."""
    # Colors to exclude (specified by user)
    excluded_colors = {
        'bacteria': ['#72c859', '#9fe18b', '#a67c28', '#bfb1d3'],  # Remove these bacteria colors
        'archaea': [],   # No archaea exclusions
        'eukaryota': ['#69c1d4', '#68536c']  # Remove light blue and gray-purple from eukaryota
    }

    # Apply exclusions
    modified_config = config.copy()
    for domain, colors_to_exclude in excluded_colors.items():
        if colors_to_exclude:
            original_colors = modified_config['color_palettes'][domain]['colors']
            filtered_colors = [color for color in original_colors if color not in colors_to_exclude]
            modified_config['color_palettes'][domain]['colors'] = filtered_colors

            print(f"🎨 {domain.upper()}: Excluded {len(colors_to_exclude)} colors")
            for color in colors_to_exclude:
                print(f"   - {color}")
            print(f"   Remaining: {len(filtered_colors)} colors")

    return modified_config

def create_color_visualization(config, output_path="color_palette_visualization.png"):
    """Create a comprehensive visualization of all color palettes."""
    
    color_palettes = config['color_palettes']
    
    # Set up the figure
    fig, axes = plt.subplots(3, 1, figsize=(16, 12))
    fig.suptitle('Comprehensive Visual Config - Color Palettes by Domain', 
                 fontsize=20, fontweight='bold', y=0.95)
    
    domains = ['bacteria', 'archaea', 'eukaryota']
    domain_titles = ['Bacteria (16 colors)', 'Archaea (5 colors)', 'Eukaryota (13 colors)']
    
    for idx, (domain, title) in enumerate(zip(domains, domain_titles)):
        ax = axes[idx]
        palette = color_palettes[domain]
        colors = palette['colors']
        
        # Create color swatches
        n_colors = len(colors)
        swatch_width = 0.8 / n_colors
        
        for i, color in enumerate(colors):
            # Create rectangle for each color
            rect = patches.Rectangle((i * swatch_width, 0.2), swatch_width * 0.9, 0.6, 
                                   facecolor=color, edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            
            # Add color code text
            ax.text(i * swatch_width + swatch_width/2, 0.1, color, 
                   ha='center', va='center', fontsize=9, rotation=45)
            
            # Add index number
            ax.text(i * swatch_width + swatch_width/2, 0.9, str(i+1), 
                   ha='center', va='center', fontsize=10, fontweight='bold')
        
        # Set up axis
        ax.set_xlim(0, 0.8)
        ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Add color group information if available
        if 'color_groups' in palette:
            group_text = "Color Groups: " + ", ".join(palette['color_groups'].keys())
            ax.text(0.4, -0.15, group_text, ha='center', va='center', 
                   fontsize=10, style='italic', transform=ax.transAxes)
        
        # Add description
        if 'description' in palette:
            ax.text(0.4, -0.25, palette['description'], ha='center', va='center', 
                   fontsize=10, transform=ax.transAxes)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, bottom=0.1)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Color palette visualization saved to: {output_path}")
    
    return fig

def print_color_summary(config):
    """Print a summary of all colors organized by domain."""
    color_palettes = config['color_palettes']
    
    print("\n" + "="*60)
    print("COLOR PALETTE SUMMARY")
    print("="*60)
    
    for domain_name, palette in color_palettes.items():
        print(f"\n{domain_name.upper()} ({len(palette['colors'])} colors):")
        print("-" * 40)
        
        for i, color in enumerate(palette['colors'], 1):
            print(f"{i:2d}. {color}")
        
        if 'color_groups' in palette:
            print(f"\nColor Groups:")
            for group_name, group_colors in palette['color_groups'].items():
                print(f"  • {group_name}: {group_colors}")
        
        if 'division_mapping' in palette:
            print(f"\nDivision Mapping:")
            for division, color in palette['division_mapping'].items():
                print(f"  • {division}: {color}")

def main():
    """Main function to run the color palette visualizer."""
    # Path to config file
    config_path = Path("config/comprehensive_visual_config.yaml")

    if not config_path.exists():
        print(f"Error: Config file not found at {config_path}")
        return

    # Load configuration
    config = load_config(config_path)

    # Apply color exclusions
    print("Applying color exclusions...")
    modified_config = apply_color_exclusions(config)

    # Print color summary to console
    print_color_summary(modified_config)

    # Create visualization with filtered colors
    create_color_visualization(modified_config, "color_palette_visualization_filtered.png")

    # Show the plot
    plt.show()

if __name__ == "__main__":
    main()
