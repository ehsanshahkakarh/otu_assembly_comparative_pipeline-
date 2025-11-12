#!/usr/bin/env python3

import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def read_color_json(json_file):
    """Read JSON color file"""
    with open(json_file, 'r') as f:
        return json.load(f)

def create_color_palette_png():
    """Create PNG visualization of all color palettes"""
    
    # Read JSON files
    bacteria_data = read_color_json('bacteria_colors.json')
    archaea_data = read_color_json('archaea_colors.json')
    eukaryota_data = read_color_json('eukaryota_colors.json')
    
    # Extract colors
    bacteria_colors = list(bacteria_data['bacteria_colors'].values())
    archaea_colors = list(archaea_data['archaea_colors'].values())
    eukaryota_colors = list(eukaryota_data['eukaryota_colors'].values())
    
    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 10))
    fig.suptitle('Color Palettes for Mega Comprehensive Stacked Visual', 
                 fontsize=16, fontweight='bold')
    
    # Plot bacteria colors
    ax1.set_title('Bacteria Colors (16)', fontsize=14, fontweight='bold')
    for i, color in enumerate(bacteria_colors):
        rect = patches.Rectangle((0, i), 1, 0.8, linewidth=2, 
                               edgecolor='white', facecolor=color)
        ax1.add_patch(rect)
        ax1.text(0.5, i + 0.4, f'color_{i+1}\n{color}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, len(bacteria_colors))
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.invert_yaxis()
    
    # Plot archaea colors
    ax2.set_title('Archaea Colors (5)', fontsize=14, fontweight='bold')
    for i, color in enumerate(archaea_colors):
        rect = patches.Rectangle((0, i), 1, 0.8, linewidth=2, 
                               edgecolor='white', facecolor=color)
        ax2.add_patch(rect)
        ax2.text(0.5, i + 0.4, f'color_{i+1}\n{color}', 
                ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, len(archaea_colors))
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.invert_yaxis()
    
    # Plot eukaryota colors
    ax3.set_title('Eukaryota Colors (13)', fontsize=14, fontweight='bold')
    for i, color in enumerate(eukaryota_colors):
        rect = patches.Rectangle((0, i), 1, 0.8, linewidth=2, 
                               edgecolor='white', facecolor=color)
        ax3.add_patch(rect)
        ax3.text(0.5, i + 0.4, f'color_{i+1}\n{color}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, len(eukaryota_colors))
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.invert_yaxis()
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('color_palette_preview.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print("Color palette preview saved to: color_palette_preview.png")
    print(f"\nColor Summary:")
    print(f"=============")
    print(f"Bacteria: {len(bacteria_colors)} colors")
    print(f"Archaea: {len(archaea_colors)} colors")
    print(f"Eukaryota: {len(eukaryota_colors)} colors")
    print(f"Total: {len(bacteria_colors) + len(archaea_colors) + len(eukaryota_colors)} colors")

if __name__ == "__main__":
    create_color_palette_png()
