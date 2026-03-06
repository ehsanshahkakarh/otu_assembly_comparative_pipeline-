# Mega Comprehensive Stacked Visual - Workshop Directory

## 📊 Purpose
This directory contains the working version of the mega comprehensive stacked visualization script.
All outputs (figures, legends, source data) will be generated in this directory.

## 🎯 What This Script Creates

### Main Output: 3-Column × 2-Row Grid
```
┌──────────────────────────────────────────────────────────────┐
│         BACTERIA      │     ARCHAEA      │    EUKARYOTA      │
├──────────────────────────────────────────────────────────────┤
│ PHYLUM  [Scatter 1]   │  [Scatter 2]     │   [Scatter 3]     │
│         (16S data)    │  (16S data)      │   (18S data)      │
├──────────────────────────────────────────────────────────────┤
│ FAMILY  [Scatter 4]   │  [Scatter 5]     │   [Scatter 6]     │
│         (16S data)    │  (16S data)      │   (18S data)      │
└──────────────────────────────────────────────────────────────┘
```

### Scatter Plot Features
- **X-axis**: Environmental OTU count (log scale)
- **Y-axis**: Genomic species count (log scale)
- **Circle size**: Cultivation success (larger = poorly cultured)
- **Circle color**: Phylum/Division-based (from shared YAML config)
- **Labels**: 
  - Top 10 novelty taxa (repel downward)
  - Top 10 overrepresentation taxa (repel upward)

## 📁 Output Files

### Main Figures
- `comprehensive_mega_stacked_visual.png` - Main figure (54" × 20" @ 300 DPI)
- `comprehensive_mega_stacked_visual.pdf` - PDF for Illustrator editing

### Legends
- `phyla_legends/combined_phyla_legend.png` - Horizontal legend strip

### Source Data
- `source_data/` - CSV files for each plot
  - `Bacteria_phylum_source_data.csv`
  - `Bacteria_family_source_data.csv`
  - `Archaea_phylum_source_data.csv`
  - `Archaea_family_source_data.csv`
  - `Eukaryota_phylum_source_data.csv`
  - `Eukaryota_family_source_data.csv`
  - `README_source_data_index.csv` - Index of all source files

## 🔧 Configuration

### Input Data Paths (Verified ✅)
- **16S data**: `../../../Eukcensus_merge/16s_merged/csv_results/`
- **18S data**: `../../../Eukcensus_merge/18s_merged/csv_results/`
- **NCBI data**: `../../../ncbi_parse/csv_ncbi/`
- **Color config**: `../../shared_config/taxonomic_color_mapping.yaml`

### Output Paths
- **All outputs**: Current directory (`mega_scrip/`)
- **Source data**: `source_data/` subdirectory

## 🚀 Usage

### Run the Script
```bash
cd parse_repaa_table/visuals/scatter_plots/mega_scrip
Rscript mega_comprehensive_stacked_visual.R
```

### Expected Runtime
- Data loading: ~10-30 seconds
- Plot generation: ~1-2 minutes
- Total: ~2-3 minutes

## 🎨 Customization Options

### Color Exclusion (lines 13-27)
```r
EXCLUDED_COLORS <- list(
  bacteria = c(),     # Add hex codes to exclude
  archaea = c(),      
  eukaryota = c()     
)
```

### Plot Configuration (lines 54-78)
```r
PLOT_CONFIG <- list(
  plot_width = 54,      # Figure width in inches
  plot_height = 20,     # Figure height in inches
  dpi = 300,            # Resolution
  top_n = 10,           # Number of top taxa to label
  text_size = 11,       # Label text size
  size_range = c(10, 22) # Circle size range
)
```

## 📋 Next Steps

1. **Run the script** to generate initial output
2. **Review the figures** in this directory
3. **Adjust parameters** as needed
4. **Iterate** until satisfied with the visualization

## ✅ Path Verification Status

All input paths have been verified and are accessible:
- ✅ 16S merged data files found
- ✅ 18S merged data files found
- ✅ NCBI reference data found
- ✅ Shared color configuration found

Ready to run!

