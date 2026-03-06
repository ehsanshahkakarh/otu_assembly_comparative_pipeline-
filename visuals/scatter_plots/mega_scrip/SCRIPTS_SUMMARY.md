# Mega Scrip Directory - Scripts Summary

## 📁 Directory Contents

This directory contains two main R scripts for generating comprehensive scatter plot visualizations and phyla legends.

---

## 🎯 Script 1: `mega_comprehensive_stacked_visual.R`

### Purpose
Generates a comprehensive 3×2 grid of scatter plots showing novelty and overrepresentation factors across all three domains of life.

### Output Files
1. **`comprehensive_mega_stacked_visual.png`** - Main figure (54" × 20" @ 300 DPI)
2. **`comprehensive_mega_stacked_visual.pdf`** - PDF version for Illustrator
3. **`phyla_legends/combined_phyla_legend.png`** - Combined phyla legend (20" × 4")
4. **`source_data/*.csv`** - Source data for each subplot (6 files)
5. **`source_data/README_source_data_index.csv`** - Index of all source files

### Grid Layout
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

### Features
- **Top 10 novelty taxa** labeled (repel downward)
- **Top 10 overrepresentation taxa** labeled (repel upward)
- **Circle size** = Cultivation success (larger = poorly cultured)
- **Circle color** = Phylum/Division (from shared YAML config)
- **Cross-domain color recycling** for unmapped taxa

### Usage
```bash
cd parse_repaa_table/visuals/scatter_plots/mega_scrip
Rscript mega_comprehensive_stacked_visual.R
```

### Runtime
~2-3 minutes

---

## 🎨 Script 2: `generate_phyla_legend.R`

### Purpose
Generates a standalone combined phyla legend showing all taxa that appear in the comprehensive scatter plots.

### Output Files
1. **`combined_phyla_legend.png`** - Horizontal legend strip (20" × 4" @ 300 DPI)
2. **`combined_phyla_legend.pdf`** - PDF version for Illustrator

### Legend Format
```
[Bacteria Taxa...] | [Archaea Taxa...] | [Eukaryota Taxa...]
```
- Colored squares with 45° angled taxon labels
- Domain separators (dotted lines)
- Only shows taxa that appear in top 10 novelty or top 10 overrepresentation

### Features
- **Identical logic** to comprehensive script (ensures consistency)
- **Cross-domain color recycling** for unmapped taxa
- **Standalone execution** (doesn't require running comprehensive script)

### Usage
```bash
cd parse_repaa_table/visuals/scatter_plots/mega_scrip
Rscript generate_phyla_legend.R
```

### Runtime
~30-60 seconds

---

## 🎨 Cross-Domain Color Recycling

Both scripts implement intelligent color recycling:

| **Domain** | **Primary Colors** | **Fallback Colors** |
|------------|-------------------|---------------------|
| **Bacteria** | `bacteria_colors` (22 phyla) | `fallback_colors.eukaryota` (6 colors) |
| **Archaea** | `archaea_colors` (5 phyla) | `fallback_colors.bacteria` + `fallback_colors.eukaryota` (18 colors) |
| **Eukaryota** | `eukaryota_colors` (13 divisions) | `fallback_colors.bacteria` (12 colors) |

**See `COLOR_RECYCLING_STRATEGY.md` for detailed explanation.**

---

## 📊 Input Data Requirements

### Required Files
1. **16S merged data** (Bacteria & Archaea):
   - `../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv`
   - `../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_family.csv`

2. **18S merged data** (Eukaryota):
   - `../../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv`
   - `../../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_family.csv`

3. **NCBI reference data**:
   - `../../../ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv`
   - `../../../ncbi_parse/csv_ncbi/ncbi_family_counts.csv`
   - `../../../ncbi_parse/csv_ncbi/ncbi_genus_counts.csv`

4. **Shared color configuration**:
   - `../../shared_config/taxonomic_color_mapping.yaml`

### Data Columns Required
- `census_otu_count` - Environmental OTU counts
- `ncbi_species_count` - Genomic species counts
- `domain` - Domain classification (Bacteria/Archaea)
- `phylum` or `division` - Phylum/Division name
- `family` - Family name (for family-level plots)

---

## 🔧 Configuration

### Plot Configuration (in both scripts)
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

### Legend Configuration (in `generate_phyla_legend.R`)
```r
LEGEND_CONFIG <- list(
  width = 20,           # Width in inches
  height = 4,           # Height in inches
  dpi = 300,            # Resolution
  output_png = "combined_phyla_legend.png",
  output_pdf = "combined_phyla_legend.pdf"
)
```

---

## 📝 Key Differences Between Scripts

| **Aspect** | **mega_comprehensive_stacked_visual.R** | **generate_phyla_legend.R** |
|------------|----------------------------------------|----------------------------|
| **Primary Output** | 6-panel scatter plot grid | Standalone phyla legend |
| **Secondary Output** | Phyla legend + source data | None |
| **Runtime** | ~2-3 minutes | ~30-60 seconds |
| **File Size** | Large (54" × 20") | Small (20" × 4") |
| **Use Case** | Main publication figure | Supplementary legend or separate use |

---

## ✅ Verification

Both scripts print detailed console output showing:
- ✅ Data loading progress
- ✅ Taxa extraction for each domain
- ✅ Color assignments (assigned vs fallback)
- ✅ Output file paths
- ✅ Final dimensions and statistics

---

## 🚀 Quick Start

### Generate Everything
```bash
cd parse_repaa_table/visuals/scatter_plots/mega_scrip

# Generate comprehensive plots (includes legend)
Rscript mega_comprehensive_stacked_visual.R
```

### Generate Only Legend
```bash
cd parse_repaa_table/visuals/scatter_plots/mega_scrip

# Generate standalone legend
Rscript generate_phyla_legend.R
```

---

## 📚 Additional Documentation

- **`README.md`** - Directory overview and setup instructions
- **`COLOR_RECYCLING_STRATEGY.md`** - Detailed color recycling explanation
- **`SCRIPTS_SUMMARY.md`** - This file

---

**Last Updated**: 2025-12-29

