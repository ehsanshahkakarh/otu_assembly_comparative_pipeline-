# Colored Flow Sankey Plots - Summary

## 🎉 What We Created

We've successfully created **standard Sankey plots with colored flows** that use your shared taxonomic color mapping configuration. These are a major improvement over the old Sankey plots!

## ✨ Key Features

### 1. **Colored Flows** 🎨
- **OLD**: Flows were gray/default color
- **NEW**: Each flow is colored to match its phylum/division color
- Makes it much easier to track specific taxa through the data pipeline

### 2. **Shared Color Mapping** 🔗
- Uses the same colors as your alluvial plots and scatter plots
- Colors defined in `shared_config/taxonomic_color_mapping.yaml`
- Ensures consistency across all visualizations

### 3. **Configuration-Driven** ⚙️
- Filtering parameters from `alluvial_filtering_config.yaml`
- Top N taxa automatically selected based on config
- Easy to adjust without editing code

### 4. **Interactive HTML** 💻
- Hover over flows to see exact counts
- Zoom and pan capabilities
- Self-contained HTML files (no external dependencies)

## 📊 Generated Plots

### 1. Bacteria Sankey (16S)
- **File**: `sankey_16s_bacteria_colored_flows.html` (4.5 MB)
- **Top phyla**: 15 + Other
- **Total genomes**: 2,559,773
- **Total sequences**: 1,055,730
- **Total OTUs**: 229,790
- **Total species**: 84,476

**Top 5 Phyla by Representation:**
1. Pseudomonadota - #FF6B6B (red)
2. Bacillota - #4ECDC4 (cyan)
3. Bacteroidota - #45B7D1 (blue)
4. Actinomycetota - #80456e (plum)
5. Campylobacterota - #FFEAA7 (yellow)

### 2. Archaea Sankey (16S)
- **File**: `sankey_16s_archaea_colored_flows.html` (4.5 MB)
- **Top phyla**: 4 + Other
- **Total genomes**: 14,849
- **Total sequences**: 41,120
- **Total OTUs**: 8,805
- **Total species**: 2,287

**All Archaea Phyla:**
1. Euryarchaeota - #f51b7f (bright pink)
2. Nitrososphaerota - #ff3f4d (red)
3. Nanoarchaeota - #8c2a50 (dark pink)
4. Thermoproteota - #d19386 (light brown)

### 3. Eukaryota Sankey (18S)
- **File**: `sankey_18s_eukaryota_colored_flows.html` (4.5 MB)
- **Top divisions**: 12 + Other
- **Total genomes**: 49,999
- **Total sequences**: 342,447
- **Total OTUs**: 56,162
- **Total species**: 20,712

**Top 5 Divisions by Representation:**
1. Opisthokonta - #FF6B6B (red)
2. Alveolata - #c73de4 (magenta)
3. Stramenopiles - #475093 (dark purple)
4. Rhizaria - #55d0ba (teal)
5. Discoba - #65417a (purple)

## 🚀 How to Use

### Quick Start
```bash
cd visuals/alluvial_plots/sankey_plots

# Generate all three plots at once
./generate_all_sankey_plots.sh

# Or generate individually
python colored_sankey_16s.py    # Bacteria (default)
python colored_sankey_18s.py    # Eukaryota
```

### View the Plots
1. Transfer HTML files to your local machine
2. Open in any web browser (Chrome, Firefox, Safari, etc.)
3. Hover over flows to see details
4. Zoom and pan to explore

## 📁 Files Created

### Python Scripts (Recommended)
- `colored_sankey_16s.py` - 16S Bacteria/Archaea Sankey generator
- `colored_sankey_18s.py` - 18S Eukaryota Sankey generator
- `generate_all_sankey_plots.sh` - Batch generator for all plots

### R Scripts (Alternative)
- `colored_sankey_16s.R` - R version (may have library issues)
- `colored_sankey_18s.R` - R version (may have library issues)

### Output Files
- `sankey_16s_bacteria_colored_flows.html` - Interactive Bacteria Sankey
- `sankey_16s_archaea_colored_flows.html` - Interactive Archaea Sankey
- `sankey_18s_eukaryota_colored_flows.html` - Interactive Eukaryota Sankey

### Documentation
- `README.md` - Detailed documentation
- `SUMMARY.md` - This file

## 🎨 Color Consistency

All plots use colors from `shared_config/taxonomic_color_mapping.yaml`:

**Bacteria Examples:**
- Planctomycetota: #ff7200 (Orange)
- Chloroflexota: #32CD32 (Lime green)
- Actinomycetota: #80456e (Plum)

**Archaea Examples:**
- Euryarchaeota: #f51b7f (Bright pink)
- Nitrososphaerota: #ff3f4d (Red)

**Eukaryota Examples:**
- Tubulinea: #416b7d (Dark blue-gray)
- Alveolata: #c73de4 (Magenta)
- Rhizaria: #55d0ba (Teal)

## 🆚 Comparison with Old Sankey Plots

| Feature | Old Sankey | New Colored Flow Sankey |
|---------|-----------|------------------------|
| Flow colors | Gray/default | Colored by taxon |
| Color mapping | Random/hardcoded | Shared taxonomic mapping |
| Configuration | Hardcoded in script | YAML-based |
| Domains | Bacteria only | Bacteria, Archaea, Eukaryota |
| Consistency | Standalone | Matches alluvial/scatter plots |
| Technology | R (networkD3) | Python (Plotly) + R |
| Status | ⚠️ Limited | ✅ Fully functional |

## 💡 Next Steps

1. **View the plots**: Transfer HTML files to your local machine and open in browser
2. **Customize colors**: Edit `shared_config/taxonomic_color_mapping.yaml` if needed
3. **Adjust filtering**: Edit `alluvial_filtering_config.yaml` to change top N taxa
4. **Export images**: Install `kaleido` for PNG export: `pip install kaleido`

## 📧 Notes

- Python version is recommended (tested and working)
- R version may require library updates
- HTML files are self-contained (no internet required to view)
- Flow opacity is set to 0.4 for better visibility of overlapping flows

---

**Created**: 2026-01-28  
**Status**: ✅ Fully functional and tested  
**Technology**: Python 3 + Plotly + YAML configuration

