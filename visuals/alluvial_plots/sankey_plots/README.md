# Sankey Plots with Colored Flows

This directory contains Sankey diagram generators that create interactive visualizations with **colored flows** matching the shared taxonomic color mapping.

## 🎨 New Colored Flow Sankey Plots

### Features
- **Colored flows**: Each flow is colored to match its phylum/division color (not just gray!)
- **Shared color mapping**: Uses the same colors as alluvial plots and scatter plots
- **Interactive HTML**: Hover over flows and nodes for detailed information
- **Static PNG export**: Optional PNG export using webshot package
- **Configuration-driven**: Uses `alluvial_filtering_config.yaml` for filtering parameters

### Scripts

#### 1. `colored_sankey_16s.py` - 16S Prokaryotic Sankey (Bacteria/Archaea) ⭐ RECOMMENDED
Creates Sankey diagrams for prokaryotic data with colored flows using Python/Plotly.

**Usage:**
```bash
cd visuals/alluvial_plots/sankey_plots
python colored_sankey_16s.py
```

**Configuration:**
- Edit line 24 to change domain: `PROCESS_DOMAIN = "Bacteria"` or `"Archaea"`
- Top N phyla controlled by `alluvial_filtering_config.yaml`
- Colors from `shared_config/taxonomic_color_mapping.yaml`

**Output:**
- `sankey_16s_bacteria_colored_flows.html` - Interactive HTML (4.5 MB)
- `sankey_16s_bacteria_colored_flows.png` - Static PNG (if kaleido installed)

**Status:** ✅ **WORKING** - Tested and verified!

#### 2. `colored_sankey_18s.py` - 18S Eukaryotic Sankey ⭐ RECOMMENDED
Creates Sankey diagrams for eukaryotic data with colored flows using Python/Plotly.

**Usage:**
```bash
cd visuals/alluvial_plots/sankey_plots
python colored_sankey_18s.py
```

**Configuration:**
- Top N divisions controlled by `alluvial_filtering_config.yaml`
- Colors from `shared_config/taxonomic_color_mapping.yaml`

**Output:**
- `sankey_18s_eukaryota_colored_flows.html` - Interactive HTML (4.5 MB)
- `sankey_18s_eukaryota_colored_flows.png` - Static PNG (if kaleido installed)

**Status:** ✅ **WORKING** - Tested and verified!

#### 3. `colored_sankey_16s.R` - 16S Prokaryotic Sankey (R version)
R version of the 16S Sankey plot (may have library compatibility issues).

**Usage:**
```bash
cd visuals/alluvial_plots/sankey_plots
Rscript colored_sankey_16s.R
```

**Status:** ⚠️ May require R library updates

#### 4. `colored_sankey_18s.R` - 18S Eukaryotic Sankey (R version)
R version of the 18S Sankey plot (may have library compatibility issues).

**Usage:**
```bash
cd visuals/alluvial_plots/sankey_plots
Rscript colored_sankey_18s.R
```

**Status:** ⚠️ May require R library updates

## 📊 Data Flow

All Sankey plots show the 4-stage data flow:

```
NCBI Total Genomes → Census Sequences → Census OTUs → NCBI Total Species
```

- **Node 1**: NCBI genome counts
- **Node 2**: EukCensus sequence counts (16S or 18S)
- **Node 3**: EukCensus OTU counts
- **Node 4**: NCBI species counts

## 🎨 Color Mapping

Colors are assigned from the shared taxonomic color mapping:

### Bacteria (examples)
- Planctomycetota: `#ff7200` (Orange)
- Acidobacteriota: `#ffb44c` (Light orange)
- Chloroflexota: `#32CD32` (Lime green)
- Verrucomicrobiota: `#46bda3` (Aqua)

### Archaea (examples)
- Euryarchaeota: `#f51b7f` (Bright pink)
- Nitrososphaerota: `#ff3f4d` (Red)
- Thermoproteota: `#d19386` (Light brown)

### Eukaryota (examples)
- Tubulinea: `#416b7d` (Dark blue-gray)
- Rhizaria: `#55d0ba` (Teal)
- Evosea: `#003ce1` (Bright blue)
- Alveolata: `#c73de4` (Magenta)

### Special
- Other: `#808080` (Gray)

## 📁 File Structure

```
sankey_plots/
├── README.md                                    # This file
├── colored_sankey_16s.R                         # NEW: 16S colored flow Sankey
├── colored_sankey_18s.R                         # NEW: 18S colored flow Sankey
├── sankey_16s_plot.R                            # OLD: Basic 16S Sankey
├── static_sankey_plot.R                         # OLD: Static stacked bars
├── sankey_16s_bacteria_colored_flows.html       # Output: Bacteria interactive
├── sankey_16s_archaea_colored_flows.html        # Output: Archaea interactive
├── sankey_18s_eukaryota_colored_flows.html      # Output: Eukaryota interactive
└── *.png                                        # Static PNG exports
```

## 🔧 Dependencies

### Python (Recommended)
Required Python packages:
- `pandas` - Data manipulation
- `plotly` - Interactive Sankey diagram creation
- `pyyaml` - Configuration loading
- `kaleido` (optional) - PNG export

Install missing packages:
```bash
pip install pandas plotly pyyaml kaleido
```

### R (Alternative)
Required R packages:
- `dplyr` - Data manipulation
- `networkD3` - Sankey diagram creation
- `htmlwidgets` - HTML export
- `scales` - Number formatting
- `yaml` - Configuration loading
- `webshot` (optional) - PNG export

Install missing packages:
```r
install.packages(c("dplyr", "networkD3", "htmlwidgets", "scales", "yaml", "webshot"))
```

## 🚀 Quick Start

Generate all Sankey plots using Python (recommended):

```bash
cd visuals/alluvial_plots/sankey_plots

# Bacteria (default)
python colored_sankey_16s.py

# Archaea (edit script first to change PROCESS_DOMAIN)
# Edit line 24: PROCESS_DOMAIN = "Archaea"
python colored_sankey_16s.py

# Eukaryota
python colored_sankey_18s.py
```

Or using R (if libraries are compatible):

```bash
cd visuals/alluvial_plots/sankey_plots

# Bacteria
Rscript colored_sankey_16s.R

# Archaea (edit script first)
Rscript colored_sankey_16s.R

# Eukaryota
Rscript colored_sankey_18s.R
```

## 📝 Notes

- **Colored flows**: The main improvement over old Sankey plots is that flows are now colored to match their phylum/division
- **Interactive**: HTML files can be opened in any web browser for interactive exploration
- **Consistent colors**: All visualizations use the same color scheme for consistency
- **Configuration**: Filtering parameters are centralized in `../config/alluvial_filtering_config.yaml`

## 🆚 Comparison with Old Sankey Plots

| Feature | Old Sankey | New Colored Flow Sankey |
|---------|-----------|------------------------|
| Flow colors | Gray/default | Colored by phylum |
| Color mapping | Random/hardcoded | Shared taxonomic mapping |
| Configuration | Hardcoded | YAML-based |
| Domains | Bacteria only | Bacteria, Archaea, Eukaryota |
| Consistency | Standalone | Matches alluvial/scatter plots |

## 📧 Author

Enhanced Sankey Team  
Date: 2026-01-28

