# 16S Alluvial Plot Visualization Directory

## Overview

This directory contains R scripts and output files for generating **alluvial (Sankey) plots** that visualize the flow of microbial diversity from **environmental 16S rRNA census data** to **genomic databases (NCBI)**. The plots reveal how well environmental microbial diversity (Bacteria and Archaea) is represented in cultivated genome collections.

### Research Question
> **"What percentage of environmental bacterial and archaeal diversity (detected via 16S rRNA surveys) is represented in genomic databases?"**

---

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    UPSTREAM DATA PROCESSING PIPELINE                         │
└─────────────────────────────────────────────────────────────────────────────┘

    ┌──────────────────────┐         ┌──────────────────────┐
    │  16S Census Data     │         │   NCBI Genome Data   │
    │  (Environmental)     │         │   (Cultivated)       │
    │                      │         │                      │
    │ ../../16S_censusparse│         │ ../../ncbi_parse     │
    │  /csv_16S/           │         │  /csv_ncbi/          │
    │                      │         │                      │
    │ eukcensus16S_by_     │         │ Genome & Species     │
    │   division.csv       │         │   counts by phylum   │
    └──────────┬───────────┘         └──────────┬───────────┘
               │                                │
               │                                │
               └────────────┬───────────────────┘
                            │
                            ▼
                ┌───────────────────────┐
                │  Data Merger          │
                │                       │
                │ ../../Eukcensus_merge │
                │  /16s_merged/         │
                │                       │
                │ 16s_ncbi_merged_      │
                │   clean_phylum.csv    │
                └───────────┬───────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    VISUALIZATION LAYER (THIS DIRECTORY)                      │
└─────────────────────────────────────────────────────────────────────────────┘

                ┌───────────────────────┐
                │  Shared Config        │
                │                       │
                │ ../../shared_config/  │
                │                       │
                │ • taxonomic_color_    │
                │   mapping.yaml        │
                │ • color_mapping_      │
                │   functions.R         │
                └───────────┬───────────┘
                            │
                            ▼
        ┌───────────────────────────────────────────┐
        │         16S ALLUVIAL SCRIPTS              │
        │                                           │
        │  ┌─────────────────────────────────────┐ │
        │  │  BACTERIA PLOTS (Domain-Specific)   │ │
        │  │  • alluvial_16s_bacteria_pct_*.R   │ │
        │  │  • alluvial_16s_bacteria_abs_*.R   │ │
        │  └─────────────────────────────────────┘ │
        │                                           │
        │  ┌─────────────────────────────────────┐ │
        │  │  ARCHAEA PLOTS (Domain-Specific)    │ │
        │  │  • alluvial_16s_archaea_pct_*.R    │ │
        │  │  • alluvial_16s_archaea_abs_*.R    │ │
        │  └─────────────────────────────────────┘ │
        │                                           │
        │  ┌─────────────────────────────────────┐ │
        │  │  LEGACY (Configurable)              │ │
        │  │  • alluvial_16s_pct_values_only.R  │ │
        │  │  • alluvial_16s_abs_values_only.R  │ │
        │  └─────────────────────────────────────┘ │
        └───────────────┬───────────────────────────┘
                        │
                        ▼
        ┌───────────────────────────────────────────┐
        │         OUTPUT DIRECTORIES                │
        │                                           │
        │  📁 figures/                              │
        │     • PNG/PDF plots                       │
        │                                           │
        │  📁 annotations/                          │
        │     • TSV flow annotations                │
        │     • TSV node descriptions               │
        │     • TSV summary statistics              │
        └───────────────────────────────────────────┘
```

---

## Alluvial Plot Structure

Each plot visualizes a **4-node flow** representing the data journey:

```
Node 1              Node 2              Node 3           Node 4
┌─────────────┐    ┌─────────────┐    ┌──────────┐    ┌─────────────┐
│  Genbank    │    │    IMG      │    │  16S     │    │  Genbank    │
│  Genomes    │───▶│  Sequences  │───▶│  OTUs    │───▶│  Species    │
│             │    │             │    │          │    │             │
│ (isolated)  │    │(Cultivated) │    │(Environ.)│    │(Cultivated) │
└─────────────┘    └─────────────┘    └──────────┘    └─────────────┘
```

### Flow Interpretation
- **Wide flows** = High representation (phylum well-represented in databases)
- **Narrow flows** = Low representation (cultivation gaps)
- **Color consistency** = Same phylum tracked across all nodes
- **Flow breaks** = Phyla present in environment but absent from databases

---

## R Scripts

### Main Production Scripts (Domain-Specific) - RECOMMENDED ⭐

#### 1. **Bacteria Plots**
- **`alluvial_16s_bacteria_pct_values_only.R`** ⭐ NEW - Bacteria percentage plot (dedicated)
  - Shows relative proportions (0-100%)
  - Top 12 bacterial phyla by total representation
  - Bacteria-only percentage calculations
  - Uses shared color configuration from YAML
  - Output: `figures/alluvial_16s_bacteria_pct_values_only.png/pdf`
  - Annotations: `annotations/alluvial_16s_bacteria_pct_*.tsv` (3 files)

- **`alluvial_16s_abs_values_only.R`** - Absolute count bacterial alluvial plot (configurable)
  - Shows actual genome/OTU counts
  - Same 12 phyla as percentage plot
  - Set `PROCESS_DOMAIN = "Bacteria"` (line 31)
  - Output: `alluvial_16s_abs_values_only.png/pdf`

#### 2. **Archaea Plots**
- **`alluvial_16s_archaea_pct_values_only.R`** - Archaea percentage plot (dedicated)
  - Archaea-specific percentage calculations
  - Top 8 archaeal phyla
  - Output: `figures/alluvial_16s_archaea_pct_values_only.png/pdf`
  - Annotations: `annotations/alluvial_16s_archaea_pct_*.tsv` (3 files)

- **`alluvial_16s_archaea_abs_values_only.R`** - Archaea absolute count plot (dedicated)
  - Archaea-specific absolute counts
  - Output: `figures/alluvial_16s_archaea_abs_values_only.png/pdf`
  - Annotations: `annotations/alluvial_16s_archaea_abs_*.tsv` (3 files)

### Legacy Scripts (Configurable)
- **`alluvial_16s_pct_values_only.R`** - Configurable percentage plot
  - Set `PROCESS_DOMAIN` variable (line 28) to "Bacteria" or "Archaea"
  - Use dedicated scripts above for cleaner workflow

### Debug/Development Scripts
- `check_plot_data.R` - Data validation
- `debug_*.R` - Troubleshooting scripts
- `test_*.R` - Testing scripts
- `old/` - Archived previous versions

---

## Output Files

All output files are organized into subdirectories for easy reference:

```
16s_alluvial/
├── figures/              # 📊 All PNG and PDF plots
│   ├── alluvial_16s_bacteria_pct_values_only.png/pdf
│   ├── alluvial_16s_archaea_pct_values_only.png/pdf
│   ├── alluvial_16s_archaea_abs_values_only.png/pdf
│   └── ...
│
└── annotations/          # 📋 All TSV annotation files
    ├── alluvial_16s_bacteria_pct_flow_annotations.tsv
    ├── alluvial_16s_bacteria_pct_node_descriptions.tsv
    ├── alluvial_16s_bacteria_pct_summary.tsv
    ├── alluvial_16s_archaea_pct_*.tsv
    ├── alluvial_16s_archaea_abs_*.tsv
    └── ...
```

### 1. **Visual Outputs** (`figures/` directory)

Each R script generates publication-quality plots in both PNG and PDF formats:

#### Bacteria Plots
- **`figures/alluvial_16s_bacteria_pct_values_only.png/pdf`** ⭐ NEW
  - Percentage-based flow diagram
  - Y-axis: 0-100% scale
  - Node labels showing total counts
  - 12 colored phyla + "Other" category
  - Bacteria-specific calculations

- **`alluvial_16s_abs_values_only.png/pdf`** (legacy location)
  - Absolute count flow diagram
  - Y-axis: Actual counts (formatted with commas)
  - Same phyla as percentage plot

#### Archaea Plots
- **`figures/alluvial_16s_archaea_pct_values_only.png/pdf`**
  - Archaea-specific percentage flows
  - 8 archaeal phyla + "Other"

- **`figures/alluvial_16s_archaea_abs_values_only.png/pdf`**
  - Archaea-specific absolute counts

### 2. **Annotation Files** (`annotations/` directory)

These TSV files provide detailed data for each plot:

#### Flow Annotations (`annotations/*_flow_annotations.tsv`)
**Example:** `annotations/alluvial_16s_bacteria_pct_flow_annotations.tsv`

Contains flow data for each taxon at each node:
```
Taxon                    Node                  Absolute_Count  Percentage  Flow_Width
Pseudomonadota          Genbank_Genome_%      1,480,382       57.83       57.83
Pseudomonadota          IMG_Genome_%          308,380         28.68       28.68
Pseudomonadota          16S_OTU_%             53,653          22.93       22.93
Pseudomonadota          Genbank_Species_%     33,394          39.53       39.53
```

**Columns:**
- `Taxon` - Phylum name
- `Node` - Which of the 4 nodes (Genbank_Genome, IMG_Genome, 16S_OTU, Genbank_Species)
- `Node_Order` - Position (1-4)
- `Absolute_Count` - Raw count at this node
- `Percentage` - Percentage of total at this node
- `Flow_Width` - Visual width in the plot

#### Node Descriptions (`annotations/*_node_descriptions.tsv`)
**Example:** `annotations/alluvial_16s_bacteria_pct_node_descriptions.tsv`

Metadata for each node:
```
Node                  Description                           Total_Count  Data_Type
Genbank_Genome_%      Genbank Total Genomes (Bacteria)     2,559,773    Percentage
IMG_Genome_%          IMG Genome Count (16S sequences)     1,075,280    Percentage
16S_OTU_%             16S OTU Count                        234,011      Percentage
Genbank_Species_%     Genbank Total Species (Bacteria)     84,476       Percentage
```

#### Summary Statistics (`annotations/*_summary.tsv`)
**Example:** `annotations/alluvial_16s_bacteria_pct_summary.tsv`

Overall plot metadata:
```
Metric                      Value
Total_Taxa_Shown            13
Top_Taxa_Count              12
Other_Category_Included     Yes
Total_Bacteria_Genomes      2,559,773
Total_16S_Sequences         1,075,280
Total_16S_OTUs              234,011
Total_Bacteria_Species      84,476
Filtering_Method            Top 12 by total representation
Color_System                Bacteria color mapping
```

### 3. **Data Tables (CSV)**
- **`16s_bacterial_alluvial_data_table.csv`** - Raw data used for bacterial plots

---

## Color Configuration System

### Shared Color Mapping

All plots use **centralized color configuration** from:
```
../../shared_config/taxonomic_color_mapping.yaml
```

This ensures **consistency across all visualizations** (alluvial plots, scatter plots, etc.).

### Color Assignment Strategy

#### 1. **Bacteria Colors** (`bacteria_colors` section in YAML)

**Shared Taxa** (appear in both scatter and alluvial plots):
- Use colors from comprehensive scatter plot analysis
- Examples:
  - `Chloroflexota: #32CD32` (Lime green)
  - `Acidobacteriota: #ffb44c` (Light orange)
  - `Verrucomicrobiota: #46bda3` (Aqua)
  - `Planctomycetota: #ff7200` (Orange)
  - `Cyanobacteriota: #bfb1d3` (Lavender)
  - `Actinomycetota: #80456e` (Plum)

**Alluvial-Only Taxa** (major phyla in 16S data):
- `Pseudomonadota: #1f77b4` (Blue) - Most abundant bacteria
- `Bacillota: #228B22` (Forest green) - Gram-positive
- `Campylobacterota: #2ca02c` (Green) - Pathogenic
- `Bacteroidota: #d62728` (Red) - Major gut bacteria
- `Thermodesulfobacteriota: #01414d` (Dark teal)

#### 2. **Archaea Colors** (`archaea_colors` section in YAML)

All 4 major archaea phyla use scatter plot colors:
- `Euryarchaeota: #f51b7f` (Bright pink) - Most common
- `Nitrososphaerota: #ff3f4d` (Red) - Ammonia-oxidizing
- `Thermoproteota: #d19386` (Light brown) - Hyperthermophilic
- `Nanoarchaeota: #8c2a50` (Dark pink) - Highest novelty

#### 3. **"Other" Category**
- Always assigned `#CCCCCC` (Light gray)
- Represents aggregated low-abundance phyla

### Color Loading Functions

Scripts use `color_mapping_functions.R`:
```r
source("../../shared_config/color_mapping_functions.R")
color_config <- load_taxonomic_colors("../../shared_config/taxonomic_color_mapping.yaml")
```

---

## Phylum Filtering Logic

### Bacteria Filtering
1. **Filter by domain**: `domain == "Bacteria"`
2. **Remove invalid entries**: NA, empty, or "N/A" phyla
3. **Calculate total representation**: Sum of percentages across all 4 nodes
4. **Select top 12 phyla** by total representation
5. **Aggregate remaining phyla** into "Other" category
6. **Include unclassified**: `Bacteria.U.phylum` (environmental sequences without database matches)

### Archaea Filtering
1. **Filter by domain**: `domain == "Archaea"`
2. **Select top 8 phyla** (fewer archaea phyla overall)
3. **Same aggregation logic** for "Other"
4. **Include unclassified**: `Archaea.U.phylum`

### ".U." Entries (Unclassified)
- **Source**: `eukcensus16S_by_division.csv` (census-only data)
- **Meaning**: Environmental sequences that don't match any database entries
- **Representation**:
  - 0% at Genbank nodes (no genomes/species)
  - >0% at 16S nodes (environmental detection only)
- **Visual**: Minimal flow width (0.1% minimum for visibility)

---

## Plot Aesthetics

### Design Principles
1. **Thin nodes** (width = 0.02) - Flows meet directly at axis lines
2. **Emphasized flows** (alpha = 0.85, width = 0.18) - Clear visual tracking
3. **Forward guidance** (`lode.guidance = "forward"`) - Reduces flow crossings
4. **Early knot positioning** (`knot.pos = 0.35`) - Tighter curves
5. **Decreasing order** (`decreasing = FALSE`) - Largest flows on top
6. **Node annotations** - Labels showing node name + total count

### Key Visual Features
- **No x-axis labels** - Clean minimal appearance
- **Bold y-axis** - Large (size 20), bold percentage/count labels
- **Legend on right** - Bold phylum names (size 16)
- **White node borders** - Clear separation between strata
- **Extended y-axis** - 110% limit to accommodate node labels

---

## Running the Scripts

### Prerequisites
```r
# Required R packages
install.packages(c("ggplot2", "ggalluvial", "dplyr", "scales", "tidyr", "yaml"))
```

### Directory Structure
The scripts will automatically create output directories if they don't exist:
- `figures/` - All PNG and PDF plot files
- `annotations/` - All TSV annotation files

### Execution

#### Recommended: Domain-Specific Scripts ⭐
```bash
cd parse_repaa_table/visuals/alluvial_plots/16s_alluvial

# Bacteria plots
Rscript alluvial_16s_bacteria_pct_values_only.R     # NEW: Bacteria percentage (dedicated)

# Archaea plots
Rscript alluvial_16s_archaea_pct_values_only.R      # Archaea percentage
Rscript alluvial_16s_archaea_abs_values_only.R      # Archaea absolute
```

#### Legacy: Configurable Scripts
```bash
# Set PROCESS_DOMAIN variable inside scripts before running
Rscript alluvial_16s_pct_values_only.R              # Configurable percentage
Rscript alluvial_16s_abs_values_only.R              # Configurable absolute
```

### Expected Runtime
- Each script: ~10-30 seconds
- Total: ~1-2 minutes for all plots

### Output Organization
After running scripts, your directory will look like:
```
16s_alluvial/
├── figures/
│   ├── alluvial_16s_bacteria_pct_values_only.png
│   ├── alluvial_16s_bacteria_pct_values_only.pdf
│   ├── alluvial_16s_archaea_pct_values_only.png
│   ├── alluvial_16s_archaea_pct_values_only.pdf
│   ├── alluvial_16s_archaea_abs_values_only.png
│   └── alluvial_16s_archaea_abs_values_only.pdf
│
└── annotations/
    ├── alluvial_16s_bacteria_pct_flow_annotations.tsv
    ├── alluvial_16s_bacteria_pct_node_descriptions.tsv
    ├── alluvial_16s_bacteria_pct_summary.tsv
    ├── alluvial_16s_archaea_pct_*.tsv (3 files)
    └── alluvial_16s_archaea_abs_*.tsv (3 files)
```

---

## Key Implementation Details

### Critical Bug Fixes

#### 1. **Alluvium ID Preservation** (Fixed Nov 2024)
**Problem**: Flows were breaking between nodes (e.g., Pseudomonadota not connecting to NCBI nodes)

**Root Cause**: Line reassigning alluvium IDs with `cur_group_id()` after initial assignment
```r
# WRONG - breaks flow continuity
long_data_f <- long_data %>%
  group_by(phylum) %>%
  mutate(alluvium = cur_group_id()) %>%  # Reassigns IDs!
  ungroup()
```

**Fix**: Keep original alluvium IDs intact
```r
# CORRECT - preserves flow continuity
long_data_f <- long_data  # Keep original alluvium IDs
```

#### 2. **Factor Level Matching for Node Labels**
**Problem**: Node annotations not appearing on plots

**Root Cause**: `node_labels` data frame x-column not matching plot factor levels

**Fix**: Explicitly set factor levels
```r
node_labels <- data.frame(
  x = factor(c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
             levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%")),
  label = c(...),
  y = 102
)
```

### Data Processing Pipeline

1. **Load merged data**: `16s_ncbi_merged_clean_phylum.csv`
2. **Filter by domain**: Bacteria or Archaea
3. **Calculate totals**: Genome, species, OTU, sequence counts
4. **Select top phyla**: By total representation across all nodes
5. **Create long format**: One row per phylum-node combination
6. **Assign alluvium IDs**: Stable ID per phylum (critical for flow continuity)
7. **Apply color mapping**: From shared YAML configuration
8. **Generate plot**: ggplot2 + ggalluvial
9. **Export outputs**: PNG, PDF, and TSV annotation files

---

## Interpretation Guide

### What the Plots Show

#### High Flow Width = Good Representation
- **Pseudomonadota** (blue): 57.83% → 28.68% → 22.93% → 39.53%
  - Well-represented across all databases
  - Strong cultivation success

#### Low Flow Width = Cultivation Gap
- **Chloroflexota** (lime green): 0.57% → 5.47% → 8.09% → 0.31%
  - High environmental abundance (8.09% of OTUs)
  - Low genome representation (0.57%)
  - **Cultivation gap identified**

#### Flow Breaks = Database Absence
- **Bacteria.U.phylum**: 0% → 0% → X% → 0%
  - Present in environment (16S OTUs)
  - Completely absent from databases
  - **Novel/uncultivated diversity**

### Percentage vs Absolute Plots

**Percentage plots** - Best for:
- Comparing relative representation
- Identifying proportional gaps
- Publication figures

**Absolute plots** - Best for:
- Understanding scale (millions of genomes vs thousands of OTUs)
- Quantifying actual database content
- Technical reports

---

## Related Directories

### Upstream Data Sources
- **`../../16S_censusparse/`** - Environmental 16S rRNA census data parsing
- **`../../ncbi_parse/`** - NCBI genome and species data extraction
- **`../../Eukcensus_merge/`** - Merges census + database data

### Parallel Visualizations
- **`../18s_alluvial/`** - Eukaryotic 18S rRNA alluvial plots (same structure)
- **`../../scatter_plots/`** - Novelty factor scatter plots (uses same colors)

### Configuration
- **`../../shared_config/`** - Centralized color and taxonomy configuration

---

## Troubleshooting

### Common Issues

**1. "Cannot find merged data file"**
- Ensure you're running from the script directory
- Check that `../../Eukcensus_merge/16s_merged/csv_results/` exists

**2. "Package 'ggalluvial' not found"**
```r
install.packages("ggalluvial")
```

**3. "Flows not connecting properly"**
- Check that alluvium IDs are NOT reassigned after initial creation
- Verify factor levels match across data frames

**4. "Node labels not appearing"**
- Ensure `node_labels$x` is a factor with correct levels
- Check y-position is within plot limits

---

## Citation & Contact

**Project**: Environmental-Genomic Comparative Pipeline
**Purpose**: Quantify representation of environmental microbial diversity in genomic databases
**Data**: 16S rRNA environmental surveys vs NCBI/GTDB genome databases

For questions about this visualization pipeline, refer to the main project README:
```
../../README.md
```

---

## Version History

- **Nov 2024**: Fixed Pseudomonadota flow continuity bug
- **Nov 2024**: Added node annotations with total counts
- **Nov 2024**: Implemented factor-level matching for labels
- **Nov 2024**: Separated bacteria and archaea into dedicated scripts
- **Earlier**: Initial development with combined bacteria/archaea plots

---

## File Inventory

### Production Scripts (4)
- `alluvial_16s_pct_values_only.R`
- `alluvial_16s_abs_values_only.R`
- `alluvial_16s_archaea_pct_values_only.R`
- `alluvial_16s_archaea_abs_values_only.R`

### Output Plots (8 files: 4 PNG + 4 PDF)
- `alluvial_16s_bacteria_pct_values_only.{png,pdf}`
- `alluvial_16s_abs_values_only.{png,pdf}`
- `alluvial_16s_archaea_pct_values_only.{png,pdf}`
- `alluvial_16s_archaea_abs_values_only.{png,pdf}`

### Annotation Files (12 TSV files: 3 per plot type)
- `alluvial_16s_bacteria_pct_{flow_annotations,node_descriptions,summary}.tsv`
- `alluvial_16s_bacteria_abs_{flow_annotations,node_descriptions,summary}.tsv`
- `alluvial_16s_archaea_pct_{flow_annotations,node_descriptions,summary}.tsv`
- `alluvial_16s_archaea_abs_{flow_annotations,node_descriptions,summary}.tsv`

### Debug/Development (~10 files)
- `check_plot_data.R`, `debug_*.R`, `test_*.R`, etc.

### Archive
- `old/` - Previous versions and deprecated scripts

