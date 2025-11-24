# Shared Taxa Analysis Between Scatter Plots and Alluvial Plots

**Analysis Date:** 2025-11-17  
**Purpose:** Identify shared taxa between comprehensive scatter plots and alluvial plots to ensure consistent color mapping

## Data Sources

### Scatter Plot Source Data
- `visuals/scatter_plots/source_data/{domain}_phylum_source_data.csv`
- Contains taxa with novelty and overrepresentation factors
- Colors from `mega_comprehensive_stacked_visual.R`

### Alluvial Plot Flow Annotations  
- `visuals/alluvial_plots/16s_alluvial/alluvial_16s_{domain}_pct_flow_annotations.tsv`
- Contains taxa that appear in alluvial flow diagrams
- Shows flow data across 4 nodes: Genbank Genomes → IMG Genomes → OTUs → Genbank Species

## Shared Taxa Analysis

### ARCHAEA (100% overlap - 4/4 taxa shared)
**All archaea phyla appear in both visualization types:**

| Taxon | Scatter Plot Color | Alluvial Status | Novelty Factor | Notes |
|-------|-------------------|-----------------|----------------|-------|
| Euryarchaeota | #f51b7f (bright pink) | ✅ Present | 2.05× | Most common archaea |
| Nitrososphaerota | #ff3f4d (red) | ✅ Present | 5.79× | Ammonia-oxidizing |
| Thermoproteota | #d19386 (light brown) | ✅ Present | 2.28× | Hyperthermophilic |
| Nanoarchaeota | #8c2a50 (dark pink) | ✅ Present | 117.70× | Highest novelty |

### BACTERIA (Partial overlap - 6+ shared taxa)
**Shared taxa between both visualization types:**

| Taxon | Scatter Plot Color | Alluvial Status | Novelty Factor | Notes |
|-------|-------------------|-----------------|----------------|-------|
| Chloroflexota | #548877 (teal-green) | ✅ Present | 18.45× | Green non-sulfur bacteria |
| Acidobacteriota | #ffb44c (light orange) | ✅ Present | 24.66× | Important soil bacteria |
| Verrucomicrobiota | #46bda3 (aqua) | ✅ Present | 16.85× | Environmental bacteria |
| Planctomycetota | #ff7200 (orange) | ✅ Present | 26.34× | Unique cell wall-less |
| Cyanobacteriota | #bfb1d3 (lavender) | ✅ Present | 0.63× | Photosynthetic bacteria |
| Actinomycetota | #80456e (plum) | ✅ Present | - | Antibiotic producers |

**Scatter plot only (high novelty taxa):**
- Abditibacteriota (#4c9b34) - 36.71× novelty
- Calditrichota (#cfe99d) - 35.00× novelty  
- Thermomicrobiota (#d58a2f) - 25.55× novelty
- Gemmatimonadota (#f5d24f) - 19.97× novelty
- Armatimonadota (#55e3ff) - 15.34× novelty
- Myxococcota (#7ac7da) - 14.31× novelty

**Alluvial plot only (major database taxa):**
- Pseudomonadota (#1f77b4) - Most abundant bacteria
- Bacillota (#ff7f0e) - Gram-positive bacteria
- Campylobacterota (#2ca02c) - Pathogenic bacteria
- Bacteroidota (#d62728) - Major gut bacteria
- Thermodesulfobacteriota (#01414d) - Sulfate reducers

### EUKARYOTA (High overlap - 9/12 shared taxa)
**Shared taxa between both visualization types:**

| Taxon | Scatter Plot Color | Alluvial Status | Novelty Factor | Notes |
|-------|-------------------|-----------------|----------------|-------|
| Tubulinea | #416b7d (dark blue-gray) | ✅ Present | 1323× | Highest novelty |
| Rhizaria | #55d0ba (teal) | ✅ Present | 341.73× | Amoeboid protists |
| Evosea | #003ce1 (bright blue) | ✅ Present | 78.42× | Amoeboid |
| Alveolata | #c73de4 (magenta) | ✅ Present | 51.67× | Diverse protists |
| Discoba | #65417a (purple) | ✅ Present | 45.42× | Flagellates |
| Metamonada | #cf8ac6 (pink) | ✅ Present | 26.05× | Anaerobic flagellates |
| Discosea | #d24390 (hot pink) | ✅ Present | 15.35× | Amoeboid |
| Stramenopiles | #475093 (dark purple) | ✅ Present | 11.12× | Diverse protists |
| Opisthokonta | #7a9dcd (light blue) | ✅ Present | - | Animals and fungi |

**Scatter plot only:**
- Chlorophyta (#663be6) - 6.75× novelty - Green algae
- Rhodophyta (#69c1d4) - Red algae
- Streptophyta (#68536c) - Land plants

## Color Assignment Strategy

1. **Shared taxa**: Use comprehensive scatter plot colors as primary source
2. **Scatter-only taxa**: Maintain scatter plot colors for consistency
3. **Alluvial-only taxa**: Use existing alluvial plot colors or assign new colors
4. **Special categories**: 
   - "Other" category: #808080 (gray)
   - ".U." entries: Light colors (#ffcc99, etc.)

## Implementation

Updated `visuals/shared_config/taxonomic_color_mapping.yaml` with:
- Comprehensive documentation of shared taxa analysis
- Priority given to scatter plot colors for shared taxa
- Clear categorization of taxa by visualization type
- Consistent color assignments across all visualization types

This ensures that when the same taxon appears in both scatter plots and alluvial plots, it will have the same color, creating visual consistency across the entire project.

## Validation Results

**Total Shared Taxa: 19 across all domains**
- ✅ **Archaea: 4/4 taxa shared (100% overlap)** - Perfect consistency
- ✅ **Bacteria: 6 taxa shared** - Major phyla covered
- ✅ **Eukaryota: 9 taxa shared** - High overlap achieved

All shared taxa have been successfully assigned consistent colors from the comprehensive scatter plot palette, ensuring visual harmony between visualization types.
