# Stacked Alluvial Plots - README

## Overview

Two new R scripts have been created to generate **stacked alluvial plots with single shared legends**:

1. **18S Eukaryotic Stacked Alluvial** - `18s_alluvial/scripts/alluvial_18s_stacked_pct.R`
2. **16S Prokaryotic Stacked Alluvial** - `16s_alluvial/scripts/alluvial_16s_stacked_pct.R`

---

## 1. 18S Eukaryotic Stacked Alluvial

### Script Location
```
visuals/alluvial_plots/18s_alluvial/scripts/alluvial_18s_stacked_pct.R
```

### Description
Creates a stacked alluvial plot for 18S eukaryotic data with a single shared legend.

### Features
- **Data Flow**: NCBI Eukaryota Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Eukaryota Species
- **Top Taxa**: Top 8 eukaryotic divisions by total representation
- **Color Scheme**: Professional eukaryotic colors from shared configuration
- **Percentage-based**: All values normalized to 0-100%
- **Single Legend**: Shared legend for all stacked plots

### Output Files
```
visuals/alluvial_plots/18s_alluvial/figures/alluvial_18s_stacked_pct.png
visuals/alluvial_plots/18s_alluvial/figures/alluvial_18s_stacked_pct.pdf
```

### How to Run
```bash
cd visuals/alluvial_plots/18s_alluvial/scripts
Rscript alluvial_18s_stacked_pct.R
```

---

## 2. 16S Prokaryotic Stacked Alluvial (Bacteria + Archaea)

### Script Location
```
visuals/alluvial_plots/16s_alluvial/scripts/alluvial_16s_stacked_pct.R
```

### Description
Creates a stacked alluvial plot combining **Bacteria** and **Archaea** vertically with a single shared legend.

### Features
- **Data Flow**: Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species
- **Bacteria**: Top 12 bacterial phyla by total representation
- **Archaea**: Top 6-8 archaeal phyla by total representation
- **Color Scheme**: Professional bacteria and archaea colors from shared configuration
- **Percentage-based**: All values normalized to 0-100% (domain-specific calculations)
- **Single Legend**: Shared legend combining both bacteria and archaea phyla
- **Vertical Stacking**: Bacteria plot on top, Archaea plot on bottom

### Output Files
```
visuals/alluvial_plots/16s_alluvial/figures/alluvial_16s_stacked_pct.png
visuals/alluvial_plots/16s_alluvial/figures/alluvial_16s_stacked_pct.pdf
```

### How to Run
```bash
cd visuals/alluvial_plots/16s_alluvial/scripts
Rscript alluvial_16s_stacked_pct.R
```

---

## Technical Details

### Dependencies
Both scripts require the following R packages:
- `ggplot2` - Base plotting
- `ggalluvial` - Alluvial/Sankey diagrams
- `dplyr` - Data manipulation
- `scales` - Number formatting
- `tidyr` - Data tidying
- `yaml` - Configuration loading
- `patchwork` - Plot composition and stacking

### Plot Composition
The scripts use the `patchwork` package to:
1. Create individual alluvial plots
2. Combine them vertically using the `/` operator
3. Extract and share a single legend using `plot_layout(guides = "collect")`
4. Add overall titles using `plot_annotation()`

### Color Configuration
Both scripts use the shared color configuration system:
- **18S**: `visuals/shared_config/taxonomic_color_mapping.yaml` → `eukaryota_colors`
- **16S Bacteria**: `visuals/shared_config/taxonomic_color_mapping.yaml` → `bacteria_colors`
- **16S Archaea**: `visuals/shared_config/taxonomic_color_mapping.yaml` → `archaea_colors`

### Data Sources
- **18S Merged Data**: `Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv`
- **18S Census Data**: `18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv`
- **16S Merged Data**: `Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv`
- **16S Census Data**: `16S_censusparse/csv_16S/eukcensus16S_by_division.csv`

---

## Current Status

### ⚠️ R Environment Issue

There is currently a **library compatibility issue** with the R environment:

```
Error: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
```

This affects **all R scripts** in the project, not just the new stacked alluvial scripts.

### Resolution Options

1. **Update system libraries**: Install/update `libstdc++` to include GLIBCXX_3.4.29
2. **Use different R installation**: Switch to an R installation with compatible libraries
3. **Rebuild R packages**: Reinstall R packages with current system libraries
4. **Use conda environment**: Create a conda environment with compatible R and libraries

### Scripts Are Ready

The stacked alluvial scripts have been **successfully created** and are ready to run once the R environment issue is resolved. The scripts follow the same structure and patterns as the existing working alluvial scripts.

---

## Comparison with Existing Scripts

### 18S Alluvials
- **Existing**: `alluvial_18s_pct_values_only.R` - Single 18S eukaryotic plot
- **New**: `alluvial_18s_stacked_pct.R` - Stacked 18S eukaryotic plots with shared legend

### 16S Alluvials
- **Existing**: 
  - `alluvial_16s_bacteria_pct_values_only.R` - Bacteria only
  - `alluvial_16s_archaea_pct_values_only.R` - Archaea only
- **New**: `alluvial_16s_stacked_pct.R` - Bacteria + Archaea stacked with shared legend

---

## Next Steps

1. **Fix R environment** - Resolve the GLIBCXX library compatibility issue
2. **Run scripts** - Execute both stacked alluvial scripts
3. **Verify outputs** - Check that PNG and PDF files are generated correctly
4. **Review visualizations** - Ensure the stacked plots and shared legends look correct
5. **Adjust if needed** - Fine-tune colors, spacing, or layout as desired

---

## Questions or Issues?

If you encounter any problems or need modifications to the scripts, please let me know!

