# Cross-Domain Color Recycling Strategy

## 🎨 Overview

Both the comprehensive scatter plot script and the standalone phyla legend generator now implement **cross-domain color recycling** to maximize color diversity and visual distinction between domains.

## 🔄 Recycling Rules

### **Bacteria** → Uses **Eukaryota** Fallback Colors
- **Primary**: Assigned colors from `bacteria_colors` in YAML config
- **Fallback**: When a bacterial phylum is NOT in the config, it uses colors from `fallback_colors.eukaryota`
- **Rationale**: Eukaryotic color palette is visually distinct from bacterial assigned colors

### **Eukaryota** → Uses **Bacteria** Fallback Colors
- **Primary**: Assigned colors from `eukaryota_colors` in YAML config
- **Fallback**: When a eukaryotic division is NOT in the config, it uses colors from `fallback_colors.bacteria`
- **Rationale**: Bacterial color palette is visually distinct from eukaryotic assigned colors

### **Archaea** → Uses **Combined Bacteria + Eukaryota** Fallback Colors
- **Primary**: Assigned colors from `archaea_colors` in YAML config
- **Fallback**: When an archaeal phylum is NOT in the config, it uses combined colors from both `fallback_colors.bacteria` AND `fallback_colors.eukaryota`
- **Rationale**: Archaea typically has fewer taxa, so it benefits from the largest fallback pool

## 📊 Implementation Details

### In `generate_phyla_legend.R` (Lines 272-324)

```r
# Bacteria uses Eukaryota fallback
if (phylum %in% names(color_config$bacteria_colors)) {
  return(color_config$bacteria_colors[[phylum]])
} else {
  shared_pool <- unique(unlist(color_config$fallback_colors$eukaryota))
  pool_index <- ((match(phylum, bacteria_phyla) - 1) %% length(shared_pool)) + 1
  return(shared_pool[pool_index])
}

# Eukaryota uses Bacteria fallback
if (div %in% names(color_config$eukaryota_colors)) {
  return(color_config$eukaryota_colors[[div]])
} else {
  shared_pool <- unique(unlist(color_config$fallback_colors$bacteria))
  pool_index <- ((match(div, eukaryota_divisions) - 1) %% length(shared_pool)) + 1
  return(shared_pool[pool_index])
}

# Archaea uses combined Bacteria + Eukaryota fallback
if (phylum %in% names(color_config$archaea_colors)) {
  return(color_config$archaea_colors[[phylum]])
} else {
  shared_pool <- unique(c(
    unlist(color_config$fallback_colors$bacteria),
    unlist(color_config$fallback_colors$eukaryota)
  ))
  pool_index <- ((match(phylum, archaea_phyla) - 1) %% length(shared_pool)) + 1
  return(shared_pool[pool_index])
}
```

### In `mega_comprehensive_stacked_visual.R` (Lines 745-815 & 239-315)

Same logic applied in:
1. **Scatter plot color assignment** (lines 745-815)
2. **Legend generation** (lines 239-315)

## 🎯 Benefits

### 1. **Maximum Color Diversity**
- Each domain gets access to a different fallback palette
- Reduces color repetition across domains
- Better visual separation between Bacteria, Archaea, and Eukaryota

### 2. **Consistent Cross-Domain Aesthetics**
- Bacteria and Eukaryota "share" each other's fallback palettes
- Creates visual harmony while maintaining distinction
- Archaea gets the best of both worlds

### 3. **Automatic Conflict Avoidance**
- Assigned colors are still removed from fallback pools
- No risk of duplicate colors within a domain
- Cross-domain recycling only affects unmapped taxa

## 📋 Fallback Color Pools (from YAML)

### `fallback_colors.bacteria` (12 colors)
```yaml
["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
 "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA"]
```
→ **Used by**: Eukaryota (for unmapped divisions)

### `fallback_colors.eukaryota` (6 colors)
```yaml
["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD"]
```
→ **Used by**: Bacteria (for unmapped phyla)

### `fallback_colors.archaea` (6 colors)
```yaml
["#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222", "#FF4500"]
```
→ **Not currently used** (Archaea uses Bacteria + Eukaryota instead)

## 🔍 Example Scenario

### Bacteria Domain
- **Abditibacteriota**: `#4c9b34` (assigned in config) ✅
- **Planctomycetota**: `#ff7200` (assigned in config) ✅
- **UnmappedPhylum1**: `#FF6B6B` (from eukaryota fallback) 🔄
- **UnmappedPhylum2**: `#4ECDC4` (from eukaryota fallback) 🔄

### Eukaryota Domain
- **Tubulinea**: `#416b7d` (assigned in config) ✅
- **Rhizaria**: `#55d0ba` (assigned in config) ✅
- **UnmappedDivision1**: `#FF6B6B` (from bacteria fallback) 🔄
- **UnmappedDivision2**: `#4ECDC4` (from bacteria fallback) 🔄

### Archaea Domain
- **Euryarchaeota**: `#f51b7f` (assigned in config) ✅
- **Nanoarchaeota**: `#8c2a50` (assigned in config) ✅
- **UnmappedPhylum1**: `#FF6B6B` (from bacteria fallback) 🔄
- **UnmappedPhylum2**: `#4ECDC4` (from eukaryota fallback) 🔄

## ✅ Scripts Updated

1. ✅ `generate_phyla_legend.R` - Standalone legend generator
2. ✅ `mega_comprehensive_stacked_visual.R` - Comprehensive scatter plots

Both scripts now use identical cross-domain color recycling logic for consistency.

## 🚀 Usage

No changes needed to run the scripts - the cross-domain recycling happens automatically:

```bash
# Generate standalone legend
Rscript generate_phyla_legend.R

# Generate comprehensive plots (includes legend)
Rscript mega_comprehensive_stacked_visual.R
```

The color assignments will be printed to console showing which colors are "assigned" vs "from shared pool".

