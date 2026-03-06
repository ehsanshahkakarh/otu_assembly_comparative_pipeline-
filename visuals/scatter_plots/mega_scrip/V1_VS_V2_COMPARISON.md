# V1 vs V2 Comparison

## 📋 **Overview**

Both versions use the **EXACT SAME** color assignment logic (lines 42-93).
The only differences are the **output file names** to allow side-by-side comparison.

---

## 🔧 **Color Assignment Logic (IDENTICAL)**

Both versions use this simple system:

```r
# Global map: taxon → color
TAXON_COLORS <- list()

assign_taxon_color <- function(taxon, domain, color_config) {
  # 1. Already seen? Return it
  if (taxon %in% names(TAXON_COLORS)) return(TAXON_COLORS[[taxon]])
  
  # 2. Hardcoded in YAML? Use it
  if (domain == "Bacteria" && taxon %in% names(color_config$bacteria_colors)) {
    color <- color_config$bacteria_colors[[taxon]]
  } # ... (same for Archaea, Eukaryota)
  
  # 3. Not hardcoded? Use cross-domain fallback (round-robin)
  if (is.null(color)) {
    pool <- get_fallback_pool(domain, color_config)
    index <- (length(TAXON_COLORS) %% length(pool)) + 1
    color <- pool[index]
  }
  
  # 4. Save it
  TAXON_COLORS[[taxon]] <<- color
  return(color)
}
```

---

## 📁 **Output File Differences**

| **File Type** | **V1 (Original)** | **V2 (Comparison)** |
|---------------|-------------------|---------------------|
| **Main PNG** | `comprehensive_mega_stacked_visual.png` | `comprehensive_mega_stacked_visual_v2.png` |
| **Main PDF** | `comprehensive_mega_stacked_visual.pdf` | `comprehensive_mega_stacked_visual_v2.pdf` |
| **Legend** | `phyla_legends/combined_phyla_legend.png` | `phyla_legends/combined_phyla_legend_v2.png` |
| **Source Data** | `source_data/` | `source_data_v2/` |

---

## 🎯 **Purpose**

The V2 version exists **only** to allow you to:
1. Run both scripts
2. Compare the outputs side-by-side
3. Verify they produce identical results
4. Choose which one to keep

---

## ✅ **Expected Result**

Both versions should produce **IDENTICAL** visualizations because they use the **SAME** color assignment logic.

---

## 🚀 **How to Compare**

### **Run V1:**
```r
source('mega_comprehensive_stacked_visual.R')
main()
```

### **Run V2:**
```r
source('mega_comprehensive_stacked_visual_v2.R')
main()
```

### **Compare Outputs:**
```bash
# Compare PNG files
ls -lh comprehensive_mega_stacked_visual*.png

# Compare legends
ls -lh phyla_legends/combined_phyla_legend*.png

# Compare source data
diff -r source_data/ source_data_v2/
```

---

## 📊 **Console Output Differences**

### **V1:**
```
Comprehensive Mega Stacked Visual Creation
==========================================
✅ Initialized color assignment (YAML colors + cross-domain recycling)
```

### **V2:**
```
Comprehensive Mega Stacked Visual Creation (V2 - SIMPLIFIED)
=============================================================
✅ Initialized color assignment (V2: Simple round-robin)
```

---

## 🎨 **Color Assignment (IDENTICAL)**

Both versions:
- ✅ Use YAML hardcoded colors for 21 taxa
- ✅ Use cross-domain recycling for non-hardcoded taxa
- ✅ Ensure consistency (same taxon = same color)
- ✅ Use simple round-robin assignment

---

## 🗑️ **After Comparison**

Once you verify they're identical, you can:
1. Keep V1 and delete V2
2. Keep V2 and delete V1
3. Keep both (if you want)

---

**Created:** 2026-01-05  
**Purpose:** Side-by-side comparison of identical color assignment logic

