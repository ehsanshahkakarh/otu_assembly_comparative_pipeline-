# Color Assignment System - Complete Review

## 📋 **System Overview**

**File:** `mega_comprehensive_stacked_visual.R`  
**Lines:** 42-93 (51 lines total)  
**Purpose:** Assign colors to taxa consistently across all plots using YAML config + cross-domain recycling

---

## ✅ **Current Implementation**

### **Core Components:**

1. **Global Storage** (Line 47)
   ```r
   TAXON_COLORS <- list()  # Maps taxon name → color
   ```

2. **Assignment Function** (Lines 50-88)
   ```r
   assign_taxon_color(taxon, domain, color_config)
   ```

3. **Initialization** (Lines 91-93)
   ```r
   init_color_registry()  # Clears TAXON_COLORS
   ```

---

## 🔄 **Assignment Logic**

### **Step 1: Check if Already Assigned**
```r
if (taxon %in% names(TAXON_COLORS)) {
  return(TAXON_COLORS[[taxon]])  # CONSISTENCY!
}
```

### **Step 2: Check YAML Hardcoded Colors**
```r
if (domain == "Bacteria" && taxon %in% names(color_config$bacteria_colors)) {
  color <- color_config$bacteria_colors[[taxon]]
}
# ... same for Archaea, Eukaryota
```

### **Step 3: Use Cross-Domain Fallback (Round-Robin)**
```r
if (is.null(color)) {
  # Get fallback pool
  if (domain == "Bacteria") {
    pool <- unlist(color_config$fallback_colors$eukaryota)  # 6 colors
  } else if (domain == "Eukaryota") {
    pool <- unlist(color_config$fallback_colors$bacteria)   # 12 colors
  } else {  # Archaea
    pool <- c(
      unlist(color_config$fallback_colors$bacteria),
      unlist(color_config$fallback_colors$eukaryota)        # 18 colors
    )
  }
  
  # Round-robin assignment
  index <- (length(TAXON_COLORS) %% length(pool)) + 1
  color <- pool[index]
}
```

### **Step 4: Store Mapping**
```r
TAXON_COLORS[[taxon]] <<- color
return(color)
```

---

## 📊 **Cross-Domain Recycling**

| **Domain** | **Hardcoded (YAML)** | **Fallback Pool** | **Pool Size** |
|------------|---------------------|-------------------|---------------|
| **Bacteria** | 7 phyla | Eukaryota fallback | 6 colors |
| **Archaea** | 4 phyla | Bacteria + Eukaryota | 18 colors |
| **Eukaryota** | 10 divisions | Bacteria fallback | 12 colors |

---

## 🎯 **Usage in Script**

### **1. Scatter Plot Creation** (Line 870)
```r
group_colors <- sapply(plot_groups, function(taxon) {
  assign_taxon_color(taxon, domain, color_config)
})
```

### **2. Legend Creation** (Lines 310, 320, 330)
```r
bacteria_colors <- sapply(bacteria_phyla, function(phylum) {
  assign_taxon_color(phylum, "Bacteria", color_config)
})
```

### **3. Main Function** (Line 1029)
```r
init_color_registry()
```

---

## ✅ **Key Features**

1. ✅ **Consistency** - Same taxon always gets same color
2. ✅ **YAML-based** - All colors from config file
3. ✅ **Cross-domain** - Bacteria uses Eukaryota fallback, etc.
4. ✅ **Simple** - Just 51 lines, easy to understand
5. ✅ **Round-robin** - Cycles through fallback pool predictably

---

## 📝 **Example Flow**

### **Bacteria Phylum Plot:**
```
1. Planctomycetota → Check YAML → "#ff7200" (hardcoded) ✅
   Store: TAXON_COLORS["Planctomycetota"] = "#ff7200"

2. UnmappedPhylum1 → Not in YAML → Use eukaryota fallback
   index = (1 %% 6) + 1 = 2
   color = eukaryota_fallback[2] = "#4ECDC4"
   Store: TAXON_COLORS["UnmappedPhylum1"] = "#4ECDC4"

3. UnmappedPhylum2 → Not in YAML → Use eukaryota fallback
   index = (2 %% 6) + 1 = 3
   color = eukaryota_fallback[3] = "#45B7D1"
   Store: TAXON_COLORS["UnmappedPhylum2"] = "#45B7D1"
```

### **Bacteria Family Plot (Later):**
```
1. Planctomycetota → Already in TAXON_COLORS!
   Return: "#ff7200" (CONSISTENT!) ✅

2. UnmappedFamily1 → Not in YAML → Use eukaryota fallback
   index = (3 %% 6) + 1 = 4
   color = eukaryota_fallback[4] = "#96CEB4"
   Store: TAXON_COLORS["UnmappedFamily1"] = "#96CEB4"
```

---

## 🚀 **Status: READY TO USE**

The system is:
- ✅ **Implemented** (lines 42-93)
- ✅ **Integrated** (scatter plots + legends)
- ✅ **Tested** (ran successfully in terminal)
- ✅ **Simple** (just 51 lines)
- ✅ **Documented** (this file + GLOBAL_COLOR_REGISTRY.md)

---

**Last Updated:** 2026-01-05

