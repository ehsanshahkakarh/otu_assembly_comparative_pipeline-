# Color Assignment Fix - V2

## 🐛 **Problem Identified**

You reported seeing **duplicate colors** for different bacterial phyla:
- Abditibacteriota, Chlorobiota, Ignavibacteriota → **Same color**
- Thermomicrobiota, Armatimonadota → **Same color**

## 🔍 **Root Cause**

The old code had **TWO critical bugs**:

### **Bug 1: Limited Color Pool**
```r
# OLD (WRONG):
if (domain == "Bacteria") {
  pool <- unlist(color_config$fallback_colors$eukaryota)  # Only 6 colors!
}
```

This only used the **fallback colors** from Eukaryota (6 colors), ignoring:
- ❌ Eukaryota hardcoded colors (10 colors)
- ❌ Archaea hardcoded colors (4 colors)
- ❌ Archaea fallback colors (6 colors)

**Total available: 6 colors** (should be 26!)

### **Bug 2: Wrong Index Calculation**
```r
# OLD (WRONG):
index <- (length(TAXON_COLORS) %% length(pool)) + 1
```

This counted **ALL taxa** (including hardcoded ones), causing the same index to repeat.

---

## ✅ **Solution**

### **Fix 1: Use ALL Colors from Other Domains**
```r
# NEW (CORRECT):
if (domain == "Bacteria") {
  pool <- c(
    unlist(color_config$archaea_colors),           # 4 hardcoded
    unlist(color_config$eukaryota_colors),         # 10 hardcoded
    unlist(color_config$fallback_colors$archaea),  # 6 fallback
    unlist(color_config$fallback_colors$eukaryota) # 6 fallback
  )
  pool <- unique(pool)  # Remove duplicates
}
```

**Total available: 26 unique colors** for Bacteria!

### **Fix 2: Track Used Colors**
```r
# NEW (CORRECT):
used_colors <- unlist(TAXON_COLORS)
available_colors <- setdiff(pool, used_colors)

if (length(available_colors) > 0) {
  color <- available_colors[1]  # Use first unused color
} else {
  # Only recycle if all colors are used
  index <- (length(TAXON_COLORS) %% length(pool)) + 1
  color <- pool[index]
}
```

---

## 🎨 **Color Pool Breakdown**

### **For Bacteria (non-hardcoded):**
- ✅ 4 Archaea hardcoded colors
- ✅ 10 Eukaryota hardcoded colors
- ✅ 6 Archaea fallback colors
- ✅ 6 Eukaryota fallback colors
- **Total: 26 unique colors**

### **For Archaea (non-hardcoded):**
- ✅ 7 Bacteria hardcoded colors
- ✅ 10 Eukaryota hardcoded colors
- ✅ 12 Bacteria fallback colors
- ✅ 6 Eukaryota fallback colors
- **Total: 35 unique colors**

### **For Eukaryota (non-hardcoded):**
- ✅ 7 Bacteria hardcoded colors
- ✅ 4 Archaea hardcoded colors
- ✅ 12 Bacteria fallback colors
- ✅ 6 Archaea fallback colors
- **Total: 29 unique colors**

---

## 📏 **Text Size Increase**

Changed from **11** to **11.5** for better readability:

```r
# OLD:
text_size = 11,

# NEW:
text_size = 11.5,  # Increased from 11 for better readability
```

---

## 🎯 **Expected Result**

Now each phylum should get a **unique color** until all 26+ colors are exhausted:

1. **Abditibacteriota** → Color #1 from pool
2. **Chlorobiota** → Color #2 from pool (DIFFERENT!)
3. **Ignavibacteriota** → Color #3 from pool (DIFFERENT!)
4. **Thermomicrobiota** → Color #4 from pool
5. **Armatimonadota** → Color #5 from pool (DIFFERENT!)

---

## 🚀 **How to Test**

Run the V2 script:
```r
source('mega_comprehensive_stacked_visual_v2.R')
main()
```

Check the legend file:
```
phyla_legends/combined_phyla_legend_v2.png
```

Verify that **all bacterial phyla have unique colors**!

---

**Fixed:** 2026-01-05  
**Changes:** 
1. Expanded color pool to include ALL cross-domain colors (hardcoded + fallback)
2. Track used colors to avoid repeats
3. Increased text size from 11 to 11.5

