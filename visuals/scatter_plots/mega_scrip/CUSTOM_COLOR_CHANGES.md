# Custom Color Changes - Chrysiogenota & Rhodophyta

## 🎨 **Changes Made**

### **Problem:**
- **Chrysiogenota** (Bacteria) and **Verrucomicrobiota** (Bacteria) had similar colors
- Needed to distinguish them clearly

### **Solution:**
Changed both **Chrysiogenota** (Bacteria) and **Rhodophyta** (Eukaryota) to **#008f00** (dark green)

---

## 📝 **Updated YAML Configuration**

### **File:** `parse_repaa_table/visuals/shared_config/taxonomic_color_mapping.yaml`

#### **Bacteria Colors (line 33):**
```yaml
bacteria_colors:
  Planctomycetota: "#ff7200"         # Orange
  Acidobacteriota: "#ffb44c"         # Light orange
  Chloroflexota: "#32CD32"           # Lime green
  Verrucomicrobiota: "#46bda3"       # Aqua (UNCHANGED)
  Myxococcota: "#7ac7da"             # Light blue
  Cyanobacteriota: "#bfb1d3"         # Lavender
  Actinomycetota: "#80456e"          # Plum
  Chrysiogenota: "#008f00"           # Dark green ← NEW!
```

#### **Eukaryota Colors (line 58):**
```yaml
eukaryota_colors:
  Tubulinea: "#416b7d"               # Dark blue-gray
  Rhizaria: "#55d0ba"                # Teal
  Evosea: "#003ce1"                  # Bright blue
  Alveolata: "#c73de4"               # Magenta
  Discoba: "#65417a"                 # Purple
  Metamonada: "#cf8ac6"              # Pink
  Discosea: "#d24390"                # Hot pink
  Stramenopiles: "#475093"           # Dark purple
  Chlorophyta: "#663be6"             # Violet
  Streptophyta: "#68536c"            # Gray-purple
  Rhodophyta: "#008f00"              # Dark green ← NEW!
```

---

## 🔄 **How Family-Level Colors Work**

The script automatically handles family-level color inheritance:

### **For Phylum-Level Plots:**
- Colors are assigned via `assign_taxon_color(taxon, domain, color_config)`
- Checks YAML first for hardcoded colors
- Falls back to cross-domain pool if not found

### **For Family-Level Plots:**
- Family taxa inherit their parent phylum/division color
- This happens automatically via the `color_col` variable:
  ```r
  if (domain == "Eukaryota") {
    color_col <- "Division"  # Uses Division column
  } else {
    color_col <- "Phylum"    # Uses Phylum column
  }
  ```

### **Example:**
If a family belongs to **Chrysiogenota** phylum:
1. The family data has `Phylum = "Chrysiogenota"`
2. The plot uses `aes(fill = Phylum)` or `aes(color = Phylum)`
3. The color scale maps `Chrysiogenota → #008f00`
4. **All families under Chrysiogenota automatically get #008f00**

Same logic applies to **Rhodophyta** families!

---

## ✅ **What This Affects**

### **Phylum-Level Plots:**
- ✅ Chrysiogenota → #008f00 (dark green)
- ✅ Rhodophyta → #008f00 (dark green)

### **Family-Level Plots:**
- ✅ All families under Chrysiogenota → #008f00
- ✅ All families under Rhodophyta → #008f00

### **Legends:**
- ✅ Phylum legend shows Chrysiogenota in dark green
- ✅ Phylum legend shows Rhodophyta in dark green
- ✅ Family legend inherits these colors automatically

---

## 🎯 **No Code Changes Needed!**

The existing script already handles this correctly:

1. **`assign_taxon_color()`** checks YAML for hardcoded colors
2. **Family plots** use the Phylum/Division column for coloring
3. **Legends** are built dynamically from the data

**Everything works automatically!** 🎉

---

## 🧪 **Testing**

Run the V2 script:
```r
source('mega_comprehensive_stacked_visual_v2.R')
main()
```

Check:
1. **Phylum plots:** Chrysiogenota and Rhodophyta should be dark green (#008f00)
2. **Family plots:** Families under these phyla should also be dark green
3. **Legends:** Should show the correct dark green color

---

## 📊 **Color Comparison**

| Taxon | Old Color | New Color | Visual |
|-------|-----------|-----------|--------|
| **Verrucomicrobiota** | #46bda3 (Aqua) | #46bda3 (Aqua) | 🟦 UNCHANGED |
| **Chrysiogenota** | (from pool) | **#008f00 (Dark Green)** | 🟩 **NEW** |
| **Rhodophyta** | (from pool) | **#008f00 (Dark Green)** | 🟩 **NEW** |

---

**Updated:** 2026-01-05  
**Files Modified:** `taxonomic_color_mapping.yaml`  
**Script Compatibility:** ✅ No changes needed - works automatically!

