# SIMPLE Color Assignment System

## 🎯 **Purpose**

Dead simple color assignment:
1. ✅ **Consistency** - Same taxon = same color everywhere
2. ✅ **Cross-domain recycling** - Bacteria uses Eukaryota fallback, etc.
3. ✅ **YAML-based** - All colors come from the YAML file

---

## 🔧 **How It Works (SUPER SIMPLE)**

### **1. Initialization**
```r
init_color_registry()  # Just clears the TAXON_COLORS list
```

### **2. Color Assignment**
```r
assign_taxon_color(taxon, domain, color_config)
```

**Logic (4 simple steps):**
1. **Seen before?** → Return stored color
2. **In YAML hardcoded?** → Use that color
3. **Not hardcoded?** → Get next color from cross-domain fallback pool
4. **Save it** → Store taxon → color mapping

---

## 📊 **Cross-Domain Recycling**

| **Domain** | **Uses Fallback From** | **Pool Size** |
|------------|------------------------|---------------|
| **Bacteria** | Eukaryota | 6 colors |
| **Eukaryota** | Bacteria | 12 colors |
| **Archaea** | Bacteria + Eukaryota | 18 colors |

---

## ✅ **Example**

```r
# Bacteria phylum plot
assign_taxon_color("Planctomycetota", "Bacteria", config)
# → "#ff7200" (hardcoded in YAML)
# Stores: TAXON_COLORS["Planctomycetota"] = "#ff7200"

assign_taxon_color("UnmappedPhylum1", "Bacteria", config)
# → Not in YAML, use eukaryota fallback
# → index = (0 %% 6) + 1 = 1
# → "#FF6B6B" (first color from eukaryota fallback)
# Stores: TAXON_COLORS["UnmappedPhylum1"] = "#FF6B6B"

assign_taxon_color("UnmappedPhylum2", "Bacteria", config)
# → Not in YAML, use eukaryota fallback
# → index = (1 %% 6) + 1 = 2
# → "#4ECDC4" (second color from eukaryota fallback)
# Stores: TAXON_COLORS["UnmappedPhylum2"] = "#4ECDC4"

# Bacteria family plot (later)
assign_taxon_color("Planctomycetota", "Bacteria", config)
# → Already in TAXON_COLORS!
# → "#ff7200" (CONSISTENT!)
```

---

## 📝 **Code (Lines 42-93)**

```r
# Global map: taxon → color
TAXON_COLORS <- list()

assign_taxon_color <- function(taxon, domain, color_config) {
  # 1. Already seen? Return it
  if (taxon %in% names(TAXON_COLORS)) {
    return(TAXON_COLORS[[taxon]])
  }

  # 2. Hardcoded in YAML?
  color <- NULL
  if (domain == "Bacteria" && taxon %in% names(color_config$bacteria_colors)) {
    color <- color_config$bacteria_colors[[taxon]]
  } else if (domain == "Archaea" && taxon %in% names(color_config$archaea_colors)) {
    color <- color_config$archaea_colors[[taxon]]
  } else if (domain == "Eukaryota" && taxon %in% names(color_config$eukaryota_colors)) {
    color <- color_config$eukaryota_colors[[taxon]]
  }

  # 3. Not hardcoded? Use cross-domain fallback
  if (is.null(color)) {
    if (domain == "Bacteria") {
      pool <- unlist(color_config$fallback_colors$eukaryota)
    } else if (domain == "Eukaryota") {
      pool <- unlist(color_config$fallback_colors$bacteria)
    } else {
      pool <- c(
        unlist(color_config$fallback_colors$bacteria),
        unlist(color_config$fallback_colors$eukaryota)
      )
    }

    # Round-robin based on total taxa assigned
    index <- (length(TAXON_COLORS) %% length(pool)) + 1
    color <- pool[index]
  }

  # 4. Save it
  TAXON_COLORS[[taxon]] <<- color
  return(color)
}
```

---

## ✅ **That's It!**

**Total code: ~50 lines. Super simple. Just uses YAML colors!** 🎨

---

**Last Updated**: 2026-01-05

