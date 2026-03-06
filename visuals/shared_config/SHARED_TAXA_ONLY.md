# Shared Taxa Color Configuration

## 🎯 **Philosophy**

This YAML configuration **ONLY hardcodes colors for taxa that appear in BOTH scatter plots AND alluvial plots**.

All other taxa use the **cross-domain fallback color pools** for maximum color diversity.

---

## ✅ **Hardcoded Shared Taxa (21 total)**

### 🦠 **Bacteria (7 shared taxa)**
| **Phylum** | **Color** | **Novelty Factor** |
|------------|-----------|-------------------|
| Planctomycetota | `#ff7200` (Orange) | 26.34× |
| Acidobacteriota | `#ffb44c` (Light orange) | 24.66× |
| Chloroflexota | `#32CD32` (Lime green) | 18.45× |
| Verrucomicrobiota | `#46bda3` (Aqua) | 16.85× |
| Myxococcota | `#7ac7da` (Light blue) | 14.36× |
| Cyanobacteriota | `#bfb1d3` (Lavender) | - |
| Actinomycetota | `#80456e` (Plum) | - |

### 🔥 **Archaea (4 shared taxa)**
| **Phylum** | **Color** | **Novelty Factor** |
|------------|-----------|-------------------|
| Euryarchaeota | `#f51b7f` (Bright pink) | 2.05× |
| Nitrososphaerota | `#ff3f4d` (Red) | 5.79× |
| Thermoproteota | `#d19386` (Light brown) | 2.28× |
| Nanoarchaeota | `#8c2a50` (Dark pink) | 117.70× |

### 🌿 **Eukaryota (10 shared taxa)**
| **Division** | **Color** | **Novelty Factor** |
|------------|-----------|-------------------|
| Tubulinea | `#416b7d` (Dark blue-gray) | 1323× |
| Rhizaria | `#55d0ba` (Teal) | 341.73× |
| Evosea | `#003ce1` (Bright blue) | 78.42× |
| Alveolata | `#c73de4` (Magenta) | 51.67× |
| Discoba | `#65417a` (Purple) | 45.42× |
| Metamonada | `#cf8ac6` (Pink) | 26.05× |
| Discosea | `#d24390` (Hot pink) | 15.35× |
| Stramenopiles | `#475093` (Dark purple) | 11.12× |
| Chlorophyta | `#663be6` (Violet) | 6.75× |
| Streptophyta | `#68536c` (Gray-purple) | (pct alluvial) |

---

## 🔄 **Cross-Domain Fallback Pools**

All **non-shared taxa** use cross-domain recycling:

- **Bacteria** → Uses **Eukaryota** fallback colors (6 colors)
- **Eukaryota** → Uses **Bacteria** fallback colors (12 colors)
- **Archaea** → Uses **Bacteria + Eukaryota** fallback colors (18 colors)

---

## 📊 **Examples**

### **Scatter Plot Only Taxa** (use fallback pools):
- **Bacteria**: Abditibacteriota, Calditrichota, Thermomicrobiota, Gemmatimonadota, Armatimonadota, Chlorobiota, Chlamydiota, Ignavibacteriota, etc.
- **Eukaryota**: Opisthokonta, Rhodophyta

### **Alluvial Plot Only Taxa** (use fallback pools):
- **Bacteria**: Pseudomonadota, Bacillota, Campylobacterota, Bacteroidota, Thermodesulfobacteriota

---

## ✅ **Benefits**

1. ✅ **Consistency** - Shared taxa have identical colors across all visualizations
2. ✅ **Simplicity** - Only 21 hardcoded colors instead of 40+
3. ✅ **Flexibility** - Non-shared taxa get diverse colors from cross-domain pools
4. ✅ **Maintainability** - Easy to update when new shared taxa are identified

---

## 📝 **Last Updated**
2026-01-04

**Shared Taxa Counts:**
- Bacteria: 7
- Archaea: 4
- Eukaryota: 10
- **Total: 21 hardcoded colors**

