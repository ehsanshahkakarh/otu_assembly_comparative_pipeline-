# iTOL Tree Upload Guide

## 🌳 Taxonomic Trees Generated

Four hierarchical taxonomic trees have been generated from the species-grouped dataset:

### Files Generated:

| Domain | Tree File | Metadata File | Species | Genomes |
|--------|-----------|---------------|---------|---------|
| **Bacteria** | `bacteria_tree_*.nwk` | `bacteria_tree_metadata_*.txt` | 92,757 | 2,619,208 |
| **Archaea** | `archaea_tree_*.nwk` | `archaea_tree_metadata_*.txt` | 3,319 | 29,744 |
| **Eukaryota** | `eukaryota_tree_*.nwk` | `eukaryota_tree_metadata_*.txt` | 20,945 | 50,636 |
| **Viruses** | `viruses_tree_*.nwk` | `viruses_tree_metadata_*.txt` | 32,099 | 220,797 |

---

## 📤 How to Upload to iTOL

### Step 1: Upload Tree File

1. Go to iTOL: https://itol.embl.de/upload.cgi
2. Click "Choose File" and select one of the `.nwk` files (e.g., `bacteria_tree_20260129_182506.nwk`)
3. Click "Upload"
4. Wait for the tree to load and render

### Step 2: Add Metadata (Optional)

1. Once the tree is displayed, click "Datasets" in the control panel
2. Click "Upload annotation files"
3. Select the corresponding metadata file (e.g., `bacteria_tree_metadata_20260129_182506.txt`)
4. The genome counts will be displayed as bar charts on the tree

### Step 3: Customize Visualization

iTOL provides many customization options:
- **Display mode**: Normal, Circular, Unrooted
- **Branch lengths**: Ignore or use
- **Labels**: Show/hide, font size, rotation
- **Colors**: Customize by clade
- **Export**: PDF, SVG, PNG formats

---

## 🔍 Tree Structure

Each tree is built hierarchically from the taxonomic lineage strings:

```
Domain
  └─ Phylum
      └─ Class
          └─ Order
              └─ Family
                  └─ Genus
                      └─ Species
```

### Node Naming Convention:

Each node is labeled as: `TaxonName_Rank_TaxID`

Example: `Escherichia_genus_561`

This ensures unique identifiers even when taxon names are reused across different lineages.

---

## 📊 Tree Statistics

### Bacteria Tree:
- **Size**: 4.4 MB (4,557,030 characters)
- **Species**: 92,757
- **Genomes**: 2,619,208
- **Isolate**: 79.4%
- **Uncultured**: 20.6%

### Archaea Tree:
- **Size**: 196 KB (199,781 characters)
- **Species**: 3,319
- **Genomes**: 29,744
- **Isolate**: 9.0%
- **Uncultured**: 91.0% (mostly MAGs!)

### Eukaryota Tree:
- **Size**: 1.2 MB (1,228,669 characters)
- **Species**: 20,945
- **Genomes**: 50,636
- **Isolate**: 96.5%
- **Uncultured**: 3.5%

### Viruses Tree:
- **Size**: 1.5 MB (1,520,133 characters)
- **Species**: 32,099
- **Genomes**: 220,797
- **Isolate**: 94.0%
- **Uncultured**: 6.0%

---

## ⚠️ Important Notes

1. **Large Trees**: The Bacteria tree is very large (92,757 species). iTOL may take some time to render it.

2. **Browser Performance**: For best performance with large trees, use:
   - Chrome or Firefox (latest versions)
   - Close other browser tabs
   - Increase browser memory if needed

3. **Metadata**: The metadata files contain genome counts for each species, which will be displayed as bar charts.

4. **Newick Format**: All trees are in standard Newick format, compatible with any phylogenetic tree viewer.

---

## 🎨 Suggested iTOL Settings

For best visualization:

- **Display mode**: Circular (for large trees)
- **Branch lengths**: Ignore (since this is taxonomic, not phylogenetic)
- **Label display**: Show only at certain zoom levels
- **Leaf sorting**: By name or by dataset value
- **Color ranges**: Use gradient for genome counts

---

## 📁 File Locations

All tree files are in:
```
ncbi_parse/metadata/nev_parse_meth/sanity_check/output/
```

---

## 🔗 Useful Links

- **iTOL**: https://itol.embl.de/
- **iTOL Help**: https://itol.embl.de/help.cgi
- **iTOL Datasets**: https://itol.embl.de/help.cgi#datasets

---

Generated: 2026-01-29

