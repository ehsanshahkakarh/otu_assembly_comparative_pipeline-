# 18S Systematic Resolver - User Guide

**What it does:** Resolves taxonomic names that are NOT in NCBI taxonomy by mapping them to known parent taxids.

---

## 📊 Overview

The systematic resolver handles **environmental clades** and **non-NCBI taxonomic names** that appear in 18S census data but don't exist in the NCBI taxonomy database.

### **The Problem:**
- 18S census data contains **122 unmapped taxa** (14.8% of total)
- These include environmental clades like "Dino-Group-II.U.family" and "MAST-12"
- Taxonkit can't find them because they're not in NCBI

### **The Solution:**
- Curated database of **80 known parent mappings**
- Maps unmapped taxa to their closest NCBI parent
- Generates complete lineages by appending to parent lineage

---

## 🗂️ Resolver Components

### **1. `known_parents.py` (293 lines)**
**The curated database of parent taxid mappings**

**Contains:**
- 44 family-level mappings
- 36 genus-level mappings
- 31 unique parent taxa

**Top categories:**
- **Dinoflagellates**: 28 taxa → Dinophyceae (taxid: 2864)
- **Apicomplexans**: 5 taxa → Apicomplexa (taxid: 5794)
- **Stramenopiles**: 4 taxa → Stramenopiles (taxid: 33634)
- **SAR clades**: 4 taxa → Sar (taxid: 2698737)

**Example entries:**
```python
"Dino-Group-II.U.family": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family")
"MAST-12": ("33634", "Stramenopiles", "Marine Stramenopile clade", "family")
"Vermamoebidae": ("554915", "Echinamoebida", "Amoebozoa family", "family")
```

**Functions:**
- `get_parent_info(taxon_name)` - Get parent taxid for a taxon
- `get_statistics()` - Get database statistics

---

### **2. `resolution_builder.py` (167 lines)**
**Builds systematic resolutions from unmapped taxa**

**What it does:**
1. Reads unmapped taxa from taxonkit parser log
2. Looks up each taxon in `KNOWN_PARENTS` database
3. Uses taxonkit to get parent's full lineage
4. Appends unmapped taxon to parent lineage
5. Saves resolutions to JSON file

**Key functions:**
- `build_resolution(taxon_name, parent_taxid, parent_name, rank, env)` - Build single resolution
- `build_all_resolutions(unmapped_log)` - Process all unmapped taxa
- `save_resolutions(resolutions, output_file)` - Save to JSON

**Example resolution:**
```json
{
  "Dino-Group-II.U.family": {
    "lineage": "Eukaryota;Sar;Alveolata;Dinoflagellata;Dinophyceae;Dino-Group-II.U.family",
    "lineage_ranks": "superkingdom;clade;clade;phylum;class;family",
    "lineage_taxids": "2759;2698737;33630;2864;Dino-Group-II.U.family",
    "parent_taxid": "2864",
    "parent_name": "Dinophyceae"
  }
}
```

---

### **3. `resolution_applier.py` (212 lines)**
**Applies resolutions to CSV files**

**What it does:**
1. Loads resolutions from JSON file
2. Reads CSV files from taxonkit parser
3. Finds rows with unmapped taxa
4. Updates lineage fields with resolved values
5. Writes updated CSV files

**Key functions:**
- `load_resolutions(resolution_file)` - Load JSON resolutions
- `apply_resolutions_to_csv(input_csv, output_csv, resolutions, level_name)` - Apply to one CSV
- `create_final_unmapped_log(division_csv, family_csv, genus_csv, output_log)` - Create final unmapped log

**CSV fields updated:**
- `lineage` - Full taxonomic lineage
- `lineage_ranks` - Ranks for each level
- `lineage_taxids` - Taxids for each level

---

## 🔄 Workflow

### **Step 1: Taxonkit Parser**
```bash
python3 run_taxonkit_parser.py
```
**Output:**
- `eukcensus_taxonkit_only_by_division.csv`
- `eukcensus_taxonkit_only_by_family.csv`
- `eukcensus_taxonkit_only_by_genus.csv`
- `eukcensus_taxonkit_only_unmapped.log` ← **122 unmapped taxa**

### **Step 2: Systematic Resolver**
```bash
python3 run_systematic_resolver.py
```
**Process:**
1. Read unmapped log (122 taxa)
2. Build resolutions using `KNOWN_PARENTS` (80 taxa resolved)
3. Apply resolutions to CSV files
4. Create final unmapped log (42 taxa still unmapped)

**Output:**
- `eukcensus_18S_by_division.csv` ← **Final output**
- `eukcensus_18S_by_family.csv` ← **Final output**
- `eukcensus_18S_by_genus.csv` ← **Final output**
- `eukcensus_18S_unmapped_final.log` ← **42 still unmapped**

---

## 📈 Results

### **Before Resolution:**
- Total unmapped: **122 taxa (14.8%)**
  - Division: 2 unmapped
  - Family: 64 unmapped
  - Genus: 56 unmapped

### **After Resolution:**
- Resolved: **80 taxa (65.6% of unmapped)**
- Still unmapped: **42 taxa (5.1% of total)**
  - Family: 20 unmapped
  - Genus: 20 unmapped

---

## 🔍 What You Can Look At

### **1. View the database:**
```bash
cd py_18S
python3 -c "from src.known_parents import KNOWN_PARENTS; print(KNOWN_PARENTS)"
```

### **2. Get statistics:**
```bash
python3 -c "from src.known_parents import get_statistics; import json; print(json.dumps(get_statistics(), indent=2))"
```

### **3. Check unmapped taxa:**
```bash
cat logs/eukcensus_taxonkit_only_unmapped.log | grep "^FAMILY"
cat logs/eukcensus_18S_unmapped_final.log | head -50
```

### **4. View resolutions:**
```bash
# Find the resolutions file (if it exists)
find . -name "systematic_resolutions.json"
```

### **5. Run the resolver:**
```bash
python3 run_18S_pipeline.py --all
```

---

## 🎯 Key Insights

1. **Dinoflagellates dominate** - 28 of 80 resolutions are dinoflagellate clades
2. **Environmental clades** - Most unmapped taxa are environmental/uncultured groups
3. **High success rate** - 65.6% of unmapped taxa can be resolved
4. **Remaining unmapped** - 42 taxa have no known parent (need manual curation)

---

**Want to add more resolutions?** Edit `src/known_parents.py` and add entries to the `KNOWN_PARENTS` dictionary!

