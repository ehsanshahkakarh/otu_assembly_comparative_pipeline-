# 16S Resolution System - Technical Documentation

## Overview

The 16S census parser includes a **systematic resolution system** for prokaryotic taxa (bacteria and archaea) that are not in NCBI taxonomy. This system maps environmental and GTDB-specific names to their NCBI parent taxids.

---

## How It Works

### The Resolution Pipeline

```
Unmapped Taxon
    ↓
Check known_parents.py
    ↓
Found? → Get Parent Taxid → Taxonkit Lineage → Append Taxon → Complete!
    ↓
Not Found? → Log as unmapped → AI-assisted resolution (optional)
```

### Example Resolution

**Input:** `"Lokiarchaeia"` (Asgard archaea, not in NCBI)

**Resolution:**
1. Check `known_parents.py` → Found!
2. Parent taxid: `1935183` (Asgardarchaeota)
3. Taxonkit lookup: Get full lineage for 1935183
4. Append "Lokiarchaeia" to lineage
5. **Result:** Complete lineage with environmental name preserved

---

## The Resolution Database

### Location

<augment_code_snippet path="00_gaps_taxonomic/00parse_database/16S_censusparse/py_16S/census_parser/known_parents.py" mode="EXCERPT">
````python
# Database of taxonomic names mapped to their parent taxids
# Format: {taxon_name: (parent_taxid, parent_name, notes, rank)}
KNOWN_PARENTS = {
    # Asgard Archaea
    "Lokiarchaeia": ("1935183", "Asgardarchaeota", "Asgard archaea phylum", "phylum"),
    
    # CPR (Candidate Phyla Radiation)
    "Microgenomatia": ("1783273", "Patescibacteria", "CPR group", "family"),
}
````
</augment_code_snippet>

### Database Statistics

**Total mappings:** 273 entries

**By rank:**
- **Phylum level:** ~20 mappings
  - Asgard archaea (Lokiarchaeia, Thorarchaeia, etc.)
  - Modern bacterial phyla (Pseudomonadota, Bacillota)
  - Candidate Phyla Radiation (CPR)
  
- **Family level:** ~150 mappings
  - Environmental clades (ABY, WWE, etc.)
  - GTDB-specific families
  - Uncultured lineages
  
- **Genus level:** ~100 mappings
  - Candidatus genera
  - Uncultured genera
  - GTDB-specific names

### Taxonomic Focus

**Prokaryotic diversity:**
- **Bacteria:** ~85% of mappings
  - Pseudomonadota (formerly Proteobacteria)
  - Bacillota (formerly Firmicutes)
  - Bacteroidota (formerly Bacteroidetes)
  - CPR (Candidate Phyla Radiation)
  
- **Archaea:** ~15% of mappings
  - Asgard archaea (Lokiarchaeia, Thorarchaeia, etc.)
  - DPANN superphylum
  - Euryarchaeota variants

---

## Why This System Exists

### Problem: NCBI vs GTDB Taxonomy

**NCBI Taxonomy:**
- Conservative, slow to update
- Uses older phylum names (Proteobacteria, Firmicutes)
- Missing many environmental clades

**GTDB Taxonomy:**
- Modern, frequently updated
- Uses new phylum names (Pseudomonadota, Bacillota)
- Includes environmental clades

**Environmental Data:**
- Uses GTDB-like names
- Includes uncultured lineages
- Contains candidate divisions

**Result:** Many environmental names don't map to NCBI taxids.

### Solution: Curated Mapping Database

The `known_parents.py` database bridges this gap by mapping:
- GTDB names → NCBI parent taxids
- Environmental clades → NCBI parent taxids
- Uncultured lineages → NCBI parent taxids

---

## Resolution Examples

### Example 1: Modern Phylum Name

**Input:** `"Pseudomonadota"` (modern GTDB name)

**Resolution:**
```python
"Pseudomonadota" → Parent: 1224 (Proteobacteria)
```

**Result:**
```
Lineage: Bacteria;Proteobacteria;Pseudomonadota
Ranks:   superkingdom;phylum;phylum
Taxids:  2;1224;1224
```

### Example 2: Asgard Archaea

**Input:** `"Lokiarchaeia"` (Asgard archaea class)

**Resolution:**
```python
"Lokiarchaeia" → Parent: 1935183 (Asgardarchaeota)
```

**Result:**
```
Lineage: Archaea;Asgardarchaeota;Lokiarchaeia
Ranks:   superkingdom;phylum;phylum
Taxids:  2157;1935183;1935183
```

### Example 3: CPR Group

**Input:** `"Microgenomatia"` (CPR family)

**Resolution:**
```python
"Microgenomatia" → Parent: 1783273 (Patescibacteria)
```

**Result:**
```
Lineage: Bacteria;Patescibacteria;Microgenomatia
Ranks:   superkingdom;phylum;family
Taxids:  2;1783273;1783273
```

---

## Reproducibility

### Is the resolution deterministic?

**YES!** ✅

The resolution database is:
- **Part of the code** (known_parents.py)
- **Version controlled** (in git)
- **Curated manually** (no external API calls)

**Result:** Same resolutions on any computer, any time.

### Dependencies

1. **known_parents.py** - ✅ In code
2. **taxonkit** - ✅ In environment.yml
3. **NCBI taxonomy database** - ✅ Auto-downloaded

**No internet required** after initial setup!

---

## Maintenance

### Adding New Resolutions

To add a new resolution:

1. **Identify the taxon** that needs resolution
2. **Find the NCBI parent** using NCBI taxonomy browser
3. **Add to known_parents.py:**

```python
"New_Taxon_Name": ("parent_taxid", "Parent Name", "Description", "rank"),
```

4. **Test the resolution:**

```bash
python run_16S_parser.py
# Check logs for successful resolution
```

### Validation

The system includes validation:
- Checks if parent taxid exists in NCBI
- Verifies lineage can be retrieved
- Logs any resolution failures

---

## Comparison to 18S System

| Feature | 18S (Eukaryotes) | 16S (Prokaryotes) |
|---------|------------------|-------------------|
| **Focus** | Dinoflagellates, stramenopiles | Bacteria, archaea |
| **Mappings** | 80 entries | 273 entries |
| **Main groups** | Environmental eukaryotes | CPR, Asgard, GTDB names |
| **Organelles** | Extensive (chloroplast, mito) | Minimal |
| **Database** | known_parents.py | known_parents.py |
| **Deterministic** | ✅ Yes | ✅ Yes |

---

## Future Enhancements

### Potential Improvements

1. **Automated GTDB sync:**
   - Automatically map GTDB names to NCBI
   - Update database quarterly

2. **Validation suite:**
   - Unit tests for all resolutions
   - Verify parent taxids still exist

3. **Web interface:**
   - Browse resolution database
   - Submit new resolutions

4. **Integration with GTDB:**
   - Direct GTDB taxonomy lookups
   - Hybrid NCBI/GTDB lineages

---

## Summary

The 16S resolution system is:
- ✅ **Deterministic** - Same results every time
- ✅ **Reproducible** - Works on any computer
- ✅ **Comprehensive** - 273 prokaryotic mappings
- ✅ **Maintainable** - Easy to add new resolutions
- ✅ **Well-documented** - Clear examples and usage

**Bottom line:** Bridges the gap between environmental/GTDB names and NCBI taxonomy!

