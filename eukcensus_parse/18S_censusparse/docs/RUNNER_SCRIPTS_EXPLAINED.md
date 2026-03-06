# Runner Scripts Explained

## Overview
The 18S parser uses a **3-stage pipeline** to handle the complex naming conventions in EukCensus metadata. Each runner script represents one stage.

---

## 🔷 Script 1: `run_taxonkit_parser.py`

### Purpose
**Pure NCBI taxonomy lookup** - No custom resolutions, just what taxonkit can find.

### What It Does

#### Step 1: Aggregate Data
```python
# Groups clusters by division/family/genus
# Counts: OTU occurrences + sequence counts
process_taxonomic_level(df, 'division', division_data)
process_taxonomic_level(df, 'family', family_data)
process_taxonomic_level(df, 'genus', genus_data)
```

#### Step 2: Get Taxids (4-Tier Fallback)
```python
# Tier 1: Direct lookup - "Dinophyceae" → taxid 2864
# Tier 2: Genus fallback - "Alexandrium tamarense" → "Alexandrium"
# Tier 3: Number stripping - "Dino-Group-II-Clade-7" → "Dino Group II Clade"
# Tier 4: Hyphenated extraction - "MAST-12" → "MAST"
all_taxid_results, all_failed_names = get_taxids_for_names(all_names)
```

#### Step 3: Get Lineages
```python
# For each valid taxid, get full lineage from taxonkit
taxid_to_lineage = get_lineages_for_taxids(all_taxids, env)
```

#### Step 4: Write CSV Files
```python
# Creates 3 files with NCBI-only data:
# - eukcensus_taxonkit_only_by_division.csv
# - eukcensus_taxonkit_only_by_family.csv
# - eukcensus_taxonkit_only_by_genus.csv
```

#### Step 5: Log Unmapped Taxa
```python
# Creates unmapped log for next stage:
# - eukcensus_taxonkit_only_unmapped.log
```

### Output Example
```csv
Name_to_use,taxid,otu_count,size_count,lineage,lineage_ranks,lineage_taxids
Dinophyceae,2864,1234,5678,cellular organisms;Eukaryota;...;Dinophyceae,cellular root;domain;...;class,131567;2759;...;2864
Dino-Group-II.U.family,NA,456,789,,,
```

### Key Point
**Environmental clades get `taxid=NA`** because they're not in NCBI taxonomy.

---

## 🔷 Script 2: `run_systematic_resolver.py`

### Purpose
**Apply curated resolutions** for environmental clades using the `known_parents.py` database.

### What It Does

#### Step 1: Build Resolutions
```python
# Reads unmapped log from Script 1
# For each unmapped name, checks known_parents database
# If found, queries parent taxid to build full lineage

# Example:
# "Dino-Group-II.U.family" → parent_taxid=2864 (Dinophyceae)
# Query taxonkit for lineage of 2864
# Append "Dino-Group-II.U.family" to end of lineage
```

#### Step 2: Apply Resolutions to CSV
```python
# Reads: eukcensus_taxonkit_only_by_*.csv
# Updates entries with systematic resolutions
# Writes: eukcensus_18S_by_*.csv (final output)
```

#### Step 3: Create Final Unmapped Log
```python
# Logs taxa that are STILL unmapped after systematic resolution
# These need division context (Script 3)
```

#### Step 4: Clean Up Intermediate Files
```python
# Deletes eukcensus_taxonkit_only_* files
# Keeps only final eukcensus_18S_* files
```

### Output Example
```csv
Name_to_use,taxid,otu_count,size_count,lineage,lineage_ranks,lineage_taxids
Dino-Group-II.U.family,NA,456,789,cellular organisms;Eukaryota;...;Dinophyceae;Dino-Group-II.U.family,cellular root;domain;...;class;family,131567;2759;...;2864;NA
```

### Key Point
**Environmental clades now have lineages** built from their known parent taxids.

---

## 🔷 Script 3: `run_division_context_adder.py`

### Purpose
**Add minimal taxonomic context** for entries that are STILL unmapped after systematic resolution.

### What It Does

#### Step 1: Load Division Mapping
```python
# Reads raw census file: eukcensus_18S.clusters.97.tsv
# Builds mapping: family/genus → division
# Example: "WIM80-lineage" → "Evosea"
```

#### Step 2: Process CSV Files
```python
# For each entry with taxid=NA and no lineage:
#   1. Look up division from raw metadata
#   2. Add division as parent in lineage
#   3. Update lineage_ranks and lineage_taxids

# Before:
# WIM80-lineage | WIM80-lineage | original_name | NA

# After:
# WIM80-lineage | Evosea;WIM80-lineage | division;family | NA;NA
```

#### Step 3: Write Updated CSV
```python
# Creates files with '_with_division_context' suffix:
# - eukcensus_18S_by_family_with_division_context.csv
# - eukcensus_18S_by_genus_with_division_context.csv
# - eukcensus_18S_by_division_with_division_context.csv
```

### Output Example
```csv
Name_to_use,taxid,otu_count,size_count,lineage,lineage_ranks,lineage_taxids
WIM80-lineage,NA,12,34,Evosea;WIM80-lineage,division;family,NA;NA
```

### Key Point
**Provides SOME context** even when we can't build a full lineage.

---

## 📊 Complete Pipeline Flow

```
INPUT: eukcensus_18S.clusters.97.tsv (raw metadata)
  ↓
┌─────────────────────────────────────────────────────────┐
│ Script 1: run_taxonkit_parser.py                       │
│ - Pure NCBI taxonomy lookup                            │
│ - 4-tier fallback system                               │
│ OUTPUT: eukcensus_taxonkit_only_by_*.csv (NCBI only)   │
│         eukcensus_taxonkit_only_unmapped.log           │
└─────────────────────────────────────────────────────────┘
  ↓
┌─────────────────────────────────────────────────────────┐
│ Script 2: run_systematic_resolver.py                   │
│ - Apply known_parents database                         │
│ - Build lineages from parent taxids                    │
│ OUTPUT: eukcensus_18S_by_*.csv (with resolutions)      │
│         eukcensus_18S_unmapped_final.log               │
└─────────────────────────────────────────────────────────┘
  ↓
┌─────────────────────────────────────────────────────────┐
│ Script 3: run_division_context_adder.py                │
│ - Add division context from raw metadata               │
│ - Minimal context for remaining unmapped               │
│ OUTPUT: eukcensus_18S_by_*_with_division_context.csv   │
└─────────────────────────────────────────────────────────┘
  ↓
FINAL OUTPUT: CSV files with maximum possible lineage coverage
```

---

## 🎯 Success Metrics

After all 3 scripts:
- **~95%** of entries have full NCBI lineages (Script 1)
- **~3-4%** have systematic resolutions (Script 2)
- **~1-2%** have division context only (Script 3)
- **<1%** remain completely unmapped

---

## 💡 Key Insight

The 3-script design reflects the **3 types of naming conventions** in EukCensus:
1. **Standard NCBI names** → Script 1 handles
2. **Environmental clades with known parents** → Script 2 handles
3. **Completely novel names** → Script 3 provides minimal context

This is why reorganizing `src/` into `core/`, `processing/`, and `resolution/` makes sense!

