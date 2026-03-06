# Source Code Organization Analysis

## Overview
The `src/` directory contains 12 modules organized into functional groups. This analysis examines the current organization and identifies opportunities for improvement.

---

## Current Module Inventory

### 1. **Core Infrastructure** (3 modules)
- `config.py` - Directory paths, environment setup, logging configuration
- `taxonkit_utils.py` - Taxonkit command wrappers, 4-tier fallback system
- `unmapped_logger.py` - Logging for unmapped taxa

### 2. **Data Processing** (2 modules)
- `level_processor.py` - Taxonomic level aggregation, CSV writing
- `lineage_processor.py` - Lineage manipulation, CSV field cleaning

### 3. **Name Cleaning & Validation** (2 modules)
- `taxon_cleaner.py` - Name cleaning (underscores, numbers, special chars)
- `taxon_validator.py` - Filtering logic (empty, null, invalid entries)

### 4. **Resolution Systems** (5 modules)
These handle the **Dinoflagellate problem** and other naming convention issues:

#### a) **Systematic Resolution** (3 modules)
- `known_parents.py` - **DATABASE**: Curated mappings of non-NCBI names → parent taxids
- `resolution_builder.py` - Builds resolutions from unmapped log using known_parents
- `resolution_applier.py` - Applies resolutions to CSV files

#### b) **Division-Based Resolution** (2 modules)
- `division_lineage_inferrer.py` - Infers lineage from division field in raw metadata
- `division_context_adder.py` - Adds division as minimal context for unmapped entries

---

## The Dinoflagellate Problem & Resolution Strategy

### Problem
Original EukCensus metadata used **environmental clade names** not in NCBI taxonomy:
- `Dino-Group-II.U.family`
- `Dino-Group-I-Clade-5`
- `MAST-12` (Marine Stramenopiles)
- `WIM80-lineage` (Amoebozoa)

These names have **NO taxids** in NCBI, so taxonkit returns `NA`.

### Current Solution (3-Layer Approach)

#### **Layer 1: Taxonkit Parser** (`run_taxonkit_parser.py`)
- Uses `taxonkit_utils.py` with 4-tier fallback:
  1. Direct lookup
  2. Genus fallback
  3. Number stripping
  4. Hyphenated pattern extraction
- **Result**: Pure NCBI taxonomy, many environmental clades unmapped

#### **Layer 2: Systematic Resolver** (`run_systematic_resolver.py`)
- Uses `known_parents.py` database (294 entries!)
- Maps environmental clades to their **known parent taxids**
- Example: `Dino-Group-II.U.family` → parent taxid `2864` (Dinophyceae)
- Builds full lineage by querying parent taxid with taxonkit
- **Result**: Recovers most environmental clades with proper lineages

#### **Layer 3: Division Context Adder** (`run_division_context_adder.py`)
- For remaining unmapped entries
- Looks up division field from raw metadata
- Adds division as minimal taxonomic context
- Example: `WIM80-lineage` → adds `Evosea` as parent
- **Result**: Provides some context for otherwise isolated entries

---

## Critical Observation: Missing Integration

### **The Gap You Identified**
The taxid fallback resolver we created for the NCBI species grouper (which searches the assembly file for alternative taxids) **should have been integrated into the 18S parser workflow**.

### Where It Should Go
**Option 1**: New module in `src/`
```
src/metadata_lineage_resolver.py
```
- Searches raw census metadata for lineage information
- Used when taxonkit returns NA but metadata has division/family/genus context
- Similar to `taxid_fallback_resolver.py` from NCBI species grouper

**Option 2**: Extend existing modules
- Add to `division_lineage_inferrer.py` (already does metadata lookup)
- Or add to `taxonkit_utils.py` as Tier 5 fallback

### Why This Matters
Currently, the workflow is:
1. Taxonkit lookup → NA
2. Check known_parents database → Found or Not Found
3. If not found, add division context (minimal)

**Better workflow** would be:
1. Taxonkit lookup → NA
2. **Search raw metadata for any lineage info** ← MISSING!
3. Check known_parents database
4. Add division context as last resort

---

## Reorganization Recommendations

### **Proposed Structure**

```
src/
├── core/
│   ├── config.py
│   ├── taxonkit_utils.py
│   └── unmapped_logger.py
│
├── processing/
│   ├── level_processor.py
│   ├── lineage_processor.py
│   ├── taxon_cleaner.py
│   └── taxon_validator.py
│
└── resolution/
    ├── known_parents.py              # Database
    ├── systematic_resolver.py         # Builder + Applier combined
    ├── division_resolver.py           # Inferrer + Adder combined
    └── metadata_lineage_resolver.py   # NEW - searches raw metadata
```

### Benefits
1. **Clear separation**: Core, Processing, Resolution
2. **Easier to find**: Related modules grouped together
3. **Better for new features**: Clear where to add metadata resolver
4. **Matches workflow**: Resolution modules match the 3-layer approach

---

## Next Steps

1. **Create `metadata_lineage_resolver.py`**
   - Search raw census file for lineage information
   - Integrate into taxonkit_utils as Tier 5 fallback
   - Or call separately between taxonkit and systematic resolver

2. **Reorganize into subdirectories**
   - `core/`, `processing/`, `resolution/`
   - Update imports in runner scripts

3. **Consolidate resolution modules**
   - Combine builder + applier
   - Combine inferrer + adder
   - Reduce from 5 modules to 3

4. **Document the workflow**
   - Create flowchart showing all resolution layers
   - Document when each layer is used

