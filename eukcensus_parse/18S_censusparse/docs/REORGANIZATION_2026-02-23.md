# 18S Census Parser - Reorganization Summary

**Date**: 2026-02-23  
**Status**: ✅ **COMPLETE**

---

## 🎯 Goal

Consolidate all source modules into a single `src/` directory, following the clean modular pattern from `new_merger`, `nev_parse_meth`, and other parsers in the codebase.

---

## 📁 New Directory Structure

```
py_18S/
├── src/                          # ✅ ALL source modules consolidated here
│   ├── __init__.py
│   ├── config.py                 # Configuration and path setup
│   ├── level_processor.py        # Process taxonomic levels
│   ├── lineage_processor.py      # Process lineage strings
│   ├── taxon_cleaner.py          # Clean taxon names
│   ├── taxon_validator.py        # Validate taxon names
│   ├── taxonkit_utils.py         # Taxonkit wrapper functions
│   ├── unmapped_logger.py        # Log unmapped taxa
│   ├── known_parents.py          # Curated parent taxid database
│   ├── resolution_applier.py     # Apply resolutions to CSV
│   └── resolution_builder.py     # Build resolutions from unmapped log
│
├── run_18S_pipeline.py           # Main entry point (orchestrator)
├── run_taxonkit_parser.py        # Taxonkit parser entry point
├── run_systematic_resolver.py    # Systematic resolver entry point
│
├── sanity_checks/                # Verification scripts (unchanged)
├── docs/                         # Documentation (unchanged)
├── logs/                         # Log files (unchanged)
├── archive/                      # Old code (unchanged)
│   ├── legacy_parsers/
│   ├── experimental/
│   └── resolution_tools/         # ✅ Moved here (was not used by pipeline)
└── archive_current_src/          # Previous src backup (unchanged)
```

---

## 🔄 Changes Made

### **1. Created `src/` directory**
- Consolidated all source modules from `taxonkit_parser/` and `systematic_resolver/`
- Created `src/__init__.py` with package metadata

### **2. Moved modules to `src/`**

**From `taxonkit_parser/`:**
- `config.py`
- `level_processor.py`
- `lineage_processor.py`
- `taxon_cleaner.py`
- `taxon_validator.py`
- `taxonkit_utils.py`
- `unmapped_logger.py`

**From `systematic_resolver/`:**
- `known_parents.py`
- `resolution_applier.py`
- `resolution_builder.py`

### **3. Created new entry point scripts**

**`run_taxonkit_parser.py`** - Runs the NCBI taxonkit parser
- Imports from `src.config`, `src.level_processor`, etc.
- Replaces `taxonkit_parser/run_taxonkit_parser.py`

**`run_systematic_resolver.py`** - Runs the systematic family/genus resolver
- Imports from `src.resolution_builder`, `src.resolution_applier`, etc.
- Replaces `systematic_resolver/run_systematic_resolver.py`

### **4. Updated `run_18S_pipeline.py`**
- Changed imports from `taxonkit_parser.run_taxonkit_parser` → `run_taxonkit_parser`
- Changed imports from `systematic_resolver.run_systematic_resolver` → `run_systematic_resolver`

### **5. Archived unused code**
- Moved `resolution_tools/` to `archive/resolution_tools/`
- This directory contained old resolver scripts that were NOT used by the pipeline

### **6. Removed old directories**
- Deleted `taxonkit_parser/` (modules moved to `src/`)
- Deleted `systematic_resolver/` (modules moved to `src/`)

---

## ✅ Verification

### **Import tests:**
```bash
cd py_18S
python3 -c "from src import config; print('✅ src.config works')"
python3 -c "from src.known_parents import get_statistics; print('✅ known_parents works')"
```

### **Pipeline help:**
```bash
python3 run_18S_pipeline.py --help
```

**Output:**
```
usage: run_18S_pipeline.py [-h] [--step {taxonkit,resolve}] [--all]

18S Census Parser Pipeline

options:
  -h, --help            show this help message and exit
  --step {taxonkit,resolve}
                        Run specific pipeline step(s). Can be specified multiple
                        times.
  --all                 Run all pipeline steps in sequence
```

---

## 🚀 Usage

### **Run complete pipeline:**
```bash
python3 run_18S_pipeline.py --all
```

### **Run individual steps:**
```bash
# Step 1: Taxonkit parser (NCBI lookups only)
python3 run_18S_pipeline.py --step taxonkit

# Step 2: Systematic resolver (apply custom resolutions)
python3 run_18S_pipeline.py --step resolve
```

### **Or run scripts directly:**
```bash
python3 run_taxonkit_parser.py
python3 run_systematic_resolver.py
```

---

## 📊 Benefits

1. **✅ Clean structure** - All source code in one place (`src/`)
2. **✅ Consistent with other parsers** - Follows pattern from `new_merger`, `nev_parse_meth`
3. **✅ Easy to navigate** - No confusion between `resolution_tools/` and `systematic_resolver/`
4. **✅ Clear entry points** - Three main scripts at root level
5. **✅ Archives preserved** - Old code safely stored in `archive/`

---

## 🗑️ What Was Removed

- `taxonkit_parser/` directory (modules moved to `src/`)
- `systematic_resolver/` directory (modules moved to `src/`)
- `resolution_tools/` directory (moved to `archive/` - was not used)

---

**Reorganization complete!** The 18S census parser now has a clean, modular structure. 🎉

