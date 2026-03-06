# 16S Census Parser Reorganization Summary

**Date:** 2026-03-06  
**Status:** ✅ COMPLETE

---

## Overview

The 16S census parser has been reorganized to match the clean, modular structure of the 18S parser. This improves maintainability, reproducibility, and consistency across the project.

---

## New Directory Structure

```
16S_censusparse/
├── run_16S_pipeline.py          # Main entry point (NEW)
├── src/                          # Core modules (NEW)
│   ├── __init__.py
│   ├── ai_resolution_cache.py
│   ├── cached_ai_resolver.py
│   ├── config.py
│   ├── known_parents.py
│   ├── level_processor.py
│   ├── lineage_processor.py
│   ├── organelle_handler.py
│   ├── resolution_applier.py
│   ├── run_census_parser.py
│   ├── taxon_cleaner.py
│   ├── taxonkit_utils.py
│   └── unmapped_logger.py
├── output/                       # CSV outputs (NEW)
│   ├── eukcensus16S_by_division.csv
│   ├── eukcensus16S_by_family.csv
│   └── eukcensus16S_by_genus.csv
├── logs/                         # Processing logs (NEW)
│   ├── eukcensus16S_processing.log
│   └── eukcensus16S_comprehensive_unmapped.log
├── cache/                        # AI resolution cache (NEW)
│   └── ai_resolutions_16S.json
├── docs/                         # Documentation (NEW)
│   ├── ENHANCED_PARSER_IMPROVEMENTS.md
│   ├── ENHANCED_PARSER_README.md
│   └── RESOLUTIONS_ADDED.md
├── old_pipeline/                 # Deprecated scripts (NEW)
│   ├── README.md
│   ├── 16S_eukcensus_parser.py
│   ├── 16S_eukcensus_parser_01.py
│   └── validate_data_integrity.py
├── metadata/                     # Input data
│   └── eukcensus_16S.clusters.97.tsv
├── requirements.txt
├── environment.yml
├── .gitignore                    # Updated
└── README.md
```

---

## Changes Made

### 1. Created New Directories
- ✅ `src/` - Modular components (moved from `py_16S/census_parser/`)
- ✅ `output/` - CSV outputs (replaces `csv_16S/`)
- ✅ `logs/` - Processing logs (replaces `py_16S/logs/`)
- ✅ `cache/` - AI resolution cache
- ✅ `docs/` - Documentation files
- ✅ `old_pipeline/` - Deprecated scripts

### 2. Moved Files
- ✅ `py_16S/census_parser/*` → `src/`
- ✅ `py_16S/16S_eukcensus_parser.py` → `old_pipeline/`
- ✅ `py_16S/16S_eukcensus_parser_01.py` → `old_pipeline/`
- ✅ `py_16S/validate_data_integrity.py` → `old_pipeline/`
- ✅ `py_16S/ENHANCED_PARSER_*.md` → `docs/`
- ✅ `csv_16S/*.csv` → `output/`
- ✅ `py_16S/logs/*.log` → `logs/`

### 3. Created New Entry Point
- ✅ `run_16S_pipeline.py` - Top-level runner script (replaces `py_16S/run_16S_parser.py`)

### 4. Updated Configuration
- ✅ `src/config.py` - Updated paths to use new directory structure
  - `csv_output_dir`: `csv_16S/` → `output/`
  - `log_dir`: `py_16S/logs/` → `logs/`
  - Added `cache_dir`: `cache/`

### 5. Updated .gitignore
- ✅ Added `py_16S/` (deprecated)
- ✅ Added `csv_16S/` (deprecated)
- ✅ Added `logs/*.log`
- ✅ Added `cache/`
- ✅ Added `__pycache__/`

---

## Usage

### Old Way (Deprecated)
```bash
cd py_16S
python run_16S_parser.py
```

### New Way ✅
```bash
python run_16S_pipeline.py
```

---

## Benefits

1. **Cleaner Structure:** Matches 18S parser organization
2. **Better Separation:** Code, data, logs, and docs in separate directories
3. **Easier Navigation:** Flat structure at root level
4. **Consistent Naming:** Both pipelines use `run_*_pipeline.py`
5. **Deprecated Code Isolated:** Old scripts in `old_pipeline/` with README
6. **Git-Friendly:** Proper .gitignore for temporary files

---

## Backward Compatibility

The old `py_16S/` directory structure is **deprecated** but still present for reference. All functionality has been migrated to the new structure.

**Do not use the old scripts** - they are not maintained and may not work with the new directory structure.

---

## Testing

The reorganized pipeline has been tested and confirmed working:

```bash
$ python run_16S_pipeline.py
================================================================================
16S CENSUS PARSER
================================================================================

✅ Successfully processed 287,468 entries
✅ Generated 3 CSV output files
✅ Created comprehensive unmapped log
✅ Applied 10 systematic resolutions
```

---

## Next Steps

1. ✅ Test the new structure
2. ✅ Verify all outputs are correct
3. ⏳ Commit changes to git
4. ⏳ Update main README to reference new structure
5. ⏳ Remove old `py_16S/` directory after confirmation

---

## Comparison with 18S Structure

Both pipelines now have identical organization:

| Component | 16S | 18S |
|-----------|-----|-----|
| Entry point | `run_16S_pipeline.py` | `run_18S_pipeline.py` |
| Core modules | `src/` | `src/` |
| Outputs | `output/` | `output/` |
| Logs | `logs/` | `logs/` |
| Cache | `cache/` | `cache/` |
| Docs | `docs/` | `docs/` |
| Old code | `old_pipeline/` | N/A |

---

## Summary

✅ **Reorganization complete!**  
✅ **Pipeline tested and working!**  
✅ **Structure matches 18S parser!**  
✅ **Old code preserved in old_pipeline/!**  
✅ **Ready for git commit!**

