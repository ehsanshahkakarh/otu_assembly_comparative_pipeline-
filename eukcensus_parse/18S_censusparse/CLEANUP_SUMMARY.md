# 18S Census Parse Cleanup Summary
**Date:** 2026-03-03  
**Purpose:** Reorganize to match clean ncbi_parse structure

## Quick Summary

**Before:** Nested, cluttered structure with 18+ markdown files, multiple archives, redundant scripts  
**After:** Flat, clean structure with single entry point, organized directories

## Execution

```bash
cd 00_gaps_taxonomic/parse_repaa_table/18S_censusparse
bash cleanup_18S.sh
mv 18S_censusparse_old ../old_pipeline/
```

## Before → After

### Directory Structure

**BEFORE:**
```
18S_censusparse/
├── csv_outputs/                    # Outputs
│   ├── archive/
│   ├── archive_before_integration/
│   └── 3 CSV files
├── metadata/                       # Input data
├── py_18S/                         # ← NESTED STRUCTURE
│   ├── AI_CACHE_IMPLEMENTATION_SUMMARY.md
│   ├── AI_RESOLUTION_WORKFLOW.md
│   ├── PIPELINE_STRUCTURE.md
│   ├── QUICKSTART.md
│   ├── README.md
│   ├── archive/                    # Old code
│   ├── archive_current_src/        # More old code
│   ├── cache/                      # + 8 old log files
│   ├── docs/                       # 13 MD files!
│   ├── logs/
│   ├── sanity_checks/
│   ├── src/                        # Actual code
│   ├── run_18S_pipeline.py
│   ├── run_complete_pipeline.py    # Redundant
│   ├── run_division_context_adder.py  # Redundant
│   ├── run_systematic_resolver.py  # Redundant
│   ├── run_taxonkit_parser.py      # Redundant
│   └── fix_corrupted_csv.py        # Redundant
└── sanity_check_parser.py
```

**AFTER:**
```
18S_censusparse/
├── README.md                       # Main overview
├── run_18S_pipeline.py             # Single entry point
├── example_ai_resolutions.json     # Example data
├── src/                            # All source code (17 modules)
│   ├── ai_cache_cli.py
│   ├── ai_resolution_cache.py
│   ├── cached_ai_resolver.py
│   ├── config.py
│   ├── known_parents.py
│   └── ... (12 more modules)
├── cache/                          # Clean cache
│   └── ai_resolutions.json
├── logs/                           # Pipeline logs (5 files)
│   ├── eukcensus_18S_unmapped_final.log
│   ├── eukcensus_taxonkit_only_unmapped_from_taxonkit.log
│   └── systematic_resolver.log
├── output/                         # Final outputs (renamed from csv_outputs)
│   ├── eukcensus_18S_by_division.csv
│   ├── eukcensus_18S_by_family.csv
│   └── eukcensus_18S_by_genus.csv
├── metadata/                       # Input data
│   └── eukcensus_18S.clusters.97.tsv
└── docs/                           # Consolidated docs (~5 key files)
    ├── README.md
    ├── QUICKSTART.md
    ├── AI_RESOLUTION_WORKFLOW.md
    ├── ARCHITECTURE_README.md
    └── CLEANUP_PLAN.md

old_pipeline/18S_censusparse_old/   # Archived
├── archive/
├── archive_current_src/
├── sanity_checks/
└── README.md
```

## Changes Made

### 1. Flattened Structure
- ❌ Removed nested `py_18S/` directory
- ✅ Moved all contents to `18S_censusparse/` root

### 2. Consolidated Documentation
- ❌ 18 markdown files scattered across root and docs/
- ✅ ~5 key docs in `docs/` directory
- Moved: AI_CACHE_IMPLEMENTATION_SUMMARY.md, AI_RESOLUTION_WORKFLOW.md, PIPELINE_STRUCTURE.md, QUICKSTART.md

### 3. Removed Redundant Scripts
Deleted:
- `run_complete_pipeline.py`
- `run_division_context_adder.py`
- `run_systematic_resolver.py`
- `run_taxonkit_parser.py`
- `fix_corrupted_csv.py`

Kept:
- `run_18S_pipeline.py` (single entry point)

### 4. Cleaned Cache Directory
- ❌ 8 old log files (ai_cache_run*.log, ai_resolution_review_*.txt)
- ✅ Only `ai_resolutions.json`

### 5. Renamed Outputs
- ❌ `csv_outputs/` (unclear name)
- ✅ `output/` (clear, matches ncbi_parse)

### 6. Archived Old Code
Moved to `18S_censusparse_old/`:
- `archive/`
- `archive_current_src/`
- `sanity_checks/`
- `__pycache__/`

## File Count Comparison

| Category | Before | After | Change |
|----------|--------|-------|--------|
| Root MD files | 5 | 1 | -4 |
| Docs MD files | 13 | 5 | -8 |
| Run scripts | 6 | 1 | -5 |
| Cache files | 9 | 1 | -8 |
| Archive dirs | 2 | 0 | -2 |
| **Total clutter** | **35** | **8** | **-27** |

## Benefits

1. ✅ **Flat structure**: No confusing nested directories
2. ✅ **Single entry point**: `run_18S_pipeline.py`
3. ✅ **Clear organization**: src/, cache/, logs/, output/, docs/
4. ✅ **Matches ncbi_parse**: Consistent structure across projects
5. ✅ **Less clutter**: 27 fewer files/directories
6. ✅ **Easier navigation**: Everything in logical places

## Usage After Cleanup

```bash
# Main pipeline
cd 18S_censusparse
python run_18S_pipeline.py --all

# AI cache management
python -m src.ai_cache_cli stats
python -m src.ai_cache_cli add
python -m src.ai_cache_cli import example_ai_resolutions.json

# View outputs
ls output/

# View logs
ls logs/

# Read docs
ls docs/
```

## Verification

After running cleanup, verify:
```bash
# Check structure
ls -la 18S_censusparse/

# Verify pipeline still works
cd 18S_censusparse
python run_18S_pipeline.py --step ai-cache

# Check archive was created
ls -la ../old_pipeline/18S_censusparse_old/
```

## Rollback (if needed)

If something goes wrong:
```bash
# The script doesn't delete anything, just moves
# You can restore from 18S_censusparse_old/
```

## Status

- [x] Cleanup script created
- [ ] Script executed
- [ ] Archive moved to old_pipeline/
- [ ] Verified pipeline works
- [ ] Updated main README

## Next Steps

1. Run `bash cleanup_18S.sh`
2. Verify pipeline works
3. Move archive: `mv 18S_censusparse_old ../old_pipeline/`
4. Update main README.md
5. Commit changes

