# 18S Cleanup - Quick Start

## TL;DR

```bash
cd 00_gaps_taxonomic/parse_repaa_table/18S_censusparse
bash cleanup_18S.sh
mv 18S_censusparse_old ../old_pipeline/
```

## What It Does

1. **Flattens structure**: Removes nested `py_18S/` directory
2. **Consolidates docs**: Moves all MD files to `docs/`
3. **Removes clutter**: Deletes redundant scripts and old logs
4. **Renames outputs**: `csv_outputs/` → `output/`
5. **Archives old code**: Creates `18S_censusparse_old/`

## Before vs After

| Before | After |
|--------|-------|
| `18S_censusparse/py_18S/src/` | `18S_censusparse/src/` |
| `18S_censusparse/csv_outputs/` | `18S_censusparse/output/` |
| 6 run_* scripts | 1 run_18S_pipeline.py |
| 18 MD files | 5 key docs |
| 2 archive dirs | 0 (moved to old_pipeline) |

## After Cleanup

```bash
# Run pipeline
cd 18S_censusparse
python run_18S_pipeline.py --all

# AI cache
python -m src.ai_cache_cli stats

# View outputs
ls output/
```

## Safety

- ✅ Nothing is deleted, only moved
- ✅ Can rollback from `18S_censusparse_old/`
- ✅ Script exits on any error (`set -e`)

## Verification

```bash
# Check structure
ls -la 18S_censusparse/

# Test pipeline
cd 18S_censusparse
python run_18S_pipeline.py --step ai-cache

# Verify archive
ls -la old_pipeline/18S_censusparse_old/
```

## Files Removed

- `run_complete_pipeline.py` (redundant)
- `run_division_context_adder.py` (redundant)
- `run_systematic_resolver.py` (redundant)
- `run_taxonkit_parser.py` (redundant)
- `fix_corrupted_csv.py` (redundant)
- `cache/ai_cache_run*.log` (old logs)
- `cache/ai_resolution_review_*.txt` (old logs)

## Files Moved to Archive

- `archive/`
- `archive_current_src/`
- `sanity_checks/`
- `__pycache__/`

## Result

**Clean, flat structure matching ncbi_parse:**
```
18S_censusparse/
├── run_18S_pipeline.py    # Single entry point
├── src/                   # All code
├── cache/                 # AI cache
├── logs/                  # Pipeline logs
├── output/                # CSV outputs
└── docs/                  # Documentation
```

## Questions?

See `CLEANUP_SUMMARY.md` for full details.

