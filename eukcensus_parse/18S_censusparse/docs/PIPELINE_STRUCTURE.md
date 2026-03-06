# 18S Census Parser - Pipeline Structure

## Overview

The 18S census parser has **one main runner script** and **three individual step scripts**.

## Main Runner (Use This!)

**`run_complete_pipeline.py`** - Runs all 3 steps in sequence
- ✅ **Recommended for normal use**
- Runs: Taxonkit Parser → Systematic Resolver → Division Context Adder
- Single command: `python run_complete_pipeline.py`

## Individual Step Scripts (Advanced Use Only)

These are called by `run_complete_pipeline.py` but can be run individually if needed:

### 1. `run_taxonkit_parser.py`
- **What it does:** NCBI taxonomy lookups using taxonkit
- **Input:** `metadata/eukcensus_18S.clusters.97.tsv`
- **Output:** 
  - `csv_outputs/eukcensus_taxonkit_only_by_family.csv`
  - `csv_outputs/eukcensus_taxonkit_only_by_genus.csv`
  - `csv_outputs/eukcensus_taxonkit_only_by_division.csv`
  - `logs/eukcensus_taxonkit_only_unmapped.log`

### 2. `run_systematic_resolver.py`
- **What it does:** Applies custom resolutions for 21 families not in NCBI
- **Input:** Output from Step 1
- **Output:**
  - `csv_outputs/eukcensus_18S_by_family.csv`
  - `csv_outputs/eukcensus_18S_by_genus.csv`
  - `csv_outputs/eukcensus_18S_by_division.csv`
  - `logs/eukcensus_18S_unmapped_final.log`

### 3. `run_division_context_adder.py`
- **What it does:** Adds division context to unmapped entries with no parent lineage
- **Input:** Output from Step 2
- **Output:**
  - `csv_outputs/eukcensus_18S_by_family_with_division_context.csv`
  - `csv_outputs/eukcensus_18S_by_genus_with_division_context.csv`
  - `csv_outputs/eukcensus_18S_by_division_with_division_context.csv`

## Legacy Scripts

### `run_18S_pipeline.py`
- **Status:** Legacy orchestrator (runs Steps 1+2 only, not Step 3)
- **Use:** Only if you want to skip the division context adder
- **Recommendation:** Use `run_complete_pipeline.py` instead

## Why Multiple Scripts?

**Historical reasons:** The pipeline was built incrementally:
1. First: Taxonkit parser was created
2. Then: Systematic resolver was added
3. Recently: Division context adder was added

**Current recommendation:** Just use `run_complete_pipeline.py` for everything.

## File Organization

```
py_18S/
├── run_complete_pipeline.py      ← USE THIS (runs all 3 steps)
├── run_taxonkit_parser.py         ← Step 1 (can run individually)
├── run_systematic_resolver.py     ← Step 2 (can run individually)
├── run_division_context_adder.py  ← Step 3 (can run individually)
├── run_18S_pipeline.py            ← Legacy (Steps 1+2 only)
├── src/                           ← Core modules used by all scripts
│   ├── taxonkit_utils.py          ← Taxonkit integration
│   ├── level_processor.py         ← CSV writing
│   ├── resolution_applier.py      ← Systematic resolutions
│   ├── division_context_adder.py  ← Division context logic
│   └── ...
└── csv_outputs/                   ← Final output files
```

## Quick Reference

| Task | Command |
|------|---------|
| Run complete pipeline | `python run_complete_pipeline.py` |
| Run only taxonkit parser | `python run_taxonkit_parser.py` |
| Run only systematic resolver | `python run_systematic_resolver.py` |
| Run only division context adder | `python run_division_context_adder.py` |
| Run Steps 1+2 only (legacy) | `python run_18S_pipeline.py --all` |

## Recent Fix (2026-03-02)

**Issue:** Taxonkit was returning multiple taxids for some names (e.g., "Craniata" → 89593 and 115366), causing corrupted CSV output with embedded newlines in the taxid field.

**Fix:** Updated `src/taxonkit_utils.py` to detect and skip multi-match results, preventing newlines from being embedded in taxid fields.

**Files modified:**
- `src/taxonkit_utils.py` - Added newline detection in taxid parsing

