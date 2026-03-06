# 18S Census Parse Cleanup Plan
## Date: 2026-03-03

## Current Issues
- Too many markdown files in root (5 MD files)
- Too many docs files (13 MD files in docs/)
- Multiple archive directories (archive/, archive_current_src/)
- Redundant run_* scripts
- Unclear structure

## Proposed Clean Structure

```
18S_censusparse/
├── README.md                          # Main overview
├── run_18S_pipeline.py                # Single entry point
├── example_ai_resolutions.json        # Example data
├── src/                               # All source modules
│   ├── __init__.py
│   ├── ai_cache_cli.py
│   ├── ai_resolution_cache.py
│   ├── cached_ai_resolver.py
│   ├── config.py
│   ├── known_parents.py
│   ├── pipeline_*.py (3 files)
│   └── ... (other modules)
├── cache/                             # AI cache data
│   └── ai_resolutions.json
├── logs/                              # Pipeline logs
│   ├── eukcensus_18S_unmapped_final.log
│   ├── eukcensus_taxonkit_only_unmapped_from_taxonkit.log
│   └── systematic_resolver.log
├── output/                            # Final CSV outputs
│   ├── eukcensus_18S_by_division.csv
│   ├── eukcensus_18S_by_family.csv
│   └── eukcensus_18S_by_genus.csv
├── metadata/                          # Input data
│   └── eukcensus_18S.clusters.97.tsv
└── docs/                              # All documentation
    ├── README.md
    ├── QUICKSTART.md
    ├── AI_RESOLUTION_WORKFLOW.md
    └── ARCHITECTURE.md

18S_censusparse_old/                   # Archive (move to old_pipeline/)
├── archive/
├── archive_current_src/
├── sanity_checks/
├── old_docs/
└── README.md
```

## Cleanup Steps

### 1. Create Archive Directory
```bash
cd 00_gaps_taxonomic/parse_repaa_table/18S_censusparse
mkdir -p 18S_censusparse_old
```

### 2. Move Archives
```bash
cd py_18S
mv archive ../18S_censusparse_old/
mv archive_current_src ../18S_censusparse_old/
mv sanity_checks ../18S_censusparse_old/
```

### 3. Consolidate Documentation
```bash
# Move all docs to docs/ directory
mv AI_CACHE_IMPLEMENTATION_SUMMARY.md docs/
mv AI_RESOLUTION_WORKFLOW.md docs/
mv PIPELINE_STRUCTURE.md docs/
mv QUICKSTART.md docs/

# Keep only main README in root
```

### 4. Clean Up Redundant Scripts
```bash
# Remove redundant runner scripts (keep only run_18S_pipeline.py)
rm run_complete_pipeline.py
rm run_division_context_adder.py
rm run_systematic_resolver.py
rm run_taxonkit_parser.py
rm fix_corrupted_csv.py
```

### 5. Reorganize Outputs
```bash
# Move csv_outputs to parent and rename
cd ..
mv csv_outputs output
```

### 6. Clean Cache Directory
```bash
cd py_18S/cache
# Keep only ai_resolutions.json, remove old logs
rm ai_cache_run*.log
rm ai_resolution_review_*.txt
```

### 7. Move py_18S contents up
```bash
cd ..
# Move everything from py_18S/ to 18S_censusparse/
mv py_18S/src .
mv py_18S/cache .
mv py_18S/logs .
mv py_18S/docs .
mv py_18S/run_18S_pipeline.py .
mv py_18S/example_ai_resolutions.json .
mv py_18S/README.md docs/PIPELINE_README.md
mv py_18S/__pycache__ 18S_censusparse_old/
rmdir py_18S
```

### 8. Final Archive Move
```bash
cd ..
mv 18S_censusparse_old ../old_pipeline/
```

## Final Structure

```
18S_censusparse/
├── README.md                    # Main overview
├── run_18S_pipeline.py          # Single entry point  
├── example_ai_resolutions.json  # Example data
├── src/                         # 17 Python modules
├── cache/                       # ai_resolutions.json only
├── logs/                        # 5 log files
├── output/                      # 3 CSV files
├── metadata/                    # Input TSV
└── docs/                        # ~5 key docs (consolidated)

old_pipeline/18S_censusparse_old/
├── archive/
├── archive_current_src/
├── sanity_checks/
└── old_docs/
```

## Benefits

1. **Flat structure**: No nested py_18S directory
2. **Clear separation**: src/, cache/, logs/, output/, docs/
3. **Single entry point**: run_18S_pipeline.py
4. **Clean docs**: 5 key docs instead of 18
5. **Archived clutter**: All old code in old_pipeline/

## Execution

Run the cleanup script:
```bash
bash cleanup_18S.sh
```

Or execute steps manually as listed above.

