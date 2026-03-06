# 18S Parser Reorganization Plan

## Goal
Separate the parser into two independent systems:
1. **Clean Taxonkit Parser** - Pure NCBI taxonomy lookup
2. **Systematic Resolution System** - Handles unmapped families with custom resolutions

## New Directory Structure

```
py_18S/
├── archive_current_src/          # Backup of mixed implementation
│   └── census_parser_18S/
├── taxonkit_parser/              # NEW: Clean NCBI-only parser
│   ├── __init__.py
│   ├── config.py                 # Directory paths, logging
│   ├── level_processor.py        # CSV writing (NO resolution integration)
│   ├── lineage_processor.py      # Lineage manipulation
│   ├── taxon_cleaner.py          # Name cleaning
│   ├── taxon_validator.py        # Validation
│   ├── taxonkit_utils.py         # Taxonkit subprocess calls
│   ├── unmapped_logger.py        # Simple unmapped log (NO resolution filtering)
│   └── run_taxonkit_parser.py    # Main entry point
├── systematic_resolver/          # NEW: Resolution system
│   ├── __init__.py
│   ├── resolution_builder.py     # Build resolutions from known_parents
│   ├── resolution_applier.py     # Apply resolutions to unmapped families
│   ├── known_parents.py          # Database of parent taxids
│   ├── outputs/
│   │   └── systematic_resolutions.json
│   └── run_systematic_resolver.py
├── run_18S_pipeline.py           # NEW: Orchestrates both systems
└── resolution_tools/             # Keep existing research tools
    ├── systematic_family_resolver.py
    ├── analyze_unmapped_families.py
    └── outputs/
```

## Workflow

### Step 1: Clean Taxonkit Parse
```bash
python run_18S_pipeline.py --step taxonkit
```
- Reads eukcensus_18S.clusters.97.tsv
- Uses ONLY taxonkit for lookups
- Writes CSV files with NA for unmapped families
- Creates `unmapped_families.log` (simple list)

### Step 2: Systematic Resolution
```bash
python run_18S_pipeline.py --step resolve
```
- Reads `unmapped_families.log`
- Applies systematic resolutions from known_parents database
- Creates `systematic_resolutions.json`
- Creates `resolved_families.csv` and `still_unmapped.log`

### Step 3: Merge Results
```bash
python run_18S_pipeline.py --step merge
```
- Reads original CSV files
- Applies systematic resolutions
- Writes final CSV files with resolved lineages
- Creates final unmapped log

### Or run all steps:
```bash
python run_18S_pipeline.py --all
```

## Benefits

1. **Clean Separation**: Taxonkit parser has NO knowledge of systematic resolutions
2. **Reproducible**: Each step can be run independently
3. **Debuggable**: Easy to see what taxonkit found vs what we resolved
4. **Maintainable**: Clear boundaries between NCBI data and custom resolutions
5. **Testable**: Each system can be tested independently

## Migration Steps

1. ✅ Archive current mixed implementation
2. Create `taxonkit_parser/` with clean NCBI-only code
3. Create `systematic_resolver/` with resolution logic
4. Create `run_18S_pipeline.py` orchestrator
5. Test each step independently
6. Run full pipeline and verify results

