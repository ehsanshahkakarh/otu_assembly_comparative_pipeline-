# 18S Census Parser - Clean Architecture

## Overview

This is a **clean separation** of the 18S census parser into two independent systems:

1. **Taxonkit Parser** - Pure NCBI taxonomy lookups (no custom resolutions)
2. **Systematic Resolver** - Custom resolutions for families not in NCBI

This architecture eliminates the issues with mixed implementations and Python bytecode caching.

## Directory Structure

```
py_18S/
├── run_18S_pipeline.py           # Main orchestrator - run this!
│
├── taxonkit_parser/               # Clean NCBI-only parser
│   ├── __init__.py
│   ├── config.py                  # Paths and environment setup
│   ├── taxon_cleaner.py           # Name cleaning utilities
│   ├── taxon_validator.py         # Validation rules
│   ├── taxonkit_utils.py          # Taxonkit subprocess calls
│   ├── lineage_processor.py       # Lineage manipulation
│   ├── level_processor.py         # CSV writing (NO resolutions)
│   ├── unmapped_logger.py         # Simple unmapped log
│   └── run_taxonkit_parser.py     # Entry point
│
├── systematic_resolver/           # Custom resolution system
│   ├── __init__.py
│   ├── known_parents.py           # Database: 21 families → parent taxids
│   ├── resolution_builder.py      # Build resolutions from database
│   ├── resolution_applier.py      # Apply resolutions to CSV files
│   ├── run_systematic_resolver.py # Entry point
│   └── outputs/
│       └── systematic_resolutions.json
│
├── archive_current_src/           # Backup of old mixed implementation
│   └── census_parser_18S/
│
└── resolution_tools/              # Research and analysis tools
    ├── systematic_family_resolver.py
    └── analyze_unmapped_families.py
```

## How It Works

### Step 1: Clean Taxonkit Parser

**What it does:**
- Reads `eukcensus_18S.clusters.97.tsv`
- Uses **ONLY** taxonkit for NCBI taxonomy lookups
- Writes CSV files with `NA` for unmapped families
- Creates unmapped log showing what NCBI doesn't have

**Output files:**
- `csv_outputs/eukcensus_taxonkit_only_by_division.csv`
- `csv_outputs/eukcensus_taxonkit_only_by_family.csv`
- `csv_outputs/eukcensus_taxonkit_only_by_genus.csv`
- `logs/eukcensus_taxonkit_only_unmapped_from_taxonkit.log`

**Run it:**
```bash
cd 18S_censusparse/py_18S
python run_18S_pipeline.py --step taxonkit
```

### Step 2: Systematic Resolver

**What it does:**
- Reads unmapped families from taxonkit log
- Looks up parent taxids in `known_parents.py` database
- Uses taxonkit to get parent lineages
- Appends family names to parent lineages
- Applies resolutions to CSV files

**Known Parents Database (21 families):**
- **13 Dinoflagellate families** → Dinophyceae (taxid: 2864)
- **2 MAST clades** → Stramenopiles (taxid: 33634)
- **6 other families** → Various parent taxa

**Output files:**
- `systematic_resolver/outputs/systematic_resolutions.json`
- `csv_outputs/eukcensus_18S_by_family.csv` (with resolutions)
- `logs/eukcensus_18S_unmapped_final.log` (only still-unmapped families)

**Run it:**
```bash
cd 18S_censusparse/py_18S
python run_18S_pipeline.py --step resolve
```

## Quick Start

### Run Complete Pipeline (Recommended)

```bash
cd 18S_censusparse/py_18S
python run_18S_pipeline.py --all
```

This runs both steps in sequence:
1. Clean taxonkit parser
2. Systematic resolver

### Run Individual Steps

```bash
# Step 1 only
python run_18S_pipeline.py --step taxonkit

# Step 2 only (requires Step 1 to be run first)
python run_18S_pipeline.py --step resolve
```

## Expected Results

### Before Systematic Resolution
- **64 families unmapped** (not in NCBI taxonomy)
- Includes all 13 dinoflagellate families
- CSV files show `NA` for lineages

### After Systematic Resolution
- **21 families resolved** (13 dinoflagellates + 8 others)
- **43 families still unmapped** (64 - 21)
- CSV files show complete lineages for resolved families
- **3,259 OTUs** (4.6% of dataset) gain taxonomic context

## Key Benefits

✅ **Clean Separation**: Taxonkit parser has NO knowledge of resolutions
✅ **No Caching Issues**: Separate modules = no bytecode conflicts
✅ **Reproducible**: Each step can be run independently
✅ **Debuggable**: Easy to see what taxonkit found vs what we resolved
✅ **Maintainable**: Clear boundaries between NCBI data and custom resolutions
✅ **No Duplicate Names**: Resolution already includes family name, no double-append

## Adding New Resolutions

To add more families to the systematic resolver:

1. Edit `systematic_resolver/known_parents.py`
2. Add entry to `KNOWN_PARENTS` dictionary:
   ```python
   "Your-Family-Name": ("parent_taxid", "Parent Name", "Notes")
   ```
3. Run the resolver:
   ```bash
   python run_18S_pipeline.py --step resolve
   ```

## Troubleshooting

**Problem:** "Input CSV not found"
**Solution:** Run taxonkit parser first: `python run_18S_pipeline.py --step taxonkit`

**Problem:** "No resolutions built"
**Solution:** Check that unmapped log exists and contains families in known_parents database

**Problem:** Resolutions not appearing in CSV
**Solution:** Make sure you're looking at the final output files (eukcensus_18S_*.csv), not the taxonkit-only files

## File Naming Convention

- `eukcensus_taxonkit_only_*` = Output from clean taxonkit parser (Step 1)
- `eukcensus_18S_*` = Final output with resolutions applied (Step 2)

