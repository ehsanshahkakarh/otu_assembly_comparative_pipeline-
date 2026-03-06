# 18S Census Parser

Three-step pipeline for processing 18S rRNA environmental sequencing data.

## Quick Start

```bash
# Run complete pipeline (all 3 steps)
python run_complete_pipeline.py
```

## Pipeline Steps

### Step 1: Taxonkit Parser
**Purpose:** Pure NCBI taxonomy lookups using taxonkit
**Output:** CSV files with `NA` for unmapped families + unmapped log

### Step 2: Systematic Resolver
**Purpose:** Custom resolutions for 21 families not in NCBI
**Output:** CSV files with resolutions applied + final unmapped log

### Step 3: Division Context Adder
**Purpose:** Add division context to remaining unmapped entries
**Output:** Final CSV files with division context added

## Usage

### Run Complete Pipeline (Recommended)
```bash
# Single command to run all 3 steps
python run_complete_pipeline.py
```

### Run Individual Steps (Advanced)
```bash
# Step 1: Taxonkit parser only
python run_taxonkit_parser.py

# Step 2: Systematic resolver only (requires Step 1 output)
python run_systematic_resolver.py

# Step 3: Division context adder only (requires Step 2 output)
python run_division_context_adder.py

# Or use the old orchestrator for Steps 1+2 only
python run_18S_pipeline.py --all
```

## Output Files

**Final CSV files:**
- `csv_outputs/eukcensus_18S_by_division.csv`
- `csv_outputs/eukcensus_18S_by_family.csv`
- `csv_outputs/eukcensus_18S_by_genus.csv`

**Logs:**
- `logs/eukcensus_18S_unmapped_final.log` - Families still unmapped after resolution

**Note:** The pipeline automatically deletes intermediate `eukcensus_taxonkit_only_*.csv` files after the systematic resolver completes. Only the final `eukcensus_18S_by_*.csv` files should exist in `csv_outputs/`.

## Documentation

See `docs/` directory for detailed documentation:
- `ARCHITECTURE_README.md` - Complete architecture guide
- `REORGANIZATION_COMPLETE.md` - Summary of changes
- `REORGANIZATION_PLAN.md` - Design plan

## Directory Structure

```
py_18S/
├── run_18S_pipeline.py          # Main entry point
├── taxonkit_parser/             # Clean NCBI-only parser
├── systematic_resolver/         # Custom resolution system
├── resolution_tools/            # Research and analysis tools
├── sanity_checks/               # Validation scripts
├── logs/                        # Output logs
├── docs/                        # Documentation
└── archive_current_src/         # Backup of old implementation
```

## Adding New Resolutions

1. Edit `systematic_resolver/known_parents.py`
2. Add family to `KNOWN_PARENTS` dictionary:
   ```python
   "Family-Name": ("parent_taxid", "Parent Name", "Notes")
   ```
3. Run resolver: `python run_18S_pipeline.py --step resolve`

## Results

- **21 families resolved** (13 dinoflagellates + 8 others)
- **~3,259 OTUs** (4.6% of dataset) gain complete taxonomic lineages
- **43 families still unmapped** (down from 64)

