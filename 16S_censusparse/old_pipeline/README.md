# Old Pipeline Scripts

This directory contains legacy parser scripts that have been replaced by the modular `src/` architecture.

## Deprecated Scripts

### 16S_eukcensus_parser.py
- **Status:** Deprecated
- **Replaced by:** `src/run_census_parser.py` + modular components
- **Reason:** Monolithic script, hard to maintain
- **Date deprecated:** 2026-03-06

### 16S_eukcensus_parser_01.py
- **Status:** Deprecated
- **Replaced by:** `src/run_census_parser.py` + modular components
- **Reason:** Variant of the main parser, functionality merged into modular system
- **Date deprecated:** 2026-03-06

### validate_data_integrity.py
- **Status:** Deprecated
- **Replaced by:** Integrated validation in `src/level_processor.py`
- **Reason:** Validation now built into the main pipeline
- **Date deprecated:** 2026-03-06

## Current Pipeline

Use the new modular pipeline instead:

```bash
# Run the new pipeline
python run_16S_pipeline.py
```

The new pipeline uses:
- `src/` - Modular components
- `run_16S_pipeline.py` - Main entry point
- `output/` - CSV outputs
- `logs/` - Processing logs
- `cache/` - AI resolution cache

## Migration Notes

The old scripts are kept here for reference only. All functionality has been migrated to the new modular architecture with improvements:

1. **Better organization:** Separate modules for each function
2. **Easier testing:** Each module can be tested independently
3. **Better logging:** Comprehensive logging system
4. **AI resolution:** Integrated AI-assisted taxonomic resolution
5. **Systematic resolution:** Manual database for known taxa
6. **Reproducibility:** Full environment specification

## Do Not Use

These scripts are **not maintained** and may not work with the current directory structure.

