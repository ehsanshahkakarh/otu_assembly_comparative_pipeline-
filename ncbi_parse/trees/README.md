# Sanity Check Pipeline

This folder contains analysis scripts to validate and analyze the species-grouped output dataset.

## Structure

```
sanity_check/
├── README.md                    # This file
├── domain_parser.py             # Parse and analyze by domain (Bacteria, Archaea, Eukaryota, Viruses)
├── genus_extractor.py           # Extract and analyze genus-level statistics
├── taxonomic_coverage.py        # Analyze taxonomic coverage at all levels
└── run_all_checks.py            # Run all sanity checks
```

## Usage

### Run individual checks:

```bash
# Parse by domain
python sanity_check/domain_parser.py

# Extract genus statistics
python sanity_check/genus_extractor.py

# Check taxonomic coverage
python sanity_check/taxonomic_coverage.py
```

### Run all checks:

```bash
python sanity_check/run_all_checks.py
```

## Output

All analysis results are saved to `sanity_check/output/` with timestamped filenames.

