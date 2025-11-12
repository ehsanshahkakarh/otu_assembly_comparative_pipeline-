# EukProt Taxonomy Sanity Checks

This directory contains scripts and resources for validating and analyzing EukProt taxonomic data.

## Workflows

### 1. Validation Workflow

```
┌─────────────────────────┐
│ INPUT                   │
│ Eukprot_included_datasets.txt │
└───────────┬─────────────┘
            │
            ▼
┌───────────────────────────┐
│ VALIDATION                │
│ validate_eukprot_source.py│
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ OUTPUT                    │
│ eukprot_validated_with_taxids.csv │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ LOGS                      │
│ eukprot_validation.log    │
│ fishy_matches_manual_review.log │
└───────────────────────────┘
```

This workflow validates organism names in the EukProt dataset against NCBI taxonomy, applying various matching strategies and corrections for known misspellings.

### 2. Taxonomic Mapping Workflow

```
┌─────────────────────────┐
│ INPUT                   │
│ Eukprot_included_datasets.txt │
└───────────┬─────────────┘
            │
            ▼
┌───────────────────────────┐
│ EXTRACT NAMES             │
│ extract_names.py          │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ OUTPUT                    │
│ eukprot_names_to_use.csv  │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ MAP TO TAXIDS             │
│ map_to_taxids.py          │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ OUTPUT                    │
│ eukprot_with_taxids.csv   │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ GENERATE LINEAGES         │
│ generate_lineages.py      │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ OUTPUT                    │
│ eukprot_with_lineages.csv │
└───────────────────────────┘
```

This workflow maps EukProt organism names to NCBI taxids and generates complete taxonomic lineages.

### 3. Analysis Workflow

```
┌─────────────────────────┐
│ INPUT                   │
│ eukprot_with_lineages.csv │
└───────────┬─────────────┘
            │
            ▼
┌───────────────────────────┐
│ ANALYZE TAXONOMIC GROUPS  │
│ analyze_taxonomic_groups.py │
└───────────┬───────────────┘
            │
            ▼
┌───────────────────────────┐
│ OUTPUT                    │
│ taxonomic_analysis.md     │
│ taxonomic_visualizations/ │
└───────────────────────────┘
```

This workflow analyzes the taxonomic distribution of EukProt organisms and generates reports and visualizations.

## Script Descriptions

### Validation Scripts

- **`validate_eukprot_source.py`**: Validates organism names against NCBI taxonomy
  - Input: `Eukprot_included_datasets.txt`
  - Output: `eukprot_validated_with_taxids.csv`
  - Logs: `eukprot_validation.log`, `fishy_matches_manual_review.log`

### Mapping Scripts

- **`extract_names.py`**: Extracts organism names from the source file
  - Input: `Eukprot_included_datasets.txt`
  - Output: `eukprot_names_to_use.csv`

- **`map_to_taxids.py`**: Maps organism names to NCBI taxids
  - Input: `eukprot_names_to_use.csv`
  - Output: `eukprot_with_taxids.csv`

- **`generate_lineages.py`**: Generates taxonomic lineages for each taxid
  - Input: `eukprot_with_taxids.csv`
  - Output: `eukprot_with_lineages.csv`

### Analysis Scripts

- **`analyze_taxonomic_groups.py`**: Analyzes taxonomic groups in EukProt
  - Input: `eukprot_with_lineages.csv`
  - Output: `taxonomic_analysis.md`, visualizations

## Migration Guide

To migrate from the current structure to the proposed structure:

1. Create the new directory structure
2. Move files to their appropriate locations:
   - Move `validate_eukprot_source.py` to `validation/`
   - Move `new_thing/extract_name_to_use.py` to `mapping/extract_names.py`
   - Move `new_thing/map_names_to_taxids_parallel.py` to `mapping/parallel/map_to_taxids_parallel.py`
   - Move `new_thing/generate_taxonomic_lineages.py` to `mapping/generate_lineages.py`
   - Move `Taxgroup_UniEuk_parse.py` to `analysis/analyze_taxonomic_groups.py`
   - Move `old/` contents to `archive/old/`
3. Update import statements in each script to reflect the new structure
4. Create utility modules by extracting common functionality
5. Update file paths in scripts to use the new directory structure

## Best Practices

1. **Use descriptive names**: Name directories and files based on their purpose
2. **Separate concerns**: Keep validation, mapping, and analysis scripts separate
3. **Centralize common functionality**: Extract shared code into utility modules
4. **Document workflows**: Create README files for each major component
5. **Version control**: Use git tags to mark stable versions
6. **Archive, don't delete**: Move obsolete scripts to the archive directory
7. **Standardize interfaces**: Use consistent input/output formats
8. **Log everything**: Implement comprehensive logging
9. **Parallelize where possible**: Use parallel processing for performance-critical tasks
10. **Validate inputs and outputs**: Add validation checks at each step

## Proposed Directory Structure

```
sanity_checks/
├── README.md                           # This file
├── validation/                         # Validation scripts
│   ├── validate_eukprot_source.py      # Main validation script
│   ├── verify_taxonomic_ranks.py       # Rank verification
│   └── logs/                           # Validation logs
│       ├── eukprot_validation.log      # Main validation log
│       └── fishy_matches_manual_review.log  # Suspicious matches
├── mapping/                            # Taxonomic mapping scripts
│   ├── extract_names.py                # Extract organism names
│   ├── map_to_taxids.py                # Map names to taxids
│   ├── generate_lineages.py            # Generate taxonomic lineages
│   └── parallel/                       # Parallel processing versions
│       ├── map_to_taxids_parallel.py   # Parallel taxid mapping
│       └── generate_lineages_parallel.py  # Parallel lineage generation
├── analysis/                           # Analysis scripts
│   ├── analyze_taxonomic_groups.py     # Analyze taxonomic groups
│   ├── compare_databases.py            # Compare with other databases
│   └── visualize_taxonomy.py           # Generate taxonomic visualizations
├── data/                               # Input/output data
│   ├── input/                          # Input files
│   │   └── eukprot_names_to_use.csv    # Extracted names
│   ├── output/                         # Output files
│   │   ├── eukprot_with_taxids.csv     # Names with taxids
│   │   └── eukprot_with_lineages.csv   # Complete lineages
│   └── reference/                      # Reference data
│       └── name_corrections.csv        # Known name corrections
├── utils/                              # Utility functions
│   ├── taxon_utils.py                  # Taxonomy utilities
│   └── file_utils.py                   # File handling utilities
└── archive/                            # Archived scripts
    ├── old/                            # Old versions
    │   └── phylum_eukprot_parser.py    # Old parser
    └── experimental/                   # Experimental scripts
        └── new_approach.py             # Experimental approaches
```