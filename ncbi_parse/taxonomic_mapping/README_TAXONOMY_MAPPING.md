# NCBI Taxonomy Mapping Scripts

This directory contains scripts for mapping NCBI taxonomy IDs to their respective phylum, family, and genus levels, along with domain information.

## Scripts Overview

1. `phylum_taxid.py`: Maps taxids to phylum level
2. `family_taxid.py`: Maps taxids to family level
3. `genus_taxid.py`: Maps taxids to genus level

## Required Files

Place the following files in the appropriate directories:

```
ncbi_parse_scripts/taxonomic_mapping/
├── taxdump_ncbi/
│   ├── nodes.dmp
│   └── names.dmp
├── 00assembly_summary_genbank.txt
├── phylum_taxid.py
├── family_taxid.py
└── genus_taxid.py
```

## Output Files

Each script generates a CSV file in the same directory:
- `taxid_to_phylum.csv`
- `taxid_to_family.csv`
- `taxid_to_genus.csv`

## Features

### Common Features Across All Scripts

1. **Path Handling**
   - Uses `pathlib` for cross-platform path handling
   - Automatically locates required files in the correct directories
   - Validates file existence before processing

2. **Taxonomy Processing**
   - Loads and processes NCBI taxonomy files (nodes.dmp and names.dmp)
   - Builds optimized taxonomy maps for efficient lookups
   - Handles missing or invalid taxonomy data gracefully

3. **Domain Assignment**
   - Assigns domains (Eukaryota, Bacteria, Archaea, Viruses)
   - Infers domains from taxonomic names when direct assignment fails
   - Filters out invalid domain assignments

4. **Candidatus Handling**
   - Preserves Candidatus taxa in the final output for complete taxonomic coverage
   - Tracks and reports the number of Candidatus entries found

5. **Special Cases Handling**
   - Tracks special cases for debugging
   - Implements domain inference based on taxonomic name patterns
   - Handles ambiguous domain assignments

### Script-Specific Features

#### phylum_taxid.py
- Handles phylum-level taxonomic assignments
- Implements special domain inference for phylum names:
  - Names ending in "ota" or "mycota" → Eukaryota
  - Names containing "bacter" → Bacteria
  - Names containing "archae" → Archaea

#### family_taxid.py
- Handles family-level taxonomic assignments
- Implements special domain inference for family names:
  - Names ending in "aceae" → Bacteria
  - Names ending in "idae" → Eukaryota

#### genus_taxid.py
- Handles genus-level taxonomic assignments
- Implements special domain inference for genus names:
  - Names ending in "aceae" → Bacteria
  - Names ending in "idae" → Eukaryota

## Output Format

Each script produces a CSV file with the following columns:
- `taxid`: NCBI taxonomy ID
- `[rank]`: The taxonomic rank (phylum/family/genus)
- `domain`: The assigned domain
- `depth`: The depth of traversal in the taxonomy tree

## Terminal Output

All scripts provide consistent, clean terminal output:

```
Loading taxonomy files...
Building taxonomy maps...
Processing assembly file...
Adding domain information...

Processing Statistics:
Total taxids processed: <number>
Successfully mapped: <number>
Candidatus taxa found and preserved: <number>
Output file: <path>
Done.
```

## Error Handling

- Validates input files before processing
- Handles missing taxonomy data gracefully
- Provides clear error messages for missing files
- Implements cycle detection in taxonomy traversal
- Limits traversal depth to prevent infinite loops

## Performance Considerations

- Uses pandas for efficient data processing
- Implements optimized taxonomy maps for fast lookups
- Processes data in chunks to manage memory usage
- Uses numpy for efficient numerical operations

## Recent Updates

1. Standardized terminal output across all scripts
2. Removed hardcoded special cases in favor of pattern-based inference
3. Improved domain inference logic
4. Enhanced error handling and validation
5. Optimized file path handling
6. Added comprehensive statistics reporting
