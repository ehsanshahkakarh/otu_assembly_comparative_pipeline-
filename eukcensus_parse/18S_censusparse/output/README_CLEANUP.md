# CSV Files Cleanup - 2026-03-02

## What Was Done

The CSV output directory has been cleaned up and reorganized for the pipeline merger.

### Files Removed
The following redundant files (without division context) were removed:
- ~~`eukcensus_18S_by_division.csv`~~ (old version)
- ~~`eukcensus_18S_by_family.csv`~~ (old version)
- ~~`eukcensus_18S_by_genus.csv`~~ (old version)

### Files Renamed
The `_with_division_context` files were renamed to the standard names:
- `eukcensus_18S_by_division_with_division_context.csv` → **`eukcensus_18S_by_division.csv`**
- `eukcensus_18S_by_family_with_division_context.csv` → **`eukcensus_18S_by_family.csv`**
- `eukcensus_18S_by_genus_with_division_context.csv` → **`eukcensus_18S_by_genus.csv`**

## Current Files (Ready for Pipeline)

### Main Output Files (CLEANED & READY)
1. **`eukcensus_18S_by_division.csv`** - Division-level aggregation with division context
2. **`eukcensus_18S_by_family.csv`** - Family-level aggregation with division context
3. **`eukcensus_18S_by_genus.csv`** - Genus-level aggregation with division context

### Taxonkit-Only Files (Reference)
These contain ONLY entries that were successfully mapped via taxonkit (no NA taxids):
- `eukcensus_taxonkit_only_by_division.csv`
- `eukcensus_taxonkit_only_by_family.csv`
- `eukcensus_taxonkit_only_by_genus.csv`

### Backup Files
All original files before corruption fix are saved with `.backup` extension:
- `eukcensus_18S_by_division_with_division_context.csv.backup`
- `eukcensus_18S_by_family_with_division_context.csv.backup`
- `eukcensus_18S_by_genus_with_division_context.csv.backup`

## File Format

All main CSV files have the following structure:

```csv
Name_to_use,taxid,otu_count,size_count,lineage,lineage_ranks,lineage_taxids
```

### Columns:
- **Name_to_use**: Original taxonomic name from EukCensus
- **taxid**: NCBI taxonomy ID (or "NA" if not found)
- **otu_count**: Number of OTU clusters for this taxon
- **size_count**: Total sequence count (sum of cluster sizes)
- **lineage**: Full taxonomic lineage (semicolon-separated)
- **lineage_ranks**: Ranks for each lineage level (semicolon-separated)
- **lineage_taxids**: Taxids for each lineage level (semicolon-separated)

## Key Features

### Division Context
The main files include division context for unmapped entries. For example:

**Before (no context):**
```csv
WIM80-lineage,NA,12,34,WIM80-lineage,original_name,NA
```

**After (with division context):**
```csv
WIM80-lineage,NA,12,34,Evosea;WIM80-lineage,division;family,NA;NA
```

This provides minimal taxonomic context even for entries that couldn't be mapped to NCBI taxonomy.

### Corruption Fixed
All files have been cleaned of corrupted taxid fields that contained embedded newlines.
See `../CORRUPTION_FIX_REPORT.md` for details.

## Pipeline Compatibility

These files are now ready to be used by the merger pipeline with the standard naming convention:
- `eukcensus_18S_by_division.csv`
- `eukcensus_18S_by_family.csv`
- `eukcensus_18S_by_genus.csv`

The pipeline should be able to process these files without any issues.

## Statistics

- **Division file**: 22 entries
- **Family file**: 314 entries  
- **Genus file**: 491 entries

All entries include proper lineage information and cleaned taxid fields.

