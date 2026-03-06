# Enhanced EukCensus Parser with Organelle Handling

## Overview

The enhanced EukCensus parser addresses key issues with organellar sequence entries and taxonomic rank filtering. This parser provides improved handling of entries like `Vitis_vinifera:plas.Chloroplast` and ensures appropriate taxonomic rank filtering.

## Key Improvements

### 1. Enhanced Organelle Detection and Cleaning

The parser now properly detects and handles organellar sequences:

**Supported Organelle Types:**
- **Chloroplast**: `.Chloroplast`, `:plas.Chloroplast`, `.plas.Chloroplast`
- **Mitochondria**: `.Mitochondria`, `:mito.Mitochondria`, `.mito.Mitochondria`
- **Plastid**: `.Plastid`, `:plas.Plastid`, `.plas.Plastid`
- **Apicoplast**: `.Apicoplast`, `:api.Apicoplast`, `.api.Apicoplast`

**Examples:**
```
Vitis_vinifera:plas.Chloroplast → Vitis vinifera
uncultured_bacterium.Mitochondria → uncultured bacterium
Arabidopsis_thaliana.Plastid → Arabidopsis thaliana
Plasmodium_falciparum.Apicoplast → Plasmodium falciparum
```

### 2. Taxonomic Rank Filtering

The parser now filters entries based on their appropriateness for the target taxonomic rank:

**Genus-level parsing:**
- ✅ Species entries → Extract genus (e.g., `Vitis vinifera` → `Vitis`)
- ✅ Genus entries → Keep as-is
- ❌ Family/higher entries → Filter out

**Family-level parsing:**
- ❌ Species/genus entries → Filter out (too specific)
- ✅ Family entries → Keep as-is
- ✅ Higher-level entries → Keep as-is

**Phylum-level parsing:**
- ❌ Species/genus/family entries → Filter out (too specific)
- ✅ Phylum entries → Keep as-is
- ✅ Higher-level entries → Keep as-is

### 3. Improved Name Cleaning

**Features:**
- Removes trailing numbers (e.g., `Theileria1` → `Theileria`)
- Handles Candidatus prefixes (e.g., `Candidatus Pelagibacter` → `Pelagibacter`)
- Converts underscores to spaces for species names
- Filters out unidentified/unknown entries with `.U.` patterns

## Files

### Core Scripts

1. **`enhanced_eukcensus_parser.py`** - Main enhanced parser
2. **`test_enhanced_parser.py`** - Test suite demonstrating functionality
3. **`optimized_eukcensus_16s_parser.py`** - Original parser (for comparison)

### Output Files

The enhanced parser generates:
- `eukcensus_enhanced_by_phylum.csv`
- `eukcensus_enhanced_by_family.csv`
- `eukcensus_enhanced_by_genus.csv`
- `eukcensus_enhanced_processing.log`

### Enhanced Output Format

The enhanced parser includes an additional column `original_names` that tracks all original taxon names that were consolidated:

```csv
Name_to_use,taxid,size_count,count,lineage,lineage_ranks,lineage_taxids,original_names
Vitis,12345,1000,50,cellular organisms;Eukaryota;...,cellular root;domain;...,131567;2759;...,Vitis_vinifera:plas.Chloroplast;Vitis_riparia.Mitochondria
```

## Usage

### Basic Usage

```bash
# Run with default input file (eukcensus_16S.clusters.97.tsv)
python enhanced_eukcensus_parser.py

# Run with custom input file
python enhanced_eukcensus_parser.py custom_input.tsv

# Run with custom input and output prefix
python enhanced_eukcensus_parser.py custom_input.tsv custom_output
```

### Testing

```bash
# Run the test suite to see functionality demonstrations
python test_enhanced_parser.py
```

## Key Functions

### `detect_organelle_type(taxon_name)`
Detects organelle type and extracts base organism name.

### `clean_organelle_taxon_name(taxon_name)`
Comprehensive cleaning including organelle removal, underscore conversion, and number stripping.

### `extract_appropriate_rank_name(taxon_name, target_rank)`
Extracts the appropriate taxonomic rank name and filters inappropriate entries.

### `validate_rank_appropriateness(taxon_name, target_rank, lineage_info)`
Uses lineage information to validate rank appropriateness.

## Example Workflow

For an entry like `Vitis_vinifera:plas.Chloroplast`:

1. **Organelle Detection**: Detects chloroplast organelle
2. **Name Cleaning**: `Vitis_vinifera:plas.Chloroplast` → `Vitis vinifera`
3. **Rank Extraction**:
   - Genus level: `Vitis` ✅
   - Family level: Filtered out ❌
   - Phylum level: Filtered out ❌
4. **Taxonkit Lookup**: Uses cleaned name `Vitis vinifera` for taxid lookup
5. **Lineage Retrieval**: Gets full taxonomic lineage
6. **Output**: Records `Vitis` with original name tracking

## Benefits

1. **Accurate Taxonomic Assignment**: Organellar sequences are properly attributed to their host organisms
2. **Rank-Appropriate Filtering**: Prevents species-level entries from contaminating higher-rank analyses
3. **Improved Data Quality**: Better cleaning and validation of taxonomic names
4. **Traceability**: Original names are preserved for verification
5. **Enhanced Logging**: Detailed processing logs for troubleshooting

## Comparison with Original Parser

| Feature | Original Parser | Enhanced Parser |
|---------|----------------|-----------------|
| Organelle Detection | Basic (some patterns) | Comprehensive (4 types, multiple patterns) |
| Rank Filtering | None | Intelligent filtering by target rank |
| Species→Genus | Limited | Robust extraction |
| Name Cleaning | Basic | Enhanced (numbers, Candidatus, etc.) |
| Traceability | None | Original names tracked |
| Logging | Basic | Detailed with organelle detection |

## Requirements

- Python 3.6+
- pandas
- taxonkit (in `../../taxonkit` relative path)
- NCBI taxonomy database (for taxonkit)

## Notes

- The parser expects the same input format as the original parser
- Taxonkit must be properly configured with NCBI taxonomy database
- The enhanced parser is backward compatible with existing workflows
- Processing time may be slightly longer due to enhanced validation