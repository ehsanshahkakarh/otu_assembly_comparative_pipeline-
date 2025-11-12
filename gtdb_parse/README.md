# GTDB Taxonomy Parser

This directory contains Python scripts for parsing GTDB (Genome Taxonomy Database) taxonomy files and extracting taxonomic information at different hierarchical levels (phylum, family, genus).

## Directory Structure

```
gtdb_parse/
├── py_gtdb/                    # Python parsing scripts (in metadata_proj)
│   ├── phylum_gtdb_parse.py   # Phylum-level parser
│   ├── family_gtdb_parse.py   # Family-level parser
│   ├── genus_gtdb_parse.py    # Genus-level parser
│   ├── 00bac120_taxonomy.tsv  # GTDB bacterial taxonomy input
│   ├── 00ar53_taxonomy.tsv    # GTDB archaeal taxonomy input
│   └── taxid.map              # Accession to taxid mapping
├── csv_gtdb/                  # Output CSV files (in metadata_proj)
│   ├── gtdb_phylum_counts.csv
│   ├── gtdb_phylum_with_accessions.csv
│   ├── gtdb_family_counts.csv
│   ├── gtdb_family_with_accessions.csv
│   ├── gtdb_genus_counts.csv
│   ├── gtdb_genus_with_accessions.csv
│   └── logs/                  # Execution logs
└── taxdump_gtdp/              # GTDB taxdump database
    └── gtdb-taxdump-R226/     # GTDB taxonomy database files
        ├── names.dmp
        ├── nodes.dmp
        ├── merged.dmp
        └── delnodes.dmp
```

**Note**: The actual parsing scripts are located in `metadata_proj/parse_repaa_table/gtdb_parse/py_gtdb/` while this directory contains the GTDB taxdump database.

## Input Files

### 1. GTDB Taxonomy Files
- **`00bac120_taxonomy.tsv`**: Bacterial taxonomy data from GTDB
- **`00ar53_taxonomy.tsv`**: Archaeal taxonomy data from GTDB
- **Format**: `accession\ttaxonomy_string`
- **Example**: `GB_GCA_000005825.2\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...`

### 2. Taxid Mapping File
- **`taxid.map`**: Maps GTDB accessions to taxids
- **Format**: `accession\ttaxid`
- **Example**: `GB_GCA_000005825.2\t511145`

### 3. GTDB Taxdump Database
- **Location**: `taxdump_gtdp/gtdb-taxdump-R226/`
- **Files**: Standard taxdump format (names.dmp, nodes.dmp, etc.)
- **Purpose**: Provides taxid-to-lineage mapping for GTDB taxonomy

## Output Files

### Summary Files (Counts)
Contains aggregated species counts by taxonomic group with lineage information.

**Files**: `gtdb_{level}_counts.csv`

**Columns**:
- `{level}`: Taxonomic name (e.g., phylum, family, genus)
- `domain`: Bacteria or Archaea
- `{level}_species_count`: Number of unique species in this group
- `taxid`: GTDB or NCBI taxid
- `lineage`: Full taxonomic lineage (semicolon-separated)
- `lineage_ranks`: Taxonomic ranks (semicolon-separated)
- `lineage_taxids`: Lineage taxids (semicolon-separated)

**Example** (`gtdb_family_counts.csv`):
```csv
family,domain,family_species_count,taxid,lineage,lineage_ranks,lineage_taxids
Enterobacteriaceae,Bacteria,150,543,Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae,superkingdom;phylum;class;order;family,2;1224;1236;91347;543
SM23-39,Bacteria,30,589911515,Bacteria;Chlamydiota;Chlamydiia;Chlamydiales;SM23-39,superkingdom;phylum;class;order;family,81602897;1417588462;1218949594;1568705753;164208376
```

### Detailed Files (Accessions)
Contains individual accession mappings to taxonomic classifications.

**Files**: `gtdb_{level}_with_accessions.csv`

**Columns**:
- `accession_clean`: Cleaned accession ID (prefixes removed)
- `{level}`: Taxonomic name
- `domain`: Bacteria or Archaea
- `taxid`: GTDB or NCBI taxid (if available)

**Example** (`gtdb_family_with_accessions.csv`):
```csv
accession_clean,family,domain,taxid
GCA_000005825.2,Enterobacteriaceae,Bacteria,543
GCA_000006745.1,Enterobacteriaceae,Bacteria,543
```

## Comprehensive Fallback System

The parser implements a sophisticated dual-database approach for maximum taxonomic coverage:

### 1. Primary GTDB Lookup
```
Taxonomic Name → GTDB Taxdump → GTDB Taxid → GTDB Lineage
```
- **Success Rate**: 60-80%
- **Database**: `taxdump_gtdp/gtdb-taxdump-R226/`
- **Advantages**: GTDB-specific taxa, updated taxonomy

### 2. NCBI Fallback
```
Failed GTDB Names → NCBI Taxdump → NCBI Taxid → NCBI Lineage
```
- **Additional Success**: +10-20%
- **Database**: Standard NCBI taxdump
- **Advantages**: Broader coverage, established taxa

### 3. Source-Aware Lineage Processing
```
GTDB Taxids → GTDB Lineage Lookup
NCBI Taxids → NCBI Lineage Lookup
Combined Results → Unified Output
```

## Parser Scripts

### Core Features
All parser scripts share the same architecture and features:

1. **Dual Database Support**: GTDB taxdump + NCBI fallback
2. **Species Deduplication**: One genome per species (keep='first')
3. **Name Standardization**: Removes suffixes like `_A`, `_B`
4. **Comprehensive Logging**: Success rates, failed matches, performance metrics
5. **Lineage Cleanup**: Truncates lineages to target taxonomic level

### Script-Specific Details

#### `phylum_gtdb_parse.py`
- **Target Level**: Phylum
- **Output**: `gtdb_phylum_counts.csv`, `gtdb_phylum_with_accessions.csv`
- **Lineage Cleanup**: Keeps domain + phylum levels only

#### `family_gtdb_parse.py`
- **Target Level**: Family
- **Output**: `gtdb_family_counts.csv`, `gtdb_family_with_accessions.csv`
- **Lineage Cleanup**: Keeps domain through family levels

#### `genus_gtdb_parse.py`
- **Target Level**: Genus
- **Output**: `gtdb_genus_counts.csv`, `gtdb_genus_with_accessions.csv`
- **Lineage Cleanup**: Keeps domain through genus levels

## Usage

### Basic Usage
```bash
cd metadata_proj/parse_repaa_table/gtdb_parse/py_gtdb/
python3 phylum_gtdb_parse.py
python3 family_gtdb_parse.py
python3 genus_gtdb_parse.py
```

### Advanced Options
```bash
python3 phylum_gtdb_parse.py \
    --input-dir /path/to/input \
    --output-dir /path/to/output \
    --log-level DEBUG
```

### Command Line Arguments
- `--input-dir`: Directory containing input files (default: script directory)
- `--output-dir`: Directory for output files (default: ../csv_gtdb)
- `--log-level`: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- `--verbose`: Enable verbose output

## Workflow Overview

```
Input Data → Load & Combine → Extract Ranks → Clean & Dedupe → Count by Level
     ↓
Taxid Resolution (GTDB → NCBI Fallback) → Lineage Resolution (Source-Aware)
     ↓
Lineage Cleanup → Generate Outputs (Summary + Detailed)
```

### Key Processing Steps:

1. **Data Loading**: Combine bacterial and archaeal taxonomy files
2. **Rank Extraction**: Parse taxonomy strings into structured ranks
3. **Data Cleaning**: Standardize names, remove duplicates
4. **Taxonomic Counting**: Group and count by target level
5. **Taxid Resolution**: GTDB primary + NCBI fallback
6. **Lineage Resolution**: Source-aware lineage lookup
7. **Output Generation**: Summary counts + detailed accessions

## Dependencies

### Required Python Packages
```bash
pip install pandas tqdm pathlib
```

### External Tools
- **taxonkit**: Required for taxid and lineage lookups
- **Installation**: `conda install -c bioconda taxonkit`

## Performance

### Typical Processing Times
- **Phylum level**: ~2-5 minutes
- **Family level**: ~5-10 minutes  
- **Genus level**: ~10-15 minutes

### Memory Usage
- **Peak RAM**: ~2-4 GB (depends on input size)
- **Chunked processing**: Handles large files efficiently

## Data Quality

### Success Rates (Typical)
- **GTDB taxid lookup**: 60-80%
- **NCBI fallback**: +10-20%
- **Combined success**: 70-90%
- **Lineage retrieval**: 95-99% of found taxids

### Validation
- Taxonomic hierarchy consistency checks
- Missing value detection
- Accession format validation
- Species count verification

## Troubleshooting

### Common Issues

1. **"taxonomy data not found"**
   - Check GTDB taxdump path in script
   - Ensure taxdump files exist in `taxdump_gtdp/gtdb-taxdump-R226/`

2. **Low taxid success rates**
   - Normal for GTDB-specific taxa
   - NCBI fallback provides additional coverage

3. **Memory errors**
   - Reduce `CHUNK_SIZE` in script
   - Ensure sufficient RAM available

### Log Files
- **Location**: `csv_gtdb/logs/`
- **Format**: `gtdb_{level}_parse_YYYYMMDD_HHMMSS.log`
- **Content**: Detailed execution logs, error messages, performance metrics

## Important Notes

### Taxid Interpretation
- **GTDB taxids**: May not exist in NCBI (e.g., 589911515 for SM23-39)
- **NCBI taxids**: Standard NCBI taxonomy identifiers
- **Lineage accuracy**: Manually verified to be correct regardless of source

### File Locations
- **Scripts**: `metadata_proj/parse_repaa_table/gtdb_parse/py_gtdb/`
- **Outputs**: `metadata_proj/parse_repaa_table/gtdb_parse/csv_gtdb/`
- **GTDB Database**: `gtdb/parse_repaa_table/gtdb_parse/taxdump_gtdp/`

For detailed workflow information, see `WORKFLOW.md` in the metadata_proj directory.