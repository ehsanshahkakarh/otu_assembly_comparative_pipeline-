# Modular NCBI Merger Architecture
## New Merger - Modular Source Code Structure

This document describes the modular architecture of the NCBI Census Merger.

---

## 📁 Directory Structure

```
new_merger/
├── config.py                      # Centralized configuration
├── run_18s_ncbi_merger.py        # 18S entry point (Eukaryotes)
├── run_16s_ncbi_merger.py        # 16S entry point (Prokaryotes)
├── outputs/                       # Merged CSV outputs
├── logs/                          # Processing logs
├── src/                           # SHARED source modules
│   ├── __init__.py
│   ├── data_loader.py             # Load census + NCBI data
│   ├── domain_filter.py           # Domain-specific filtering
│   ├── lineage_matcher.py         # Vectorized lineage matching
│   ├── isolate_tracker.py         # Isolate genome analysis
│   ├── metrics_calculator.py      # Novelty, overrepresentation, etc.
│   └── output_writer.py           # Save merged CSVs
└── py_merger/                     # OLD monolithic scripts (deprecated)
    ├── 18s_ncbi_merger.py
    ├── 16s_ncbi_merger.py
    ├── merge_division.py
    ├── merge_family.py
    └── merge_genus.py
```

---

## 🏗️ Architecture Overview

### **Design Pattern: Modular Pipeline**

The merger follows a **modular pipeline architecture** similar to:
- `ncbi_parse/metadata/nev_parse_meth/` (species parser)
- `18S_censusparse/py_18S/` (census parser)
- `database_merger/py_dtb/` (triple-anchor merger)

### **Key Principles:**

1. **Single Responsibility**: Each module has one clear purpose
2. **Shared Code**: One `src/` directory serves both 18S and 16S
3. **Configuration-Driven**: Domain filtering and paths via `config.py`
4. **Logging**: Comprehensive logging at each step
5. **Testability**: Each module can be tested independently

---

## 📦 Module Descriptions

### **`config.py`** - Configuration Module

**Purpose**: Centralized configuration for paths, domains, and taxonomic levels

**Key Classes**:
- `PathConfig`: File paths for inputs/outputs
- `DomainConfig`: Domain filtering (Eukaryota vs Bacteria+Archaea)
- `TaxonomicConfig`: Taxonomic level mappings
- `Config`: Main configuration container

**Functions**:
- `setup_logging()`: Configure logging with file and console handlers

**Usage**:
```python
from config import Config
config = Config()
census_file = config.paths.get_census_file('18S', 'phylum')
domains = config.domains.get_domains('18S')  # ['Eukaryota']
```

---

### **`src/data_loader.py`** - Data Loading Module

**Purpose**: Load census and NCBI data files with validation

**Functions**:
- `load_census_data(census_file)`: Load 18S/16S census CSV
- `load_ncbi_data(ncbi_file, level)`: Load NCBI counts CSV
- `load_isolate_data(isolate_file, level)`: Load isolate genome data (optional)

**Returns**: Validated pandas DataFrames

**Error Handling**: Raises `FileNotFoundError` or `ValueError` for missing files/columns

---

### **`src/domain_filter.py`** - Domain Filtering Module

**Purpose**: Filter NCBI data to specific domains

**Functions**:
- `filter_by_domain(ncbi_df, domains, isolate_df)`: Filter to Eukaryota or Bacteria+Archaea

**Logic**:
- 18S: Filters to `['Eukaryota']`
- 16S: Filters to `['Bacteria', 'Archaea']`

**Returns**: Tuple of (filtered_ncbi_df, filtered_isolate_df)

---

### **`src/lineage_matcher.py`** - Lineage Matching Module

**Purpose**: Vectorized lineage-based matching between census and NCBI

**Functions**:
- `match_taxa_by_lineage(census_df, ncbi_df, level, census_level)`: Match taxa using regex

**Matching Strategy**:
```python
# For each census taxon, search NCBI lineages using regex:
pattern = f';{taxon};|^{taxon};|;{taxon}$|^{taxon}$'
matches = ncbi_df['lineage'].str.contains(pattern, regex=True)
```

**Aggregation**: Sums genome/species counts from all matching NCBI entries

**Returns**: DataFrame with matched results including:
- Census data (OTU counts, percentages)
- NCBI data (genome/species counts, percentages)
- Match statistics (direct, lineage, total matches)

---

### **`src/isolate_tracker.py`** - Isolate Tracking Module

**Purpose**: Add isolate genome information to matched results

**Functions**:
- `add_isolate_information(matched_df, isolate_df, level, census_level)`: Add isolate counts

**Logic**:
- Looks up isolate counts for each matched taxon
- Calculates isolate percentage: `isolate_count / total_genomes * 100`

**Returns**: DataFrame with added `isolate_count` and `isolate_percentage` columns

---

### **`src/metrics_calculator.py`** - Metrics Calculation Module

**Purpose**: Calculate novelty, overrepresentation, and coverage metrics

**Functions**:
- `calculate_metrics(matched_df)`: Calculate all metrics

**Metrics**:
1. **Novelty Factor**: `census_otu_count / ncbi_species_count`
   - Higher = more environmental diversity than genomic representation
   - >10 = massive cultivation gaps
   - <2 = well-represented taxa

2. **Overrepresentation Factor**: `ncbi_species_count / census_otu_count`
   - Higher = database bias toward cultured taxa
   - >2 = overrepresented in databases
   - <0.5 = environmental diversity dominates

3. **Coverage Percentage**: Placeholder for future implementation

**Returns**: DataFrame with added metric columns

---

### **`src/output_writer.py`** - Output Writing Module

**Purpose**: Save merged results to CSV files with proper sorting

**Functions**:
- `save_merged_output(merged_df, output_file, census_level)`: Save sorted CSV

**Sorting Logic**:
1. Matched taxa first (sorted by genome count descending)
2. Census-only taxa second (sorted by OTU count descending)

**Returns**: Path to saved output file

---

## 🚀 Entry Point Scripts

### **`run_18s_ncbi_merger.py`** - 18S Eukaryotic Merger

**Purpose**: Merge 18S eukcensus with NCBI eukaryotic genomes

**Pipeline Steps**:
1. Load data (census, NCBI, isolate)
2. Filter by domain (Eukaryota only)
3. Match taxa by lineage
4. Add isolate information
5. Calculate metrics
6. Save output

**Processes**: Division, Family, Genus levels

**Usage**:
```bash
python run_18s_ncbi_merger.py
```

---

### **`run_16s_ncbi_merger.py`** - 16S Prokaryotic Merger

**Purpose**: Merge 16S eukcensus with NCBI prokaryotic genomes

**Pipeline Steps**: Same as 18S merger

**Domain Filter**: Bacteria + Archaea only

**Processes**: Division, Family, Genus levels

**Usage**:
```bash
python run_16s_ncbi_merger.py
```

---

## 🔄 Data Flow

```
[Census CSV] ──┐
               ├──> [Load Data] ──> [Filter Domain] ──> [Match Lineage] ──┐
[NCBI CSV]  ───┤                                                           │
               │                                                           ├──> [Add Isolates] ──> [Calculate Metrics] ──> [Save Output]
[Isolate CSV]──┘                                                           │
                                                                           │
                                                                    [Matched DataFrame]
```

---

## 📊 Output Schema

Each output CSV contains:

| Column | Description |
|--------|-------------|
| `division/family/genus` | Taxon name |
| `census_otu_count` | Environmental OTU count |
| `census_size_count` | Sequence abundance |
| `otu_percentage` | % of total census OTUs |
| `size_percentage` | % of total census size |
| `ncbi_genome_count` | Total genomes from NCBI |
| `ncbi_species_count` | Total species from NCBI |
| `isolate_count` | Number of isolate genomes |
| `genome_pct_db` | % of NCBI database |
| `species_pct` | % of NCBI species |
| `isolate_percentage` | % of genomes that are isolates |
| `novelty_factor` | census_otu / ncbi_species |
| `overrepresentation_factor` | ncbi_species / census_otu |
| `direct_matches` | Exact name matches |
| `lineage_matches` | Hierarchical matches |
| `total_matches` | Total NCBI entries matched |
| `match_status` | 'matched' or 'census_only' |
| `domain` | Most common domain |


---

## 🆚 Comparison: Modular vs Monolithic

| Aspect | Modular (NEW) | Monolithic (OLD) |
|--------|---------------|------------------|
| **Files** | 9 modules + 2 entry points | 5 standalone scripts |
| **Code Reuse** | 100% shared between 18S/16S | ~95% duplicated |
| **Maintainability** | High (single source of truth) | Low (must update 2 files) |
| **Testability** | High (test each module) | Low (test entire script) |
| **Readability** | High (clear separation) | Medium (all in one file) |
| **Logging** | Comprehensive (per module) | Basic (console only) |
| **Configuration** | Centralized (config.py) | Hardcoded paths |

---

## 🧪 Testing

Each module can be tested independently:

```python
# Test data loader
from src.data_loader import load_census_data
df = load_census_data(Path('test_census.csv'))

# Test domain filter
from src.domain_filter import filter_by_domain
filtered_df, _ = filter_by_domain(ncbi_df, ['Eukaryota'])

# Test lineage matcher
from src.lineage_matcher import match_taxa_by_lineage
matched_df = match_taxa_by_lineage(census_df, ncbi_df, 'phylum', 'division')
```

---

**Created**: 2026-02-19  
**Author**: Ehsan Shah Kakar (with Augment AI)  
**Architecture Pattern**: Modular Pipeline (following nev_parse_meth pattern)

