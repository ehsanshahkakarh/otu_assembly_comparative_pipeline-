# New Merger: Eukcensus Census + NCBI Data

## Overview

This directory contains scripts and outputs for merging eukcensus census data (both 16S and 18S) with NCBI genome/species data. The mergers use the **eukcensus census files as the backbone** and find matching NCBI entries based on lineage-based matching with comprehensive metrics.

## 📚 Documentation

For detailed documentation, see the `docs/` folder:
- **[QUICK_START.md](docs/QUICK_START.md)** - Quick start guide for running mergers
- **[MERGER_COMPARISON.md](docs/MERGER_COMPARISON.md)** - Comparison of different merger approaches
- **[MODULAR_ARCHITECTURE.md](docs/MODULAR_ARCHITECTURE.md)** - Technical architecture details

## 📁 Directory Structure

```
new_merger/
├── README.md                      # This file
├── docs/                          # Documentation
│   ├── QUICK_START.md            # Quick start guide
│   ├── MERGER_COMPARISON.md      # Comparison of merger approaches
│   └── MODULAR_ARCHITECTURE.md   # Architecture documentation
├── config.py                      # Configuration
├── run_18s_ncbi_merger.py        # 18S eukaryotic merger (MAIN)
├── run_16s_ncbi_merger.py        # 16S prokaryotic merger (MAIN)
├── src/                           # Core source modules
│   ├── data_loader.py
│   ├── domain_filter.py
│   ├── lineage_matcher.py
│   ├── isolate_tracker.py
│   ├── metrics_calculator.py
│   └── output_writer.py
├── py_merger/                     # Legacy monolithic scripts
│   ├── 18s_ncbi_merger.py
│   ├── 16s_ncbi_merger.py
│   └── merge_*.py
├── outputs/                       # Merged CSV outputs (MAIN)
├── logs/                          # Processing logs
├── utils/                         # Utility scripts
│   ├── analyze_hierarchical_matching.py
│   ├── sanity_check_merger.py
│   └── sort_by_novelty.py
├── test_outputs/                  # Test/experimental outputs
├── test_src/                      # Test implementations
└── optimized_src/                 # Optimized implementations
```

## Methodology

### Merging Strategy

1. **Backbone**: Eukcensus 18S files serve as the reference/backbone
2. **Matching**: For each eukcensus entry, we find all NCBI entries that share any taxid in their lineage
3. **Aggregation**: When multiple NCBI entries match, we aggregate their data:
   - Sum genome counts
   - Sum species counts
   - Concatenate matched names
   - Track number of matches

### Merger Scripts

Five merger scripts are available:

1. **18s_ncbi_merger.py** - 18S eukcensus ↔ NCBI eukaryotic genomes (vectorized, with isolate tracking)
2. **16s_ncbi_merger.py** - 16S eukcensus ↔ NCBI prokaryotic genomes (vectorized, with isolate tracking)
3. **merge_division.py** - 18S division-level merger (original simple version)
4. **merge_family.py** - 18S family-level merger (original simple version)
5. **merge_genus.py** - 18S genus-level merger (original simple version)

### Taxonomic Levels

Three taxonomic levels are processed:

1. **Division/Phylum** (eukcensus) ↔ **Phylum** (NCBI)
2. **Family** (eukcensus) ↔ **Family** (NCBI)
3. **Genus** (eukcensus) ↔ **Genus** (NCBI)

### Input Files

**Eukcensus 18S (Eukaryotic):**
- `18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv`
- `18S_censusparse/csv_outputs/eukcensus_18S_by_family.csv`
- `18S_censusparse/csv_outputs/eukcensus_18S_by_genus.csv`

**Eukcensus 16S (Prokaryotic):**
- `16S_censusparse/csv_16S/eukcensus16S_by_division.csv`
- `16S_censusparse/csv_16S/eukcensus16S_by_family.csv`
- `16S_censusparse/csv_16S/eukcensus16S_by_genus.csv`

**NCBI:**
- `ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv`
- `ncbi_parse/csv_ncbi/ncbi_family_counts.csv`
- `ncbi_parse/csv_ncbi/ncbi_genus_counts.csv`
- `ncbi_parse/csv_ncbi/ncbi_*_with_accessions.csv` (for isolate analysis)

### Output Files

All outputs are saved to `new_merger/outputs/`:

**18S NCBI Merger:**
- `18s_ncbi_merged_division.csv`
- `18s_ncbi_merged_family.csv`
- `18s_ncbi_merged_genus.csv`

**16S NCBI Merger:**
- `16s_ncbi_merged_division.csv`
- `16s_ncbi_merged_family.csv`
- `16s_ncbi_merged_genus.csv`

**Original Simple Mergers:**
- `merged_division.csv`
- `merged_family.csv`
- `merged_genus.csv`

### Output Schema

#### **New NCBI Merger Scripts (18s_ncbi_merger.py, 16s_ncbi_merger.py)**

Each merged file contains comprehensive metrics:

**From Census:**
- `division/family/genus`: Taxon name
- `census_otu_count`: Environmental OTU count
- `census_size_count`: Sequence abundance
- `otu_percentage`: Percentage of total census OTUs
- `size_percentage`: Percentage of total census size

**From NCBI (aggregated from all matches):**
- `ncbi_genome_count`: Total genome count from matched entries
- `ncbi_species_count`: Total species count from matched entries
- `isolate_count`: Number of isolate genomes
- `genome_pct_db`: Percentage of NCBI database (domain-specific)
- `species_pct`: Percentage of NCBI species (domain-specific)
- `isolate_percentage`: Percentage of genomes that are isolates

**Comparative Metrics:**
- `novelty_factor`: census_otu_count / ncbi_species_count (higher = more novel)
- `overrepresentation_factor`: ncbi_species_count / census_otu_count (higher = database bias)

**Match Statistics:**
- `direct_matches`: Exact name matches
- `lineage_matches`: Hierarchical lineage matches
- `total_matches`: Total NCBI entries matched
- `match_status`: 'matched' or 'census_only'
- `domain`: Most common domain from matched entries

#### **Original Simple Merger Scripts (merge_*.py)**

Each merged file contains:

**From Eukcensus:**
- `Division/Family/Genus`: Name to use
- `Eukcensus_Taxid`: Taxonomic ID
- `Eukcensus_OTU_Count`: OTU count
- `Eukcensus_Size_Count`: Size count
- `Eukcensus_Lineage`: Full lineage string
- `Eukcensus_Lineage_Ranks`: Rank names
- `Eukcensus_Lineage_Taxids`: Taxid lineage

**From NCBI (aggregated):**
- `NCBI_Domain`: Domain(s) from matched entries
- `NCBI_Genome_Count`: Total genome count
- `NCBI_Species_Count`: Total species count
- `NCBI_Match_Count`: Number of NCBI entries matched
- `NCBI_Matched_Names`: Names of matched NCBI entries

**Genus level only:**
- `NCBI_Isolate_Genome_Count`: Total isolate genome count

## Usage

### **Recommended: New NCBI Merger Scripts (Vectorized, Production-Ready)**

Run the comprehensive NCBI mergers with isolate tracking and advanced metrics:

```bash
# 18S Eukaryotic merger (all three levels at once)
python new_merger/py_merger/18s_ncbi_merger.py

# 16S Prokaryotic merger (all three levels at once)
python new_merger/py_merger/16s_ncbi_merger.py
```

**Features:**
- ✅ Vectorized pandas operations (4x faster)
- ✅ Domain-specific filtering (Eukaryota or Bacteria+Archaea)
- ✅ Isolate genome tracking
- ✅ Comprehensive metrics (novelty, overrepresentation, coverage)
- ✅ Lineage-based matching with aggregation
- ✅ Progress indicators and summary statistics

### **Alternative: Original Simple Merger Scripts**

Run individual merger scripts for 18S data only:

```bash
# Division level
python new_merger/py_merger/merge_division.py

# Family level
python new_merger/py_merger/merge_family.py

# Genus level
python new_merger/py_merger/merge_genus.py
```

Or run all at once:

```bash
cd new_merger/py_merger
for script in merge_*.py; do
    echo "Running $script..."
    python "$script"
done
```

## Key Features & Differences

### **New NCBI Merger Scripts vs Original Scripts**

| Feature | New NCBI Mergers | Original Mergers |
|---------|------------------|------------------|
| **Performance** | Vectorized (4x faster) | Iterative loops |
| **Matching Method** | Lineage-based regex | Simple name matching |
| **Isolate Tracking** | ✅ Full isolate analysis | ❌ Not available |
| **Domain Filtering** | ✅ Eukaryota or Bacteria+Archaea | ❌ No filtering |
| **Metrics** | Comprehensive (novelty, overrepresentation, coverage) | Basic counts only |
| **Aggregation** | Sums from multiple matches | Simple 1:1 mapping |
| **Output Format** | Standardized, sorted by match status | Basic CSV |
| **Progress Indicators** | ✅ Detailed progress output | ❌ Minimal output |
| **Data Sources** | Both 16S and 18S | 18S only |

### **New NCBI Mergers vs Eukcensus_merge Scripts**

The new NCBI merger scripts are **simplified, streamlined versions** of the production Eukcensus_merge scripts:

**Similarities:**
- ✅ Vectorized pandas operations
- ✅ Lineage-based matching with regex
- ✅ Isolate tracking and analysis
- ✅ Comprehensive metrics calculation
- ✅ Domain-specific filtering

**Differences:**
- ❌ No taxid-based matching (lineage-only)
- ❌ No detailed match logging to separate files
- ❌ No class-based architecture
- ❌ Simpler output structure (no NCBI-only taxa)
- ✅ Single-file, easier to understand
- ✅ All levels processed in one run

## Notes

- **New NCBI mergers**: Use lineage-based matching for maximum compatibility
- **Original mergers**: Simple name matching, may miss hierarchical relationships
- NCBI entries can match multiple census entries through lineage matching
- The `total_matches` field indicates how many NCBI entries contributed to aggregated data
- Isolate percentages show the proportion of cultured vs MAG/SAG genomes

