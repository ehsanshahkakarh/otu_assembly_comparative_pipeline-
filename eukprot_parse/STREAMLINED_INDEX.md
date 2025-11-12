# 🧬 EukProt Parse - Streamlined Directory Index

## 📋 Overview
Streamlined workflow for processing eukaryotic taxonomic data from EukProt database with NCBI taxonomy integration.

## 🗂️ Reorganized Directory Structure

```
eukprot_parse/
├── 📁 py_eukprot/                    # Core EukProt lineage generation
├── 📁 merge/                         # Cross-database integration & analysis
├── 📁 csv_output/                    # Final processed outputs
├── 📁 metadata/                      # Input datasets
└── 📁 visuals/                       # Analysis visualizations
```

## 🔄 Streamlined Workflow

### 1. **Core Processing** (`py_eukprot/`)
**Primary Script**: `improv_eukprot_lineage.py` *(Recently optimized - reduced from 2691 to 2302 lines)*
- **Input**: `metadata/Eukprot_included_datasets.txt`
- **Output**: `csv_output/eukprot_new_lineages.csv`
- **Features**:
  - Parallel processing with vectorized operations
  - Multi-stage taxonomic name mapping
  - Complete NCBI lineage generation
  - Optimized path handling for new directory structure

### 2. **Analysis & Integration** (`merge/`)
**Analysis Scripts**:
- `lineage_merger_div.py` - Division-level cross-database analysis
- `lineage_merger_family.py` - Family-level comparisons  
- `lineage_merger_genus.py` - Genus-level species matching
- `create_bias_visualizations.py` - Generate analysis charts
- `run_bias_analysis.py` - Comprehensive bias analysis

### 3. **Visualization** (`visuals/`)
- Taxonomic hierarchy visualizations
- Species overlap analysis
- Summary statistics charts

## 🎯 Recent Optimizations (v2.0)

### Code Reduction Achievements
- **Lines reduced**: 374 lines total (13.9% reduction)
  - High priority removals: 227 lines (deprecated functions, shell scripts)
  - Non-functional optimizations: 147 lines (documentation, logging)
- **Functions streamlined**: 8+ functions simplified
- **Path handling**: Updated for reorganized directory structure
- **Logging centralized**: All logs directed to `py_eukprot/log/` directory

### Performance Improvements
- ✅ Removed deprecated fallback functions
- ✅ Simplified shell script generation → direct subprocess calls
- ✅ Consolidated environment setup into single function
- ✅ Streamlined error handling and logging
- ✅ Optimized progress tracking
- ✅ Centralized log file management

### Maintainability Enhancements
- ✅ Single source of truth for environment setup
- ✅ Cleaner function signatures
- ✅ Reduced code duplication
- ✅ Simplified documentation
- ✅ Better path handling for reorganized structure
- ✅ Centralized log directory management with helper function

## 📊 Key Output Files

### Primary Outputs
1. **`csv_output/eukprot_new_lineages.csv`** - Complete EukProt lineages
2. **`merge/results/division_analysis_summary.csv`** - Cross-database analysis
3. **`visuals/*.png`** - Analysis visualizations

## 🚀 Quick Start

### Generate EukProt Lineages
```bash
cd py_eukprot/
python improv_eukprot_lineage.py
```

### Run Analysis
```bash
cd merge/
python lineage_merger_div.py
python create_bias_visualizations.py
```

## 🔧 Dependencies
- **Python 3.7+**: pandas, numpy, tqdm, multiprocessing
- **taxonkit**: NCBI taxonomy command-line tool
- **NCBI taxdump**: Taxonomy database files

## 📈 Directory Benefits

### Streamlined Structure
- **Centralized outputs**: All results in `csv_output/`
- **Clear separation**: Input (`metadata/`) vs Output (`csv_output/`)
- **Focused workflow**: Removed redundant 16S/18S processing directories
- **Better organization**: Analysis tools consolidated in `merge/`

### Improved Maintainability
- **Reduced complexity**: Fewer directories to manage
- **Clear data flow**: metadata → py_eukprot → csv_output → merge → visuals
- **Optimized scripts**: Cleaner, more efficient code
- **Better documentation**: Focused on essential information

This streamlined structure provides the same powerful taxonomic analysis capabilities with improved organization and significantly optimized code.
