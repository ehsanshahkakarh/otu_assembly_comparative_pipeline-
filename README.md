# Environmental-Genomic Comparative Pipeline

## 🌍 **Project Overview**

This comprehensive pipeline quantifies **how well genomic databases represent environmental microbial diversity** by comparing environmental census data (16S/18S rRNA surveys) with major genomic databases (NCBI, GTDB, EukProt). The pipeline reveals cultivation gaps, database biases, and novel environmental diversity awaiting discovery.

### **Core Research Question**
> **"What percentage of environmental microbial diversity is represented in our genomic databases, and where are the biggest gaps?"**

## 🎯 **Pipeline Architecture**

```
                    ENVIRONMENTAL BASELINES
                           ↓
        ┌─────────────────────────────────────────────────────┐
        │  16S_censusparse     │     18S_censusparse          │
        │  (Prokaryotic)       │     (Eukaryotic)             │
        │  287K OTUs           │     70K OTUs                 │
        │  4,578 genera        │     847 genera               │
        └─────────────────────────────────────────────────────┘
                           ↓
                    GENOMIC DATABASES
                           ↓
        ┌─────────────────────────────────────────────────────┐
        │  ncbi_parse    │  gtdb_parse    │  eukprot_parse    │
        │  (Universal)   │  (Prokaryotic) │  (Eukaryotic)     │
        │  2.9M genomes  │  400K genomes  │  Protein DB       │
        │  142K species  │  196 phyla     │  Modern taxonomy  │
        └─────────────────────────────────────────────────────┘
                           ↓
                    COMPARATIVE ANALYSIS
                           ↓
        ┌─────────────────────────────────────────────────────┐
        │                Eukcensus_merge                      │
        │                                                     │
        │  • 16S ↔ GTDB/NCBI (Prokaryotic comparison)        │
        │  • 18S ↔ EukProt/NCBI (Eukaryotic comparison)      │
        │  • Coverage analysis & novelty quantification      │
        │  • Cultivation gap identification                  │
        └─────────────────────────────────────────────────────┘
                           ↓
                    VISUALIZATION & INSIGHTS
                           ↓
        ┌─────────────────────────────────────────────────────┐
        │                    visuals/                         │
        │                                                     │
        │  • Scatter plots (Coverage vs Novelty)             │
        │  • Alluvial plots (Database flow analysis)         │
        │  • Domain-specific visualizations                  │
        │  • Publication-ready figures                       │
        └─────────────────────────────────────────────────────┘
```

## 📊 **Key Findings & Impact**

### **Prokaryotic Diversity Gaps (16S Analysis)**
- **90%+ environmental genera** lack genomic representation
- **Massive cultivation gaps** in environmental prokaryotes
- **Database bias** toward easily cultured organisms
- **Novel diversity** in environmental samples vastly exceeds genomic databases

### **Eukaryotic Diversity Gaps (18S Analysis)**
- **94.4% division-level coverage** but **36.4% genus-level coverage**
- **Extreme novelty factors** (>100x) for major eukaryotic groups
- **Protist cultivation crisis**: Vast environmental diversity uncultured
- **Fungal gaps**: Environmental fungi poorly represented in databases

### **Cross-Database Validation**
- **GTDB vs NCBI**: Modern phylogenetic taxonomy vs traditional systems
- **Taxonomic consistency**: Resolving naming conflicts and polyphyletic groups
- **Database completeness**: Quantifying representation across domains of life

## 📁 **Directory Structure & Pipeline Flow**

```
parse_repaa_table/
├── 16S_censusparse/                      # 🦠 Prokaryotic Environmental Baseline
│   ├── py_16S/                           # Enhanced 16S parser with organelle handling
│   ├── csv_16S/                          # 287K OTUs across 4,578 genera
│   └── metadata/                         # Raw EukCensus 16S cluster data
│
├── 18S_censusparse/                      # 🧬 Eukaryotic Environmental Baseline
│   ├── py_18S/                           # Advanced 18S parser with pattern recognition
│   ├── csv_outputs/                      # 70K OTUs across 847 genera
│   └── metadata/                         # Raw EukCensus 18S cluster data
│
├── ncbi_parse/                           # 🗄️ Universal Genomic Database Reference
│   ├── py_ncbi/                          # Unified parser for 2.9M genomes
│   ├── csv_ncbi/                         # 142K species across all domains
│   ├── taxonomic_mapping/                # Taxid-to-name conversion pipeline
│   └── metadata/                         # NCBI assembly summary data
│
├── gtdb_parse/                           # 🔬 Modern Prokaryotic Taxonomy Reference
│   ├── py_gtdb/                          # GTDB R226 parser (400K genomes)
│   ├── csv_gtdb/                         # Phylogenetically-consistent taxonomy
│   ├── taxdump_gtdb/                     # GTDB taxonomy database
│   └── metadata/                         # GTDB bacterial/archaeal taxonomy
│
├── eukprot_parse/                        # 🧪 Eukaryotic Protein Database Reference
│   ├── py_eukprot/                       # EukProt parser with lineage integration
│   ├── csv_output/                       # Structured eukaryotic protein data
│   └── metadata/                         # Raw EukProt database files
│
├── Eukcensus_merge/                      # ⚡ Comparative Analysis Engine
│   ├── py_mergers/                       # Advanced merger scripts with vectorized processing
│   ├── 16s_merged/                       # Prokaryotic comparison results
│   ├── 18s_merged/                       # Eukaryotic comparison results
│   └── analysis_summary/                 # Coverage metrics and novelty analysis
│
├── visuals/                              # 📈 Visualization & Publication Pipeline
│   ├── scatter_plots/                    # Coverage vs novelty analysis
│   ├── alluvial_plots/                   # Database flow visualizations
│   ├── new_visualizations/               # Domain-specific plots
│   └── final_visualizations/             # Publication-ready figures
│
└── README.md                             # This comprehensive overview
```

## 🔬 **Technical Innovation & Methodology**

### **Advanced Data Processing**
- **Vectorized merging algorithms**: 4x performance improvement over traditional approaches
- **Multi-stream matching**: Direct, taxid, lineage, and pattern-based taxonomic matching
- **Enhanced name resolution**: Handles complex naming conventions and organelle sequences
- **Quality control systems**: Comprehensive validation and error checking throughout

### **Taxonomic Integration**
- **NCBI Taxonomy Integration**: Full lineage resolution using taxonkit
- **GTDB Modern Taxonomy**: Phylogenetically-consistent prokaryotic classification
- **Cross-database validation**: Resolving naming conflicts between taxonomic systems
- **Hierarchical processing**: Species → Genus → Family → Phylum aggregation

### **Coverage & Novelty Metrics**
```
Coverage = (Genomic_Species_Count / Environmental_OTU_Count) × 100
Novelty_Factor = Environmental_OTU_Count / Genomic_Species_Count
Overrepresentation_Factor = Genomic_Species_Count / Environmental_OTU_Count
```

## 🚀 **Quick Start Guide**

### **1. Environmental Baseline Generation**
```bash
# Process 16S prokaryotic environmental data
cd 16S_censusparse/py_16S/
python 16S_eukcensus_parser.py

# Process 18S eukaryotic environmental data
cd ../../18S_censusparse/py_18S/
python 18S_eukcensus_parser.py
```

### **2. Genomic Database Processing**
```bash
# Process NCBI genomic database
cd ../../ncbi_parse/py_ncbi/
python ncbi_parser_clean.py --level phylum
python ncbi_parser_clean.py --level family
python ncbi_parser_clean.py --level genus

# Process GTDB modern taxonomy
cd ../../gtdb_parse/py_gtdb/
python phylum_gtdb_parse.py
python family_gtdb_parse.py
python genus_gtdb_parse.py

# Process EukProt protein database
cd ../../eukprot_parse/py_eukprot/
python eukprot_parser.py
```

### **3. Comparative Analysis**
```bash
# Run environmental-genomic comparisons
cd ../../Eukcensus_merge/py_mergers/

# Prokaryotic comparisons
python 16s_gtdb_merger.py    # 16S vs GTDB
python 16s_ncbi_merger.py    # 16S vs NCBI

# Eukaryotic comparisons
python 18s_eukprot_merger.py # 18S vs EukProt
python 18s_ncbi_merger.py    # 18S vs NCBI
```

### **4. Visualization Generation**
```bash
# Generate publication-ready visualizations
cd ../../visuals/scatter_plots/
Rscript mega_comprehensive_stacked_visual.R

cd ../alluvial_plots/
Rscript create_alluvial_plots.R
```

## 📈 **Research Applications & Impact**

### **Microbial Ecology**
- **Cultivation Gap Analysis**: Quantifies uncultured microbial diversity
- **Environmental Importance**: Links environmental abundance to genomic representation
- **Taxonomic Priorities**: Identifies lineages requiring cultivation efforts
- **Biodiversity Assessment**: Measures true vs represented microbial diversity

### **Database Development**
- **Coverage Assessment**: Evaluates genomic database completeness
- **Bias Identification**: Reveals systematic underrepresentation of taxa
- **Quality Metrics**: Provides quantitative database evaluation frameworks
- **Improvement Guidance**: Directs future sequencing and cultivation efforts

### **Comparative Genomics**
- **Cross-database validation**: GTDB vs NCBI taxonomic consistency
- **Modern taxonomy adoption**: Phylogenetic vs traditional classification systems
- **Standardization efforts**: Unified taxonomic nomenclature across databases
- **Quality control**: Identifies taxonomic inconsistencies and errors

## 🎯 **Key Performance Metrics**

### **Processing Scale**
- **Environmental Data**: 357K total OTUs (287K prokaryotic + 70K eukaryotic)
- **Genomic Data**: 3.3M+ genome assemblies across all databases
- **Taxonomic Resolution**: Species through phylum levels with full lineages
- **Cross-database Integration**: 6 major database comparisons

### **Analysis Efficiency**
- **Vectorized Processing**: 4x speed improvement over legacy methods
- **Memory Optimization**: Handles large datasets efficiently
- **Parallel Processing**: Concurrent operations where applicable
- **Quality Control**: Comprehensive validation throughout pipeline

### **Research Output**
- **Coverage Quantification**: Precise measurement of database representation
- **Novelty Discovery**: Identification of uncultured environmental diversity
- **Publication-ready Figures**: Professional visualizations for research dissemination
- **Reproducible Workflows**: Documented, standardized processing pipelines

## 🔧 **Dependencies & Requirements**

### **Core Dependencies**
- **Python 3.7+**: pandas, numpy, tqdm, pathlib
- **R 4.0+**: ggplot2, dplyr, ggrepel, viridis, alluvial
- **taxonkit**: NCBI taxonomy operations and lineage resolution
- **System**: Linux/Unix environment with sufficient RAM (16GB+ recommended)

### **Database Requirements**
- **NCBI Taxonomy**: taxdump files for taxonomic resolution
- **GTDB Database**: R226 or later for modern prokaryotic taxonomy
- **EukCensus Data**: 16S/18S environmental cluster data
- **Storage**: ~50GB for complete pipeline data and outputs

---

**This pipeline represents a comprehensive framework for quantifying the relationship between environmental microbial diversity and genomic database representation, revealing the vast "microbial dark matter" that remains uncultured and providing quantitative guidance for future cultivation and sequencing efforts.**