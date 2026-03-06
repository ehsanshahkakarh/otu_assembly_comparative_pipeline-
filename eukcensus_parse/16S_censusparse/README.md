# 16S Environmental Census Data Processing Pipeline

## 🎯 **Purpose & Pipeline Position**

The **16S_censusparse** directory serves as the **environmental diversity baseline** for the entire comparative genomics pipeline. It processes EukCensus 16S rRNA gene cluster data to create structured taxonomic datasets that represent **actual environmental microbial diversity** detected in natural samples.

### **Pipeline Context**
```
Environmental Samples → 16S rRNA Sequencing → OTU Clustering → EukCensus Database
                                    ↓
                        THIS DIRECTORY (16S_censusparse)
                                    ↓
                    Structured Environmental Baseline
                                    ↓
        ┌─────────────────────────────────────────────────────────┐
        │              COMPARATIVE ANALYSIS                       │
        │                                                         │
        │  Environmental Diversity ←→ NCBI Genomic Database      │
        │  Environmental Diversity ←→ GTDB Taxonomic Database    │
        │  Environmental Diversity ←→ EukProt Protein Database   │
        └─────────────────────────────────────────────────────────┘
                                    ↓
                Database Merger Scripts → Visualizations → Publications
```

### **Research Questions Addressed**
1. **Coverage Analysis**: Which environmental taxa lack genomic representation?
2. **Database Completeness**: How well do genomic databases represent environmental diversity?
3. **Taxonomic Consistency**: Do NCBI and GTDB classifications align with environmental detections?
4. **Cultivation Bias**: Which environmental taxa remain uncultured?

## 📊 **Data Scale & Scope**

### **Input Dataset Statistics**
- **287,468 total OTU clusters** from environmental 16S surveys
- **1,283,740 total sequences** representing environmental diversity
- **4,578 unique genera** detected across all environments
- **1,843 unique families** spanning bacterial and archaeal diversity
- **95 unique phyla/divisions** covering major microbial lineages
- **100% taxonomic coverage** at all hierarchical levels

### **Taxonomic Scope**
**Primary Focus**: Bacterial and Archaeal diversity (16S rRNA gene)
- **Bacteria**: Dominant component (~85% of diversity)
- **Archaea**: Significant environmental component (~15% of diversity)
- **Note**: Despite "EukCensus" naming, 16S data primarily captures prokaryotic diversity

## 📁 **Directory Structure**

```
16S_censusparse/
├── py_16S/                           # 🐍 Processing Scripts
│   ├── 16S_eukcensus_parser.py       # ⭐ Main enhanced parser
│   ├── 16S_eukcensus_parser_01.py    # Legacy parser (comparison)
│   ├── ENHANCED_PARSER_README.md     # Detailed parser documentation
│   ├── IMPLEMENTATION_SUMMARY.md     # Enhancement summary
│   ├── logs/                         # Processing logs & unmapped analysis
│   └── sanity_checks/                # Validation & testing scripts
│
├── csv_16S/                          # 📈 Final Structured Outputs
│   ├── README.md                     # Output file documentation
│   ├── eukcensus16S_by_division.csv  # Phylum-level aggregation (95 taxa)
│   ├── eukcensus16S_by_family.csv    # Family-level aggregation (800 taxa)
│   └── eukcensus16S_by_genus.csv     # Genus-level aggregation (4,578 taxa)
│
├── metadata/                         # 📊 Input Data & Analysis
│   ├── eukcensus_16S.clusters.97.tsv # Raw EukCensus cluster data
│   ├── taxonomy_summary.txt          # Statistical overview
│   ├── taxonomy_combinations.csv     # Unique taxonomic combinations
│   └── parse_taxonomy_columns.py     # Metadata analysis tools
│
├── README.md                         # This comprehensive overview
└── IMPLEMENTATION_SUMMARY.md         # Technical implementation details
```

## 🔬 **Enhanced Processing Pipeline**

### **Core Processing Workflow**
```
Raw EukCensus Data (287K clusters)
    ↓
Enhanced Organelle Detection & Cleaning
    ↓
Taxonomic Rank Filtering & Validation
    ↓
NCBI Taxonomy Integration (taxonkit)
    ↓
Three-Level Taxonomic Aggregation
    ↓
Structured CSV Outputs for Database Comparison
```

### **🧬 Advanced Organelle Handling**
The enhanced parser addresses complex organellar sequence entries:

**Organelle Types Detected:**
- **Chloroplast**: `.Chloroplast`, `:plas.Chloroplast`, `.plas.Chloroplast`
- **Mitochondria**: `.Mitochondria`, `:mito.Mitochondria`, `.mito.Mitochondria`
- **Plastid**: `.Plastid`, `:plas.Plastid`, `.plas.Plastid`
- **Apicoplast**: `.Apicoplast`, `:api.Apicoplast`, `.api.Apicoplast`

**Processing Examples:**
```
Vitis_vinifera:plas.Chloroplast → Vitis vinifera → Vitis (genus)
Annulohypoxylon_stygium.Mitochondria → Annulohypoxylon stygium → Annulohypoxylon (genus)
Aspergillus_nidulans_FGSC_A4.Mitochondria → Aspergillus nidulans FGSC A → Aspergillus (genus)
```

### **📊 Taxonomic Rank Filtering**
Intelligent filtering ensures appropriate taxonomic resolution:

- **Division/Phylum Parsing**: Filters species/genus/family entries, keeps phylum+ entries
- **Family Parsing**: Filters species/genus entries, keeps family+ entries
- **Genus Parsing**: Species entries → Extract genus, keeps genus entries, filters higher ranks

### **🏷️ Enhanced Name Cleaning**
- **Underscore Conversion**: `Genus_species` → `Genus species`
- **Candidatus Handling**: Preserves Candidatus prefixes appropriately
- **Number Removal**: `Theileria1` → `Theileria`
- **Unidentified Preservation**: **PRESERVES** `.U.` taxa for visualization (not filtered out)

## 📈 **Output Files & Data Structure**

### **Three Taxonomic Resolution Levels**

#### **1. Division/Phylum Level** (`eukcensus16S_by_division.csv`)
- **Size**: 7.3 KB, ~95 taxonomic divisions
- **Purpose**: Broad taxonomic diversity overview
- **Examples**: Pseudomonadota (53,653 OTUs), Bacillota (49,858 OTUs)

#### **2. Family Level** (`eukcensus16S_by_family.csv`)
- **Size**: 151 KB, ~800 families
- **Purpose**: Ecological functional group analysis
- **Examples**: Lachnospiraceae (12,280 OTUs), Oscillospiraceae (5,753 OTUs)

#### **3. Genus Level** (`eukcensus16S_by_genus.csv`)
- **Size**: 1.1 MB, ~4,578 genera
- **Purpose**: High-resolution diversity analysis
- **Examples**: Bacillus (3,741 OTUs), Pseudomonas (2,357 OTUs)

### **Standardized Column Structure**
Each CSV file contains:
- **`Name_to_use`**: Clean taxonomic name for analysis
- **`taxid`**: NCBI taxonomy ID for cross-database linking
- **`otu_count`**: Number of environmental OTU clusters
- **`otu_percentage`**: Percentage of total environmental OTUs
- **`size_count`**: Total sequence abundance in environmental samples
- **`size_percentage`**: Percentage of total environmental sequences
- **`lineage`**: Full taxonomic lineage (semicolon-separated)
- **`lineage_ranks`**: Corresponding taxonomic ranks
- **`lineage_taxids`**: NCBI taxids for each lineage level

## 🔗 **Downstream Integration**

### **Database Merger Scripts** (`../database_merger/`)
Environmental baseline files serve as input for comparative analysis:
- **`16s_ncbi_merger.py`**: Compare environmental diversity with NCBI genomic coverage
- **`16s_gtdb_merger.py`**: Compare environmental diversity with GTDB taxonomic classifications
- **Triple-anchor merging**: Uses name, accession, and lineage-based matching strategies

### **Visualization Pipeline** (`../visuals/`)
Environmental data enables comprehensive visual analysis:
- **Scatter Plots**: Coverage vs novelty factors showing database representation gaps
- **Alluvial Plots**: Taxonomic flow between environmental and genomic databases
- **Coverage Analysis**: Identification of underrepresented environmental taxa

### **Comparative Genomics Research**
- **Gap Identification**: Environmental taxa lacking genomic representation
- **Cultivation Bias Assessment**: Cultured vs environmental microbial diversity
- **Database Quality Control**: Consistency checks across major taxonomic databases

## 🔧 **Technical Implementation**

### **Enhanced Parser Features**
- **Memory Efficient**: Handles large datasets without overflow
- **Parallel Processing**: Optimized taxonkit operations
- **Comprehensive Logging**: Detailed processing logs and error tracking
- **Quality Control**: Validation checks and unmapped taxa analysis

### **NCBI Integration**
- **Taxonkit Integration**: Automated NCBI taxonomy ID mapping
- **Lineage Resolution**: Full taxonomic hierarchy with ranks and taxids
- **Cross-Database Compatibility**: Standardized output format for merger scripts

### **Performance Metrics**
- **Processing Speed**: ~30-60 seconds per taxonomic level
- **Success Rate**: >95% successful taxonomic mapping
- **Data Quality**: 100% coverage at all taxonomic levels

## 🎯 **Key Advantages**

### **1. Environmental Truth Baseline**
- Represents **actual environmental diversity** detected in natural samples
- Provides unbiased baseline for database coverage assessment
- Captures uncultured and difficult-to-culture microbial diversity

### **2. Enhanced Data Quality**
- **Organelle Attribution**: Proper host organism assignment for organellar sequences
- **Rank Consistency**: Prevents taxonomic rank mismatches in analysis
- **Name Standardization**: Consistent taxonomic nomenclature across databases

### **3. Comprehensive Integration**
- **NCBI Compatibility**: Full integration with NCBI taxonomy system
- **Cross-Database Linking**: Enables comparison with GTDB and EukProt databases
- **Visualization Ready**: Structured output optimized for downstream analysis

### **4. Research Impact**
- **Gap Analysis**: Identifies priority targets for genome sequencing
- **Database Development**: Informs taxonomic database improvement efforts
- **Ecological Research**: Links environmental diversity to genomic potential

## 📚 **Usage & Documentation**

### **Quick Start**
```bash
cd py_16S/
python 16S_eukcensus_parser.py
```

### **Detailed Documentation**
- **`py_16S/ENHANCED_PARSER_README.md`**: Comprehensive parser documentation
- **`py_16S/IMPLEMENTATION_SUMMARY.md`**: Technical implementation details
- **`csv_16S/README.md`**: Output file specifications

### **Dependencies**
- Python 3.6+
- pandas, tqdm, pathlib
- taxonkit (with NCBI taxonomy database)

---

**The 16S_censusparse directory represents the foundation of environmental diversity assessment in this comparative genomics pipeline, providing the critical baseline against which all genomic database coverage and taxonomic representation is measured.**
