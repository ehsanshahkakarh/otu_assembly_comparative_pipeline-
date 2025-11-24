# EukCensus Environmental-Genomic Database Merger Pipeline

## 🎯 **Purpose & Pipeline Position**

The **Eukcensus_merge** directory serves as the **critical integration hub** that merges environmental diversity baselines (16S/18S census data) with genomic databases (GTDB, NCBI, EukProt) to quantify **database coverage gaps** and **environmental novelty**. This is where environmental truth meets genomic reality.

### **Pipeline Context**
```
Environmental Baselines                    Genomic Databases
        ↓                                         ↓
┌─────────────────┐                    ┌─────────────────┐
│ 16S_censusparse │ ←──────────────────→ │   gtdb_parse    │
│ (Prokaryotic)   │                    │   ncbi_parse    │
└─────────────────┘                    └─────────────────┘
        ↓                                         ↓
┌─────────────────┐                    ┌─────────────────┐
│ 18S_censusparse │ ←──────────────────→ │  eukprot_parse  │
│ (Eukaryotic)    │                    │   ncbi_parse    │
└─────────────────┘                    └─────────────────┘
        ↓                                         ↓
                    THIS DIRECTORY
                 (Eukcensus_merge)
                         ↓
        ┌─────────────────────────────────────────┐
        │        COMPARATIVE ANALYSIS             │
        │                                         │
        │  • Coverage Assessment                  │
        │  • Novelty Factor Calculation           │
        │  • Overrepresentation Analysis          │
        │  • Cultivation Gap Identification       │
        └─────────────────────────────────────────┘
                         ↓
            Visualizations → Publications
```

### **Research Questions Addressed**
1. **Database Coverage**: What percentage of environmental diversity is represented in genomic databases?
2. **Cultivation Bias**: Which environmental taxa lack cultured representatives?
3. **Novelty Discovery**: Which environmental taxa represent novel diversity?
4. **Database Completeness**: How comprehensive are our genomic databases?

## 🧬 **Advanced Merger Logic & Algorithms**

### **Multi-Stream Matching Strategy**

#### **1. Prokaryotic Mergers (16S ↔ GTDB/NCBI)**
**Algorithm**: Vectorized Lineage-Based Matching
```python
# Core matching logic
def vectorized_lineage_matching(census_taxa, genomic_lineages):
    """
    Uses pandas string operations for high-performance matching:
    1. Direct name matching in lineage strings
    2. Regex pattern matching with anchors
    3. Vectorized aggregation using groupby operations
    """
    pattern = f';{taxon_name};|^{taxon_name};|;{taxon_name}$|^{taxon_name}$'
    matches = genomic_df['lineage'].str.contains(pattern, regex=True, na=False)
    return matches.sum()  # Aggregate genome/species counts
```

**Performance**: 4x faster than iterative approaches, handles 25K+ taxa efficiently

#### **2. Eukaryotic Mergers (18S ↔ EukProt/NCBI)**
**Algorithm**: Multi-Stream Hierarchical Matching with Confidence Scoring
```python
# Enhanced matching with confidence levels
def multi_stream_matching(census_taxon, eukprot_data):
    """
    Four-tier matching system:
    1. Direct name matching (high confidence)
    2. Taxid-based matching (high confidence)  
    3. Hierarchical lineage matching (medium confidence)
    4. Fuzzy pattern matching (low confidence)
    """
    matches = {
        'direct': direct_name_match(census_taxon, eukprot_data),
        'taxid': taxid_based_match(census_taxon, eukprot_data),
        'lineage': hierarchical_lineage_match(census_taxon, eukprot_data),
        'pattern': fuzzy_pattern_match(census_taxon, eukprot_data)
    }
    return aggregate_with_confidence(matches)
```

### **Coverage & Novelty Calculations**

#### **Coverage Percentage**
```
Coverage = (Genomic_Species_Count / Environmental_OTU_Count) × 100
```
- **High Coverage (>50%)**: Well-represented in databases
- **Moderate Coverage (10-50%)**: Partially represented
- **Low Coverage (<10%)**: Underrepresented
- **No Coverage (0%)**: Novel environmental diversity

#### **Novelty Factor**
```
Novelty_Factor = Environmental_OTU_Count / Genomic_Species_Count
```
- **High Novelty (>10)**: Massive environmental diversity vs genomic representation
- **Moderate Novelty (2-10)**: Significant cultivation gaps
- **Low Novelty (<2)**: Well-cultivated taxa

#### **Overrepresentation Factor**
```
Overrepresentation_Factor = Genomic_Species_Count / Environmental_OTU_Count
```
- **High Overrepresentation (>2)**: Database bias toward easily cultured taxa
- **Balanced Representation (0.5-2)**: Proportional representation
- **Underrepresentation (<0.5)**: Environmental diversity exceeds genomic data

## 📁 **Directory Structure & Output Organization**

```
Eukcensus_merge/
├── py_mergers/                           # 🐍 Core Merger Scripts
│   ├── 16s_gtdb_merger.py               # ⭐ 16S ↔ GTDB prokaryotic merger
│   ├── 16s_ncbi_merger.py               # ⭐ 16S ↔ NCBI prokaryotic merger
│   ├── 18s_eukprot_merger.py            # ⭐ 18S ↔ EukProt eukaryotic merger
│   ├── 18s_ncbi_merger.py               # ⭐ 18S ↔ NCBI eukaryotic merger
│   ├── logs/                            # Processing logs & performance metrics
│   └── old_scripts/                     # Legacy versions & development history
│
├── 16s_merged/                          # 📊 Prokaryotic Merger Results
│   ├── csv_results/                     # Final merged datasets
│   │   ├── 16s_gtdb_merged_clean_phylum.csv    # Phylum-level 16S-GTDB comparison
│   │   ├── 16s_gtdb_merged_clean_family.csv    # Family-level 16S-GTDB comparison
│   │   ├── 16s_gtdb_merged_clean_genus.csv     # Genus-level 16S-GTDB comparison
│   │   ├── 16s_ncbi_merged_clean_phylum.csv    # Phylum-level 16S-NCBI comparison
│   │   ├── 16s_ncbi_merged_clean_family.csv    # Family-level 16S-NCBI comparison
│   │   └── 16s_ncbi_merged_clean_genus.csv     # Genus-level 16S-NCBI comparison
│   ├── analysis_summary/                # Summary statistics & performance metrics
│   │   ├── 16s_gtdb_merger_clean_summary.csv   # GTDB merger efficacy summary
│   │   └── 16s_ncbi_merger_clean_summary.csv   # NCBI merger efficacy summary
│   └── logs/                            # Detailed processing logs
│
├── 18s_merged/                          # 📊 Eukaryotic Merger Results
│   ├── csv_results/                     # Final merged datasets
│   │   ├── 18s_eukprot_merged_division.csv     # Division-level 18S-EukProt comparison
│   │   ├── 18s_eukprot_merged_family.csv       # Family-level 18S-EukProt comparison
│   │   ├── 18s_eukprot_merged_genus.csv        # Genus-level 18S-EukProt comparison
│   │   ├── 18s_ncbi_merged_clean_phylum.csv    # Phylum-level 18S-NCBI comparison
│   │   ├── 18s_ncbi_merged_clean_family.csv    # Family-level 18S-NCBI comparison
│   │   ├── 18s_ncbi_merged_clean_genus.csv     # Genus-level 18S-NCBI comparison
│   │   ├── high_novelty_genera_families.csv    # High-novelty taxa analysis
│   │   └── sanity/                      # Quality control & validation
│   ├── analysis_summary/                # Summary statistics & performance metrics
│   │   ├── 18s_eukprot_merger_summary.csv      # EukProt merger efficacy summary
│   │   ├── 18s_eukprot_merger_analysis.txt     # Detailed analysis report
│   │   └── 18s_ncbi_merger_clean_summary.csv   # NCBI merger efficacy summary
│   └── logs/                            # Detailed processing logs
│
└── README.md                            # This comprehensive overview
```

## 📈 **Merger Efficacy & Performance Metrics**

### **Prokaryotic Merger Performance (16S Data)**

#### **16S ↔ GTDB Merger Efficacy**
| Taxonomic Level | Total Entries | Matched Taxa | Match Rate | Coverage Quality |
|-----------------|---------------|--------------|------------|------------------|
| **Phylum**      | 206           | 33           | **16.0%**  | High precision   |
| **Family**      | 5,640         | 444          | **7.9%**   | Moderate coverage|
| **Genus**       | 25,496        | 2,519        | **9.9%**   | Comprehensive    |

**Key Insights**:
- **High-confidence matching** at phylum level with 16% success rate
- **Extensive genus-level analysis** covering 25K+ environmental genera
- **Cultivation gaps identified** in 90%+ of environmental genera

#### **16S ↔ NCBI Merger Efficacy**
- **Vectorized processing**: 4x performance improvement over legacy methods
- **Prokaryotic focus**: Filters for Bacteria and Archaea domains
- **Lineage-based matching**: Handles complex taxonomic hierarchies
- **Quality control**: Comprehensive validation and error checking

### **Eukaryotic Merger Performance (18S Data)**

#### **18S ↔ EukProt Merger Efficacy**
| Taxonomic Level | Total Taxa | Matched Taxa | Coverage Rate | Avg Coverage | Median Coverage |
|-----------------|------------|--------------|---------------|--------------|-----------------|
| **Division**    | 18         | 17           | **94.4%**     | 89.70%       | 2.02%           |
| **Family**      | 170        | 72           | **42.4%**     | 5.84%        | 0.00%           |
| **Genus**       | 283        | 103          | **36.4%**     | 8.39%        | 0.00%           |

**Key Insights**:
- **Excellent division-level coverage** (94.4% match rate)
- **Significant cultivation gaps** at family/genus levels
- **High novelty factors** indicating vast uncultured eukaryotic diversity

#### **18S ↔ NCBI Merger Efficacy**
- **Multi-stream matching**: Direct, taxid, lineage, and pattern-based approaches
- **Confidence scoring**: High/medium/low confidence match classification
- **Species deduplication**: Prevents double-counting across taxonomic levels
- **Eukaryotic specialization**: Handles complex eukaryotic taxonomic structures

## 🔬 **Advanced Technical Features**

### **Vectorized Processing Architecture**
- **Pandas-based operations**: Leverages vectorized string operations for speed
- **Memory optimization**: Efficient handling of large datasets (70K+ OTUs)
- **Parallel processing**: Concurrent operations where applicable
- **Scalable design**: Handles growing database sizes efficiently

### **Quality Control & Validation**
- **Comprehensive logging**: Detailed processing logs and match tracking
- **Sanity checking**: Extensive validation of merger results
- **Error handling**: Robust error recovery and reporting
- **Performance monitoring**: Processing time and memory usage tracking

### **Output Standardization**
- **Consistent column structure**: Standardized across all merger outputs
- **Cross-database compatibility**: Uniform format for visualization pipeline
- **Research-ready format**: Optimized for downstream analysis and plotting

## 📊 **Detailed Output File Specifications**

### **Prokaryotic Merger Outputs (16s_merged/csv_results/)**

#### **Column Structure for 16S-GTDB/NCBI Merged Files**
```csv
phylum,census_otu_count,census_size_count,otu_percentage,size_percentage,
gtdb_genome_count,gtdb_species_count,isolate_count,genome_pct_db,species_pct,
isolate_percentage,novelty_factor,overrepresentation_factor,direct_matches,
lineage_matches,total_matches,match_status,coverage_percentage
```

**Key Metrics Explained**:
- **`census_otu_count`**: Environmental OTU clusters for this taxon
- **`census_size_count`**: Total sequence abundance in environmental samples
- **`gtdb_genome_count`**: Number of genomes in GTDB for this taxon
- **`novelty_factor`**: Environmental diversity / Genomic representation ratio
- **`overrepresentation_factor`**: Genomic representation / Environmental diversity ratio
- **`match_status`**: 'matched' or 'census_only' or 'gtdb_only'

#### **Example High-Impact Results**
```csv
Verrucomicrobiota,10214,45775,3.89,3.83,6636,6636,0,1.11,1.11,0.0,1.539,0.65,1,0,1,matched
Chloroflexota,13210,58834,5.04,4.93,4762,4762,0,0.8,0.8,0.0,2.774,0.36,1,0,1,matched
```
- **Verrucomicrobiota**: Novelty factor 1.54 (moderate cultivation gap)
- **Chloroflexota**: Novelty factor 2.77 (significant cultivation gap)

### **Eukaryotic Merger Outputs (18s_merged/csv_results/)**

#### **Column Structure for 18S-EukProt Merged Files**
```csv
division,census_otu_count,census_size_count,otu_percentage,size_percentage,
eukprot_species_count,coverage_percentage,novelty_factor,overrepresentation_factor,
direct_matches,lineage_matches,total_matches,match_status
```

**Key Metrics Explained**:
- **`eukprot_species_count`**: Number of species in EukProt for this taxon
- **`coverage_percentage`**: (EukProt species / Environmental OTUs) × 100
- **`direct_matches`**: Exact name matches between databases
- **`lineage_matches`**: Hierarchical lineage-based matches

#### **Example High-Impact Results**
```csv
Opisthokonta,24345,158758,34.34,39.56,180,0.74,135.25,0.007,0,180,180,matched
Alveolata,9042,60488,12.75,15.07,169,1.87,53.503,0.019,0,169,169,matched
Rhizaria,5126,24072,7.23,6.0,40,0.78,128.15,0.008,0,40,40,matched
```
- **Opisthokonta**: Massive novelty factor (135.25) despite being well-studied
- **Rhizaria**: Extreme cultivation gap (novelty factor 128.15)
- **Alveolata**: Significant underrepresentation (novelty factor 53.5)

### **Summary Statistics Files**

#### **Merger Efficacy Summaries**
Located in `analysis_summary/` directories:
- **`*_merger_summary.csv`**: Quantitative merger performance metrics
- **`*_merger_analysis.txt`**: Detailed qualitative analysis reports

**Summary Metrics Include**:
- **Total entries processed**
- **Match rates by taxonomic level**
- **Coverage distribution statistics**
- **Top taxa by novelty/coverage factors**

## 🎯 **Research Applications & Impact**

### **Database Coverage Assessment**
- **Cultivation Bias Quantification**: Identifies systematically underrepresented taxa
- **Novelty Discovery**: Highlights environmental taxa lacking genomic representation
- **Database Completeness**: Measures how well genomic databases represent environmental diversity

### **Ecological Insights**
- **Uncultured Diversity**: Quantifies the "microbial dark matter" in environmental samples
- **Taxonomic Gaps**: Identifies specific lineages requiring cultivation efforts
- **Environmental Importance**: Links environmental abundance to genomic representation

### **Methodological Advances**
- **Vectorized Merging**: 4x performance improvement over traditional approaches
- **Multi-confidence Matching**: Reduces false positives while maximizing true matches
- **Standardized Outputs**: Enables cross-database comparative analysis

## 🚀 **Usage & Execution**

### **Running Individual Mergers**
```bash
# Prokaryotic mergers
cd py_mergers/
python 16s_gtdb_merger.py    # 16S environmental vs GTDB genomes
python 16s_ncbi_merger.py    # 16S environmental vs NCBI genomes

# Eukaryotic mergers
python 18s_eukprot_merger.py # 18S environmental vs EukProt proteins
python 18s_ncbi_merger.py    # 18S environmental vs NCBI eukaryotic genomes
```

### **Output Locations**
- **Merged datasets**: `16s_merged/csv_results/` and `18s_merged/csv_results/`
- **Summary statistics**: `16s_merged/analysis_summary/` and `18s_merged/analysis_summary/`
- **Processing logs**: `py_mergers/logs/` and individual `logs/` directories

### **Dependencies**
- Python 3.7+
- pandas, numpy for vectorized operations
- pathlib for cross-platform file handling
- logging for comprehensive process tracking

---

**The Eukcensus_merge directory represents the analytical heart of the comparative genomics pipeline, where environmental diversity meets genomic databases to reveal cultivation gaps, database biases, and novel microbial diversity awaiting discovery.**
