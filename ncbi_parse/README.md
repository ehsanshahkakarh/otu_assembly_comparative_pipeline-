# NCBI Genomic Database Processing Pipeline

## 🎯 **Purpose & Pipeline Position**

The **ncbi_parse** directory processes **NCBI's massive genomic assembly database** (~3M genomes) into structured taxonomic datasets for comparative analysis with environmental diversity baselines. This serves as the **genomic database reference** for assessing how well cultured organisms represent environmental microbial diversity.

### **Pipeline Context**
```
NCBI GenBank Database (~3M Genomes)
                ↓
        THIS DIRECTORY (ncbi_parse)
                ↓
    Structured Genomic Database Reference
                ↓
        ┌─────────────────────────────────────────────────────────┐
        │              COMPARATIVE GENOMICS ANALYSIS              │
        │                                                         │
        │  NCBI Genomes ←→ 16S Environmental Diversity (Prokaryotic) │
        │  NCBI Genomes ←→ 18S Environmental Diversity (Eukaryotic)  │
        │  NCBI Genomes ←→ GTDB Taxonomy (Cross-validation)          │
        │  NCBI Genomes ←→ EukProt Proteins (Eukaryotic focus)       │
        └─────────────────────────────────────────────────────────┘
                ↓
        Database Merger Scripts → Coverage Analysis → Visualizations
```

### **Research Questions Addressed**
1. **Genomic Representation**: How many species/genera/families have genomic representation in NCBI?
2. **Cultivation Success**: Which taxonomic groups are well-represented vs underrepresented?
3. **Database Completeness**: What is the taxonomic breadth of NCBI's genomic collection?
4. **Cross-Database Validation**: How does NCBI taxonomy compare with GTDB and environmental surveys?

## 📊 **Data Scale & Processing Power**

### **Input Dataset Statistics**
- **2.9M+ genome assemblies** from NCBI GenBank
- **~142K unique species** across all domains of life
- **297 phyla** representing taxonomic breadth
- **Bacteria**: 1.48M genomes (Pseudomonadota dominant)
- **Viruses**: 178K genomes (emerging viral diversity)
- **Eukaryota**: 17K genomes (fungi, plants, animals)
- **Archaea**: Smaller but significant representation

### **Taxonomic Processing Scope**
- **Species-level accuracy**: Uses species_taxid for true biological diversity
- **Multi-domain coverage**: Bacteria, Archaea, Eukaryota, Viruses
- **Hierarchical mapping**: Species → Genus → Family → Phylum
- **Isolate classification**: Distinguishes cultured vs environmental samples

## 🔬 **Advanced Processing Pipeline**

### **Two-Stage Processing Architecture**

#### **Stage 1: Taxonomic Mapping** (`taxonomic_mapping/`)
```
Raw NCBI Taxids → NCBI Taxonomy Database → Taxonomic Names
```
- **Input**: 2.9M assembly records with taxids
- **Process**: NCBI taxdump parsing and hierarchical resolution
- **Output**: Taxid-to-name mapping files for phylum/family/genus levels
- **Quality Control**: Unmapped taxid analysis and error logging

#### **Stage 2: Data Aggregation** (`py_ncbi/`)
```
Assembly Data + Taxonomic Mappings → Structured Counts → Research-Ready Datasets
```
- **Unified Parser**: Single script handles all taxonomic levels
- **Dual Counting**: Genome counts AND species counts
- **Classification**: Isolate vs uncultured organism identification
- **Validation**: Cross-referencing and sanity checking

### **Species Counting Methodology**
**Core Principle**: One species_taxid = One species (avoids assembly duplication)
- **Genome Count**: Total assemblies per taxon (includes multiple strains)
- **Species Count**: Unique species per taxon (biological diversity metric)
- **Isolate Count**: Cultured organisms vs environmental samples
- **Percentage Calculations**: Relative abundance within taxonomic levels

## 📁 **Directory Structure & Outputs**

```
ncbi_parse/
├── py_ncbi/                              # 🐍 Main Processing Engine
│   ├── ncbi_parser_clean.py              # ⭐ Unified parser for all levels
│   ├── README.md                         # Detailed parser documentation
│   ├── old_parsers/                      # Legacy scripts (deprecated)
│   └── error_log/                        # Processing logs & debugging
│
├── taxonomic_mapping/                    # 🗺️ Taxid → Name Conversion
│   ├── phylum_taxid_improved.py          # Phylum-level mapping
│   ├── family_taxid_improved.py          # Family-level mapping  
│   ├── genus_taxid.py                    # Genus-level mapping
│   ├── taxdump_ncbi/                     # NCBI taxonomy database
│   ├── taxid_to_phylum.csv               # Phylum mapping results
│   ├── taxid_to_family.csv               # Family mapping results
│   ├── taxid_to_genus.csv                # Genus mapping results
│   └── error_log/                        # Unmapped taxids analysis
│
├── csv_ncbi/                             # 📈 Final Structured Outputs
│   ├── ncbi_phylum_counts.csv            # Phylum-level aggregation (297 taxa)
│   ├── ncbi_family_counts.csv            # Family-level aggregation
│   ├── ncbi_genus_counts.csv             # Genus-level aggregation
│   ├── ncbi_*_with_accessions.csv        # Detailed accession-level data
│   ├── sanity_check/                     # Quality control validation
│   └── backup/                           # Previous versions
│
├── metadata/                             # 📊 Input Data & Quality Control
│   ├── 00assembly_summary_genbank.txt    # Raw NCBI assembly data (1.3GB)
│   ├── SPECIES_COUNTING_METHOD_EXPLAINED.md # Methodology documentation
│   ├── species_count_verification.py     # Validation scripts
│   └── sanity_check/                     # Comprehensive quality analysis
│
└── README.md                             # This comprehensive overview
```

## 📈 **Output File Specifications**

### **Standardized Column Structure**
All count files contain:
- **`{level}`**: Taxonomic name (phylum/family/genus)
- **`domain`**: Biological domain (Bacteria/Archaea/Eukaryota/Viruses)
- **`{level}_genome_count`**: Total genome assemblies
- **`{level}_genome_percentage`**: Percentage of total genomes
- **`{level}_species_count`**: Unique species count
- **`{level}_species_percentage`**: Percentage of total species
- **`taxid`**: NCBI taxonomy ID
- **`lineage`**: Full taxonomic lineage
- **`lineage_ranks`**: Corresponding taxonomic ranks
- **`lineage_taxids`**: NCBI taxids for each lineage level

### **Key Output Statistics**
- **Pseudomonadota**: 1.48M genomes (51.35%), 33K species (23.55%)
- **Bacillota**: 605K genomes (21.0%), 14K species (10.16%)
- **Ascomycota**: 17K genomes (0.6%), 4.6K species (3.24%)
- **Viral Diversity**: 178K genomes across multiple phyla

## 🔗 **Downstream Integration**

### **Database Merger Scripts** (`../Eukcensus_merge/`)
NCBI structured data serves as genomic reference for:
- **16S-NCBI Merger**: Prokaryotic environmental vs genomic comparison
- **18S-NCBI Merger**: Eukaryotic environmental vs genomic comparison
- **Cross-validation**: NCBI vs GTDB taxonomic consistency checks

### **Visualization Pipeline** (`../visuals/`)
NCBI data enables comprehensive visual analysis:
- **Coverage Assessment**: Environmental diversity vs genomic representation
- **Cultivation Success**: Well-represented vs underrepresented taxa
- **Domain Comparisons**: Bacteria vs Archaea vs Eukaryota representation
- **Temporal Analysis**: Database growth and taxonomic expansion

## 🚀 **Usage & Execution**

### **Complete Pipeline Execution**
```bash
# Stage 1: Generate taxonomic mappings
cd taxonomic_mapping/
python phylum_taxid_improved.py
python family_taxid_improved.py  
python genus_taxid.py

# Stage 2: Process and aggregate data
cd ../py_ncbi/
python ncbi_parser_clean.py --level phylum
python ncbi_parser_clean.py --level family
python ncbi_parser_clean.py --level genus
```

### **Quality Control**
```bash
# Validate species counting methodology
cd metadata/
python species_count_verification.py

# Run comprehensive sanity checks
cd ../csv_ncbi/
python ncbi_sanity_check.py
```

---

**The ncbi_parse directory transforms NCBI's massive genomic database into research-ready taxonomic datasets, providing the genomic reference baseline for comparative analysis with environmental microbial diversity and enabling quantification of cultivation success across the tree of life.**
