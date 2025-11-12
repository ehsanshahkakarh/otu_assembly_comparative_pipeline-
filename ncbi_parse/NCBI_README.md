# NCBI Genomic Database Parser & Analysis Pipeline

This directory contains a comprehensive pipeline for processing, analyzing, and preparing NCBI genomic data for comparative taxonomic analysis. The pipeline transforms raw NCBI assembly data into structured taxonomic datasets compatible with downstream comparative genomics workflows.

## 🎯 **Purpose & Objectives**

### **Primary Goal**
Process NCBI's massive genomic assembly database (~3M genomes) into taxonomically organized datasets for comparison with:
- **GTDB** (Genome Taxonomy Database)
- **EukProt** (Eukaryotic protein database)
- **Environmental census data** (eukcensus, 16S/18S rRNA surveys)

### **Why These Three Taxonomic Levels?**
The pipeline focuses on **phylum**, **family**, and **genus** levels because:

1. **EukCensus Compatibility**: Environmental census datasets (eukcensus) are organized by these exact taxonomic ranks
2. **Comparative Analysis**: Enables direct comparison across databases at standardized taxonomic levels
3. **Biological Relevance**: These ranks provide meaningful ecological and evolutionary groupings
4. **Statistical Power**: Sufficient resolution for diversity analysis while maintaining adequate sample sizes

### **Why Taxonomic Mapping is Essential**
NCBI assembly data contains only **taxids** (numeric taxonomy identifiers), not taxonomic names. The mapping process:
- **Converts taxids → taxonomic names** using NCBI's taxonomy database
- **Resolves hierarchical relationships** (species → genus → family → phylum)
- **Handles complex cases** (missing ranks, Candidatus taxa, environmental samples)
- **Enables cross-database comparisons** with standardized taxonomic nomenclature

## 📁 **Directory Structure**

```
00ncbi_parse/
├── py_ncbi/                           # 🐍 Main parsing engine
│   ├── ncbi_parser_clean.py           # ⭐ Unified parser for all taxonomic levels
│   ├── README.md                      # Detailed parser documentation
│   ├── old_parsers/                   # Legacy individual parsers (deprecated)
│   ├── alt_parse/                     # Alternative implementations (polars, dask)
│   └── error_log/                     # Processing logs and debugging info
│
├── taxonomic_mapping/                 # 🗺️ Taxid → Taxonomic name conversion
│   ├── phylum_taxid_improved.py       # Phylum-level mapping script
│   ├── family_taxid_improved.py       # Family-level mapping script
│   ├── genus_taxid.py                 # Genus-level mapping script
│   ├── taxdump_ncbi/                  # NCBI taxonomy database files
│   │   ├── nodes.dmp                  # Taxonomy tree structure
│   │   ├── names.dmp                  # Scientific names database
│   │   └── ...                        # Additional NCBI taxonomy files
│   ├── taxid_to_phylum.csv            # 📋 Phylum mapping output
│   ├── taxid_to_family.csv            # 📋 Family mapping output
│   ├── taxid_to_genus.csv             # 📋 Genus mapping output
│   └── error_log/                     # 🔍 Unmapped taxids analysis & quality control
│       ├── analyze_unmapped_patterns.py
│       ├── potential_loss/             # Interesting unmapped cases requiring attention
│       └── comprehensive_unmapped_summary.py
│
├── metadata/                          # 📊 Input data & statistics
│   ├── 00assembly_summary_genbank.txt # 🗃️ NCBI assembly database (~3M genomes)
│   ├── ncbi_species_statistics.py     # Species counting methodology
│   └── SPECIES_COUNTING_METHOD_EXPLAINED.md
│
└── csv_ncbi/                          # 📈 Final processed outputs
    ├── ncbi_phylum_counts.csv         # Phylum-level genome & species counts
    ├── ncbi_phylum_with_accessions.csv # Individual genome accessions by phylum
    ├── ncbi_family_counts.csv         # Family-level genome & species counts
    ├── ncbi_family_with_accessions.csv # Individual genome accessions by family
    ├── ncbi_genus_counts.csv          # Genus-level genome & species counts
    ├── ncbi_genus_with_accessions.csv # Individual genome accessions by genus
    └── backup/                        # Previous versions & backups
```

## 🔄 **Complete Processing Workflow**

### **Step 1: Taxonomic Mapping (Prerequisites)**
```bash
cd taxonomic_mapping/
python phylum_taxid_improved.py    # Creates taxid_to_phylum.csv
python family_taxid_improved.py    # Creates taxid_to_family.csv
python genus_taxid.py              # Creates taxid_to_genus.csv
```

**What happens:**
- Processes NCBI taxonomy database (`nodes.dmp`, `names.dmp`)
- Maps ~3M taxids to taxonomic names at each level
- Handles complex cases (missing ranks, Candidatus taxa)
- Generates error logs for unmapped entries

### **Step 2: Genome Data Processing**
```bash
cd py_ncbi/
python ncbi_parser_clean.py --level phylum
python ncbi_parser_clean.py --level family
python ncbi_parser_clean.py --level genus
```

**What happens:**
- Loads NCBI assembly data (~3M genomes)
- Merges with taxonomic mappings
- Classifies genomes (isolate vs uncultured)
- Generates counts and accession files

### **Step 3: Quality Control & Analysis**
```bash
cd taxonomic_mapping/error_log/
python analyze_unmapped_patterns.py     # Categorize unmapped entries
python comprehensive_unmapped_summary.py # Generate summary reports
```

**What happens:**
- Analyzes unmapped taxids (viruses, unclassified taxa, etc.)
- Identifies true data loss vs expected unmapped cases
- Generates quality control reports

## 📊 **Why Multiple CSV Files? Understanding the Output Structure**

### **The Six Essential Output Files**

The pipeline generates **6 core CSV files** (2 per taxonomic level) because each serves a distinct purpose in downstream comparative analysis:

#### **📈 Count Files** (`*_counts.csv`)
- **Purpose**: Aggregate statistics for merger scripts and comparative analysis
- **Content**: Taxonomic groups with genome counts, species counts, percentages, and lineage data
- **Usage**: Direct input for database comparison scripts (`*_merger.py`)
- **Size**: Small (~295 phylum, ~8K family, ~35K genus entries)

#### **🗃️ Accession Files** (`*_with_accessions.csv`)
- **Purpose**: Individual genome records for detailed analysis and filtering
- **Content**: Every genome with its taxonomic classification and metadata
- **Usage**: Subset creation, isolate/uncultured analysis, quality control
- **Size**: Large (~2.9M entries each, 200-250MB files)

### **Why These Specific Taxonomic Levels?**

#### **🎯 EukCensus Alignment**
Environmental census datasets are organized by:
- **Phylum level**: Broad taxonomic diversity (e.g., Proteobacteria, Firmicutes)
- **Family level**: Ecological functional groups (e.g., Enterobacteriaceae, Rhizobiaceae)
- **Genus level**: Fine-scale diversity patterns (e.g., Escherichia, Bacillus)

#### **🔬 Comparative Analysis Requirements**
```
NCBI Phylum ←→ GTDB Phylum ←→ EukCensus Phylum
NCBI Family ←→ GTDB Family ←→ EukCensus Family
NCBI Genus  ←→ GTDB Genus  ←→ EukCensus Genus
```

#### **📊 Statistical Considerations**
- **Phylum**: Sufficient sample sizes for robust statistics (~295 groups)
- **Family**: Balance between resolution and statistical power (~8K groups)
- **Genus**: High resolution for diversity analysis (~35K groups)

## 🗂️ **Detailed File Descriptions**

### **Input Files**

#### **`metadata/00assembly_summary_genbank.txt`**
- **Source**: NCBI GenBank assembly summary (updated monthly)
- **Size**: ~3M assemblies, ~50GB uncompressed
- **Format**: Tab-separated with header `#assembly_accession`
- **Key columns**: `assembly_accession`, `taxid`, `organism_name`, `assembly_level`, `excluded_from_refseq`
- **Purpose**: Primary source of genome metadata and taxid associations

#### **`taxonomic_mapping/taxdump_ncbi/nodes.dmp`**
- **Source**: NCBI taxonomy database dump
- **Format**: `tax_id | parent_tax_id | rank | ...` (pipe-delimited)
- **Purpose**: Hierarchical taxonomy tree structure
- **Usage**: Traverse from species → genus → family → phylum

#### **`taxonomic_mapping/taxdump_ncbi/names.dmp`**
- **Source**: NCBI taxonomy database dump
- **Format**: `tax_id | name_txt | unique_name | name_class | ...`
- **Purpose**: Scientific names for taxonomy nodes
- **Usage**: Convert taxids to readable taxonomic names

### **Output Files Specifications**

#### **📈 Count Files** (`csv_ncbi/ncbi_{level}_counts.csv`)

**Columns:**
- `{level}`: Taxonomic name (e.g., "Enterobacteriaceae", "Pseudomonadota")
- `domain`: Bacteria, Archaea, Eukaryota, or Viruses
- `{level}_genome_count`: Total genomes in this taxonomic group
- `{level}_genome_percentage`: Percentage of total genomes
- `{level}_species_count`: Unique species in this taxonomic group
- `{level}_species_percentage`: Percentage of total species
- `taxid`: NCBI taxid for this taxonomic group
- `lineage`: Full taxonomic lineage (semicolon-separated)
- `lineage_ranks`: Corresponding ranks (semicolon-separated)
- `lineage_taxids`: Lineage taxids (semicolon-separated)

**Example** (`ncbi_family_counts.csv`):
```csv
family,domain,family_genome_count,family_genome_percentage,family_species_count,family_species_percentage,taxid,lineage,lineage_ranks,lineage_taxids
Enterobacteriaceae,Bacteria,485123,16.84,4021,2.84,543,cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae,no rank;superkingdom;phylum;class;order;family,131567;2;1224;1236;91347;543
```

#### **🗃️ Accession Files** (`csv_ncbi/ncbi_{level}_with_accessions.csv`)

**Columns:**
- `accession_clean`: Assembly accession (e.g., "GCA_000005825.2")
- `{level}`: Taxonomic classification
- `domain`: Bacteria, Archaea, Eukaryota, or Viruses
- `organism_name`: Full organism name from NCBI
- `assembly_level`: Complete Genome, Chromosome, Scaffold, or Contig
- `genome_source`: **isolate** or **uncultured** (automated classification)
- `taxid`: NCBI taxid for this specific organism

**Example** (`ncbi_family_with_accessions.csv`):
```csv
accession_clean,family,domain,organism_name,assembly_level,genome_source,taxid
GCA_000005825.2,Enterobacteriaceae,Bacteria,Escherichia coli str. K-12 substr. MG1655,Complete Genome,isolate,511145
GCA_000006745.1,Enterobacteriaceae,Bacteria,Escherichia coli CFT073,Complete Genome,isolate,199310
```

## 🔗 **Downstream Integration & Usage**

### **Database Merger Scripts**
The CSV outputs are specifically designed for comparative analysis:

```bash
# Example merger script usage
python 16s_ncbi_merger.py    # Uses ncbi_phylum_counts.csv + ncbi_family_counts.csv + ncbi_genus_counts.csv
python 18s_eukprot_merger.py # Uses ncbi_phylum_counts.csv for eukaryotic comparisons
```

**Integration points:**
- **Count files** → Direct input for merger scripts
- **Accession files** → Subset creation and detailed analysis
- **Species counts** → Accurate comparative metrics (not just genome counts)

### **Visualization Pipeline**
```bash
# Alluvial/Sankey plots showing taxonomic flow between databases
Rscript create_alluvial_plots.R

# Scatter plots comparing coverage and novelty factors
Rscript create_scatter_comparisons.R
```

**Data flow:**
- **NCBI counts** → **Merger analysis** → **Visualization** → **Publication figures**

### **Quality Control Integration**
```bash
# Consistency checks across taxonomic levels
python check_unmapped_in_counts.py

# Verify no data loss between processing steps
python comprehensive_unmapped_summary.py
```

## 🧬 **Genome Source Classification**

### **Isolate vs Uncultured Classification**
The parser automatically classifies each genome:

**Isolate** (~81% of genomes):
- Cultured, characterized strains
- Laboratory isolates with strain names
- Type strains and reference genomes

**Uncultured** (~19% of genomes):
- Metagenome-assembled genomes (MAGs)
- Single amplified genomes (SAGs)
- Environmental samples
- Unclassified/unknown organisms

**Classification logic:**
```python
# Patterns that trigger 'uncultured' classification
uncultured_patterns = [
    'uncultured', 'environmental', 'metagenome', 'unclassified',
    'unknown', 'unidentified', 'mixed culture', 'enrichment culture',
    'derived from metagenome', 'metagenome-assembled', 'mag',
    'single amplified genome', 'sag', 'environmental sample'
]
```

**Why this matters:**
- **Comparative analysis**: Different databases have different isolate/uncultured ratios
- **Quality assessment**: Isolates generally have higher quality assemblies
- **Ecological interpretation**: Uncultured genomes represent environmental diversity

## 🔬 **Data Processing Logic & Methodology**

### **Taxonomic Mapping Process**
```
Raw NCBI Data → Taxid Extraction → Hierarchy Traversal → Name Resolution → Domain Assignment
```

**Detailed steps:**
1. **Load NCBI taxonomy** (`nodes.dmp`, `names.dmp`) into memory-efficient dictionaries
2. **Build hierarchy maps**: parent-child relationships, ranks, scientific names
3. **For each taxid in assembly data**:
   - Traverse taxonomy tree upward until target rank found
   - Handle missing ranks (use higher-level classification)
   - Resolve domain by continuing to superkingdom level
4. **Quality control**: Filter to major domains, preserve Candidatus taxa
5. **Output mapping**: CSV files for parser consumption

### **Assembly Data Processing**
```
Assembly File → Header Parsing → Taxid Merging → Source Classification → Counting → Output Generation
```

**Detailed steps:**
1. **Parse assembly file**: Handle NCBI's complex header format
2. **Merge with taxonomy**: Join on taxid, handle data type mismatches
3. **Classify genome source**: Pattern matching for isolate/uncultured determination
4. **Extract species taxids**: Use taxonkit for species-level resolution
5. **Generate counts**: Aggregate by taxonomic groups with percentages
6. **Add lineage data**: Full taxonomic lineage using taxonkit
7. **Create outputs**: Both summary counts and detailed accession files

### **Species Counting Methodology**
**Critical distinction**: Species counts ≠ Genome counts

- **Genome counts**: Total number of assemblies per taxonomic group
- **Species counts**: Number of unique species (by species_taxid) per taxonomic group
- **Species subset**: One representative genome per species (isolate-preferential)

**Why this matters:**
- **Comparative analysis**: Different databases may have different genome:species ratios
- **Bias correction**: Prevents over-representation of well-studied species
- **Ecological relevance**: Species diversity is more meaningful than assembly abundance

## 📊 **Performance & Scale**

### **Processing Statistics**
- **Input**: 2,945,415 NCBI assemblies
- **Successfully processed**: 2,882,656 (97.9% success rate)
- **Taxonomic groups generated**:
  - Phylum: 295 groups
  - Family: ~8,000 groups
  - Genus: ~35,000 groups

### **File Sizes**
- **Count files**: 64KB - 7.4MB (small, fast loading)
- **Accession files**: 200-250MB each (detailed, comprehensive)
- **Processing time**: 30-60 seconds per taxonomic level

### **Quality Metrics**
- **Mapping success rate**: 97.9% of assemblies successfully mapped
- **Error categorization**: 93.7% of "unmapped" entries are expected (viruses, unclassified)
- **True data loss**: <4% of entries require manual attention

## 🌍 **Project Context & Role in Comparative Genomics Pipeline**

### **Position in Larger Workflow**
```
Environmental Samples → 16S/18S rRNA Sequencing → OTU Clustering → Taxonomic Assignment
                                    ↓
                            EukCensus Database
                                    ↓
                    ┌─────────────────────────────────────┐
                    │     COMPARATIVE ANALYSIS            │
                    │                                     │
NCBI Assembly → This Pipeline → NCBI Taxonomic Data ────┤
                                                         │ → Merger Scripts → Visualizations → Publications
GTDB Assembly → GTDB Pipeline → GTDB Taxonomic Data ────┤
                                                         │
EukProt Database ────────────────────────────────────────┘
```

### **Why This Pipeline Exists**
**Research Question**: How well do genomic databases (NCBI, GTDB) represent environmental microbial diversity captured in rRNA surveys?

**Key Comparisons**:
1. **NCBI vs Environmental Census**: Genome availability vs environmental detection
2. **NCBI vs GTDB**: Taxonomic classification consistency between major databases
3. **NCBI vs EukProt**: Eukaryotic representation in genomic vs protein databases
4. **Coverage Analysis**: Which environmental taxa lack genomic representation?

### **Critical Design Decisions**

#### **Why Phylum/Family/Genus Levels?**
- **Environmental data compatibility**: EukCensus organized by these exact ranks
- **Statistical robustness**: Sufficient sample sizes at each level
- **Biological relevance**: Meaningful ecological and evolutionary groupings
- **Cross-database standardization**: Enables direct comparisons

#### **Why Species Counts Matter**
- **Bias correction**: Prevents over-representation of model organisms
- **True diversity assessment**: Species richness vs assembly abundance
- **Comparative accuracy**: Different databases have different sequencing biases

#### **Why Isolate/Uncultured Classification?**
- **Data quality assessment**: Isolates typically have higher quality assemblies
- **Cultivation bias analysis**: Identifies gaps in cultured diversity
- **Environmental relevance**: Uncultured genomes represent environmental diversity

### **Downstream Applications**

#### **Scientific Publications**
- **Taxonomic coverage analysis**: Which environmental taxa lack genome representation?
- **Database comparison studies**: NCBI vs GTDB classification differences
- **Cultivation bias assessment**: Cultured vs environmental microbial diversity

#### **Database Development**
- **Gap identification**: Priority targets for genome sequencing
- **Quality control**: Consistency checks across major databases
- **Metadata enhancement**: Improved taxonomic annotations

#### **Ecological Research**
- **Environmental surveys**: Link rRNA diversity to genomic potential
- **Biogeography studies**: Geographic patterns in genome availability
- **Functional analysis**: Connect taxonomic diversity to metabolic potential

## 🔧 **Technical Implementation Notes**

### **Memory Efficiency**
- **Streaming processing**: Handle 3M+ assemblies without memory overflow
- **Chunked operations**: Process large files in manageable pieces
- **Optimized data structures**: Efficient dictionaries for taxonomy lookups

### **Error Handling**
- **Graceful degradation**: Continue processing despite individual failures
- **Comprehensive logging**: Track all unmapped entries for quality control
- **Validation checks**: Ensure data consistency across processing steps

### **Scalability**
- **Modular design**: Easy to add new taxonomic levels or databases
- **Parallel processing**: Taxonkit operations can be parallelized
- **Incremental updates**: Process only new assemblies in database updates

This pipeline represents a critical component in understanding the relationship between environmental microbial diversity and available genomic resources, enabling researchers to identify gaps in our genomic knowledge and prioritize future sequencing efforts.
