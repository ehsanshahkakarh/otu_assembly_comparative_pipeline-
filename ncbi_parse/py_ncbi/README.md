# NCBI Taxonomy Parser

This directory contains a unified Python parser for processing NCBI assembly data and extracting taxonomic information at different hierarchical levels (phylum, family, genus) with genome source classification.

## Main Parser

**`ncbi_parser_clean.py`** - Unified, clean parser for all taxonomic levels

### Usage

#### Basic Usage (Auto-detect directories)
```bash
python ncbi_parser_clean.py --level phylum
python ncbi_parser_clean.py --level family
python ncbi_parser_clean.py --level genus
```

#### Advanced Usage (Custom directories)
```bash
# Custom output directory
python ncbi_parser_clean.py --level phylum --output-dir /path/to/custom/output

# Custom input directory (for assembly file)
python ncbi_parser_clean.py --level family --input-dir /path/to/custom/input

# Both custom directories
python ncbi_parser_clean.py --level genus --input-dir /path/to/input --output-dir /path/to/output
```

#### Command-line Options
- `--level` (required): Taxonomic level to process (`phylum`, `family`, or `genus`)
- `--input-dir` (optional): Input directory containing assembly file (default: auto-detect `../metadata/`)
- `--output-dir` (optional): Output directory for results (default: `../csv_ncbi/`)

#### Examples
```bash
# Process phylum level with default directories
python ncbi_parser_clean.py --level phylum

# Process family level with custom output to current directory
python ncbi_parser_clean.py --level family --output-dir ./

# Process genus level with custom paths
python ncbi_parser_clean.py --level genus --input-dir /data/ncbi --output-dir /results/ncbi
```

### Key Features
- **Unified codebase**: Single script handles all taxonomic levels
- **Genome source classification**: Distinguishes isolates from uncultured/environmental samples
- **Species-level processing**: Extracts species taxids and creates species subsets
- **Lineage integration**: Adds full taxonomic lineage using taxonkit
- **Clean output**: Minimal terminal output with essential progress information

## Directory Structure

```
py_ncbi/
├── ncbi_parser_clean.py        # Main unified parser
├── csv_ncbi/                   # Output directory
│   ├── ncbi_phylum_counts.csv
│   ├── ncbi_phylum_with_accessions.csv
│   ├── ncbi_family_counts.csv
│   ├── ncbi_family_with_accessions.csv
│   ├── ncbi_genus_counts.csv
│   └── ncbi_genus_with_accessions.csv
├── alt_parse/                  # Legacy/alternative parsers
├── old_parsers/                # Deprecated individual parsers
└── error_log/                  # Processing logs
```

## Input Files

### 1. NCBI Assembly File
- **`../metadata/00assembly_summary_genbank.txt`**: NCBI GenBank assembly summary
- **Format**: Tab-separated values with header starting with `#assembly_accession`
- **Key Columns**: `assembly_accession`, `taxid`, `organism_name`, `assembly_level`, `isolate`, `excluded_from_refseq`
- **Size**: ~3M entries, ~50GB file
- **Example**: `GCA_000005825.2\t511145\tEscherichia coli str. K-12 substr. MG1655\tComplete Genome\tna\t`

### 2. Taxonomic Mapping File
- **`../taxonomic_mapping/taxid_to_{level}.csv`**: Maps taxids to taxonomic classifications
- **Format**: CSV with columns `taxid`, `domain`, `{level}`
- **Levels**: phylum, family, genus
- **Example**: `511145,Bacteria,Enterobacteriaceae`
- **Coverage**: ~224K taxid mappings

### 3. External Dependencies
- **taxonkit**: Command-line tool for NCBI taxonomy operations
- **NCBI taxdump**: Built-in database used by taxonkit for lineage resolution
- **Purpose**: Provides name→taxid mapping and full lineage information

## Output Files

### 1. Summary Files (Counts)
Aggregated genome and species counts by taxonomic group with full lineage information.

**Files**: `csv_ncbi/ncbi_{level}_counts.csv`

**Columns**:
- `{level}`: Taxonomic name (e.g., Enterobacteriaceae, Pseudomonadota)
- `domain`: Bacteria, Archaea, Eukaryota, or Viruses
- `{level}_genome_count`: Total number of genomes in this group
- `{level}_species_count`: Number of unique species in this group
- `taxid`: NCBI taxid for this taxonomic group
- `lineage`: Full taxonomic lineage (semicolon-separated)
- `lineage_ranks`: Taxonomic ranks (semicolon-separated)
- `lineage_taxids`: Lineage taxids (semicolon-separated)

**Example** (`ncbi_family_counts.csv`):
```csv
family,domain,family_genome_count,family_species_count,taxid,lineage,lineage_ranks,lineage_taxids
Enterobacteriaceae,Bacteria,485123,4021,543,cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae,no rank;superkingdom;phylum;class;order;family,131567;2;1224;1236;91347;543
Streptomycetaceae,Bacteria,89456,6996,2062,cellular organisms;Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae,no rank;superkingdom;phylum;class;order;family,131567;2;201174;1760;85011;2062
```

### 2. Detailed Files (Accessions)
Individual genome accessions with taxonomic classifications and genome source information.

**Files**: `csv_ncbi/ncbi_{level}_with_accessions.csv`

**Columns**:
- `accession_clean`: Assembly accession ID (e.g., GCA_000005825.2)
- `{level}`: Taxonomic name
- `domain`: Bacteria, Archaea, Eukaryota, or Viruses
- `organism_name`: Full organism name from NCBI
- `assembly_level`: Complete Genome, Chromosome, Scaffold, or Contig
- `genome_source`: **isolate** or **uncultured** (see classification below)
- `taxid`: NCBI taxid

**Example** (`ncbi_family_with_accessions.csv`):
```csv
accession_clean,family,domain,organism_name,assembly_level,genome_source,taxid
GCA_000005825.2,Enterobacteriaceae,Bacteria,Escherichia coli str. K-12 substr. MG1655,Complete Genome,isolate,511145
GCA_000006745.1,Enterobacteriaceae,Bacteria,Escherichia coli CFT073,Complete Genome,isolate,199310
GCA_000007405.1,Enterobacteriaceae,Bacteria,Salmonella enterica subsp. enterica serovar Typhi str. CT18,Complete Genome,isolate,220341
```

## Genome Source Classification

The parser automatically classifies genomes as **isolate** or **uncultured** based on organism names and exclusion notes:

### Classification Logic
1. **Default**: All genomes start as `isolate`
2. **Reclassified to `uncultured`** if organism name or excluded_from_refseq contains:
   - `uncultured`, `environmental`, `metagenome`, `unclassified`
   - `unknown`, `unidentified`, `mixed culture`, `enrichment culture`
   - `derived from metagenome`, `metagenome-assembled`, `mag`
   - `single amplified genome`, `sag`, `environmental sample`

### Special Cases
- **Enrichment cultures with strain names** remain as `isolate`
- **MAGs and SAGs** are classified as `uncultured`
- **Environmental samples** are classified as `uncultured`

### Typical Results
- **~81% isolate**: Cultured, characterized strains
- **~19% uncultured**: MAGs, SAGs, environmental samples, metagenomes

## Processing Pipeline

### 1. Data Loading & Merging
```
Assembly File (3M entries) → Parse Header → Merge with Taxonomy Mapping → 2.9M matched entries
```

### 2. Genome Source Classification
```
Organism Names + Exclusion Notes → Pattern Matching → Isolate/Uncultured Assignment
```

### 3. Species Processing
```
Taxids → Species Taxid Extraction → Species Subset Creation (Isolate-Preferential)
```

### 4. Taxonomic Counting
```
Grouped Data → Genome Counts + Species Counts → Domain Separation
```

### 5. Lineage Integration
```
Taxonomic Names → taxonkit name2taxid → taxonkit lineage → Full Lineage Data
```

## Script Functions & Architecture

### **🏗️ Class Constructor & Setup**
- **`__init__(level, input_dir, output_dir)`** - Initialize parser with taxonomic level and directories
- **`_find_input_dir()`** - Auto-detect best input directory with fallback logic
- **`_find_assembly_file()`** - Locate assembly file across multiple potential locations
- **`_find_mapping_file()`** - Find taxonomic mapping file for specified level

### **📖 Data Loading & Processing**
- **`load_data()`** - Load and merge assembly + mapping data with proper header parsing
- **`classify_genome_source(df)`** - Classify genomes as 'isolate' vs 'uncultured' using pattern matching
- **`extract_species_taxid(df)`** - Extract species-level taxids using taxonkit when not available
- **`create_species_subset(df, species_col)`** - Create species representatives with isolate preference

### **🔍 Taxonomic Information Retrieval**
- **`get_taxids_from_names(names)`** - Convert taxonomic names to taxids via `taxonkit name2taxid`
- **`get_lineages_from_taxids(taxids)`** - Get full lineages via `taxonkit lineage -R -t`

### **📊 Analysis & Counting**
- **`generate_counts(df, species_df)`** - Generate genome/species counts with percentages and lineage data

### **💾 Output Generation**
- **`save_outputs(merged_df, counts_df)`** - Save counts CSV and accessions CSV files

### **🚀 Main Workflow**
- **`run()`** - Orchestrate complete parsing workflow
- **`main()`** - Command-line interface and argument parsing

## Downstream Analysis Goals

### **🔄 Database Comparison Pipeline**
The NCBI parser outputs are designed for comparative analysis with other genomic databases:

#### **Primary Comparisons**
- **NCBI vs GTDB**: Compare taxonomic coverage and classification differences
- **NCBI vs EukProt**: Eukaryotic protein database comparison for 18S rRNA analysis
- **NCBI vs Census Data**: Compare against environmental census/survey data

#### **Merger Scripts Integration**
- **`*_merger.py` scripts**: Use count files for database intersection analysis
- **Species-level comparisons**: Leverage species counts for accurate comparative metrics
- **Taxonomic flow analysis**: Track how organisms are classified across databases

### **📈 Visualization & Analysis Outputs**
- **Alluvial/Sankey plots**: Show taxonomic flow between databases
- **Scatter plots**: Compare coverage and novelty factors across databases
- **Summary statistics**: Generate comparative coverage reports

### **🎯 Key Metrics Generated**
- **Coverage percentages**: How much of each database overlaps
- **Novelty factors**: Ratio of OTU counts to species counts
- **Taxonomic completeness**: Coverage across phylum/family/genus levels
- **Isolate vs uncultured ratios**: Genome source distribution analysis

## Taxonomic Mapping Pipeline Integration

### **📋 Prerequisites**
Before running the parser, ensure taxonomic mapping files exist:
```bash
# Required mapping files (generated by taxonomic_mapping/ scripts)
../taxonomic_mapping/taxid_to_phylum.csv
../taxonomic_mapping/taxid_to_family.csv
../taxonomic_mapping/taxid_to_genus.csv
```

### **🔗 Mapping File Generation**
```bash
# Generate mapping files using taxonomic_mapping scripts
cd ../taxonomic_mapping/
python phylum_taxid_improved.py    # Creates taxid_to_phylum.csv
python family_taxid_improved.py    # Creates taxid_to_family.csv
python genus_taxid.py              # Creates taxid_to_genus.csv
```

### **⚠️ Error Handling & Unmapped Entries**
The parser handles taxonomic mapping failures gracefully:
- **Unmapped taxids**: Logged to `error_log/unmapped_taxids_{level}.csv`
- **Coverage analysis**: ~97.9% of assemblies successfully mapped
- **Error categorization**: Viruses, unclassified taxa, environmental samples filtered
- **Data loss assessment**: True mapping issues identified in `potential_loss/` directory

### **🔍 Quality Control**
- **Consistency checks**: Verify unmapped entries don't appear in other taxonomic levels
- **Species counting validation**: Ensure species counts are reasonable vs genome counts
- **Lineage verification**: Cross-reference taxonomic names with NCBI taxonomy database

## Performance & Statistics

### Processing Time
- **Phylum**: ~30 seconds (295 groups)
- **Family**: ~45 seconds (~8,000 groups)
- **Genus**: ~60 seconds (~35,000 groups)

### Data Coverage
- **Input**: 2,945,415 total assemblies
- **Processed**: 2,882,656 entries (97.9% success rate)
- **Species Representatives**: ~142,000 (isolate-preferential selection)

### Output Sizes
- **Phylum counts**: ~295 rows, 64KB
- **Family counts**: ~8,000 rows, 1.7MB
- **Genus counts**: ~35,000 rows, 7.4MB
- **Accession files**: ~2.9M rows, 200-250MB each