# Merged Tables

This directory contains scripts for merging and analyzing taxonomic data from multiple sources, including GTDB (Genome Taxonomy Database), NCBI (National Center for Biotechnology Information), and EukProt.

## Directory Structure

```
gtdb/parse_repaa_table/merged_tables/
├── csv_output/                        # Output directory for merged taxonomic data
│   ├── merged_phylum_taxonomy.csv     # Merged phylum-level data
│   ├── merged_family_taxonomy.csv     # Merged family-level data
│   ├── merged_genus_taxonomy.csv      # Merged genus-level data
│   ├── phylum_mismatches.csv          # Records of phylum classification mismatches
│   ├── family_mismatches.csv          # Records of family classification mismatches
│   └── genus_mismatches.csv           # Records of genus classification mismatches
│
├── domain_taxonomic_tables/           # Domain-specific taxonomic tables
│   ├── bacteria_phylum.csv            # Bacteria phylum-level data
│   ├── bacteria_family.csv            # Bacteria family-level data
│   ├── bacteria_genus.csv             # Bacteria genus-level data
│   ├── archaea_phylum.csv             # Archaea phylum-level data
│   ├── archaea_family.csv             # Archaea family-level data
│   ├── archaea_genus.csv              # Archaea genus-level data
│   ├── eukaryota_phylum.csv           # Eukaryota phylum-level data
│   ├── eukaryota_family.csv           # Eukaryota family-level data
│   ├── eukaryota_genus.csv            # Eukaryota genus-level data
│   └── mismatches/                    # Domain-specific taxonomic mismatches
│       ├── bacteria_family_mismatches.csv  # Bacteria family mismatches
│       ├── archaea_family_mismatches.csv   # Archaea family mismatches
│       └── ...                        # Other domain-specific mismatch files
│
├── eukaryota_taxonomic_tables/        # Eukaryota-specific taxonomic tables
│   ├── eukaryota_phylum.csv           # Eukaryota phylum-level data
│   ├── eukaryota_family.csv           # Eukaryota family-level data
│   └── eukaryota_genus.csv            # Eukaryota genus-level data
│
├── merged_taxonomic_tables/           # Merged taxonomic tables across domains
│   ├── merged_phylum.csv              # Merged phylum data across all domains
│   ├── merged_family.csv              # Merged family data across all domains
│   └── merged_genus.csv               # Merged genus data across all domains
│
├── older/                             # Older versions of scripts and data
│   ├── py_scripts/                    # Python scripts from earlier iterations
│   │   └── mtdata_library/            # Library scripts for metadata processing
│   │       ├── ncbi_accesion_with_mapping.py  # NCBI accession mapping script
│   │       └── merged_accession_phylum_cleaned.py  # Cleaned merged accession script
│   │
│   └── raw_total_test/                # Raw test data and scripts
│       ├── initial_raw_merg.py        # Initial raw merge script
│       ├── parse_by_domain_and_rank.py  # Domain and rank parsing script
│       └── domain_taxonomic_tables/   # Domain-specific taxonomic tables from testing
│
├── merge_taxonomic_ranks.py           # Script to merge GTDB and NCBI taxonomic data
├── parse_by_domain.py                 # Script to parse data by domain
├── merge_domain_tables.py             # Script to merge domain-specific tables
└── create_merged_domain_tables.py     # Script to create merged domain tables for Eukaryota
```

## File Descriptions

### Main Scripts

1. **`merge_taxonomic_ranks.py`**:
   - Merges taxonomic data from GTDB and NCBI sources
   - Processes data at phylum, family, and genus levels
   - Identifies and records taxonomic mismatches between sources
   - Outputs merged taxonomic data to CSV files in `csv_output/`
   - Input paths:
     - NCBI: `/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/csv_ncbi`
     - GTDB: `/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/gtdb_parse/csv_gtdb`

2. **`parse_by_domain.py`**:
   - Separates merged taxonomic data by domain (Bacteria, Archaea, Eukaryota)
   - Creates domain-specific taxonomic tables for each rank
   - Outputs domain-specific data to CSV files in `domain_taxonomic_tables/`
   - Processes data from `csv_output/` directory

3. **`merge_domain_tables.py`**:
   - Combines domain-specific taxonomic tables into unified tables
   - Creates comprehensive taxonomic tables across all domains
   - Outputs merged data to CSV files in `merged_taxonomic_tables/`
   - Processes data from `domain_taxonomic_tables/` directory

4. **`create_merged_domain_tables.py`**:
   - Creates Eukaryota-specific taxonomic tables comparing NCBI and EukProt data
   - Merges Eukaryota data from different sources
   - Outputs Eukaryota-specific data to CSV files in `eukaryota_taxonomic_tables/`
   - Input paths:
     - NCBI: `/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/csv_ncbi`
     - EukProt: `/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/eukprot_parse/csv_eukprot`

### Output Files

#### CSV Output Files

1. **Merged Taxonomic Data**:
   - `merged_phylum_taxonomy.csv`, `merged_family_taxonomy.csv`, `merged_genus_taxonomy.csv`
   - Columns: taxon_name (e.g., phylum_gtdb, phylum_ncbi), domain, genome_count, accessions, data_source, taxonomy_match

2. **Taxonomic Mismatches**:
   - `phylum_mismatches.csv`, `family_mismatches.csv`, `genus_mismatches.csv`
   - Columns: {rank}_gtdb, domain_gtdb, gtdb_genome_count, accession_clean, {rank}_ncbi, domain_ncbi, ncbi_genome_count, data_source, taxonomy_match

#### Domain-Specific Tables

1. **Domain Taxonomic Tables**:
   - `bacteria_phylum.csv`, `archaea_phylum.csv`, `eukaryota_phylum.csv`, etc.
   - Columns: taxon_name, domain, genome_count, accessions, source, data_source

2. **Domain-Specific Mismatches**:
   - `bacteria_family_mismatches.csv`, `archaea_family_mismatches.csv`, etc.
   - Columns: {rank}_gtdb, domain_gtdb, gtdb_genome_count, accession_clean, {rank}_ncbi, domain_ncbi, ncbi_genome_count, data_source, taxonomy_match, domain

#### Merged Taxonomic Tables

1. **Merged Domain Tables**:
   - `merged_phylum.csv`, `merged_family.csv`, `merged_genus.csv`
   - Columns: taxon_name, domain, genome_count, accessions, source, data_source

### Older Scripts and Data

1. **`older/py_scripts/mtdata_library/ncbi_accesion_with_mapping.py`**:
   - Processes NCBI assembly summary data
   - Maps GCA genomes to their paired GCF genomes
   - Creates clean accession-to-phylum mapping table

2. **`older/py_scripts/mtdata_library/merged_accession_phylum_cleaned.py`**:
   - Cleans and processes merged accession-phylum data
   - Filters for prokaryotic data
   - Handles version numbers in accessions

3. **`older/raw_total_test/initial_raw_merg.py`**:
   - Early version of the merge script
   - Merges GTDB and NCBI metadata
   - Creates phylum-level summaries with accession lists

4. **`older/raw_total_test/parse_by_domain_and_rank.py`**:
   - Early version of domain parsing script
   - Creates domain-specific taxonomic tables

## Data Flow

1. **Input Data Sources**:
   - NCBI assembly summary and taxonomy data
   - GTDB taxonomy data
   - EukProt taxonomy data (for Eukaryota)

2. **Processing Pipeline**:
   - Merge taxonomic data from different sources (`merge_taxonomic_ranks.py`)
   - Separate data by domain (`parse_by_domain.py`)
   - Create domain-specific tables for Eukaryota (`create_merged_domain_tables.py`)
   - Merge domain-specific tables into unified tables (`merge_domain_tables.py`)

3. **Output Data**:
   - Merged taxonomic data across sources
   - Domain-specific taxonomic tables
   - Taxonomic mismatches between sources
   - Comprehensive taxonomic tables across all domains

## Usage

The scripts in this directory are designed to be run in sequence:

1. First, run `merge_taxonomic_ranks.py` to merge GTDB and NCBI data
2. Then, run `parse_by_domain.py` to separate data by domain
3. For Eukaryota-specific analysis, run `create_merged_domain_tables.py`
4. Finally, run `merge_domain_tables.py` to create unified taxonomic tables

Each script reads from the output of the previous script, creating a pipeline for taxonomic data processing and analysis.