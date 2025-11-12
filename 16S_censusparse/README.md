# EukCensus Data Processing

This directory contains scripts and data for processing EukCensus 16S cluster data, focusing on eukaryotic taxonomic information.

## Overview

The scripts in this directory process the `eukcensus_16S.clusters.97.tsv` file to generate taxonomic summaries at different levels (phylum, family, genus) and filter out non-eukaryotic entries.

## Input Files

- `eukcensus_16S.clusters.97.tsv`: Tab-separated file containing 16S cluster data with columns:
  - `centroid`: Cluster centroid identifier
  - `members`: Array of member identifiers
  - `size`: Cluster size
  - `phylum`: Taxonomic phylum
  - `familiy`: Taxonomic family (note the typo in the column name)
  - `genus`: Taxonomic genus

## Scripts

### 1. parse_eukcensus_clusters.py

This script processes the `eukcensus_16S.clusters.97.tsv` file and generates three CSV files organized by phylum, family, and genus.

**Features:**
- Groups data by taxonomic rank (phylum, family, genus)
- Calculates total member size for each group
- Collects all centroid entries for each group
- Uses taxonkit to get NCBI taxids for each taxon name
- Handles special cases:
  - Removes underscores in taxon names for better taxid matching
  - Ignores organelle information (e.g., ".Mitochondria", ".Chloroplast")
  - Falls back to genus-level matching when full name matching fails

**Output files:**
- `eukcensus_by_phylum.csv`: Grouped by phylum
- `eukcensus_by_family.csv`: Grouped by family
- `eukcensus_by_genus.csv`: Grouped by genus

Each output file contains columns:
- `taxid`: NCBI taxonomy ID
- `taxon_name`: Name of the taxonomic group
- `member_size`: Total number of members in the group
- `centroids`: Semicolon-separated list of centroid identifiers

**Usage:**
```bash
python parse_eukcensus_clusters.py
```

### 2. filter_eukaryotes.py

This script filters the `eukcensus_by_genus.csv` file to remove any non-eukaryotic entries.

**Purpose:**
- Remove bacterial and archaeal entries from the EukCensus data
- Ensure that only true eukaryotic entries remain in the dataset
- Properly handle entries with organelle information (e.g., ".Mitochondria", ".Chloroplast")

**Features:**
- Filters out entries with bacterial or archaeal terms in the taxon name
- Checks if taxids belong to Eukaryota using taxonkit
- Uses batch processing and multiprocessing for high performance
- Handles organelle information in taxon names
- Provides progress bars for each processing step
- Generates detailed statistics on the filtering process

**Key Functions:**
- `extract_organism_name()`: Extracts the organism part from taxon names with organelle information
- `get_taxid_for_name()`: Gets NCBI taxid for a taxon name with fallback to genus-level matching
- `process_taxon_batch()`: Processes a batch of taxon names in parallel
- `check_if_eukaryotic()`: Coordinates the parallel processing of taxon names

**Performance Optimizations:**
- Parallel processing with multiple CPU cores
- Batch processing of taxon names
- Progress tracking with tqdm progress bars
- Reduced verbosity for faster execution

**Output file:**
- `eukcensus_by_genus_eukaryotes_only.csv`: Contains only eukaryotic entries

**Usage:**
```bash
python filter_eukaryotes.py
```

## Taxonomic ID Mapping

The scripts use taxonkit to map taxonomic names to NCBI taxids. The taxonkit database is located at:
```
/clusterfs/jgi/scratch/science/mgs/nelli/ehsan/UNI56v2/00data/refgenomes/gtdb/parse_repaa_table/ncbi_parse_scripts/taxonomic_mapping/taxdump_ncbi
```

## Handling Special Cases

### Underscores in Names
Taxon names often contain underscores (e.g., "Genus_species"). The scripts replace these with spaces before looking up taxids.

### Organelle Information
Some taxon names include organelle information after a period (e.g., "Genus_species.Mitochondria"). The scripts ignore this part when mapping to taxids.

### Genus Fallback
If a full taxon name doesn't match any NCBI taxid, the scripts attempt to match just the genus part.

## Results

The filtering process identified and removed approximately 608 non-eukaryotic entries from the genus-level data, which represents about 13.3% of the original entries.

## Dependencies

- Python 3.6+
- pandas
- taxonkit (with NCBI taxonomy database)
