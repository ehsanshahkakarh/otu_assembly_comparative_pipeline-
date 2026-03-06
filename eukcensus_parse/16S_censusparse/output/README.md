# 16S Census Parse CSV Files

## Directory Contents

This directory contains three CSV files derived from 16S rRNA gene census data:

1. **eukcensus16S_by_division.csv** (7.3 KB)
   - Division/phylum-level taxonomic groupings
   - 53,653 OTUs for Pseudomonadota, 49,858 for Bacillota, etc.

2. **eukcensus16S_by_family.csv** (151 KB)
   - Family-level taxonomic groupings
   - Lachnospiraceae (12,280 OTUs), Oscillospiraceae (5,753 OTUs), etc.

3. **eukcensus16S_by_genus.csv** (1.1 MB)
   - Genus-level taxonomic groupings
   - Bacillus (3,741 OTUs), Pseudomonas (2,357 OTUs), etc.

## Data Structure

All files share the same column structure:
- `Name_to_use`: Taxonomic name
- `taxid`: NCBI taxonomy ID
- `otu_count`: Number of OTUs in this taxonomic group
- `lineage`: Full taxonomic lineage (semicolon-separated)
- `lineage_ranks`: Taxonomic ranks for each level
- `lineage_taxids`: NCBI taxonomy IDs for each level

## Purpose

These files provide a structured representation of 16S rRNA gene diversity across taxonomic levels, with NCBI taxonomy integration. They serve as input for:

1. Comparative analysis with genomic databases (NCBI, GTDB)
2. Coverage assessment of environmental vs. genomic diversity
3. Visualization of taxonomic distributions

## Usage

These files are used by scripts in the `OTU_97_eukcensus_merger` directory, particularly:
- `16s_ncbi_merger.py`
- `16s_gtdb_merger.py`
- `create_16s_gtdb_visualizations.py`