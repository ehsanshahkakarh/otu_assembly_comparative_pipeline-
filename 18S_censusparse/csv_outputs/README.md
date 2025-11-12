# 18S Census Parse CSV Files

## Directory Contents

This directory contains three CSV files derived from 18S rRNA gene census data:

1. **eukcensus_by_division.csv** (3.3 KB)
   - Division/clade-level taxonomic groupings of eukaryotes
   - 9,042 OTUs for Alveolata, 5,496 for Discoba, etc.

2. **eukcensus_by_family.csv** (49.4 KB)
   - Family-level taxonomic groupings of eukaryotes
   - Includes families like Acanthamoebidae (123 OTUs)

3. **eukcensus_by_genus.csv** (108.3 KB)
   - Genus-level taxonomic groupings of eukaryotes
   - Includes genera like Acanthamoeba (104 OTUs)

## Data Structure

All files share the same column structure:
- `Name_to_use`: Taxonomic name
- `taxid`: NCBI taxonomy ID
- `otu_count`: Number of OTUs (Operational Taxonomic Units) in this taxonomic group
- `lineage`: Full taxonomic lineage (semicolon-separated)
- `lineage_ranks`: Taxonomic ranks for each level
- `lineage_taxids`: NCBI taxonomy IDs for each level

## Purpose

These files provide a structured representation of 18S rRNA gene diversity across eukaryotic taxonomic levels, with NCBI taxonomy integration. They serve as input for:

1. Comparative analysis with eukaryotic genomic databases (NCBI, EukProt)
2. Coverage assessment of environmental vs. genomic diversity
3. Visualization of eukaryotic taxonomic distributions

## Usage

These files are used by scripts in the `OTU_97_eukcensus_merger` directory, particularly:
- `18s_eukprot_merger.py`
- `18s_ncbi_merger.py`
- `create_18s_eukprot_visualizations.py`

## Generation

These files were generated using the script:
- `18S_eukcensus_parser.py` located in the `../py_18S/` directory

The script processes the original EukCensus 18S cluster data, groups it by taxonomic ranks, and enriches it with NCBI taxonomy information using taxonkit.