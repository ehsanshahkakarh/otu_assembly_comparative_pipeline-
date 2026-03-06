# Centralized Metadata Directory

This directory contains all reference metadata files used across the database parsing pipelines.

## Directory Structure

```
metadata/
├── gtdb/              # GTDB taxonomy files
├── ncbi/              # NCBI assembly summary files
├── eukprot/           # EukProt dataset information
├── eukcensus_16S/     # 16S rRNA census data
└── eukcensus_18S/     # 18S rRNA census data
```

## File Inventory

### GTDB (101M)
- `ar53_taxonomy.tsv` - Archaeal taxonomy (GTDB R220)
- `bac120_taxonomy.tsv` - Bacterial taxonomy (GTDB R220)

### NCBI (1.5G)
- `assembly_summary_genbank.txt` - GenBank assembly summary

### EukProt (770K)
- `Eukprot_included_datasets.txt` - EukProt dataset metadata

### EukCensus 16S (129M)
- `eukcensus_16S.clusters.97.tsv` - 16S rRNA OTU clusters (97% similarity)

### EukCensus 18S (56M)
- `eukcensus_18S.clusters.97.tsv` - 18S rRNA OTU clusters (97% similarity)

## Usage

All parsing scripts should reference files from this centralized location:

```python
# Example path construction
from pathlib import Path

METADATA_DIR = Path(__file__).parent.parent / "metadata"
GTDB_ARCHAEA = METADATA_DIR / "gtdb" / "ar53_taxonomy.tsv"
GTDB_BACTERIA = METADATA_DIR / "gtdb" / "bac120_taxonomy.tsv"
NCBI_ASSEMBLY = METADATA_DIR / "ncbi" / "assembly_summary_genbank.txt"
EUKPROT_DATASETS = METADATA_DIR / "eukprot" / "Eukprot_included_datasets.txt"
EUKCENSUS_16S = METADATA_DIR / "eukcensus_16S" / "eukcensus_16S.clusters.97.tsv"
EUKCENSUS_18S = METADATA_DIR / "eukcensus_18S" / "eukcensus_18S.clusters.97.tsv"
```

## Notes

- Original metadata files remain in their respective pipeline directories for backward compatibility
- This centralized location is the **canonical source** for all metadata
- File sizes are approximate and may vary with database updates
- NCBI assembly file renamed from `00assembly_summary_genbank.txt` to `assembly_summary_genbank.txt`

## Last Updated
2026-03-05

