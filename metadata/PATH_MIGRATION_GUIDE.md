# Metadata Path Migration Guide

## Overview
All metadata files have been centralized to `00parse_database/metadata/` directory.

## Path Changes

### GTDB Files
**Old paths:**
- `gtdb_parse/metadata/ar53_taxonomy.tsv`
- `gtdb_parse/metadata/00ar53_taxonomy.tsv`
- `gtdb_parse/metadata/bac120_taxonomy.tsv`
- `gtdb_parse/metadata/00bac120_taxonomy.tsv`

**New paths:**
- `metadata/gtdb/ar53_taxonomy.tsv`
- `metadata/gtdb/bac120_taxonomy.tsv`

### NCBI Files
**Old paths:**
- `old_pipeline/ncbi_parse_old/metadata/00assembly_summary_genbank.txt`
- `ncbi_parse/00assembly_summary_genbank.txt`

**New path:**
- `metadata/ncbi/assembly_summary_genbank.txt`

### EukProt Files
**Old path:**
- `eukprot_parse/metadata/Eukprot_included_datasets.txt`

**New path:**
- `metadata/eukprot/Eukprot_included_datasets.txt`

### EukCensus 16S Files
**Old path:**
- `16S_censusparse/metadata/eukcensus_16S.clusters.97.tsv`

**New path:**
- `metadata/eukcensus_16S/eukcensus_16S.clusters.97.tsv`

### EukCensus 18S Files
**Old path:**
- `18S_censusparse/metadata/eukcensus_18S.clusters.97.tsv`

**New path:**
- `metadata/eukcensus_18S/eukcensus_18S.clusters.97.tsv`

## Scripts Requiring Updates

### GTDB Parse Scripts
- `gtdb_parse/py_gtdb/family_gtdb_parse.py`
- `gtdb_parse/py_gtdb/genus_gtdb_parse.py`
- `gtdb_parse/py_gtdb/phylum_gtdb_parse.py`
- `gtdb_parse/py_gtdb/*_species_subset.py`
- `gtdb_parse/py_gtdb/gtdb_parser/config.py`

### EukProt Parse Scripts
- `eukprot_parse/py_eukprot/improv_eukprot_lineage.py`

### 16S Census Parse Scripts
- `16S_censusparse/py_16S/16S_eukcensus_parser.py`
- `16S_censusparse/py_16S/census_parser/config.py`

### 18S Census Parse Scripts
- `18S_censusparse/src/config.py`
- `18S_censusparse/src/pipeline_taxonkit.py`
- `18S_censusparse/src/division_context_adder.py`

### Database Merger Scripts
- `database_merger/gtdb_ncbi_mapper/config.yaml` ✅ UPDATED

## Migration Status
- [x] Metadata directory created
- [x] Files copied to centralized location
- [x] gtdb_ncbi_mapper config updated
- [ ] GTDB parser scripts updated
- [ ] EukProt parser scripts updated
- [ ] 16S census parser scripts updated
- [ ] 18S census parser scripts updated
- [ ] NCBI parser scripts updated

## Notes
- Original files remain in place for backward compatibility
- Scripts should be updated to use centralized paths
- Use relative paths from script location to metadata directory

