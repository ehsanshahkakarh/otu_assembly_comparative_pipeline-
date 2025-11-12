# GTDB Taxonomy Parsing Workflow

This document describes the detailed workflow for parsing GTDB taxonomy files, including the comprehensive fallback system for taxid and lineage resolution.

## Overview

The GTDB parsing workflow processes raw taxonomy strings into structured taxonomic data with species counts, taxids, and full lineages using a dual-database approach (GTDB + NCBI fallback).

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           GTDB PARSING WORKFLOW                            │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────┐    ┌─────────────────┐
│ 00bac120_       │    │ 00ar53_         │
│ taxonomy.tsv    │    │ taxonomy.tsv    │
│ (Bacteria)      │    │ (Archaea)       │
└─────────┬───────┘    └─────────┬───────┘
          │                      │
          └──────────┬───────────┘
                     │
                     ▼
          ┌─────────────────────┐
          │   LOAD & COMBINE    │
          │   Taxonomy Data     │
          └─────────┬───────────┘
                    │
                    ▼
          ┌─────────────────────┐
          │  EXTRACT TAXONOMIC  │
          │      RANKS          │
          │ (domain→species)    │
          └─────────┬───────────┘
                    │
                    ▼
          ┌─────────────────────┐
          │   CLEAN & DEDUPE    │
          │ • Remove _A, _B     │
          │ • One per species   │
          └─────────┬───────────┘
                    │
                    ▼
          ┌─────────────────────┐
          │  COUNT BY TARGET    │
          │   TAXONOMIC LEVEL   │
          │ (phylum/family/     │
          │     genus)          │
          └─────────┬───────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        TAXID RESOLUTION SYSTEM                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────┐     SUCCESS     ┌─────────────────┐                   │
│  │   GTDB TAXDUMP  │ ──────────────► │  STORE TAXIDS   │                   │
│  │   name2taxid    │                 │   (Primary)     │                   │
│  │                 │                 └─────────────────┘                   │
│  └─────────┬───────┘                                                       │
│            │                                                               │
│            │ FAILED                                                        │
│            ▼                                                               │
│  ┌─────────────────┐     SUCCESS     ┌─────────────────┐                   │
│  │   NCBI FALLBACK │ ──────────────► │  STORE TAXIDS   │                   │
│  │   name2taxid    │                 │   (Fallback)    │                   │
│  └─────────────────┘                 └─────────────────┘                   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       LINEAGE RESOLUTION SYSTEM                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────┐                 ┌─────────────────┐                   │
│  │  SEPARATE BY    │                 │  SEPARATE BY    │                   │
│  │     SOURCE      │                 │     SOURCE      │                   │
│  └─────────┬───────┘                 └─────────┬───────┘                   │
│            │                                   │                           │
│            ▼                                   ▼                           │
│  ┌─────────────────┐                 ┌─────────────────┐                   │
│  │ GTDB TAXIDS     │                 │ NCBI TAXIDS     │                   │
│  │ ↓               │                 │ ↓               │                   │
│  │ GTDB lineage    │                 │ NCBI lineage    │                   │
│  │ -R -t           │                 │ -R -t           │                   │
│  └─────────┬───────┘                 └─────────┬───────┘                   │
│            │                                   │                           │
│            └─────────────┬─────────────────────┘                           │
│                          │                                                 │
│                          ▼                                                 │
│                ┌─────────────────┐                                         │
│                │ COMBINE LINEAGE │                                         │
│                │      DATA       │                                         │
│                └─────────────────┘                                         │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                    │
                    ▼
          ┌─────────────────────┐
          │  CLEANUP LINEAGES   │
          │ • Truncate to level │
          │ • Format output     │
          └─────────┬───────────┘
                    │
                    ▼
          ┌─────────────────────┐
          │   GENERATE OUTPUT   │
          │ • Summary counts    │
          │ • Detailed accessions│
          └─────────────────────┘
```

## Detailed Process Steps

### 1. Data Loading Phase

```
Input Files:
├── 00bac120_taxonomy.tsv  (Bacterial genomes)
├── 00ar53_taxonomy.tsv    (Archaeal genomes)
└── taxid.map              (Optional: accession→taxid mapping)

Processing:
┌─────────────────┐
│ Load in chunks  │ → Memory efficient processing
│ Combine datasets│ → Single unified DataFrame
│ Extract ranks   │ → Parse taxonomy strings
└─────────────────┘

Output: DataFrame with columns:
- accession_clean
- domain, phylum, class, order, family, genus, species
- taxid (if available from taxid.map)
```

### 2. Data Cleaning Phase

```
┌─────────────────────────────────────────────────────────────────┐
│                    CLEANING OPERATIONS                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ 1. Name Standardization:                                        │
│    Bacillota_A → Bacillota                                      │
│    Pseudomonas_B → Pseudomonas                                  │
│                                                                 │
│ 2. Species Deduplication:                                       │
│    Multiple genomes per species → Keep first (best quality)     │
│                                                                 │
│ 3. Accession Cleaning:                                          │
│    GB_GCA_000005825.2 → GCA_000005825.2                        │
│    RS_GCF_000001405.39 → GCF_000001405.39                      │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 3. Taxonomic Counting

```
Target Level: FAMILY (example)

┌─────────────────┐
│ Group by:       │
│ • family_clean  │ → Standardized family names
│ • domain        │ → Bacteria/Archaea
└─────────┬───────┘
          │
          ▼
┌─────────────────┐
│ Count unique    │
│ species per     │ → family_species_count
│ group           │
└─────────────────┘

Result: Aggregated counts by family+domain
```

### 4. Taxid Resolution System

#### Phase 4A: Primary GTDB Lookup

```
┌─────────────────────────────────────────────────────────────────┐
│                    GTDB TAXDUMP LOOKUP                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ Input: Unique family names                                      │
│ ┌─────────────────┐                                             │
│ │ Enterobacteriaceae │ ──┐                                      │
│ │ Pseudomonadaceae   │   │                                      │
│ │ SM23-39           │   │ Create temp file                     │
│ │ ...               │   │                                      │
│ └─────────────────────┘ ──┘                                      │
│                          │                                      │
│                          ▼                                      │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ TAXONKIT_DB=gtdb-taxdump-R226                               │ │
│ │ taxonkit name2taxid temp_file.txt                           │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                          │                                      │
│                          ▼                                      │
│ Results:                                                        │
│ Enterobacteriaceae → 543                                        │
│ Pseudomonadaceae → 135621                                       │
│ SM23-39 → 589911515                                             │
│                                                                 │
│ Success Rate: ~60-80%                                           │
└─────────────────────────────────────────────────────────────────┘
```

#### Phase 4B: NCBI Fallback

```
┌─────────────────────────────────────────────────────────────────┐
│                     NCBI FALLBACK LOOKUP                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ Input: Failed GTDB names                                        │
│ ┌─────────────────┐                                             │
│ │ UBA1234         │ ──┐                                         │
│ │ Novel_Family_X  │   │                                         │
│ │ ...             │   │ Create temp file                        │
│ └─────────────────────┘ ──┘                                      │
│                          │                                      │
│                          ▼                                      │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ taxonkit name2taxid temp_file.txt                           │ │
│ │ (uses default NCBI taxdump)                                 │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                          │                                      │
│                          ▼                                      │
│ Additional Success: +10-20%                                     │
│ Combined Success Rate: 70-90%                                   │
└─────────────────────────────────────────────────────────────────┘
```

### 5. Lineage Resolution System

#### Phase 5A: Source Separation

```
┌─────────────────────────────────────────────────────────────────┐
│                    TAXID SOURCE TRACKING                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ For each taxid, determine source:                               │
│                                                                 │
│ IF family_name in gtdb_results:                                 │
│     taxid_source = "GTDB"                                       │
│ ELSE:                                                           │
│     taxid_source = "NCBI"                                       │
│                                                                 │
│ ┌─────────────────┐    ┌─────────────────┐                     │
│ │ GTDB TAXIDS     │    │ NCBI TAXIDS     │                     │
│ │ 589911515       │    │ 543             │                     │
│ │ 1234567890      │    │ 135621          │                     │
│ │ ...             │    │ ...             │                     │
│ └─────────────────┘    └─────────────────┘                     │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

#### Phase 5B: Parallel Lineage Lookup

```
┌─────────────────────────────────────────────────────────────────┐
│                    GTDB LINEAGE LOOKUP                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ TAXONKIT_DB=gtdb-taxdump-R226                               │ │
│ │ taxonkit lineage -R -t gtdb_taxids.txt                      │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                                                                 │
│ Result:                                                         │
│ 589911515 → Bacteria;Chlamydiota;Chlamydiia;...;SM23-39        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                    NCBI LINEAGE LOOKUP                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ taxonkit lineage -R -t ncbi_taxids.txt                      │ │
│ │ (uses default NCBI taxdump)                                 │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                                                                 │
│ Result:                                                         │
│ 543 → Bacteria;Proteobacteria;Gammaproteobacteria;...          │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 6. Lineage Processing

```
┌─────────────────────────────────────────────────────────────────┐
│                    LINEAGE CLEANUP                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ Input: Full lineage                                             │
│ Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;   │
│ Enterobacteriaceae;Escherichia;Escherichia_coli                 │
│                                                                 │
│ Target Level: FAMILY                                            │
│                                                                 │
│ Process:                                                        │
│ 1. Split by semicolon                                           │
│ 2. Find "family" rank position                                  │
│ 3. Truncate after family level                                  │
│                                                                 │
│ Output: Cleaned lineage                                         │
│ Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;   │
│ Enterobacteriaceae                                              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 7. Output Generation

```
┌─────────────────────────────────────────────────────────────────┐
│                      OUTPUT FILES                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ 1. SUMMARY COUNTS (gtdb_family_counts.csv)                     │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ family,domain,family_species_count,taxid,lineage,...        │ │
│ │ Enterobacteriaceae,Bacteria,150,543,Bacteria;Proteo...      │ │
│ │ SM23-39,Bacteria,30,589911515,Bacteria;Chlamydiota...       │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                                                                 │
│ 2. DETAILED ACCESSIONS (gtdb_family_with_accessions.csv)       │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ accession_clean,family,domain,taxid                         │ │
│ │ GCA_000005825.2,Enterobacteriaceae,Bacteria,543            │ │
│ │ GCA_000006745.1,Enterobacteriaceae,Bacteria,543            │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Success Metrics

### Typical Success Rates
- **GTDB taxid lookup**: 60-80%
- **NCBI fallback**: +10-20% additional
- **Combined taxid success**: 70-90%
- **Lineage retrieval**: 95-99% of found taxids

### Quality Assurance
- Taxonomic hierarchy validation
- Missing value detection
- Performance monitoring
- Comprehensive logging

## Error Handling

### Common Failure Points
1. **Name not found in either database**
2. **Taxid exists but no lineage available**
3. **Database connection issues**
4. **Malformed taxonomy strings**

### Logging and Recovery
- Failed matches logged to `*_final_failed_taxids.log`
- Detailed execution logs in `logs/` directory
- Graceful degradation (partial results still useful)
- Performance metrics for optimization
