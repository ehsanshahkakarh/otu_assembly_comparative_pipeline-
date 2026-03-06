# GTDB-NCBI Database Merger

**High-precision genome-level mapping between GTDB R226 and NCBI GenBank**

This pipeline creates bidirectional mapping tables between GTDB taxonomy and NCBI metadata for **730,351 genomes (99.7% of GTDB)**.

---

## 🎯 Quick Start

### Use the Pre-Built Mapping Tables

The easiest way to use this pipeline is to directly access the mapping tables:

```python
import pandas as pd

# Load the simple map (RECOMMENDED)
df = pd.read_csv('accession_only_merger/mapping_tables/gtdb_to_ncbi_simple_map.csv')

# Look up GTDB taxonomy for any NCBI accession
genome = df[df['accession'] == 'GCA_000005845.2']
print(genome[['gtdb_phylum', 'gtdb_genus', 'gtdb_species']])

# Find all genomes in a GTDB phylum
pseudomonadota = df[df['gtdb_phylum'] == 'Pseudomonadota']
print(f"Found {len(pseudomonadota):,} genomes")
```

**📁 Mapping Tables Location:** `accession_only_merger/mapping_tables/`

See detailed documentation: [`accession_only_merger/mapping_tables/README.md`](accession_only_merger/mapping_tables/README.md)

---

## 📊 Pipeline Overview

### Current Production Pipeline: **Accession-Only Merger**

**Location:** `accession_only_merger/`

**What it does:**
- Matches GTDB genomes to NCBI assemblies using accession numbers (GCA/GCF)
- Creates clean mapping tables for bidirectional lookups
- Validates genome presence (not taxonomic nomenclature)

**Results:**
- ✅ **99.7% match rate** (730,351 / 732,475 genomes)
- ✅ **730,351 genomes mapped** with complete metadata
- ✅ **3 ready-to-use mapping tables** (simple, taxonomy, full)

**Key Files:**
- `mapping_tables/` - **⭐ USE THESE** - Clean CSV mapping tables
- `outputs/` - Raw merger outputs (mapped/unmapped CSVs)
- `merge_gtdb_ncbi_accessions.py` - Main merger script
- `create_gtdb_ncbi_map.py` - Mapping table generator
- `COMPARISON_ANALYSIS.md` - Comparison with old pipeline

---

## 📂 Directory Structure

```
database_merger/
├── README.md                          # This file
├── accession_only_merger/             # PRODUCTION PIPELINE ⭐
│   ├── mapping_tables/                # Clean mapping tables (USE THESE)
│   │   ├── gtdb_to_ncbi_simple_map.csv       (135 MB) - RECOMMENDED
│   │   ├── gtdb_to_ncbi_taxonomy_map.csv     (109 MB)
│   │   ├── gtdb_to_ncbi_full_map.csv         (483 MB)
│   │   └── README.md                          # Usage guide
│   ├── outputs/                       # Raw merger outputs
│   │   ├── archaea_mapped.csv         # 17,221 matched
│   │   ├── bacteria_mapped.csv        # 713,130 matched
│   │   ├── archaea_unmapped.csv       # 24 not in NCBI
│   │   ├── bacteria_unmapped.csv      # 2,100 not in NCBI
│   │   └── ncbi_unmapped.csv          # 2.5M NCBI not in GTDB
│   ├── logs/                          # Execution logs
│   ├── archive/                       # Old/intermediate files
│   ├── merge_gtdb_ncbi_accessions.py  # Main merger script
│   ├── create_gtdb_ncbi_map.py        # Mapping table generator
│   ├── check_unmapped_status.py       # NCBI status checker
│   └── COMPARISON_ANALYSIS.md         # Pipeline comparison
│
└── old_pipeline/                      # LEGACY (for reference only)
    ├── triple_anchor_merger/          # Old name-based merger
    │   ├── merge_outputs/             # Old outputs
    │   ├── logs/                      # Old logs
    │   └── *.py                       # Old scripts
    ├── visuals/                       # Old visualizations
    ├── README.md                      # Old documentation
    └── README_AI_assist_merged_tables.md
```

---

## 🔄 Running the Pipeline

### Step 1: Run the Merger (if needed)

```bash
cd accession_only_merger
python merge_gtdb_ncbi_accessions.py
```

This creates the raw outputs in `outputs/`.

### Step 2: Generate Mapping Tables

```bash
python create_gtdb_ncbi_map.py
```

This creates the clean mapping tables in `mapping_tables/`.

### Step 3: (Optional) Check Unmapped Status

```bash
# Check why genomes are unmapped
python check_unmapped_status.py archaea
python check_unmapped_status.py bacteria
```

---

## 📈 Results Summary

### Match Rates (GTDB R226 → NCBI GenBank)

| Domain | Total GTDB | Matched | Match Rate | Unmapped |
|--------|-----------|---------|------------|----------|
| **Archaea** | 17,245 | 17,221 | **99.9%** | 24 (0.1%) |
| **Bacteria** | 715,230 | 713,130 | **99.7%** | 2,100 (0.3%) |
| **TOTAL** | 732,475 | 730,351 | **99.7%** | 2,124 (0.3%) |

### Top 10 Phyla by Genome Count

1. Pseudomonadota: 266,237
2. Bacillota: 205,524
3. Bacteroidota: 88,978
4. Actinomycetota: 59,316
5. Campylobacterota: 12,657
6. Patescibacteriota: 10,828
7. Verrucomicrobiota: 8,309
8. Cyanobacteriota: 7,028
9. Chloroflexota: 6,844
10. Acidobacteriota: 5,430

---

## 💡 Common Use Cases

### 1. Convert NCBI Accession → GTDB Taxonomy
```python
df = pd.read_csv('accession_only_merger/mapping_tables/gtdb_to_ncbi_simple_map.csv')
result = df[df['accession'] == 'GCA_000005845.2']
print(result['gtdb_species'].values[0])
```

### 2. Get All Genomes in a GTDB Lineage
```python
family_genomes = df[df['gtdb_family'] == 'Enterobacteriaceae']
print(f"Found {len(family_genomes):,} genomes")
```

### 3. Map GTDB Species → NCBI TaxIDs
```python
species = df[df['gtdb_species'] == 'Escherichia coli']
print(species[['accession', 'ncbi_taxid', 'ncbi_organism_name']])
```

### 4. Filter High-Quality Genomes
```python
df_full = pd.read_csv('accession_only_merger/mapping_tables/gtdb_to_ncbi_full_map.csv')
complete = df_full[df_full['assembly_level'] == 'Complete Genome']
print(f"Found {len(complete):,} complete genomes")
```

---

## 🆚 Accession-Only vs Triple Anchor

| Feature | Accession-Only (NEW) | Triple Anchor (OLD) |
|---------|---------------------|---------------------|
| **Level** | Genome-level | Taxonomic name-level |
| **Matching** | Accession numbers | Taxonomic names |
| **Match Rate** | 99.7% of genomes | ~40% of phylum names |
| **Question** | "Are genomes in both DBs?" | "Do DBs use same names?" |
| **Use Case** | Genome presence validation | Nomenclature comparison |
| **Status** | ✅ Production | 📦 Archived |

**Conclusion:** For genome-level analysis, use the **accession-only merger**.

See detailed comparison: [`accession_only_merger/COMPARISON_ANALYSIS.md`](accession_only_merger/COMPARISON_ANALYSIS.md)

---

## 📚 Documentation

- **Main Pipeline:** [`accession_only_merger/README.md`](accession_only_merger/README.md) *(to be created)*
- **Mapping Tables:** [`accession_only_merger/mapping_tables/README.md`](accession_only_merger/mapping_tables/README.md)
- **Comparison Analysis:** [`accession_only_merger/COMPARISON_ANALYSIS.md`](accession_only_merger/COMPARISON_ANALYSIS.md)
- **Old Pipeline:** [`old_pipeline/README.md`](old_pipeline/README.md)

---

## ⚠️ Important Notes

1. **Unmapped genomes (0.3%):** These GTDB genomes don't exist in current NCBI GenBank because they were suppressed/removed after GTDB R226 snapshot.

2. **Match types:**
   - `GCA`: Matched via NCBI's `assembly_accession` column (GenBank)
   - `GCF`: Matched via NCBI's `gbrs_paired_asm` column (RefSeq)

3. **Data versions:**
   - GTDB: R226 (April 2025)
   - NCBI: GenBank Assembly Summary (March 2026)

---

**Created:** 2026-03-04  
**Pipeline:** Accession-Only Merger  
**Coverage:** 730,351 genomes (99.7% of GTDB R226)

