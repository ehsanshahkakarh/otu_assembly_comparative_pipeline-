# 18S-NCBI Merger: Aggregation & Calculations Explained

## Overview

This document explains in detail how the 18S-NCBI merger aggregates data and calculates metrics when merging environmental census data with NCBI genomic data.

## Data Flow

```
18S Census Data → Matching Algorithm → Aggregation → Calculations → Output CSV
```

---

## 1. INPUT DATA

### 18S Census Data (`eukcensus_18S_by_{level}.csv`)
- **Name_to_use**: Taxonomic name from environmental census
- **taxid**: NCBI taxonomy ID (or "NA" if not mapped)
- **otu_count**: Number of environmental OTU clusters (species-level diversity)
- **size_count**: Total sequence abundance (read count)
- **lineage**: Full taxonomic lineage
- **lineage_taxids**: Semicolon-separated taxids for lineage

### NCBI Species Data (`ncbi_{level}_counts.csv`)
- **species_taxid**: NCBI species-level taxonomy ID
- **total_genome_count**: Number of genome assemblies for this species
- **isolate_genome_count**: Number of cultured isolate genomes
- **uncultured_genome_count**: Number of MAG/SAG genomes
- **lineage**: Full taxonomic lineage
- **lineage_taxids**: Semicolon-separated taxids for lineage
- **domain**: Domain classification (Eukaryota, Bacteria, Archaea)

### NCBI Accessions (`ncbi_{level}_with_accessions.csv`)
- **taxon**: Taxonomic name
- **isolate_count**: Number of isolate genomes

---

## 2. MATCHING ALGORITHM

For each census entry, the algorithm finds matching NCBI species using a **hierarchical taxid matching** strategy:

### Priority 1: Direct Taxid Match
```python
census_taxid == NCBI species_taxid
```
Exact match between census taxid and NCBI species taxid.

### Priority 2: Hierarchical Taxid Match
```python
census_taxid in NCBI lineage_taxids
```
Census taxid appears anywhere in the NCBI species' lineage.

**Pattern matching**: `;taxid;` or `^taxid;` or `;taxid$` or `^taxid$`

### Example
Census entry: `Fungi.U.family` with taxid `4751` (Fungi)

Matches NCBI species like:
- *Saccharomyces cerevisiae* with lineage_taxids: `131567;2759;33154;4751;451864;...`
- *Aspergillus fumigatus* with lineage_taxids: `131567;2759;33154;4751;451864;...`
- Any species with `4751` in their lineage

### Deduplication
All matches are combined and deduplicated to avoid counting the same species twice.

---

## 3. AGGREGATION

Once all matching NCBI species are found, the data is aggregated:

### Species Count
```python
total_species = len(matched_species)
```
Count of unique NCBI species that matched.

### Genome Count
```python
total_genomes = Σ total_genome_count (across all matched species)
```
Sum of all genome assemblies from matched species.

### Isolate Count
```python
total_isolates = Σ isolate_genome_count (across all matched species)
```
Sum of cultured isolate genomes from matched species.

### Uncultured Count
```python
total_uncultured = Σ uncultured_genome_count (across all matched species)
```
Sum of MAG/SAG genomes from matched species.

### Domain
```python
domain = most_common(domain from matched species)
```
The most frequently occurring domain among matched species.

---

## 4. CALCULATIONS

### Census Percentages

**OTU Percentage**:
```python
otu_percentage = (census_otu_count / total_census_otus) × 100
```
What percentage of total environmental OTUs does this taxon represent?

**Size Percentage**:
```python
size_percentage = (census_size_count / total_census_size) × 100
```
What percentage of total sequence abundance does this taxon represent?

### NCBI Percentages

**Genome % of Database**:
```python
genome_pct_db = (total_genomes / total_ncbi_genomes) × 100
```
What percentage of total NCBI genomes (domain-filtered) does this taxon represent?

**Species %**:
```python
species_pct = (total_species / total_ncbi_species) × 100
```
What percentage of total NCBI species (domain-filtered) does this taxon represent?

### Isolate Metrics

**Isolate Percentage**:
```python
isolate_percentage = (total_isolates / total_genomes) × 100
```
What percentage of genomes for this taxon are cultured isolates (vs MAGs)?

### Novelty Metrics

**Novelty Factor**:
```python
novelty_factor = census_otu_count / ncbi_species_count
```
- **Higher values** (>1): More environmental diversity than genomic representation
- **Lower values** (<1): Well-represented in genomic databases
- **Infinite**: Taxon found in environment but completely missing from NCBI

**Overrepresentation Factor**:
```python
overrepresentation_factor = ncbi_species_count / census_otu_count
```
- **Higher values** (>1): Database bias toward cultured taxa
- **Lower values** (<1): Underrepresented in genomic databases

---

## 5. REAL EXAMPLE

### Example 1: Well-Represented Taxon

**Input**: `Eukaryota.U.family`
- Census OTU count: 9,333
- Census size count: 39,360
- Census taxid: 2759

**Matching**: Hierarchical taxid match found 24,211 NCBI species

**Aggregation**:
- Total genomes: 59,130 (sum across 24,211 species)
- Total species: 24,211 (count of matched species)
- Total isolates: 57,069

**Calculations**:
- OTU percentage: 13.16% (of total census)
- Genome % of DB: 100.0% (all eukaryotic genomes)
- Novelty factor: 0.385 (well-represented)
- Overrepresentation factor: 2.594 (database bias)

**Interpretation**: This taxon is well-represented in NCBI with 2.6x more sequenced species than environmental OTUs.

### Example 2: Census-Only Taxon

**Input**: `Embryophyceae_XX`
- Census OTU count: 2,248
- Census size count: 16,076
- Census taxid: N/A

**Matching**: No matches found

**Aggregation**:
- Total genomes: 0
- Total species: 0

**Calculations**:
- Novelty factor: ∞ (infinite)
- Match status: census_only

**Interpretation**: This taxon is found in the environment but completely missing from NCBI databases.

---

## Summary

The merger uses **hierarchical taxid matching** to find all NCBI species that belong to each census taxon, then **aggregates** genome counts and species counts, and finally **calculates** comparative metrics to assess database completeness and bias.

