# 18S EukCensus Comprehensive Sanity Check Report

**Generated:** 2025-10-15 11:17:05

---

## Executive Summary

- **Total unique taxonomic combinations:** 491
- **Total clusters:** 70,899
- **Total sequences:** 401,342
- **Average clusters per combination:** 144.4
- **Average sequences per combination:** 817.4

## Classification Completeness Analysis

- **Fully classified combinations:** 370 (75.4%)
- **Genus-only unclassified:** 38 (7.7%)
- **Multiple levels unclassified:** 83 (16.9%)

## Major Taxonomic Groups

### Top 10 Combinations by Cluster Count

| Rank | Taxonomic Combination | Clusters | Sequences |
|------|----------------------|----------|----------|
| 1 | Eukaryota.U.division → Eukaryota.U.family → Eukaryota.U.genus | 9,333 | 39,360 |
| 2 | Opisthokonta → Metazoa.U.family → Metazoa.U.genus | 4,480 | 25,091 |
| 3 | Amoebozoa.U.division → Amoebozoa.U.family → Amoebozoa.U.genus | 4,214 | 15,143 |
| 4 | Opisthokonta → Insecta → Insecta.U.genus | 2,513 | 16,928 |
| 5 | Evosea → Mastigamoebidae → Mastigamoeba | 1,609 | 10,879 |
| 6 | Streptophyta → Embryophyceae_XX → Embryophyceae_XX.U.genus | 1,535 | 7,457 |
| 7 | Discoba → Euglenida.U.family → Euglenida.U.genus | 1,526 | 5,350 |
| 8 | Opisthokonta → Nematoda.U.family → Nematoda.U.genus | 1,163 | 5,960 |
| 9 | Rhizaria → Cercozoa.U.family → Cercozoa.U.genus | 1,152 | 3,275 |
| 10 | Opisthokonta → Fungi.U.family → Fungi.U.genus | 1,112 | 2,638 |

### Top 10 Combinations by Sequence Count

| Rank | Taxonomic Combination | Sequences | Clusters |
|------|----------------------|-----------|----------|
| 1 | Eukaryota.U.division → Eukaryota.U.family → Eukaryota.U.genus | 39,360 | 9,333 |
| 2 | Opisthokonta → Metazoa.U.family → Metazoa.U.genus | 25,091 | 4,480 |
| 3 | Opisthokonta → Insecta → Insecta.U.genus | 16,928 | 2,513 |
| 4 | Amoebozoa.U.division → Amoebozoa.U.family → Amoebozoa.U.genus | 15,143 | 4,214 |
| 5 | Opisthokonta → Exobasidiomycetes → Malassezia | 11,722 | 298 |
| 6 | Evosea → Mastigamoebidae → Mastigamoeba | 10,879 | 1,609 |
| 7 | Alveolata → Alveolata.U.family → Alveolata.U.genus | 8,253 | 393 |
| 8 | Stramenopiles → Stramenopiles.U.family → Stramenopiles.U.genus | 7,804 | 1,038 |
| 9 | Alveolata → Oxytrichidae → Stylonychia | 7,715 | 800 |
| 10 | Streptophyta → Embryophyceae_XX → Embryophyceae_XX.U.genus | 7,457 | 1,535 |

## Division-Level Classification Rates

| Division | Total Clusters | Classification Rate (%) | Avg Size per Cluster |
|----------|----------------|------------------------|----------------------|
| Opisthokonta | 24,345 | 23.7% | 6.52 |
| Eukaryota | 9,333 | 0.0% | 4.22 |
| Alveolata | 9,042 | 55.1% | 6.69 |
| Discoba | 5,496 | 28.3% | 3.84 |
| Rhizaria | 5,126 | 34.6% | 4.70 |
| Amoebozoa | 4,214 | 0.0% | 3.59 |
| Stramenopiles | 3,668 | 31.3% | 6.76 |
| Evosea | 2,431 | 86.2% | 5.54 |
| Streptophyta | 2,252 | 31.8% | 7.15 |
| Chlorophyta | 1,505 | 32.8% | 9.45 |

## Size Distribution Analysis

### Cluster Count Distribution
- **Mean:** 144.4
- **Median:** 25.0
- **Range:** 1 - 9,333
- **Standard Deviation:** 550.7

### Sequence Count Distribution
- **Mean:** 817.4
- **Median:** 87.0
- **Range:** 1 - 39,360
- **Standard Deviation:** 2685.8

### Average Size per Cluster Distribution
- **Mean:** 5.09
- **Median:** 3.32
- **Range:** 1.00 - 74.54
- **Standard Deviation:** 6.44

## Key Findings

### 1. Data Completeness
- The dataset contains 491 unique taxonomic combinations
- 75.4% of combinations are fully classified to genus level
- 16.9% have multiple unclassified levels

### 2. Taxonomic Distribution
- Largest group: Eukaryota.U.division|Eukaryota.U.family|Eukaryota.U.genus with 9,333 clusters
- Opisthokonta dominates with 24,345 clusters
- High diversity in protist groups (Alveolata, Rhizaria, Stramenopiles)

### 3. Size Patterns
- Highly skewed distribution (mean cluster count: 144.4, median: 25.0)
- Large number of singleton or small combinations
- Few very large taxonomic groups dominate sequence counts

### 4. Classification Quality
- Best classified divisions:
  - Picozoa: 100.0% classification rate
  - Hemimastigophora: 100.0% classification rate
  - Nibbleridia: 100.0% classification rate
- Poorly classified divisions:
  - Eukaryota: 0.0% classification rate
  - Amoebozoa: 0.0% classification rate
  - TSAR: 0.0% classification rate

## Recommendations

1. **Focus curation efforts** on the largest unclassified groups:
   - Eukaryota.U.division|Eukaryota.U.family|Eukaryota.U.genus (9,333 clusters)
   - Opisthokonta|Metazoa.U.family|Metazoa.U.genus (4,480 clusters)
   - Amoebozoa.U.division|Amoebozoa.U.family|Amoebozoa.U.genus (4,214 clusters)
   - Opisthokonta|Insecta|Insecta.U.genus (2,513 clusters)
   - Streptophyta|Embryophyceae_XX|Embryophyceae_XX.U.genus (1,535 clusters)

2. **Investigate outliers** with unusual size patterns for potential data quality issues

3. **Prioritize classification** of high-abundance groups in poorly classified divisions

4. **Consider splitting** very large unclassified groups for more granular analysis

## Generated Files

This analysis generated the following files in the `sanity_check/` directory:

- `COMPREHENSIVE_SANITY_CHECK_REPORT.md`
- `avg_size_outliers.csv`
- `cluster_count_outliers.csv`
- `division_size_statistics.csv`
- `eukaryota_u_division_all_entries.csv`
- `eukaryota_u_division_grouped.csv`
- `eukaryota_u_division_summary.txt`
- `eukcensus_18S_division_counts.csv`
- `eukcensus_18S_family_counts.csv`
- `eukcensus_18S_genus_counts.csv`
- `fully_classified_combinations.csv`
- `genus_unclassified_combinations.csv`
- `major_division_classification_rates.csv`
- `multiple_unclassified_combinations.csv`
- `sequence_count_outliers.csv`
- `size_distribution_summary.txt`
- `taxonomic_combinations_detailed.csv`
- `taxonomic_combinations_summary.txt`
- `unclassified_division_lineages.csv`
- `unclassified_family_lineages.csv`
- `unclassified_genus_lineages.csv`
- `unclassified_patterns_summary.txt`
- `unclassified_summary_stats.csv`

---
*Report generated on 2025-10-15 11:17:05*
