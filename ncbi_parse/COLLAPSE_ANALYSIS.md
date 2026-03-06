# Genome to Species Collapse Analysis

**Date:** 2026-03-02  
**Analysis:** How 3.3 million genomes collapse to 162,731 species

---

## Executive Summary

✅ **3,315,821 genomes** collapse to **162,731 unique species**  
✅ **Collapse ratio: 20.4:1** (average 20.4 genomes per species)  
✅ **Highly skewed distribution:** 80.8% of species have only 1 genome  
✅ **Top 20 species account for 61.5% of all genomes**

---

## Overall Statistics

| Metric | Value |
|--------|-------|
| **Total genomes** | 3,315,821 |
| **Unique species (species_taxid)** | 162,731 |
| **Average genomes per species** | 20.4 |
| **Collapse ratio** | 20.4:1 |

---

## Distribution of Genomes per Species

| Genome Count | Species | Percentage |
|--------------|---------|------------|
| **1 genome** | 131,480 | 80.8% |
| **2-4 genomes** | 20,312 | 12.5% |
| **5-9 genomes** | 4,248 | 2.6% |
| **10-49 genomes** | 4,374 | 2.7% |
| **50-99 genomes** | 929 | 0.6% |
| **100-499 genomes** | 1,013 | 0.6% |
| **500-999 genomes** | 180 | 0.1% |
| **1K-4.9K genomes** | 140 | 0.1% |
| **5K-9.9K genomes** | 27 | 0.0% |
| **10K+ genomes** | 28 | 0.0% |

**Key Insight:** The distribution is extremely skewed - most species have very few genomes, but a small number of species have thousands of genomes.

---

## Top 20 Species by Genome Count

| Rank | Species Name | Genomes | % of Total |
|------|--------------|---------|------------|
| 1 | Salmonella enterica | 608,480 | 18.35% |
| 2 | Escherichia coli | 386,492 | 11.66% |
| 3 | Alphainfluenzavirus influenzae | 137,332 | 4.14% |
| 4 | Staphylococcus aureus | 121,568 | 3.67% |
| 5 | Klebsiella pneumoniae | 112,984 | 3.41% |
| 6 | Campylobacter jejuni | 106,942 | 3.23% |
| 7 | Streptococcus pneumoniae | 80,526 | 2.43% |
| 8 | Listeria monocytogenes | 72,151 | 2.18% |
| 9 | Pseudomonas aeruginosa | 50,457 | 1.52% |
| 10 | Neisseria gonorrhoeae | 49,834 | 1.50% |
| 11 | Streptococcus pyogenes | 49,805 | 1.50% |
| 12 | Acinetobacter baumannii | 48,074 | 1.45% |
| 13 | Campylobacter coli | 40,962 | 1.24% |
| 14 | Enterococcus faecium | 36,601 | 1.10% |
| 15 | Shigella sonnei | 26,665 | 0.80% |
| 16 | Streptococcus agalactiae | 25,210 | 0.76% |
| 17 | Shigella flexneri | 24,580 | 0.74% |
| 18 | Clostridioides difficile | 23,466 | 0.71% |
| 19 | uncultured Lachnospiraceae bacterium | 19,251 | 0.58% |
| 20 | Betainfluenzavirus influenzae | 18,491 | 0.56% |

**Total for Top 20:** 2,039,871 genomes (61.5% of all genomes)

---

## Concentration Analysis

### The "Long Tail" Distribution

**High-genome species (100+ genomes):**
- **Count:** 1,388 species (0.9% of all species)
- **Genomes:** 2,948,456 (88.9% of all genomes)
- **Interpretation:** Less than 1% of species account for nearly 90% of all genomes

**Single-genome species:**
- **Count:** 131,480 species (80.8% of all species)
- **Genomes:** 131,480 (4.0% of all genomes)
- **Interpretation:** Most species have only 1 genome sequenced

### Why This Distribution?

**High-genome species tend to be:**
1. **Clinically important pathogens** (e.g., Salmonella, E. coli, Staph aureus)
2. **Model organisms** (e.g., E. coli)
3. **Foodborne pathogens** (e.g., Listeria, Campylobacter)
4. **Antibiotic resistance concerns** (e.g., Acinetobacter, Klebsiella)
5. **Public health priorities** (e.g., Streptococcus, Neisseria)

**Single-genome species tend to be:**
1. **Rare or newly discovered species**
2. **Environmental isolates**
3. **Non-pathogenic organisms**
4. **Species with limited research interest**
5. **Recently described species**

---

## Implications for Analysis

### 1. Data Representation
- **Most species** (80.8%) have minimal genomic representation (1 genome)
- **Most genomes** (88.9%) come from a small set of well-studied species (1,388 species)

### 2. Statistical Considerations
- Simple species counts may not reflect biological importance
- Genome counts reflect research priorities and clinical relevance
- Need to consider both species diversity AND genome abundance

### 3. Downstream Analysis Recommendations

**For diversity studies:**
- Focus on species counts (162,731 unique species)
- Each species weighted equally

**For genome-based studies:**
- Focus on genome counts (3,315,821 genomes)
- Consider that some species dominate the dataset

**For balanced analysis:**
- Use species-level statistics (this pipeline's output)
- Consider both total_genome_count and species_genome_percentage
- Filter by isolate_genome_percentage for quality

---

## Visual Summary

```
3,315,821 genomes
       ↓
   Collapse by species_taxid
       ↓
162,731 unique species

Distribution:
  131,480 species (80.8%) → 1 genome each
   20,312 species (12.5%) → 2-4 genomes each
    4,248 species (2.6%)  → 5-9 genomes each
    4,374 species (2.7%)  → 10-49 genomes each
    1,388 species (0.9%)  → 100+ genomes each
    
Top 20 species → 2,039,871 genomes (61.5%)
Top 1,388 species → 2,948,456 genomes (88.9%)
```

---

## Conclusion

The collapse from **3.3 million genomes to 162,731 species** reveals:

1. ✅ **Massive redundancy** at the genome level (20.4:1 ratio)
2. ✅ **Highly skewed distribution** (Pareto principle applies)
3. ✅ **Research bias** toward clinically important species
4. ✅ **Long tail** of rare/understudied species

This species-level grouping provides:
- ✅ **Reduced redundancy** for downstream analysis
- ✅ **Species-level statistics** for better understanding
- ✅ **Balanced representation** of genomic diversity
- ✅ **Efficient data structure** for large-scale analysis

---

**Analysis Location:** `00_gaps_taxonomic/parse_repaa_table/ncbi_parse/`

