# Domain Splitter Analysis Summary

**Date:** 2026-03-02  
**Analysis:** Comparison of old vs new domain splitting results

---

## Executive Summary

✅ **MAJOR IMPROVEMENT CONFIRMED!**  
✅ **Unknown species reduced from 11,895 to 155 (98.7% reduction)**  
✅ **The improvement was due to updated taxonkit database (tarball)**

---

## Key Findings

### 1. Unknown Species Reduction

| Metric | OLD (Feb 27, 2026) | NEW (Mar 2, 2026) | Change |
|--------|-------------------|-------------------|--------|
| **Unknown species** | 11,895 | 155 | **-11,740 (-98.7%)** |
| **Total species** | 162,731 | 162,731 | 0 |
| **Success rate** | 92.7% | 99.9% | **+7.2%** |

### 2. Domain Distribution (NEW - Mar 2, 2026)

| Domain | Species | Total Genomes | Isolate % |
|--------|---------|---------------|-----------|
| **Bacteria** | 99,085 | 2,930,149 | 79.4% |
| **Viruses** | 35,847 | 256,154 | 94.4% |
| **Eukaryota** | 24,211 | 59,130 | 96.5% |
| **Archaea** | 3,433 | 34,539 | 9.3% |
| **Unknown** | 155 | 35,849 | 0.1% |

### 3. What Changed Between Runs?

**Comparison of pipeline logs:**

**Feb 23, 2026 (OLD):**
- Retrieved: 150,991/162,731 lineages (92.8%)
- Missing: 11,740 lineages

**Mar 1, 2026 (NEW):**
- Retrieved: 162,731/162,731 lineages (100.0%)
- Missing: 0 lineages

**Root Cause:** Updated taxonkit database (tarball)
- The code remained the same
- The taxonkit database was updated between Feb 23 and Mar 1
- This update included lineage information for 11,740 additional species_taxid values

### 4. Remaining 155 "Unknown" Species

These are **NOT missing lineages** - they are legitimately unclassified:

**Categories:**
- Metagenomes (wastewater, soil, marine, gut, etc.)
- Environmental samples
- Uncultured prokaryotes
- Unclassified sequences

**Examples:**
- `198431` - uncultured prokaryote (9,026 genomes)
- `527639` - wastewater metagenome (5,451 genomes)
- `408170` - human gut metagenome (3,396 genomes)
- `410658` - soil metagenome (2,783 genomes)

**These entries:**
- Have lineage information in NCBI
- Are classified as "unclassified entries" or "metagenomes"
- Don't belong to Bacteria, Archaea, Eukaryota, or Viruses domains
- Are correctly categorized as "Unknown" by the domain splitter

---

## Conclusion

### What We Did Differently

**Answer:** We updated the taxonkit database (tarball) between Feb 23 and Mar 1, 2026.

**Evidence:**
1. Same codebase (nev_parse_meth → ncbi_parse)
2. Same assembly data (3,315,821 genomes, 162,731 species)
3. Different taxonkit lineage retrieval success rate:
   - OLD: 92.8% success
   - NEW: 100.0% success
4. Difference: 11,740 additional lineages retrieved

### Impact

✅ **11,740 species** that were previously "unknown" now have complete lineage information  
✅ **99.9% of species** now have proper domain classification  
✅ **Only 155 species** remain as "Unknown" (metagenomes and environmental samples)  
✅ **Pipeline is working perfectly** with updated taxonkit database

---

## Files Generated

**Output directory:** `ncbi_parse/domain_meta/output/`

| File | Description | Size |
|------|-------------|------|
| `bacteria_species.csv` | 99,085 bacterial species | 34.31 MB |
| `viruses_species.csv` | 35,847 viral species | 11.21 MB |
| `eukaryota_species.csv` | 24,211 eukaryotic species | 14.70 MB |
| `archaea_species.csv` | 3,433 archaeal species | 1.21 MB |
| `unknown_species.csv` | 155 metagenomes/environmental | 0.03 MB |
| `domain_summary_20260302_192257.csv` | Summary statistics | - |

---

## Recommendations

### Immediate Actions
1. ✅ **Keep taxonkit database updated** to maintain high lineage retrieval success
2. ✅ **Use the current domain_meta output** (Mar 2, 2026) for downstream analysis
3. ✅ **Archive old domain_meta output** (Feb 27, 2026) to avoid confusion

### Future Maintenance
1. Periodically update taxonkit database to capture new species
2. Re-run domain splitter after taxonkit updates
3. Monitor the "Unknown" category for any unexpected increases

---

**Analysis Location:** `00_gaps_taxonomic/parse_repaa_table/ncbi_parse/domain_meta/`

