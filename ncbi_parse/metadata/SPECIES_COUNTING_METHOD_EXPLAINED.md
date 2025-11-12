# NCBI Species Counting Method - Complete Explanation

**Created:** 2025-08-18  
**Purpose:** Document and verify our approach to counting species in NCBI data

---

## 🎯 **Our Method Overview**

We use a **scientifically accurate, NCBI-standard approach** to count species that avoids double-counting and provides true biological diversity metrics.

### **Core Principle**
- **One species_taxid = One species**
- Multiple genome assemblies can belong to the same species (different strains, isolates, lab cultures)
- We count **unique species**, not total assemblies

---

## 📊 **Data Source & Structure**

### **Primary Data File**
- **File:** `00assembly_summary_genbank.txt` (1.3GB, 2.9M+ records)
- **Source:** NCBI GenBank assembly metadata
- **Updated:** Regularly by NCBI

### **Key Columns Used**
```
assembly_accession  → Unique genome assembly ID
taxid              → Taxonomic ID (strain/isolate level)
species_taxid      → Species-level taxonomic ID (OUR COUNTING UNIT)
organism_name      → Scientific name
```

### **Why species_taxid is Perfect**
- ✅ **NCBI Official:** Assigned by NCBI taxonomists
- ✅ **Species-level:** Represents true biological species
- ✅ **Consistent:** Same species_taxid across all strains/isolates
- ✅ **Complete:** 100% coverage in our dataset (2,945,415 records)

---

## 🔬 **Counting Methodology**

### **1. Total Species Count**
```python
total_species = df['species_taxid'].nunique()
# Result: 149,270 unique species
```

### **2. Taxonomic Group Species Count**
```python
# Example: Species per phylum
phylum_species = df.groupby('phylum')['species_taxid'].nunique()

# Example: Species per family within domain
family_species = df.groupby(['domain', 'family'])['species_taxid'].nunique()
```

### **3. Dual Metrics Approach**
We calculate **both** metrics for complete analysis:

**A. Genome Counts** (Assembly-level)
- Total number of genome assemblies per taxonomic group
- Shows research intensity and data availability

**B. Species Counts** (Species-level)  
- Number of unique species per taxonomic group
- Shows true biological diversity

---

## 📈 **Verification Results**

### **Data Quality (100% Perfect)**
- ✅ **Total records:** 2,945,415 assemblies
- ✅ **Species coverage:** 100.00% (no missing species_taxid)
- ✅ **Taxid coverage:** 100.00% (no missing taxid)
- ✅ **Unique species:** 149,270

### **Sample Verification (100K records)**
- **Average assemblies per species:** 9.78
- **Species coverage:** 100.00%
- **Method consistency:** ✅ Validated

### **Domain Distribution Example**
```
Domain        Assemblies    Species    Ratio (Assemblies/Species)
-----------   -----------   --------   -------------------------
Bacteria      1,314,591     11,873     110.7  (high research focus)
Viruses       211,871       34,313     6.2    (moderate coverage)
Eukaryota     69,088        9,035      7.6    (moderate coverage)
Archaea       344           111        3.1    (limited coverage)
```

---

## 🎯 **Why Our Method is Superior**

### **❌ Naive Approach (Wrong)**
```python
# This would OVERCOUNT species
species_count = len(df)  # Counts assemblies, not species
```

### **✅ Our Approach (Correct)**
```python
# This counts TRUE species diversity
species_count = df['species_taxid'].nunique()
```

### **Real Example**
- **Salmonella enterica** (species_taxid: 28901)
  - **Assemblies:** 20,000+ genome assemblies
  - **Species count:** 1 (correctly counted as one species)
  - **Why many assemblies:** Different strains, outbreaks, research studies

---

## 🔍 **Method Validation**

### **Cross-verification with Existing Reports**
Our verification script confirms:
- ✅ **Species count matches:** 149,270 unique species
- ✅ **Assembly count matches:** 2,945,415 total assemblies  
- ✅ **Domain breakdown matches:** Statistics report consistency
- ✅ **No data quality issues:** 100% coverage, no missing values

### **Biological Accuracy**
- **Escherichia coli:** 11,000+ assemblies → 1 species ✅
- **Staphylococcus aureus:** 3,800+ assemblies → 1 species ✅
- **Influenza A virus:** 3,800+ assemblies → 1 species ✅

---

## 🛠 **Implementation in Our Pipeline**

### **Parser Scripts Use This Method**
```python
# In our NCBI parsers (phylum_ncbi_parser.py, etc.)
species_counts = merged_df.groupby([level, 'domain']).agg({
    'species_taxid': 'nunique'  # Count unique species
}).rename(columns={'species_taxid': f'{level}_species_count'})

genome_counts = merged_df.groupby([level, 'domain']).agg({
    'assembly_accession': 'count'  # Count total assemblies  
}).rename(columns={'assembly_accession': f'{level}_genome_count'})
```

### **Output Files Include Both Metrics**
- `ncbi_phylum_counts.csv`: Phylum-level genome + species counts
- `ncbi_family_counts.csv`: Family-level genome + species counts  
- `ncbi_genus_counts.csv`: Genus-level genome + species counts

---

## 📋 **Quality Assurance**

### **Automated Verification**
Run our verification script:
```bash
python species_count_verification.py --sample-size 100000 --detailed
```

### **Expected Results**
- ✅ 100% species_taxid coverage
- ✅ Consistent species counts across samples
- ✅ Reasonable assembly-to-species ratios
- ✅ Domain distribution matches known biology

### **Red Flags to Watch For**
- ❌ Missing species_taxid values
- ❌ Unrealistic assembly-to-species ratios
- ❌ Inconsistent counts between runs
- ❌ Domain distributions that don't match biology

---

## 🎯 **Summary**

**Our method is scientifically rigorous and NCBI-standard:**

1. **Uses official NCBI species-level taxonomy** (species_taxid)
2. **Avoids double-counting** strains/isolates of same species  
3. **Provides dual metrics** (assemblies + species) for complete analysis
4. **100% data coverage** with no missing values
5. **Validated and verified** through multiple approaches
6. **Biologically accurate** species diversity measurements

**Result:** Reliable, accurate species counts that reflect true biological diversity in the NCBI database.

---

## 📁 **Related Files**

- `species_count_verification.py` - Verification script
- `ncbi_species_statistics.py` - Comprehensive statistics analysis
- `ncbi_species_statistics_report.txt` - Detailed statistics report
- `unified_ncbi_parser.py` - Main parsing implementation
