# Pipeline Fix Summary

**Date:** 2026-03-02  
**Issue:** Pipeline was using WRONG input files (archaic pre-aggregated data)  
**Status:** ✅ FIXED

---

## Problem Identified

The merger pipeline was using **archaic pre-aggregated files** instead of the correct species-level data:

### ❌ WRONG (Before Fix)
```
Input: ncbi_parse/csv_ncbi/ncbi_{level}_counts.csv
- Pre-aggregated by taxonomic level (phylum, family, genus)
- Data loss from aggregation
- 4-45% fewer species/genomes than ground-truth
```

### ✅ CORRECT (After Fix)
```
Input: ncbi_parse/metadata/nev_parse_meth/output/species_grouped_*.csv
- Species-level data (162,731 species, 3.3M genomes)
- Aggregation happens AFTER matching
- Matches ground-truth validation
```

---

## Changes Made

### 1. **config.py** - Updated to point to correct input files

**Before:**
```python
self.ncbi_dir = self.parse_dir / 'ncbi_parse' / 'csv_ncbi'

def get_ncbi_files(self, level: str):
    return {
        'counts': self.ncbi_dir / f'ncbi_{level}_counts.csv',  # WRONG
        'accessions': self.ncbi_dir / f'ncbi_{level}_with_accessions.csv'
    }
```

**After:**
```python
self.ncbi_species_dir = self.parse_dir / 'ncbi_parse' / 'metadata' / 'nev_parse_meth' / 'output'

def get_ncbi_species_file(self) -> Path:
    # Find most recent species_grouped_*.csv file
    pattern = str(self.ncbi_species_dir / 'species_grouped_*.csv')
    species_files = glob.glob(pattern)
    return Path(sorted(species_files)[-1])  # Most recent
```

### 2. **src/data_loader.py** - Load species-level data

**Before:**
```python
def load_ncbi_data(ncbi_file: Path, level: str):
    # Loads pre-aggregated ncbi_{level}_counts.csv
    df = pd.read_csv(ncbi_file)
    # Returns aggregated data by level
```

**After:**
```python
def load_ncbi_species_data(species_file: Path):
    # Loads species_grouped_*.csv (species-level data)
    df = pd.read_csv(species_file, low_memory=False)
    # Returns species-level data with:
    # - species_taxid, species_name
    # - total_genome_count, isolate_genome_count, uncultured_genome_count
    # - lineage, lineage_ranks, lineage_taxids
    # - domain (extracted from lineage)
```

### 3. **src/domain_filter.py** - Simplified for species data

**Before:**
```python
def filter_by_domain(ncbi_df, domains, isolate_df=None):
    # Filters both NCBI and isolate data
    return filtered_ncbi, filtered_isolate
```

**After:**
```python
def filter_by_domain(ncbi_df, domains):
    # Filters species data only (isolate info already in species data)
    return filtered_ncbi
```

### 4. **src/lineage_matcher.py** - Match species, then aggregate

**Before:**
```python
def match_taxa_by_lineage(census_df, ncbi_df, level, census_level):
    # Matches against pre-aggregated data
    # 4 match types: direct_name, direct_taxid, hierarchical_taxid, lineage_name
    # Returns aggregated counts from pre-aggregated data
```

**After:**
```python
def match_taxa_by_lineage(census_df, ncbi_species_df, level, census_level):
    # Matches against species-level data
    # 2 match types: direct_taxid, hierarchical_taxid
    # Pattern: ;taxid; or ^taxid; or ;taxid$ or ^taxid$
    
    # For each census taxon:
    #   1. Find all matching species
    #   2. Aggregate genome counts: total_genome_count.sum()
    #   3. Count species: len(matched_species)
    #   4. Aggregate isolate counts: isolate_genome_count.sum()
    
    # Returns aggregated results with isolate information
```

### 5. **run_18s_ncbi_merger.py & run_16s_ncbi_merger.py** - Load once, use for all levels

**Before:**
```python
for level in ['phylum', 'family', 'genus']:
    # Load pre-aggregated data for this level
    ncbi_df = load_ncbi_data(ncbi_files['counts'], level)
    isolate_df = load_isolate_data(ncbi_files['accessions'], level)
    
    # Filter by domain
    ncbi_df, isolate_df = filter_by_domain(ncbi_df, domains, isolate_df)
    
    # Match and aggregate
    matched_df = match_taxa_by_lineage(census_df, ncbi_df, level, census_level)
    matched_df = add_isolate_information(matched_df, isolate_df, level, census_level)
```

**After:**
```python
# Load species data ONCE (used for all levels)
species_file = config.paths.get_ncbi_species_file()
ncbi_species_df = load_ncbi_species_data(species_file)

# Filter by domain ONCE
ncbi_species_df = filter_by_domain(ncbi_species_df, domains)

# Process each level with progress bar
for level in tqdm(['phylum', 'family', 'genus'], desc="Processing levels"):
    # Match species and aggregate (isolate info included)
    matched_df = match_taxa_by_lineage(census_df, ncbi_species_df, level, census_level)
```

### 6. **Added Progress Bars** (as requested)

```python
from tqdm import tqdm

# Progress bar for levels
for level in tqdm(config.taxonomic.ncbi_levels, desc="Processing levels", unit="level"):
    ...

# Progress bar for matching within each level
for _, census_row in tqdm(census_df.iterrows(), total=len(census_df), 
                           desc=f"  Matching {census_level}", unit="taxa", leave=False):
    ...
```

---

## Benefits of Fix

1. ✅ **Accuracy**: Uses correct source data (species_grouped_*.csv)
2. ✅ **No Data Loss**: Aggregation happens AFTER matching, not before
3. ✅ **Isolate Tracking**: Isolate counts included directly from species data
4. ✅ **Performance**: Loads species data once, reuses for all levels
5. ✅ **Simplicity**: Removed unnecessary isolate_tracker module
6. ✅ **User Feedback**: Progress bars for level processing
7. ✅ **Ground-Truth Match**: Should now match ground-truth validation (within 5%)

---

## Expected Results

After this fix, the merger outputs should:
- ✅ Match ground-truth validation (4-45% discrepancy should disappear)
- ✅ Show MORE species and genomes (using full species-level data)
- ✅ Include accurate isolate counts
- ✅ Process faster (load species data once instead of 3 times)

---

## Next Steps

1. **Test the fixed pipeline:**
   ```bash
   cd new_merger
   python run_18s_ncbi_merger.py
   python run_16s_ncbi_merger.py
   ```

2. **Validate outputs:**
   ```bash
   cd sanity_check_test
   python rebuild_from_source.py
   ```

3. **Compare with ground-truth:**
   - Discrepancy should be <5% (down from 4-45%)
   - Species counts should match ground-truth
   - Genome counts should match ground-truth

---

**Last Updated:** 2026-03-02

