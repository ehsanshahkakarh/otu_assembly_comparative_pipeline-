# Pipeline Improvement Recommendations

## ✅ CRITICAL FIX IMPLEMENTED (2026-03-02)

**The pipeline has been FIXED to use the correct input files!**

### What Was Fixed

The pipeline was using **archaic pre-aggregated files** (`ncbi_{level}_counts.csv`) instead of the correct **species-level data** (`species_grouped_*.csv`). This has been corrected.

### Results

✅ **18S Merger Completed Successfully** (67 seconds)
- Division level: 22 taxa matched
- Family level: 314 taxa matched
- Genus level: 491 taxa matched
- Progress bars working correctly
- Isolate information included in outputs

### Key Changes Made

1. **config.py** - Points to `nev_parse_meth/output/` directory
2. **data_loader.py** - Loads species-level data with domain extraction
3. **domain_filter.py** - Simplified to work with species data only
4. **lineage_matcher.py** - Matches species first, then aggregates counts
5. **run_*.py scripts** - Load species data once, reuse for all levels
6. **Progress bars** - Added for level processing (division, family, genus)

### Expected Benefits

- ✅ **Accuracy**: Uses correct source data
- ✅ **No Data Loss**: Aggregation happens AFTER matching
- ✅ **Better Performance**: Loads data once instead of 3 times
- ✅ **Isolate Tracking**: Included directly from species data
- ✅ **Ground-Truth Match**: Should now match validation (within 5%)

---

# Additional Improvement Recommendations

Below are additional improvements that could be made to further enhance the pipeline:

**Date:** 2026-03-02  
**Status:** Based on validation work and codebase analysis

---

## Executive Summary

The current merger pipeline is **well-designed and validated**, but there are opportunities for improvement in:
1. **Performance** - Speed up processing of large datasets
2. **Accuracy** - Reduce the 4-45% discrepancy found in ground-truth validation
3. **Maintainability** - Simplify architecture and reduce code duplication
4. **Scalability** - Handle growing NCBI database more efficiently
5. **Usability** - Better error handling and user feedback

---

## 🚀 High-Priority Improvements

### 1. **Fix Aggregation Discrepancy** ⭐⭐⭐

**Problem:** Ground-truth validation shows 4-45% more species/genomes than merger output

**Root Cause:** The merger uses pre-aggregated data (`ncbi_phylum_counts.csv`) while the ground-truth validator works directly with species-level data (`species_grouped_*.csv`)

**Solution:**
```python
# CURRENT: Uses pre-aggregated data
ncbi_file = base_dir / "ncbi_parse/csv_ncbi" / f"ncbi_{level}_counts.csv"

# PROPOSED: Use species-level data directly
species_file = base_dir / "ncbi_parse/metadata/nev_parse_meth/output" / "species_grouped_*.csv"
# Then aggregate on-the-fly during matching
```

**Benefits:**
- ✅ More accurate counts (matches ground-truth)
- ✅ No data loss from pre-aggregation
- ✅ Better tracking of individual species
- ✅ Easier to debug matching issues

**Implementation:**
1. Modify `data_loader.py` to load species-level data
2. Update `lineage_matcher.py` to aggregate after matching
3. Add caching to avoid re-loading large files
4. Update tests and validation

**Estimated Impact:** Reduces discrepancy from 4-45% to <5%

---

### 2. **Implement Caching for Large Files** ⭐⭐⭐

**Problem:** `species_grouped_*.csv` is 61MB and loaded multiple times (once per level)

**Current Performance:**
- Loading time: ~5-10 seconds per level
- Total loads: 3 levels × 2 census types = 6 loads
- Total time wasted: ~30-60 seconds

**Solution:**
```python
# Add to config.py
class CacheConfig:
    """Caching configuration."""
    enable_cache = True
    cache_dir = Path("new_merger/.cache")
    cache_species_data = True

# Add to data_loader.py
import pickle
from functools import lru_cache

@lru_cache(maxsize=2)
def load_species_data_cached(species_file: Path) -> pd.DataFrame:
    """Load species data with caching."""
    cache_file = Config.cache_dir / f"{species_file.stem}.pkl"
    
    if cache_file.exists() and cache_file.stat().st_mtime > species_file.stat().st_mtime:
        logger.info(f"  Loading from cache: {cache_file.name}")
        return pd.read_pickle(cache_file)
    
    logger.info(f"  Loading from CSV: {species_file.name}")
    df = pd.read_csv(species_file, low_memory=False)
    
    # Save to cache
    Config.cache_dir.mkdir(exist_ok=True)
    df.to_pickle(cache_file)
    
    return df
```

**Benefits:**
- ✅ 10-20x faster subsequent loads
- ✅ Reduces total runtime by 30-60 seconds
- ✅ Minimal code changes

**Estimated Impact:** 40-60% faster pipeline execution

---

### 3. **Parallelize Level Processing** ⭐⭐

**Problem:** Division, family, and genus levels are processed sequentially

**Current:**
```python
for level in ['phylum', 'family', 'genus']:
    process_level(level)  # Sequential
```

**Proposed:**
```python
from concurrent.futures import ProcessPoolExecutor

def process_all_levels_parallel():
    """Process all levels in parallel."""
    levels = ['phylum', 'family', 'genus']
    
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(process_level, level): level for level in levels}
        
        for future in as_completed(futures):
            level = futures[future]
            try:
                result = future.result()
                logger.info(f"✅ {level} completed")
            except Exception as e:
                logger.error(f"❌ {level} failed: {e}")
```

**Benefits:**
- ✅ 3x faster (if CPU-bound)
- ✅ Better resource utilization
- ✅ Independent level processing

**Caveats:**
- ⚠️ Requires more memory (3x data in memory)
- ⚠️ May not help if I/O-bound

**Estimated Impact:** 2-3x faster on multi-core systems

---

### 4. **Add Progress Bars and Better Logging** ⭐⭐

**Problem:** No visual feedback during long-running operations

**Solution:**
```python
from tqdm import tqdm

def match_taxa_by_lineage(census_df, ncbi_df, level, census_level):
    """Match with progress bar."""
    results = []
    
    for _, row in tqdm(census_df.iterrows(), 
                       total=len(census_df),
                       desc=f"Matching {census_level}",
                       unit="taxa"):
        # Matching logic
        ...
    
    return results
```

**Benefits:**
- ✅ User knows progress
- ✅ Can estimate completion time
- ✅ Better debugging (see where it's slow)

---

### 5. **Implement Incremental Updates** ⭐⭐

**Problem:** Must re-run entire pipeline even if only NCBI data changed

**Solution:**
```python
class IncrementalMerger:
    """Merger with incremental update support."""
    
    def needs_update(self, output_file, input_files):
        """Check if output needs updating."""
        if not output_file.exists():
            return True
        
        output_mtime = output_file.stat().st_mtime
        for input_file in input_files:
            if input_file.stat().st_mtime > output_mtime:
                return True
        
        return False
    
    def run_incremental(self):
        """Run only if inputs changed."""
        for level in levels:
            if self.needs_update(output_file, [census_file, ncbi_file]):
                logger.info(f"  {level}: Updating (inputs changed)")
                process_level(level)
            else:
                logger.info(f"  {level}: Skipping (up-to-date)")
```

**Benefits:**
- ✅ Faster re-runs
- ✅ Only process what changed
- ✅ Make-like dependency tracking

---

## 📊 Medium-Priority Improvements

### 6. **Vectorize Matching Logic** ⭐⭐

**Problem:** Current matching uses row-by-row iteration

**Current:**
```python
for _, census_row in census_df.iterrows():
    # Match each row individually
    matches = find_matches(census_row, ncbi_df)
```

**Proposed:**
```python
# Fully vectorized matching
def vectorized_match(census_df, ncbi_df):
    """Vectorized matching using pandas operations."""
    # Pre-process lineage_taxids into sets
    ncbi_df['taxid_set'] = ncbi_df['lineage_taxids'].str.split(';').apply(set)
    
    # Vectorized taxid matching
    results = []
    for census_taxid in census_df['taxid'].unique():
        mask = ncbi_df['taxid_set'].apply(lambda x: str(census_taxid) in x)
        matched = ncbi_df[mask]
        results.append(matched)
    
    return pd.concat(results)
```

**Benefits:**
- ✅ 5-10x faster matching
- ✅ Better pandas utilization
- ✅ Less memory overhead

---

### 7. **Add Data Quality Checks** ⭐⭐

**Problem:** No validation of input data quality

**Solution:**
```python
class DataQualityChecker:
    """Validate input data quality."""
    
    def check_census_data(self, df):
        """Check census data quality."""
        issues = []
        
        # Check for missing values
        if df['otu_count'].isna().any():
            issues.append("Missing OTU counts")
        
        # Check for negative values
        if (df['otu_count'] < 0).any():
            issues.append("Negative OTU counts")
        
        # Check for duplicate taxa
        if df['Name_to_use'].duplicated().any():
            issues.append("Duplicate taxon names")
        
        return issues
    
    def check_ncbi_data(self, df):
        """Check NCBI data quality."""
        issues = []
        
        # Check lineage format
        if not df['lineage'].str.contains(';').all():
            issues.append("Invalid lineage format")
        
        return issues
```

**Benefits:**
- ✅ Catch data issues early
- ✅ Better error messages
- ✅ Prevent garbage in, garbage out

---

### 8. **Create Unified Entry Point** ⭐

**Problem:** Separate scripts for 18S and 16S

**Solution:**
```python
# run_merger.py
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--census-type', choices=['18S', '16S', 'both'], default='both')
    parser.add_argument('--levels', nargs='+', choices=['division', 'family', 'genus'])
    parser.add_argument('--parallel', action='store_true')
    args = parser.parse_args()
    
    if args.census_type in ['18S', 'both']:
        run_18s_merger(args)
    
    if args.census_type in ['16S', 'both']:
        run_16s_merger(args)
```

**Benefits:**
- ✅ Single command for all mergers
- ✅ Easier to use
- ✅ Consistent interface

---

## 🔧 Low-Priority Improvements

### 9. **Add Configuration File Support**

```yaml
# merger_config.yaml
census:
  18S:
    domain: Eukaryota
    levels: [division, family, genus]
  16S:
    domains: [Bacteria, Archaea]
    levels: [division, family, genus]

performance:
  parallel: true
  cache: true
  workers: 3

output:
  format: csv
  compression: none
```

### 10. **Generate Summary Statistics**

```python
def generate_summary_report(output_files):
    """Generate summary statistics."""
    report = {
        'total_taxa': 0,
        'matched_taxa': 0,
        'total_genomes': 0,
        'total_species': 0,
        'match_rate': 0.0
    }
    
    for file in output_files:
        df = pd.read_csv(file)
        report['total_taxa'] += len(df)
        report['matched_taxa'] += (df['match_status'] == 'matched').sum()
        # ... more stats
    
    return report
```

### 11. **Add Export Formats**

```python
def export_to_format(df, output_file, format='csv'):
    """Export to multiple formats."""
    if format == 'csv':
        df.to_csv(output_file, index=False)
    elif format == 'excel':
        df.to_excel(output_file.with_suffix('.xlsx'), index=False)
    elif format == 'json':
        df.to_json(output_file.with_suffix('.json'), orient='records')
    elif format == 'parquet':
        df.to_parquet(output_file.with_suffix('.parquet'), index=False)
```

---

## 📈 Performance Optimization Summary

| Improvement | Speed Gain | Effort | Priority |
|-------------|-----------|--------|----------|
| Fix aggregation discrepancy | Accuracy +40% | High | ⭐⭐⭐ |
| Implement caching | 40-60% faster | Low | ⭐⭐⭐ |
| Parallelize levels | 2-3x faster | Medium | ⭐⭐ |
| Vectorize matching | 5-10x faster | High | ⭐⭐ |
| Progress bars | UX only | Low | ⭐⭐ |
| Incremental updates | Variable | Medium | ⭐⭐ |

**Combined Impact:** Could achieve **5-10x overall speedup** with all optimizations

---

## 🎯 Recommended Implementation Order

1. **Phase 1 (Quick Wins)**
   - Add caching (1-2 hours)
   - Add progress bars (1 hour)
   - Fix aggregation to use species-level data (4-6 hours)

2. **Phase 2 (Performance)**
   - Vectorize matching logic (8-12 hours)
   - Parallelize level processing (4-6 hours)

3. **Phase 3 (Quality)**
   - Add data quality checks (4-6 hours)
   - Implement incremental updates (6-8 hours)

4. **Phase 4 (Polish)**
   - Unified entry point (2-4 hours)
   - Configuration file support (4-6 hours)
   - Summary statistics (2-4 hours)

**Total Estimated Effort:** 35-55 hours

---

## 🔍 Architecture Improvements

### Current Architecture Issues

1. **Code Duplication:** `run_18s_ncbi_merger.py` and `run_16s_ncbi_merger.py` are nearly identical
2. **Multiple Data Paths:** Uses both aggregated (`ncbi_*_counts.csv`) and species-level data
3. **No Dependency Tracking:** Can't tell if outputs are stale
4. **Limited Error Recovery:** Fails completely if one level fails

### Proposed Architecture

```
new_merger/
├── run_merger.py              # Unified entry point
├── config.yaml                # User-editable config
├── src/
│   ├── core/
│   │   ├── pipeline.py        # Main pipeline orchestrator
│   │   ├── cache_manager.py   # Caching logic
│   │   └── dependency_tracker.py
│   ├── loaders/
│   │   ├── census_loader.py
│   │   ├── species_loader.py  # Direct species-level loading
│   │   └── cache_loader.py
│   ├── matchers/
│   │   ├── vectorized_matcher.py
│   │   └── hierarchical_matcher.py
│   └── validators/
│       ├── data_quality.py
│       └── output_validator.py
└── tests/
    ├── test_pipeline.py
    ├── test_matchers.py
    └── test_validators.py
```

---

## ✅ Validation Improvements

Based on the ground-truth validation work:

1. **Automated Regression Testing**
   - Run ground-truth validation after each merger run
   - Alert if discrepancy > 10%
   - Track discrepancy trends over time

2. **Continuous Validation**
   - Add to CI/CD pipeline
   - Run on every code change
   - Prevent regressions

3. **Better Discrepancy Reporting**
   - Show which taxa have largest discrepancies
   - Identify systematic issues
   - Suggest fixes

---

**Last Updated:** 2026-03-02

