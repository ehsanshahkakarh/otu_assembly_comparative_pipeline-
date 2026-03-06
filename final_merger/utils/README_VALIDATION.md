# Validation Scripts for New Merger Pipeline

This directory contains validation scripts to ensure the quality and accuracy of merger outputs.

---

## Quick Start

```bash
# Run all validations at once
cd new_merger

# 1. Compare with old merger
python utils/sanity_check_merger.py

# 2. Validate data quality
python utils/validate_outputs.py

# 3. View summary dashboard
python utils/validation_dashboard.py
```

---

## Validation Scripts

### 1. `sanity_check_merger.py`

**Purpose:** Compare new merger outputs with old Eukcensus merger outputs

**What it checks:**
- Taxa overlap between old and new outputs
- Metric differences (genome counts, species counts, etc.)
- Identifies taxa present in old but not in new (and vice versa)
- Calculates percentage differences for all metrics

**Output:**
- `sanity_check/18s_sanity_check_YYYYMMDD_HHMMSS.txt`
- `sanity_check/16s_sanity_check_YYYYMMDD_HHMMSS.txt`

**Example output:**
```
DIVISION LEVEL COMPARISON
File Counts:
  Old file: 24 taxa
  New file: 22 taxa
  Matched: 22 taxa
  Old only: 2 taxa

Metric Comparisons:
  census_otu_count: 0.00% difference
  ncbi_genome_count: 5.69% difference
```

---

### 2. `validate_outputs.py`

**Purpose:** Comprehensive data quality validation

**What it checks:**
- **Data Integrity:** Missing values, negative numbers, invalid values
- **Domain Filtering:** Correct domain assignments (Eukaryota for 18S, Bacteria/Archaea for 16S)
- **Matching Logic:** Matched entries have NCBI data, census-only entries don't
- **Percentage Calculations:** OTU and size percentages sum to ~100%
- **Hierarchical Consistency:** Counts consistent across taxonomic levels

**Output:**
- `validation_reports/validation_report_YYYYMMDD_HHMMSS.txt`
- Console output with ✅/⚠️/❌ indicators

**Exit codes:**
- `0` = All checks passed
- `1` = Critical issues found

---

### 3. `validation_dashboard.py`

**Purpose:** Quick summary of all outputs

**What it shows:**
- File statistics (rows, file size, modification date)
- Key metrics (census OTUs, NCBI genomes/species)
- Comparison with old merger (match rates)
- Overall validation status

**Output:** Console only (no files generated)

**Example output:**
```
18S (EUKARYOTIC) MERGER OUTPUTS

📄 DIVISION
   File: 18s_ncbi_merged_division.csv (2.5 KB)
   Rows: 22 (14 matched, 8 census-only)
   Census OTUs: 70,899
   NCBI Genomes: 107,257
   Comparison: 91.7% match rate with old merger
```

---

## Validation Workflow

### After Running Merger

1. **Run sanity check** to compare with old outputs:
   ```bash
   python utils/sanity_check_merger.py
   ```

2. **Run comprehensive validation** to check data quality:
   ```bash
   python utils/validate_outputs.py
   ```

3. **View dashboard** for quick summary:
   ```bash
   python utils/validation_dashboard.py
   ```

4. **Review reports** in:
   - `sanity_check/` - Comparison reports
   - `validation_reports/` - Quality validation reports

---

## Understanding Validation Results

### ✅ All Checks Passed

Your outputs are valid and ready to use!

### ⚠️ Warnings

Minor issues that don't affect data quality:
- Percentage totals slightly off (96-104% is acceptable)
- Small variations across taxonomic levels
- Missing old merger files for comparison

### ❌ Critical Issues

Serious problems that need investigation:
- Missing values in critical columns
- Negative counts
- Invalid domain assignments
- Matched entries with zero NCBI data

---

## Expected Differences from Old Merger

### Fewer Taxa in New Merger

This is **expected and correct**:
- Stricter domain filtering
- Better lineage matching
- Removal of ambiguous entries
- Updated NCBI taxonomy

### Different NCBI Counts

This is **expected and correct**:
- NCBI database has grown
- Improved isolate counting
- Better domain filtering

### Match Rates

**Good match rates:**
- Division level: >80%
- Family level: >70%
- Genus level: >70%

**Low match rates at family/genus:**
- Old merger included many poorly classified taxa
- New merger is more selective and accurate

---

## Troubleshooting

### "No output files found"

Run the merger first:
```bash
python run_18s_ncbi_merger.py
python run_16s_ncbi_merger.py
```

### "Old file not found"

The old Eukcensus merger outputs are missing. This is OK - the sanity check will skip comparison but other validations will still run.

### "Critical issues found"

Check the validation report for details. Common issues:
- Corrupted CSV files
- Incomplete merger run
- Wrong input data

---

## Output Files

### Sanity Check Reports
- `sanity_check/18s_sanity_check_YYYYMMDD_HHMMSS.txt`
- `sanity_check/16s_sanity_check_YYYYMMDD_HHMMSS.txt`

### Validation Reports
- `validation_reports/validation_report_YYYYMMDD_HHMMSS.txt`

### Dashboard
- Console output only

---

## For Developers

### Adding New Validation Checks

Edit `validate_outputs.py` and add a new method to the `MergerValidator` class:

```python
def validate_my_check(self, df: pd.DataFrame, filename: str) -> bool:
    """Check something important."""
    # Your validation logic
    if problem_found:
        self.log_issue(f"{filename}: Problem description")
        return False
    
    self.log_pass(f"{filename}: Check passed")
    return True
```

Then add it to the `validate_file()` method.

---

**Last updated:** 2026-03-02

