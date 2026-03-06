# Sanity Check Tools for New Merger Outputs

This directory contains validation tools to ensure the quality and correctness of the new_merger outputs.

## Available Sanity Checks

### 1. `sanity_check_merger.py` - Compare Old vs New Merger
**Purpose:** Compare outputs from the old Eukcensus merger with the new modular merger to ensure consistency.

**Usage:**
```bash
python sanity_check_merger.py
```

**What it checks:**
- Taxon overlap between old and new outputs
- Metric comparisons (census counts, NCBI counts, novelty factors)
- Identifies taxa present in one output but not the other
- Generates detailed comparison reports

**Output:** `sanity_check/18s_sanity_check_TIMESTAMP.txt` and `16s_sanity_check_TIMESTAMP.txt`

---

### 2. `sanity_check_output_quality.py` - Validate Output Quality (NEW!)
**Purpose:** Comprehensive validation of new_merger outputs to ensure data integrity, mathematical consistency, and biological plausibility.

**Usage:**
```bash
# Validate 18S outputs only
python sanity_check_output_quality.py --gene 18S

# Validate 16S outputs only
python sanity_check_output_quality.py --gene 16S

# Validate both 18S and 16S
python sanity_check_output_quality.py --gene both
```

**What it checks:**

#### 1. Data Integrity
- ✅ All required columns are present
- ✅ No missing taxon names
- ✅ No negative counts (census_otu_count, ncbi_genome_count, etc.)
- ✅ Valid percentage ranges (0-100%)

#### 2. Mathematical Consistency
- ✅ Novelty factor = census_otu_count / ncbi_species_count
- ✅ Overrepresentation factor = ncbi_species_count / census_otu_count
- ✅ Handles special cases (division by zero → inf or 0)
- ✅ Tolerates rounding differences (0.5% relative or 0.001 absolute)

#### 3. Biological Plausibility
- ⚠️  Flags extremely high novelty factors (>1000) as warnings
- ⚠️  Identifies census_only taxa (expected for underrepresented groups)
- ⚠️  Identifies ncbi_only taxa (rare, may indicate issues)

#### 4. Cross-Level Consistency
- ✅ Total census counts are similar across division/family/genus levels
- ✅ Allows for ~10% variation (due to unclassified taxa)

**Output:** `sanity_check/quality_check_GENE_TIMESTAMP.txt`

**Exit codes:**
- `0` - All checks passed
- `1` - Errors found (review report)

---

## Example Workflow

### After running the merger:
```bash
# 1. Run the merger
cd /path/to/new_merger
python run_18s_ncbi_merger.py

# 2. Validate output quality
python utils/sanity_check_output_quality.py --gene 18S

# 3. Compare with old merger (optional)
python utils/sanity_check_merger.py
```

### Interpreting Results

#### ✅ All Checks Passed
```
✅ All validation checks passed!
The new_merger outputs are of high quality and ready for use.
```
Your outputs are ready to use!

#### ⚠️  Warnings Only
```
⚠️  Some warnings were found. These may be expected (e.g., census_only taxa)
but should be reviewed to ensure they are not data quality issues.
```
Review the warnings - they may be expected (e.g., taxa with no NCBI genomes).

#### ❌ Errors Found
```
❌ Validation found errors that need to be addressed.
Please review the errors above and fix the issues.
```
Review the report and fix the issues before using the outputs.

---

## Common Warnings and What They Mean

### `EXTREME_VALUES: X taxa with novelty_factor > 1000`
**Meaning:** Some taxa have very high novelty factors (census OTUs >> NCBI genomes).

**Is this OK?** Usually yes! This indicates underrepresented groups in genome databases.

**Example:** Tubulinea has 1323 census OTUs but only 1 NCBI species → novelty_factor = 1323

### `census_only: X taxa`
**Meaning:** Taxa found in census data but not in NCBI database.

**Is this OK?** Yes! This is expected for novel or underrepresented lineages.

**Example:** TSAR.U.division, Archaeplastida.U.division (unclassified divisions)

---

## Report Files

All reports are saved to: `new_merger/sanity_check/`

**File naming:**
- `quality_check_18s_TIMESTAMP.txt` - 18S quality validation
- `quality_check_16s_TIMESTAMP.txt` - 16S quality validation
- `quality_check_all_TIMESTAMP.txt` - Both 18S and 16S
- `18s_sanity_check_TIMESTAMP.txt` - Old vs new comparison (18S)
- `16s_sanity_check_TIMESTAMP.txt` - Old vs new comparison (16S)

---

## Troubleshooting

### "Missing file" errors
**Problem:** Output files not found in `outputs/` directory.

**Solution:** Run the merger first:
```bash
python run_18s_ncbi_merger.py  # for 18S
python run_16s_ncbi_merger.py  # for 16S
```

### Calculation errors
**Problem:** Novelty or overrepresentation factors don't match expected values.

**Likely cause:** Rounding differences (should be tolerated by the check).

**If persistent:** Check the merger calculation logic in `src/metrics_calculator.py`

---

## For Developers

### Adding New Checks

To add a new validation check:

1. Add a new method to `OutputQualityChecker` class:
```python
def check_my_new_validation(self, df: pd.DataFrame, level: str, gene: str):
    """Check something new."""
    # Your validation logic
    if problem_found:
        self.log_issue('ERROR', 'category', 'description')
```

2. Call it in `validate_gene()` method:
```python
self.check_my_new_validation(df, level, gene)
```

3. Test it:
```bash
python sanity_check_output_quality.py --gene 18S
```

---

## Summary

✅ **Use `sanity_check_output_quality.py`** after every merger run to validate outputs

✅ **Use `sanity_check_merger.py`** to compare with old merger (one-time validation)

✅ **Review warnings** - they're usually expected but worth checking

✅ **Fix errors** before using outputs downstream

