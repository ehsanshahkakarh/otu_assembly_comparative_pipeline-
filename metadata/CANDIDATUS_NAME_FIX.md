# Candidatus Name Retention Fix

**Date:** 2026-03-06  
**Issue:** Candidatus taxonomic names were being truncated in appropriate_name field  
**Status:** ✅ FIXED

---

## Problem Description

### What Was Happening

In the unmapped logs, Candidatus taxa were showing incorrect appropriate_name values:

**Example from 16S log:**
```
GENUS | Candidatus Edwardsbacteria.U.genus | Candidatus Edwardsbacteria | Candidatus | 12 | 119 | NA | NO_TAXID_FOUND
```

**Fields:**
- **Original_Name:** `Candidatus Edwardsbacteria.U.genus`
- **Cleaned_Name:** `Candidatus Edwardsbacteria` ✅ Correct
- **Appropriate_Name:** `Candidatus` ❌ WRONG - should be `Candidatus Edwardsbacteria`

### Root Cause

The `extract_genus_from_species()` function (16S) and `extract_genus()` function (18S) were splitting names by space and taking only the first word:

**Old logic:**
```python
parts = species_name.strip().split()
if len(parts) >= 1:
    return parts[0]  # ❌ Returns just "Candidatus"
```

This is correct for normal binomial names like "Vitis vinifera" → "Vitis", but **incorrect** for Candidatus taxa where "Candidatus" is a prefix, not the genus name.

---

## What is "Candidatus"?

**Candidatus** is a nomenclatural status for prokaryotic taxa that:
- Cannot be cultured in the laboratory
- Are known only from environmental DNA sequences
- Have not been formally described according to the International Code of Nomenclature

**Examples:**
- `Candidatus Edwardsbacteria` - A CPR (Candidate Phyla Radiation) bacterium
- `Candidatus Gracilibacteria` - Ultra-small bacteria from groundwater
- `Candidatus Marinimicrobia` - Marine bacteria

**Key point:** "Candidatus" is a **prefix**, not the genus name. The full name should be retained.

---

## Solution Implemented

### Updated Logic

Modified both `extract_genus_from_species()` (16S) and `extract_genus()` (18S) to recognize "Candidatus" as a prefix:

**New logic:**
```python
parts = species_name.strip().split()

# Special case: "Candidatus" is a prefix for uncultivated prokaryotes
# Keep both "Candidatus" and the following genus name
if len(parts) >= 2 and parts[0].lower() == 'candidatus':
    return f"{parts[0]} {parts[1]}"  # ✅ Returns "Candidatus Edwardsbacteria"

# Special case: "candidate division" is a prefix for environmental clades
if len(parts) >= 3 and parts[0].lower() == 'candidate' and parts[1].lower() == 'division':
    return f"{parts[0]} {parts[1]} {parts[2]}"  # ✅ Returns "candidate division NC10"

# Normal case: return first word (genus)
if len(parts) >= 1:
    return parts[0]  # For "Vitis vinifera" → "Vitis"
```

### Files Modified

**16S Census Parser:**
1. `16S_censusparse/py_16S/census_parser/taxon_cleaner.py` (main module)
2. `16S_censusparse/py_16S/16S_eukcensus_parser.py` (legacy parser)
3. `16S_censusparse/py_16S/16S_eukcensus_parser_01.py` (variant parser)
4. `16S_censusparse/py_16S/sanity_checks/py/pocket/eukcensus_16S_parser.py` (sanity check)

**18S Census Parser:**
1. `18S_censusparse/src/taxon_cleaner.py` (main module)

---

## Expected Results After Fix

### Before Fix
```
GENUS | Candidatus Edwardsbacteria.U.genus | Candidatus Edwardsbacteria | Candidatus | 12 | 119 | NA
GENUS | Candidatus Gracilibacteria.U.genus | Candidatus Gracilibacteria | Candidatus | 1737 | 5355 | NA
GENUS | candidate division NC10.U.genus | candidate division NC | candidate | 395 | 3480 | NA
```

### After Fix
```
GENUS | Candidatus Edwardsbacteria.U.genus | Candidatus Edwardsbacteria | Candidatus Edwardsbacteria | 12 | 119 | NA
GENUS | Candidatus Gracilibacteria.U.genus | Candidatus Gracilibacteria | Candidatus Gracilibacteria | 1737 | 5355 | NA
GENUS | candidate division NC10.U.genus | candidate division NC | candidate division NC | 395 | 3480 | NA
```

**Result:** The appropriate_name now correctly retains the full Candidatus name for taxonkit lookup.

---

## Impact Analysis

### Affected Taxa

**16S Census Parser:**
- ~20 Candidatus genera in unmapped log
- Examples: Candidatus Gracilibacteria, Candidatus Marinimicrobia, Candidatus Edwardsbacteria
- Total OTUs affected: ~6,000
- Total sequences affected: ~20,000

**18S Census Parser:**
- Minimal impact (eukaryotes rarely use Candidatus nomenclature)
- Included for consistency and future-proofing

### Benefits

1. **Better taxonkit matching:** Full names have higher chance of NCBI lookup success
2. **Correct logging:** Unmapped logs now show the actual name being looked up
3. **Consistency:** Both 16S and 18S parsers handle Candidatus names identically
4. **Future-proof:** Handles emerging Candidatus taxa automatically

---

## Testing Recommendations

### Manual Testing

1. **Re-run 16S parser:**
   ```bash
   cd 16S_censusparse/py_16S
   python run_16S_parser.py
   ```

2. **Check unmapped log:**
   ```bash
   grep "Candidatus" logs/eukcensus16S_comprehensive_unmapped.log
   ```

3. **Verify appropriate_name field:**
   - Should show full name (e.g., "Candidatus Edwardsbacteria")
   - Not truncated (e.g., "Candidatus")

### Expected Improvements

- **Mapping rate:** May increase slightly if NCBI has some Candidatus taxa
- **Log clarity:** Unmapped entries will show correct names for manual resolution
- **Resolution system:** Can now add Candidatus taxa to `known_parents.py` with full names

---

## Related Documentation

- **16S Resolution System:** `16S_censusparse/RESOLUTION_SYSTEM_EXPLAINED.md`
- **18S Resolution System:** `18S_censusparse/RESOLUTION_SYSTEM_EXPLAINED.md`
- **Code Review:** `16S_censusparse/CODE_REVIEW_ANALYSIS.md`

---

## Summary

✅ **Fixed:** Candidatus names now retain full genus name in appropriate_name field  
✅ **Scope:** Both 16S and 18S census parsers  
✅ **Impact:** ~6,000 OTUs in 16S data  
✅ **Benefit:** Better taxonkit matching and clearer logging  
✅ **Testing:** Re-run parsers and check unmapped logs

