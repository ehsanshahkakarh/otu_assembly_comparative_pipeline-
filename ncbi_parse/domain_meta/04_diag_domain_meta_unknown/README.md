# Diagnostic 04: Domain Meta Unknown Species

This diagnostic investigates the unknown species from the domain_meta pipeline output.

## Sub-Diagnostic 04.01: TaxID vs Species_TaxID Verification

### `01_check_all_unknown_species.py` ⭐ **Main Script**

Checks if `taxid == species_taxid` for **ALL** unknown species listed in `unknown_species.csv`.

**Purpose**: Verify whether all assemblies for unknown species have matching taxid and species_taxid values in the raw NCBI assembly file.

**Usage**:
```bash
python 01_check_all_unknown_species.py
```

**Output**:
- Console output with summary statistics
- **Text file only**: `01_unknown_species_taxid_check_<timestamp>.txt` (no CSV/TSV files)
- Result clearly states: "✓ CONFIRMED" or "✗ DENIED"

**Key Findings** (as of 2026-03-05):
- Total unknown species: 155
- Total assemblies checked: 35,849
- **Result**: ✓ CONFIRMED - All assemblies (100%) have `taxid == species_taxid`
- This means that for all unknown species, every assembly record has the same taxid and species_taxid value
- **Conclusion**: Cannot resolve unknown species through taxid differences

---

### `01_extract_taxid_749906.py` (Helper)

Extracts all rows with `species_taxid 749906` (gut metagenome) from the original NCBI assembly summary file.

**Usage**:
```bash
python 01_extract_taxid_749906.py
```

---

### `01_extract_species_taxid.py` (Helper)

General-purpose script to extract rows for **any** species_taxid from the NCBI assembly summary file.

**Usage**:
```bash
python 01_extract_species_taxid.py <species_taxid>
```

---

## Sub-Diagnostic 04.02: Name Matching Verification

### `02_verify_name_matching.py` ⭐ **Main Script**

Randomly samples genus names from 16S and 18S census data and compares them with the new_merger output to verify name-matching accuracy.

**Purpose**: Verify that the name matching logic between census data and NCBI taxonomy is working correctly by manually checking random samples.

**Usage**:
```bash
python 02_verify_name_matching.py
```

**Output**:
- Console output with sampled entries
- **Text file only**: `02_name_matching_verification_<timestamp>.txt`
- Lists 10 random samples from each dataset (16S and 18S)
- Includes merger summary statistics for context

**Key Information** (as of 2026-03-05):
- **16S Census**: 4,578 total entries
  - Genus-level match rate: 45.9% (2,777 matched / 6,049 total)
  - Census-only taxa: 1,798 genera
  - NCBI-only taxa: 1,474 genera

- **18S Census**: 491 total entries
  - Genus-level match rate: 1.5% (135 matched / 9,080 total)
  - Census-only taxa: 356 genera
  - NCBI-only taxa: 8,589 genera

**Sample Findings**:
- Random seed: 42 (for reproducibility)
- Sampled names include both regular genus names and `.U.genus` (unclassified) entries
- Some entries have `nan` TaxID values (e.g., organellar sequences)
- Manual verification needed to confirm TaxID accuracy

---

## Background

The unknown species in the domain_meta pipeline output are primarily metagenomes and environmental samples. Sub-Diagnostic 04.01 confirmed that all unknown species have `taxid == species_taxid`, which is expected for metagenomes (no strain-level classification). This means we cannot resolve these entries through taxid differences alone.

