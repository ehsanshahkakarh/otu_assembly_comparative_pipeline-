# Enhanced 16S Census Parser Improvements

## Issues Identified and Fixed

### 1. **Taxonkit Version Problem** ⭐ **MAJOR FIX**
- **Problem**: Parser was using local `./taxonkit` (v0.13.0) with incomplete taxonomy database
- **Solution**: Updated to use system `taxonkit` (v0.20.0) with complete NCBI taxonomy
- **Impact**: Resolves taxa like Procabacter (167946), Candidatus Sumerlaea (2315418), etc.

### 2. **Candidatus Taxa Handling** 
- **Problem**: Automatically stripping "Candidatus " prefix caused failures
- **Solution**: Enhanced fallback logic:
  1. Try original "Candidatus [name]" first
  2. Only try stripped version if original fails
  3. Both often resolve to same taxid (e.g., Candidatus Cardinium = Cardinium = 273135)

### 3. **Complex Uncultured Names**
- **Problem**: Names like "uncultured_Alphaproteobacteria_bacterium" reduced to just "uncultured"
- **Solution**: Extract meaningful taxonomic parts:
  - `uncultured_Alphaproteobacteria_bacterium` → `Alphaproteobacteria` (taxid 28211)
  - `uncultured_Rickettsia_sp` → `Rickettsia` (taxid 780)
  - Skip purely environmental terms (metagenome, environmental, etc.)

### 4. **Enhanced Fallback Strategy System**
Multiple strategies tried in order:
1. **Original name** (preserves Candidatus, handles exact matches)
2. **Cleaned name** (organelle removal, underscore handling)  
3. **Candidatus stripped** (for Candidatus taxa only)
4. **Meaningful part extraction** (for complex uncultured names)
5. **Number removal** (for names with trailing numbers)

## Test Results

### Before (Local taxonkit v0.13.0):
```
Procabacter                           → NO TAXID
Candidatus Sumerlaea                  → NO TAXID  
uncultured_Alphaproteobacteria_bacterium → NO TAXID
```

### After (System taxonkit v0.20.0 + Enhanced Logic):
```
Procabacter                           → 167946 ✅
Candidatus Sumerlaea                  → 2315418 ✅
Sumerlaea                            → 2315418 ✅ (fallback)
uncultured_Alphaproteobacteria_bacterium → 28211 ✅ (via Alphaproteobacteria)
uncultured_Rickettsia_sp             → 780 ✅ (via Rickettsia)
```

## Key Improvements in Code

### 1. **System Taxonkit Usage**
```python
# OLD: ["../../taxonkit", "name2taxid"]
# NEW: ["taxonkit", "name2taxid"]
```

### 2. **Enhanced Meaningful Part Extraction**
```python
def extract_meaningful_taxonomic_part(taxon_name):
    """Extract Alphaproteobacteria from uncultured_Alphaproteobacteria_bacterium"""
    if 'uncultured' in taxon_name.lower():
        parts = taxon_name.replace('_', ' ').split()
        generic_terms = ['uncultured', 'bacterium', 'organism', 'eukaryote', 'sp']
        # Find meaningful taxonomic terms
        meaningful_parts = [p for p in parts if p.lower() not in generic_terms]
        return meaningful_parts[-1] if meaningful_parts else None
```

### 3. **Multi-Strategy Fallback**
```python
def get_single_taxid_with_fallbacks(taxon_name, env):
    strategies = [
        ("original", taxon_name),
        ("cleaned", clean_organelle_taxon_name(taxon_name)),
        ("candidatus_stripped", taxon_name[11:] if taxon_name.startswith("Candidatus ") else None),
        ("meaningful_part", extract_meaningful_taxonomic_part(taxon_name)),
        ("no_numbers", strip_trailing_numbers(taxon_name))
    ]
    # Try each strategy until one succeeds
```

## Expected Impact

### Mapping Success Rate Improvements:
- **Procabacter**: 0% → 100% (8 occurrences)
- **Candidatus taxa**: ~20% → ~90% (59 names, 20,542 occurrences)
- **Complex uncultured**: ~5% → ~60% (22 names, 2,816 occurrences)

### Overall Expected Improvement:
- **Before**: ~70% mapping success
- **After**: ~85-90% mapping success

## Files Modified
- `enhanced_eukcensus_parser.py` - Main parser with all improvements
- `test_enhanced_taxonkit.py` - Test script for verification

## Usage
```bash
python enhanced_eukcensus_parser.py eukcensus_16S.clusters.97.tsv eukcensus_enhanced_v2
```

This will generate improved output files with significantly better taxonomic mapping coverage.
