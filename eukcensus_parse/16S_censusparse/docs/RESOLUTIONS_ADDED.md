# 16S Census Parser - New Resolutions Added

**Date**: 2026-03-04  
**Purpose**: Resolve unmapped prokaryotic taxa (excluding Candidatus and organelles)

## Summary

Added **4 new resolutions** to `census_parser/known_parents.py` to resolve non-Candidatus, non-organellar unmapped taxa.

## Resolutions Added

### Family Level (1 resolution)

#### 1. Procabacteriaceae
- **Parent TaxID**: 809
- **Parent Name**: Chlamydiaceae
- **Rank**: family
- **OTU Count**: 8
- **Notes**: Family for Candidatus Procabacter, amoeba endosymbionts
- **Research**: Related to Chlamydiae, amoeba-associated bacteria

---

### Genus Level (3 resolutions)

#### 1. Coenonia
- **Parent TaxID**: 78328
- **Parent Name**: Allocoenonia
- **Rank**: genus
- **OTU Count**: 3
- **Notes**: Reclassified to Allocoenonia in 2024 (Oren & Molinari Novoa)
- **Research**: Taxonomic reclassification, Flavobacteriaceae family
- **Reference**: Oren A & Molinari Novoa EA (2024), doi: 10.1099/ijsem.0.006547

#### 2. Yangia
- **Parent TaxID**: 31989
- **Parent Name**: Paracoccaceae
- **Rank**: genus
- **OTU Count**: 3
- **Notes**: Marine Alphaproteobacteria, formerly Rhodobacteraceae
- **Research**: Rhodobacteraceae was reclassified to Paracoccaceae in NCBI

#### 3. Procabacter
- **Parent TaxID**: 809
- **Parent Name**: Chlamydiaceae
- **Rank**: genus
- **OTU Count**: 8
- **Notes**: Candidatus Procabacter, amoeba endosymbiont, Chlamydiae-related
- **Research**: Obligate endosymbiont of Acanthamoeba

---

## Impact

### Before
- **Total unmapped**: 216 entries (4.0%)
- **Resolvable non-Candidatus/organelle taxa**: 4 entries
- **Mapping success**: 96.0%

### After (Expected)
- **Total unmapped**: 212 entries (3.9%)
- **Newly resolved**: 4 entries (22 OTU occurrences total)
- **Mapping success**: 96.1%

### OTU Coverage
- **Procabacteriaceae**: 8 OTUs
- **Procabacter**: 8 OTUs
- **Coenonia**: 3 OTUs
- **Yangia**: 3 OTUs
- **Total**: 22 OTUs resolved

---

## Verification

All taxids verified with taxonkit:

```bash
echo "78328
31989
809" | taxonkit lineage -i 1 -r -n -L
```

Results:
- ✅ 78328 → Allocoenonia (genus, Flavobacteriaceae)
- ✅ 31989 → Paracoccaceae (family, Alphaproteobacteria)
- ✅ 809 → Chlamydiaceae (family, Chlamydiota)

---

## Next Steps

1. ✅ Added resolutions to `known_parents.py`
2. ⏳ Rerun 16S parser to apply resolutions
3. ⏳ Verify updated unmapped log
4. ⏳ Check CSV outputs for new lineages

---

## Files Modified

- `census_parser/known_parents.py` - Added 4 new entries (now 39 total)

---

## Notes

- **Candidatus taxa**: Intentionally left unmapped (acceptable)
- **Organellar sequences**: Intentionally left unmapped (acceptable)
- **Remaining unmapped**: Mostly environmental/uncultured clades or .U. entries

