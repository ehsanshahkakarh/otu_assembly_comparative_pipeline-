# Archive: CSV Outputs Before Integration

This directory contains the original CSV output files from the 18S census parser **before** the systematic resolution integration was applied.

## Files

### `eukcensus_18S_by_division.csv`
- **Date**: January 30, 2026
- **Entries**: 22 divisions
- **Status**: Baseline before integration

### `eukcensus_18S_by_family.csv`
- **Date**: January 30, 2026
- **Entries**: 314 families
- **Unmapped**: 64 families with `NA` taxid (20.4%)
- **Status**: Baseline before integration

### `eukcensus_18S_by_genus.csv`
- **Date**: January 30, 2026
- **Entries**: 491 genera
- **Status**: Baseline before integration

## Purpose

These files are archived to enable comparison with the new outputs after systematic resolution integration:

1. **Verify resolution success**: Compare unmapped counts (should decrease from 64 to 56)
2. **Validate lineage completeness**: Check that 8 families now have complete lineages
3. **Quality control**: Ensure no data was lost during integration
4. **Historical reference**: Preserve original parser output

## Expected Changes After Integration

### Families Expected to Gain Lineages (8 total):
1. Maxillopoda → Crustacea lineage
2. Vermamoebidae → Tubulinea lineage
3. Neobodonidae → Kinetoplastea lineage
4. Tholoniidae → Ciliophora lineage
5. MAST-12 → Stramenopiles lineage
6. MAST-3 → Stramenopiles lineage
7. Ophryoglenida → Oligohymenophorea lineage
8. Haliphthorales → Peronosporomycetes lineage

### Metrics:
- **Before**: 64 unmapped families (20.4%)
- **After**: 56 unmapped families (17.8%)
- **Improvement**: 8 families resolved (12.5% of unmapped)

## Comparison Tool

To compare before/after results, run:

```bash
cd ../../py_18S
python compare_before_after_integration.py
```

This will generate a detailed comparison report showing:
- Number of families resolved
- Which families gained lineages
- Lineage completeness for each resolved family
- Differences in unmapped logs

## Retention

These files should be retained for:
- Quality assurance
- Troubleshooting integration issues
- Historical reference
- Rollback if needed

**Do not delete these files.**

