# Subfigure Text Library — Index

*Per-figure subfigure text has been split into individual files co-located with each figure. This page is an index plus the shared cross-figure cutoff scheme and supplementary plan. Numbers were pulled from the current `final_merger/outputs/` tables (alluvial label threshold = 3% of column total; `.U.` catch-alls zeroed in the merged tables).*

## Index

| Figure | Subfigure text file |
|---|---|
| Mega-comprehensive scatter (novelty × over-representation) — Bacteria / Archaea / Eukaryota at phylum + family | [`scatter_plots/subfigure_scatter_text.md`](scatter_plots/subfigure_scatter_text.md) |
| 16S Bacteria stacked alluvial (absolute + percentage) | [`alluvial_plots/16s_alluvial/subfigure_16s_bacteria_text.md`](alluvial_plots/16s_alluvial/subfigure_16s_bacteria_text.md) |
| 16S Archaea stacked alluvial (absolute + percentage) | [`alluvial_plots/16s_alluvial/subfigure_16s_archaea_text.md`](alluvial_plots/16s_alluvial/subfigure_16s_archaea_text.md) |
| 18S Eukaryota stacked alluvial (absolute + percentage) | [`alluvial_plots/18s_alluvial/subfigure_18s_eukaryota_text.md`](alluvial_plots/18s_alluvial/subfigure_18s_eukaryota_text.md) |

---

## Cross-figure: unified cutoff scheme and supplementary layout

### A. Tiered visibility scheme (recommended for methods)

| Tier | Rule | Treatment in main figure | Destination if excluded |
|---|---|---|---|
| **1 — Featured** | Per-column ≥ 3% OR top-8 by novelty/over-rep with ≥ 50 census obs | Full color, labelled (taxon name, optional value) | – |
| **2 — Visible context** | 1% ≤ per-column < 3% OR census obs ≥ 10 | Plotted, unlabelled | Supplementary Table |
| **3 — Long-tail** | Per-column < 1% AND census obs < 10 | Pooled into "Other" | Supplementary Table only |
| **4 — Artifact** | `.U.<rank>` or `_XXXX` / `_XX` Pr2 padding | Excluded entirely from main figure | Methods caveat + Supplementary Table |

Current alluvial threshold is 3% (Tier 1 vs 2 boundary). Currently Tier 3 is **not** pooled into "Other" — strata below threshold are still drawn individually, just unlabelled. Pooling them would further reduce ribbon clutter at the cost of one extra preprocessing step in the alluvial scripts.

### B. Single-sentence methods statement (drop-in)

> *"Taxa contributing < 3% of any survey axis (alluvial figures) or fewer than 50 census observations and below the top-eight novelty/over-representation rank within a panel (scatter figure) were aggregated as 'Other' in main figures and reported individually in Supplementary Table SX. Placeholder lineage strings (`*.U.<rank>`, Pr2 `_XXXX` / `_XX` rank padding) were excluded from the labelled set; their per-column counts are reported in Supplementary Table SY."*

### C. Suggested supplementary figures / tables

| Tag | Description | Source |
|---|---|---|
| Suppl. Fig S1 | Absolute-value panels of the three alluvials stacked vertically | Already saved as the upper half of the existing PDFs — split out as supplement-only |
| Suppl. Fig S2 | Full long-tail alluvial: every phylum/division (no 3% threshold) | Re-run with `LABEL_MIN_PCT <- 0` |
| Suppl. Fig S3 | Family-level alluvial for top 3 phyla per domain | New script: same logic, family ranks |
| Suppl. Fig S4 | Faint-point scatter (unlabelled markers only) for each panel | Re-run scatter with point color = grey for non-top-8 |
| Suppl. Table SX | Per-taxon counts (genome, IMG, OTU, species), percentages, novelty factor, over-rep factor, label-visible flag | Auto-export from final_merger outputs |
| Suppl. Table SY | `.U.` and `_XXXX` placeholder strings: which were zeroed, which retained, per-column counts | New script tied to `final_merger` filters |

### D. Things to keep in the main figure (even though they look "messy")

- The full four-axis flow for the alluvials (genome → IMG → OTU → species). Dropping one axis collapses the narrative.
- The "Other" stratum visible in IMG / OTU columns (≈ 10%) — its size is itself a finding.
- The unlabelled background points in the scatter — they give the reader the sense of how heavy the long tail is. Removing them entirely risks making the figure look like there are only ~ 16 taxa.

### E. Open questions for the author

1. **Do you want the absolute panel in the main figure or in the supplement?** Removing it would let the percentage panel breathe (taller bars, larger labels, more white space).
2. **Pool sub-3% strata into "Other" visually?** Currently they are drawn as separate thin ribbons; pooling them would make ribbons cleaner but hides ordering information.
3. **Drop the *named species* column in the eukaryota alluvial?** It is near-identical to the genome column and could go to supplement.
4. **Combine the two 16S alluvials (Bacteria, Archaea) into a single Figure 2 with panels A/B?** Or keep separate for breathing room.
5. **For the scatter mega-panel: are six panels in 2×3 the final layout?** Splitting into Figure 1 (phylum-level, 3 panels) + Figure 1S/3 (family-level, 3 panels) is an option.

---

*Updated: alluvial threshold = 3.0% (set in `LABEL_MIN_PCT` at the top of both alluvial scripts). All numbers in the per-figure files were extracted from the current `final_merger/outputs/` tables and the live alluvial data frames.*
