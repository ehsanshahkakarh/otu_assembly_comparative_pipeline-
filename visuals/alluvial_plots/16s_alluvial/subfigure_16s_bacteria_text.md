# Subfigure Text — 16S Bacteria stacked alluvial (absolute + percentage)

**Source script:** `visuals/alluvial_plots/16s_alluvial/scripts/alluvial_16s_stacked_pct.R`
**Outputs:** `visuals/alluvial_plots/16s_alluvial/figures/alluvial_16s_bacteria_stacked_abs_pct.{png,pdf}`
**Axis order (left → right):** Genbank Genomes → IMG/M Sequences → 16S OTUs → Genbank Species
**Domain totals:** 2,880,467 reference genomes · 1,093,021 IMG/M sequences · 239,785 16S OTUs · 96,304 named species
**Stratum label threshold:** 3% of column total (lower strata visible but unlabelled).

## 2.1 Compact caption (≈ 60 words)

> Flow of bacterial phylum-level composition across four reference frames: NCBI genome counts, IMG/M survey sequences (proxy for environmental abundance), 16S rRNA OTUs (independent amplicon-based diversity estimator), and NCBI named species. Ribbon width is proportional to relative abundance within each column; phyla are independently sorted (largest on top) per column to make magnitude comparisons direct.

## 2.2 Extended caption (≈ 180 words)

> Bacterial phylum-level composition is shown across four complementary lenses on microbial diversity. The leftmost axis is the relative abundance of each phylum among NCBI bacterial reference genomes (n = 2.88 M); the next axis is its relative abundance among IMG/M-deposited 16S/rRNA gene sequences (n = 1.09 M); the third axis is its relative abundance among 16S OTUs (n = 239,785); the rightmost axis is its relative abundance among named NCBI species (n = 96,304). Ribbons connecting axes track the same phylum across the four lenses and are proportional to the local column percentage. Phyla are independently re-ordered top-to-bottom by per-column magnitude (largest on top), so the visual order of strata always matches the absolute size of the underlying bar. Strata contributing < 3% of a given column are drawn but unlabelled to reduce crowding; values are tabulated in Supplementary Table SX. The top panel shows absolute counts on a domain-total scale; the bottom panel shows the same data normalised to 100% per column. A single shared legend identifies the phyla.

## 2.3 Key findings (current numbers, threshold-3%)

- **Pseudomonadota dominance collapses across surveys.** 57.7% of NCBI genomes (1.66 M) → 28.2% of IMG/M sequences (308 k) → 22.4% of OTUs (54 k) → 39.2% of named species (37.7 k). A ≈ 35-percentage-point drop from the genome reference to the OTU survey is the single most important story this figure tells.
- **Bacillota mirrors the pattern.** 23.5% → 15.8% → 20.8% → 16.8%. Persistent across all four lenses.
- **Bacteroidota grows in surveys.** 4.5% (genomes) → 12.2% (IMG) → 13.0% (OTUs) → 8.0% (species) — the canonical "well-amplified, under-sequenced" pattern.
- **Actinomycetota rebounds at species level.** 3.8% genomes → 10.0% IMG → 6.6% OTU → **20.8%** species — the species:genome inflation flips around (heavy historical culture base, less MAG-rich).
- **Chloroflexota / Planctomycetota / Verrucomicrobiota / Acidobacteriota** all become 3–5% strata in IMG/OTU views but are invisible (< 1%) in NCBI genome view — the "rare biosphere becomes visible" story.
- **"Other" balloons to 10–11% in IMG and OTU columns** vs ~ 1.6% in the genome column — a quantitative measure of long-tail diversity inaccessible to current reference catalogues.

## 2.4 Labels currently shown vs. hidden (3% threshold)

| Column | Visible (≥ 3%) | Hidden (< 3%, plotted) |
|---|---|---|
| Genbank Genome | Pseudomonadota (57.7), Bacillota (23.5), Campylobacterota (5.8), Bacteroidota (4.5), Actinomycetota (3.8) | Other, Chloroflexota, Acidobacteriota, Verrucomicrobiota, Planctomycetota, Patescibacteria, Thermodesulfobacteriota, Cyanobacteriota |
| IMG Genome | Pseudomon. (28.2), Bacillota (15.8), Bacteroidota (12.2), Actinomyc. (10.0), Other (10.0), Chloroflexota (5.4), Verrucomic. (4.2), Planctomyc. (4.1), Acidobact. (4.1) | Thermodesulfo. (2.4), Patescibact. (1.6), Cyano. (1.2), Campylobact. (0.9) |
| 16S OTU | Pseudomon. (22.4), Bacillota (20.8), Bacteroidota (13.0), Other (11.2), Actinomyc. (6.6), Chloroflex. (5.5), Planctomyc. (5.5), Verrucomic. (4.3), Acidobact. (3.8) | Thermodesulfo. (2.6), Patescibact. (2.3), Cyano. (1.0), Campylobact. (1.0) |
| Genbank Species | Pseudomon. (39.2), Actinomyc. (20.8), Bacillota (16.8), Bacteroidota (8.0), Cyano. (4.2), Other (4.0) | Patescibact. (3.0), Thermodesulfo. (1.0), Campylobact. (0.8), Chloroflex. (0.8), Verrucomic. (0.6), Planctomyc. (0.6), Acidobact. (0.4) |

## 2.5 Methods note (≈ 60 words)

> Stratum labels report taxon name only; per-column totals appear in the axis subtitle. Strata < 3% of their column total are drawn but not labelled (see Supplementary Table SX). The relative width of ribbons across axes is purely descriptive — no statistical comparison is implied — and stratum stacking is independently sorted by per-column magnitude so that visual rank matches numerical rank within each axis.

## 2.6 Discussion talking points

- The Pseudomonadota story (58 → 22 → 39%) is the cleanest illustration in the manuscript of *isolation bias × primer bias × species concept* interacting.
- Bacteroidota and Verrucomicrobiota are the strongest "surveys see them, references don't" cases at phylum level.
- "Other" reaching 10% in IMG/OTU columns motivates the family/genus-level follow-up in Figures 1B and S(family-alluvial).
- The four axes are *not* independent measurements of the same thing — frame them as four different observational lenses (genome catalogue, community deposit, primer-based survey, taxonomic literature). This avoids over-claiming.

## 2.7 Suggested cuts → supplement

| Item | Suggestion |
|---|---|
| Absolute-values top panel | Could move to Supplementary Figure SX; main figure keeps percentage panel + column totals in axis labels. |
| Strata < 3% (currently shown unlabelled) | Already invisible in label sense; consider also pooling visually into "Other" bin for cleaner ribbons. |
| Per-stratum numeric values (currently in label) | Already removed at 3% threshold; for shown strata, consider stripping the numeric line entirely from main figure and giving full table in supplement. |
| Side-by-side bacteria + archaea panels | Currently saved as separate figures; could combine as Figure 2A (bacteria) + 2B (archaea) using patchwork. |
