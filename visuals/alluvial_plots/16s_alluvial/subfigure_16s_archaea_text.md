# Subfigure Text — 16S Archaea stacked alluvial (absolute + percentage)

**Source script:** `visuals/alluvial_plots/16s_alluvial/scripts/alluvial_16s_stacked_pct.R` (same script as the bacterial alluvial; archaea subset)
**Outputs:** `visuals/alluvial_plots/16s_alluvial/figures/alluvial_16s_archaea_stacked_abs_pct.{png,pdf}`
**Axis order (left → right):** Genbank Genomes → IMG/M Sequences → 16S OTUs → Genbank Species
**Domain totals:** 18,661 reference genomes · 43,411 IMG/M sequences · 9,449 16S OTUs · 2,521 named species
**Stratum label threshold:** 3% of column total.

## 3.1 Compact caption (≈ 55 words)

> Flow of archaeal phylum-level composition across NCBI reference genomes, IMG/M sequences, 16S OTUs, and NCBI named species. Note the order-of-magnitude smaller absolute scale vs. bacteria; the surveys recover relatively more *Nanoarchaeota* and *Nitrososphaerota* than the named-species reference would predict, consistent with extensive recent MAG-driven expansion in these clades.

## 3.2 Extended caption (≈ 170 words)

> Archaeal phylum-level composition is shown across four reference frames (NCBI genomes, IMG/M sequences, 16S OTUs, NCBI species), with totals 18,661 / 43,411 / 9,449 / 2,521 respectively. Ribbon width is proportional to relative abundance per column; phyla are sorted top-to-bottom by per-column magnitude. Five phyla account for > 95% of each column: *Euryarchaeota*, *Nitrososphaerota*, *Thermoproteota*, *Lokiarchaeia*, and *Nanoarchaeota*. The absolute panel (top) emphasises that the IMG/M sequence pool is by far the deepest sampling of archaeal diversity (~ 2.3× the NCBI genome pool, ~ 17× the named-species count). The percentage panel (bottom) reveals the qualitative redistribution: *Nanoarchaeota* and *Lokiarchaeia* expand markedly in OTU-based surveys, while *Thermoproteota* contracts. Strata < 3% per column are drawn but unlabelled; the seventh phylum (*Heimdallarchaeia*) and the *Archaea.U.phylum* placeholder fall below this cutoff throughout and are tabulated in Supplementary Table SX.

## 3.3 Key findings (current numbers, threshold-3%)

- **Euryarchaeota dominates references but contracts in surveys.** 46.2% genomes → 43.4% IMG → 38.1% OTU → **72.6%** named species — the species:genome ratio is heavily biased toward Euryarchaeota historically (1,831 species).
- **Nanoarchaeota explodes in surveys.** 4.5% genomes → 20.9% IMG → **33.6%** OTU → < 1% species — the most extreme rare-biosphere → primer-amplified pattern in the dataset; 102.5× novelty factor confirms heavy MAG accumulation against essentially zero named species.
- **Nitrososphaerota** is consistently 15–25% across genome/IMG/OTU axes (22.7 / 24.6 / 15.4) but only 10.7% of named species — AOA expansion exceeds species formalisation.
- **Thermoproteota** contracts from 17.6% genomes → 5.9% IMG → 6.0% OTU but rebounds to 11.5% species — historically named hyperthermophile contribution.
- **Lokiarchaeia** stays consistently ~ 3–9% across all four axes — the most "balanced" Asgard signal.

## 3.4 Labels currently shown vs. hidden (3% threshold)

| Column | Visible (≥ 3%) | Hidden (< 3%) |
|---|---|---|
| Genbank Genome | Euryarchaeota (46.2), Nitrososphaerota (22.7), Thermoproteota (17.6), Lokiarchaeia (8.7), Nanoarchaeota (4.5) | Heimdallarchaeia (0.24), Odinarchaeia (0.05), `Archaea.U.phylum` (0) |
| IMG Genome | Euryarchaeota (43.4), Nitrososphaerota (24.6), Nanoarchaeota (20.9), Thermoproteota (5.9), Lokiarchaeia (3.1) | Odinarchaeia (1.2), Heimdallarchaeia (0.8), `Archaea.U.phylum` (0.2) |
| 16S OTU | Euryarchaeota (38.1), Nanoarchaeota (33.6), Nitrososphaerota (15.4), Thermoproteota (6.0), Lokiarchaeia (3.8) | Heimdallarchaeia (0.9), `Archaea.U.phylum` (0.3) |
| Genbank Species | Euryarchaeota (72.6), Thermoproteota (11.5), Nitrososphaerota (10.7), Lokiarchaeia (3.6) | Nanoarchaeota (1.2), Heimdallarchaeia (0.16), Odinarchaeia (0.12) |

## 3.5 Discussion talking points

- The *Nanoarchaeota* trajectory (4.5 → 33.6%) is the single most quantitatively striking shift in the entire taxonomic figure set and should be called out by name in the abstract / discussion.
- The fact that Asgard archaea (*Lokiarchaeia*, *Odinarchaeia*, *Heimdallarchaeia*) all sit below 5% across every axis is itself a finding — these are not yet quantitatively major contributors to environmental archaeal diversity as captured here.
- The species column is dominated by *Euryarchaeota* (73%) — a useful reminder of nomenclatural lag, since the genome catalogue already shows *Euryarchaeota* dropping to 46%.

## 3.6 Suggested cuts → supplement

| Item | Suggestion |
|---|---|
| Heimdallarchaeia, Odinarchaeia, `Archaea.U.phylum` (always < 1% per column) | Pool into "Other" for the main figure; tabulate per-column percentages in Supplementary Table SX. |
| Absolute-values panel | Same suggestion as the bacterial alluvial: move to supplement; place totals in axis label of percentage panel. |
| Combine with 16S Bacteria | The bacteria + archaea alluvials could become a single Figure 2 with panels (A) and (B); this would also justify a shared methods caption. |
