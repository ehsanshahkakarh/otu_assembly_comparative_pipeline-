# Subfigure Text — 18S Eukaryota stacked alluvial (absolute + percentage)

**Source script:** `visuals/alluvial_plots/18s_alluvial/scripts/alluvial_18s_stacked_pct.R`
**Outputs:** `visuals/alluvial_plots/18s_alluvial/figures/alluvial_18s_stacked_abs_pct.{png,pdf}`
**Axis order:** NCBI Eukaryota Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Eukaryota Species
**Domain totals:** 58,826 NCBI eukaryotic genomes · 400,105 18S sequences · 70,560 18S OTUs · 24,172 named eukaryotic species
**Stratum label threshold:** 3% of column total.

## 4.1 Compact caption (≈ 60 words)

> Flow of eukaryotic division-level composition across NCBI reference genomes, 18S EukCensus sequences, 18S OTUs, and NCBI named species. Ribbon width is proportional to per-column relative abundance; divisions are sorted by per-column magnitude. The NCBI reference is overwhelmingly *Opisthokonta* (animals + fungi); the 18S survey lenses recover a much broader protistan diversity.

## 4.2 Extended caption (≈ 200 words)

> Eukaryotic division-level composition is shown across four reference frames. The leftmost axis is the relative abundance of each division among NCBI eukaryotic reference genomes (n = 58,826); the next axis is among 18S EukCensus sequences (n = 400,105); the third axis is among 18S OTUs (n = 70,560); the rightmost axis is among NCBI named eukaryotic species (n = 24,172). The strong dominance of *Opisthokonta* (animals + fungi + close relatives) in the NCBI reference catalogue (82%) collapses to ~ 35% in the OTU survey — a quantitative measure of the well-known protistan "dark matter" gap. The IMG/EukCensus column reveals substantial *Alveolata*, *Stramenopiles*, *Rhizaria*, *Discoba*, and *Amoebozoa* contributions that the NCBI genome catalogue largely misses (each < 2% in the genome column). Pr2 lineage strings that resolve only to the *.U.division* level (`Eukaryota.U.division`, `Amoebozoa.U.division`) appear as separate strata in the survey columns because the rank-mapping cannot place them deeper; they are kept visible to faithfully represent the survey but should be interpreted as binned residuals rather than discrete divisions.

## 4.3 Key findings (current numbers, threshold-3%)

- **Opisthokonta drops from 82% (genomes) to 35% (OTUs) and back to 82% (species).** The symmetric genome/species values reflect the animal- and fungal-centric history of taxonomy; the survey-axis collapse to ~ 35% is the eukaryotic equivalent of the bacterial *Pseudomonadota* story.
- **Streptophyta (plants) holds 13.3% genomes, 12.6% species, but only 3.2% OTUs.** Confirms 18S primer under-recovery of plants (expected; standard 18S primers underrepresent Embryophyta).
- **Alveolata expands in surveys.** < 2% genomes → 15.1% IMG → 12.8% OTU → < 1% species. Protistan dark-matter case #1.
- **Stramenopiles**, **Rhizaria**, **Discoba** each move from < 2% in the genome / species columns to 5–7% in the OTU column.
- **`Eukaryota.U.division`** is 9.8% of IMG and 13.2% of OTU — a substantial "unresolvable at division level" bin that motivates the family-level follow-up.

## 4.4 Labels currently shown vs. hidden (3% threshold)

| Column | Visible (≥ 3%) | Hidden (< 3%) |
|---|---|---|
| Genbank Genome | Opisthokonta (81.8), Streptophyta (13.3) | Stramenopiles (1.8), Other (1.3), Alveolata (1.1), Discoba (0.6), Rhizaria (0.1), Eukaryota.U.div (0), Amoebozoa.U.div (0) |
| IMG Genome | Opisthokonta (39.7), Alveolata (15.1), Other (10.1), Eukaryota.U.div (9.8), Stramenopiles (6.2), Rhizaria (6.0), Discoba (5.3), Streptophyta (4.0), Amoebozoa.U.div (3.8) | (none) |
| 18S OTU | Opisthokonta (34.5), Eukaryota.U.div (13.2), Alveolata (12.8), Other (10.0), Discoba (7.8), Rhizaria (7.3), Amoebozoa.U.div (6.0), Stramenopiles (5.2), Streptophyta (3.2) | (none) |
| Genbank Species | Opisthokonta (82.3), Streptophyta (12.6) | Stramenopiles (2.0), Other (1.6), Alveolata (0.9), Discoba (0.5), Rhizaria (0.07), .U. placeholders (0) |

## 4.5 Discussion talking points

- The Opisthokonta 82 → 35 → 82% pattern is the single most quotable eukaryotic finding. Use it to motivate why genome-based diversity estimates substantially under-represent protistan lineages.
- *Streptophyta* in surveys (3.2% OTUs vs 13.3% genomes) is a known 18S primer artifact — flag in discussion to avoid mis-interpretation as a biological gap.
- The presence of two `.U.division` strata at 9–13% in the survey columns is honest reporting of taxonomic uncertainty and motivates either (a) deeper Pr2 lineage refinement or (b) reporting at supergroup level.
- The combination of high *Alveolata*/*Discoba*/*Rhizaria* in surveys with extreme novelty factors (Panel E of the scatter figure) gives the figure pair a self-reinforcing narrative.

## 4.6 Suggested cuts → supplement

| Item | Suggestion |
|---|---|
| `Eukaryota.U.division` and `Amoebozoa.U.division` strata | Either pool into "Other (unresolved)" in main figure or annotate explicitly. Discussion should acknowledge them. |
| Genbank Species column being a near-mirror of Genbank Genome column | Optional: drop the named-species axis if redundant. Saves horizontal space. |
| Absolute-values panel | Same as the 16S alluvials: move to supplement; report totals in axis label. |
| Streptophyta over-recovery in genomes vs OTUs | Annotate in caption ("18S primer bias against Embryophyta") rather than letting reader misinterpret as biology. |
