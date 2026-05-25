# Subfigure Text — Mega-comprehensive scatter (novelty × over-representation)

*Working drafts of captions, key findings, methods notes, and discussion points for every panel of the mega-comprehensive scatter figure. Numbers are pulled from `final_merger/outputs/` (`.U.` catch-alls zeroed; top-8 labels per axis with ≥ 50 census observations).*

**Source script:** `visuals/scatter_plots/mega_comprehensive_stacked_visual.R`
**Outputs:** `visuals/scatter_plots/final_visualizations/comprehensive_mega_stacked_visual.{png,pdf}`
**Panel layout (2 rows × 3 cols):**

| Position | Domain | Rank | Source file |
|---|---|---|---|
| A | Bacteria | phylum (division) | `final_merger/outputs/16s_ncbi_merged_division.csv` |
| B | Bacteria | family | `final_merger/outputs/16s_ncbi_merged_family.csv` |
| C | Archaea | phylum (division) | `final_merger/outputs/16s_ncbi_merged_division.csv` |
| D | Archaea | family | `final_merger/outputs/16s_ncbi_merged_family.csv` |
| E | Eukaryota | division | `final_merger/outputs/18s_ncbi_merged_division.csv` |
| F | Eukaryota | family | `final_merger/outputs/18s_ncbi_merged_family.csv` |

## 1.1 Compact caption (≈ 50 words, journal-style)

> Cross-database comparison of taxonomic representation between NCBI reference genomes/species and culture-independent surveys (IMG/M sequences and OTUs). Each point is a taxon; x-axis is the over-representation factor (census ÷ NCBI relative abundance), y-axis is the novelty factor (NCBI genomes ÷ NCBI species). Quadrant **I** highlights well-sampled novel clades, **II** under-sampled lineages disproportionately recovered by primers.

## 1.2 Extended caption (≈ 200 words)

> Each marker represents one taxon at the indicated taxonomic rank, with the marker's horizontal position quantifying the *over-representation factor* — the ratio of its relative abundance in the IMG/M sequence census to its relative abundance in NCBI's reference genomes — and the marker's vertical position quantifying the *novelty factor* — the ratio of NCBI genome counts to NCBI species counts (high values flag clades whose reference genomes derive from a small number of species, i.e. recently described or genome-only lineages). Markers are colored by domain (Bacteria, Archaea, Eukaryota) using the project's shared palette (`visuals/shared_config/taxonomic_color_mapping.yaml`). Taxa flagged as `.U.` (unclassified) catch-alls inflate parent ranks by definition and have been zeroed across the merged tables to prevent label collisions in the upper-right quadrant. Labeled markers represent the top eight taxa per panel ranked by novelty factor *or* by over-representation factor, restricted to taxa with ≥ 50 census observations to suppress small-count artifacts. All other taxa are plotted but unlabeled. Axis scales are panel-specific because novelty ratios vary across two orders of magnitude between bacterial and eukaryotic ranks; the dashed line at *y* = 1 marks parity between NCBI genome and species counts.

## 1.3 Panel A — Bacteria, phylum/division

**Key findings to highlight in text:**
- **Abditibacteriota** (novelty 36.7; 50 genomes from 7 species) — top novel bacterial phylum, despite low absolute census representation.
- **Calditrichota** (novelty 28.4; 372 / 16) and **Thermomicrobiota** (25.6; 314 / 20) — extremophile lineages where reference genomes outpace described species ≈ 25-fold.
- **Planctomycetota** (24.4; 11,876 / 538) and **Acidobacteriota** (24.2; 14,863 / 377) — both substantial in absolute terms; high novelty driven by genome-resolved metagenomics initiatives.
- **Over-represented in census surveys**: *Chlorobiota* (19.3×), *Ignavibacteriota* (10.9×), *Chlamydiota* (9.7×) — primer-amplifiable lineages whose reference catalogue lags behind environmental detection.
- **Cyanobacteriota** is mildly over-represented (1.64×) and has the lowest novelty among major phyla (0.61), reflecting decades of targeted culturing.

**Discussion / framing points:**
- The upper-right quadrant (high novelty, high over-rep) is essentially empty in Bacteria phyla — these dimensions are anti-correlated at this rank.
- *Acidobacteriota* and *Planctomycetota* are the canonical "rare biosphere → MAG-rich" success stories visible in this panel.
- Pelagic/marine clades (*Chlorobiota*, *Pelagibacterales* at family level) dominate the over-representation axis.

## 1.4 Panel B — Bacteria, family

**Key findings:**
- **Thermoanaerobaculaceae** (177.0 novelty; 52/4), **Blastocatellaceae** (143.3; 65/3), **Kofleriaceae** (124.8; 204/6) — all *Acidobacteriota* or *Myxococcota* radiations now MAG-rich relative to type-strain descriptions.
- **Pyrinomonadaceae** (109.8; 814 / 6 species) — extreme novelty; ~800 MAGs from only 6 named species.
- **Pelagibacteraceae** (over-rep 36.7×; 4,339/220) — the canonical "SAR11" pattern: massive in surveys, comparatively under-represented in NCBI relative abundance terms.
- **Hydrogenimonaceae**, **Bartonellaceae**, **Breoghaniaceae** all 10–20× over-represented.

**Discussion:**
- The novelty rank is dominated by Acidobacteriota subclades — useful for the gap-analysis story.
- Pelagibacteraceae is the single most over-represented family — call out as the marine-microbiology example.

## 1.5 Panel C — Archaea, phylum/division

**Key findings:**
- **Nanoarchaeota** (102.5 novelty; 849/31) — extreme genome accumulation, mainly DPANN MAGs.
- **Odinarchaeia** (56.3; 10/3) and **Heimdallarchaeia** (20.8; 44/4) — Asgard radiations; low absolute counts but very high genome:species ratio.
- **Nitrososphaerota** (5.4; 4,238/270) and **Lokiarchaeia** (3.9; 1,624/92) — ammonia-oxidizing and Asgard lineages respectively.
- *Euryarchaeota* is the only archaeal phylum with novelty ≈ 2 (most "balanced" between genomes and species).
- **Top over-rep (Archaea)**: *Thermoproteota* (0.51) and *Euryarchaeota* (0.51) — note the entire archaeal domain is *under-represented* in surveys (no factor > 1).

**Discussion:**
- The archaeal panel demonstrates that 16S primer biases under-sample archaea relative to genome-resolved approaches; every phylum sits left of *x* = 1.
- DPANN (*Nanoarchaeota*) is the most striking case for "MAGs without isolates."

## 1.6 Panel D — Archaea, family

**Key findings:**
- **Fervidicoccaceae** (9.0; 68 / 3), **Nitrosopumilaceae** (9.0; 1,035 / 59), **Nitrososphaeraceae** (8.8; 897 / 31) — all AOA-related families with strong MAG accumulation.
- **Methanotrichaceae** (6.6; 400 / 37), **Methanospirillaceae** (6.4; 66 / 9), **Thermoplasmataceae** (6.2; 64 / 6) — methanogens with moderate novelty.
- **Over-rep top**: *Picrophilaceae* (2.0×), *Thermococcaceae* (1.7×), *Methanopyraceae* (1.7×), *Methanomassiliicoccaceae* (1.6×) — only marginal over-representation; archaeal families never reach the extreme over-rep values seen in bacteria.

**Discussion:**
- Family-level archaea show a much narrower dynamic range than bacteria on both axes — readers should be told this is real and not a clipping artifact.
- Highlight the AOA cluster (Nitros\*) as a single biological story.

## 1.7 Panel E — Eukaryota, division (phylum)

**Key findings:**
- **Tubulinea** (novelty 1,323; 2 / 1) — extreme outlier driven by tiny species count; flag as a long-tail artifact and discuss in caption.
- **Rhizaria** (284.8; 72 / 18), **Evosea** (69.5; 56 / 35), **Discoba** (41.6; 354 / 132), **Alveolata** (41.3; 658 / 219) — substantial protistan novelty.
- **Over-rep top**: *Rhodophyta* (26.5×), *Streptophyta* (1.35×), *Opisthokonta* (0.82×) — *Rhodophyta* is a striking over-rep case (78 genomes / 53 species, but disproportionately recovered by 18S primers).
- *Opisthokonta* (animals + fungi + sister groups) dominates absolute counts (48,096 genomes / 19,894 species) but sits near parity on both axes.

**Discussion:**
- The eukaryotic panel is the most informative for cultivation gaps — protistan supergroups (Rhizaria, Discoba, Alveolata) all show very high novelty *and* under-representation, the textbook "dark matter" pattern in protistology.
- Tubulinea should probably be either annotated as "n=1 species (long-tail artifact)" or moved to supplement (see §1.10).

## 1.8 Panel F — Eukaryota, family

**Key findings:**
- **Mastigamoebidae** (1,609 novelty; 2 / 1), **Euglenaceae** (470; 3 / 1), **Neovahlkampfiidae** (211; 1 / 1), **Digenea** (190; 1 / 1) — extreme single-species novelty outliers; small-count artifacts.
- **Oxytrichidae** (160; 9 / 5), **Hexamitidae** (143; 41 / 4) — biologically meaningful novel ciliate / diplomonad clades.
- **Over-rep top**: *Streptophyta_XXXX* (760×), *Arthropoda_XX* (344×), *Fungi_XXX* (37×) — note the `_XXXX` / `_XX` suffix indicates Pr2 rank-padding; these are residual "supergroup catch-all" labels that bypassed the `.U.` filter and should likely be moved to supplement (see §1.10).

**Discussion:**
- The over-representation axis here is dominated by Pr2 rank-padding artifacts (`_XXXX`, `_XX`); biologically interpretable over-rep is led by *Lepidosauria* (79×) and *Amphibia* (41×) — animals captured by 18S despite low NCBI genome representation.
- Novelty ranking is dominated by single-species families; the biologically interesting bin is novelty 50–200.

## 1.9 Methods text (≈ 80 words)

> For each taxon at the indicated rank, the novelty factor was computed as (NCBI genome count) ÷ (NCBI species count); the over-representation factor as (census relative abundance) ÷ (NCBI genome relative abundance) within domain. Taxa whose names ended in `.U.<rank>` (unclassified parent placeholders) were zeroed across all count columns to avoid double-counting upstream ranks. Markers were labelled when ranked in the top eight by either axis within each panel; all other taxa were plotted but unlabelled. Panels A–B and C–D share x-axes within domain.

## 1.10 Suggested cuts → supplement

**Strong candidates to move out of the main figure:**

| Item | Where it currently appears | Suggested destination | Why |
|---|---|---|---|
| *Tubulinea* (Euk phylum, novelty 1,323 from 2/1) | Panel E label | Supplementary table SX | Single-species artifact dominates y-axis |
| *Mastigamoebidae* / *Euglenaceae* / *Neovahlkampfiidae* / *Digenea* / *Oxytrichidae* (≤ 9 genomes each) | Panel F labels | Supplementary table SX | Small-count outliers; biologically minor |
| *Streptophyta_XXXX*, *Arthropoda_XX*, *Fungi_XXX*, *Stramenopiles_XXXX* | Panel F over-rep labels | Supplementary table + methods caveat | Pr2 rank-padding artifacts, not real families |
| Unlabeled markers below census-count = 50 | All panels (currently plotted in faint color) | Supplementary scatter | Reduces ink, preserves story |
| Numerical novelty + over-rep values | Currently inline beside each label | Supplementary table SX | Frees label real-estate for taxon name only |

**Recommended cutoff statement for methods:**
> *"Taxa with fewer than 50 census observations or with `.U.`/`_XXXX` placeholder lineage strings were excluded from the labelled set in Figures 1A–F and reported individually in Supplementary Table SX."*
