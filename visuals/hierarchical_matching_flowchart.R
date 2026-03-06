#!/usr/bin/env Rscript
# Hierarchical Matching Algorithm Flowchart
# Created: 2026-01-05
# Purpose: Publication-quality flowchart using DiagrammeR
# Colors matched to eukaryotic scatter plot color scheme

# Install DiagrammeR if needed (uncomment if not installed)
# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")
# install.packages("rsvg")

suppressPackageStartupMessages({
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

cat("Creating hierarchical matching algorithm flowchart...\n")
cat("Using color scheme from eukaryotic scatter plots...\n")

# ============================================================================
# COLOR SCHEME MAPPING
# ============================================================================
# Colors sourced from: visuals/shared_config/taxonomic_color_mapping.yaml
# Section: eukaryota_colors (used in scatter plots)
#
# Flowchart Element          → Eukaryotic Phylum    → Hex Color
# -------------------------------------------------------------------------
# START (Input)              → Tubulinea            → #416b7d (Dark blue-gray)
# TIER 1 (PRIMARY, 85-90%)   → Rhizaria             → #55d0ba (Teal - bright!)
# TIER 2 (FALLBACK, 10-15%)  → Evosea               → #003ce1 (Bright blue)
# TIER 3 (VALIDATION, <5%)   → Alveolata            → #c73de4 (Magenta)
# SUCCESS (Matched)          → Chlorophyta          → #663be6 (Violet)
# FAIL (Census-Only)         → Discoba              → #65417a (Purple)
# AGGREGATION (Tier 4)       → Stramenopiles        → #475093 (Dark purple)
# METRICS (Coined)           → Discosea             → #d24390 (Hot pink)
# VALIDATION (GTDB-NCBI)     → Metamonada           → #cf8ac6 (Pink)
#
# This ensures visual consistency with the eukaryotic scatter plots
# ============================================================================

# Create the flowchart using Graphviz DOT language
flowchart <- grViz("
digraph hierarchical_matching {

  # Graph attributes
  graph [layout = dot,
         rankdir = TB,
         bgcolor = white,
         fontname = 'Arial',
         fontsize = 14,
         nodesep = 0.8,
         ranksep = 1.2]

  # Node defaults
  node [shape = box,
        style = 'rounded,filled',
        fontname = 'Arial',
        fontsize = 12,
        fontcolor = white,
        penwidth = 2.5]

  # Edge defaults
  edge [fontname = 'Arial',
        fontsize = 11,
        penwidth = 2,
        color = '#333333']

  # Define nodes with eukaryotic scatter plot colors
  START [label = 'Environmental Census Data\n(16S/18S OTUs with taxids)',
         fillcolor = '#416b7d',
         color = '#2c4a55',
         width = 4]

  T1 [label = '🥇 TIER 1: Taxid Matching\nPRIMARY METHOD\n85-90% match rate\ncensus_taxid == ncbi_taxid',
      fillcolor = '#55d0ba',
      color = '#3a9080',
      fontcolor = black,
      width = 4]

  T2 [label = '🥈 TIER 2: Lineage Matching\nFALLBACK METHOD\n10-15% match rate\nHandles misnaming',
      fillcolor = '#003ce1',
      color = '#002a9f',
      width = 4]

  T3 [label = '🥉 TIER 3: Direct Name\nVALIDATION METHOD\n<5% match rate\nEdge cases',
      fillcolor = '#c73de4',
      color = '#8f2ba0',
      width = 4]

  SUCCESS1 [label = '✅ Matched\nHighest Confidence',
            fillcolor = '#663be6',
            color = '#4a29a3',
            width = 2.5]

  SUCCESS2 [label = '✅ Matched\nHigh Confidence',
            fillcolor = '#663be6',
            color = '#4a29a3',
            width = 2.5]

  SUCCESS3 [label = '✅ Matched\nModerate Confidence',
            fillcolor = '#663be6',
            color = '#4a29a3',
            width = 2.5]

  FAIL [label = '❌ Census-Only\nNo genomic data\nnovelty_factor = ∞',
        fillcolor = '#65417a',
        color = '#472d54',
        width = 2.5]

  AGG [label = 'TIER 4: Hierarchical Aggregation\nSum all descendant genomes',
       fillcolor = '#475093',
       color = '#323867',
       width = 4]

  METRICS [label = 'Calculate COINED METRICS\n🔴 novelty_factor\n🔵 overrepresentation_factor',
           fillcolor = '#d24390',
           color = '#932f65',
           width = 4]

  VALIDATE [label = 'GTDB-NCBI Validation Pipeline\nProves complete overlap\nValidates algorithm effectiveness',
            fillcolor = '#cf8ac6',
            color = '#91608a',
            width = 3.5]
  
  # Define edges (arrows) - using eukaryotic colors
  START -> T1

  T1 -> SUCCESS1 [label = 'Match Found', color = '#663be6', penwidth = 2.5]
  T1 -> T2 [label = 'No Match']

  T2 -> SUCCESS2 [label = 'Match Found', color = '#663be6', penwidth = 2.5]
  T2 -> T3 [label = 'No Match']

  T3 -> SUCCESS3 [label = 'Match Found', color = '#663be6', penwidth = 2.5]
  T3 -> FAIL [label = 'No Match', color = '#65417a', penwidth = 2.5]

  SUCCESS1 -> AGG
  SUCCESS2 -> AGG
  SUCCESS3 -> AGG

  AGG -> METRICS

  # Validation arrows (dashed) - using Metamonada pink
  VALIDATE -> T1 [label = 'Perfected', style = dashed, color = '#cf8ac6', penwidth = 2.5]
  VALIDATE -> T2 [label = 'Perfected', style = dashed, color = '#cf8ac6', penwidth = 2.5]
  VALIDATE -> T3 [label = 'Perfected', style = dashed, color = '#cf8ac6', penwidth = 2.5]
  
  # Rank constraints to control layout
  {rank = same; T1; VALIDATE}
  {rank = same; SUCCESS1; SUCCESS2; SUCCESS3; FAIL}
}
")

# Display the flowchart
print(flowchart)

# Export as PNG
cat("Saving as PNG...\n")
flowchart %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_png("hierarchical_matching_flowchart.png", width = 3000, height = 4000)

cat("✅ Saved: hierarchical_matching_flowchart.png (3000x4000 px)\n")

# Export as PDF
cat("Saving as PDF...\n")
flowchart %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_pdf("hierarchical_matching_flowchart.pdf", width = 10, height = 13)

cat("✅ Saved: hierarchical_matching_flowchart.pdf\n")

# Export as SVG (vector format - best for papers)
cat("Saving as SVG...\n")
flowchart %>%
  export_svg() %>%
  writeLines("hierarchical_matching_flowchart.svg")

cat("✅ Saved: hierarchical_matching_flowchart.svg (vector format)\n")

cat("\n🎉 Flowchart generation complete!\n")
cat("Files saved in: visuals/\n")
cat("  - PNG: hierarchical_matching_flowchart.png (high resolution)\n")
cat("  - PDF: hierarchical_matching_flowchart.pdf (publication ready)\n")
cat("  - SVG: hierarchical_matching_flowchart.svg (vector, scalable)\n")

