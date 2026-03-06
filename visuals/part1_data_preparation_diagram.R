#!/usr/bin/env Rscript
# Part 1: Initial Data Preparation - IMG/M Census vs NCBI Assembly
# Created: 2026-01-09
# Purpose: Show how NCBI data was parsed to match IMG/M census format

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

cat("Creating Part 1: Data Preparation Diagram...\n")

# Color scheme from eukaryotic scatter plots
COLORS <- list(
  imgm_primary = "#55d0ba",      # Rhizaria - Teal (IMG/M census)
  ncbi_primary = "#003ce1",      # Evosea - Bright blue (NCBI)
  process = "#c73de4",           # Alveolata - Magenta (processing)
  output = "#663be6",            # Chlorophyta - Violet (output)
  text_dark = "#2c2c2c",
  arrow = "#333333"
)

create_part1_diagram <- function() {

  # Create a blank canvas (extended to accommodate bottom section)
  p <- ggplot() +
    xlim(0, 20) +
    ylim(-3, 24) +
    theme_void() +
    theme(plot.margin = margin(20, 20, 20, 20))
  
  # ============================================================================
  # TITLE
  # ============================================================================
  p <- p + 
    annotate("text", x = 10, y = 23, 
             label = "Part 1: Initial Data Preparation & Parsing Strategy",
             size = 7, fontface = "bold", color = COLORS$text_dark)
  
  # ============================================================================
  # LEFT SIDE: IMG/M CENSUS DATA (16S/18S)
  # ============================================================================
  
  # Main box - IMG/M Census
  p <- p +
    annotate("rect", xmin = 1, xmax = 8, ymin = 18, ymax = 21,
             fill = COLORS$imgm_primary, color = "white", linewidth = 2, alpha = 0.9) +
    annotate("text", x = 4.5, y = 20.2,
             label = "IMG/M Environmental Census", 
             size = 5, fontface = "bold", color = "white") +
    annotate("text", x = 4.5, y = 19.5,
             label = "16S rRNA (Prokaryotes)", 
             size = 4, color = "white") +
    annotate("text", x = 4.5, y = 18.8,
             label = "18S rRNA (Eukaryotes)", 
             size = 4, color = "white") +
    annotate("text", x = 4.5, y = 18.3,
             label = "Source: Dr. Miguel Romero", 
             size = 3, color = "white", fontface = "italic")
  
  # IMG/M Data characteristics
  p <- p +
    annotate("rect", xmin = 1, xmax = 8, ymin = 14.5, ymax = 17.5,
             fill = "#7de0d0", color = "white", linewidth = 1.5, alpha = 0.8) +
    annotate("text", x = 4.5, y = 17,
             label = "Data Structure:", 
             size = 4, fontface = "bold", color = COLORS$text_dark) +
    annotate("text", x = 4.5, y = 16.4, hjust = 0.5,
             label = "• OTU clusters (97% similarity)\n• Taxonomic assignments\n• Phylum/Family/Genus levels\n• OTU counts + sequence sizes\n• Environmental diversity baseline",
             size = 3.5, color = COLORS$text_dark, lineheight = 1.1)
  
  # Arrow down
  p <- p +
    annotate("segment", x = 4.5, xend = 4.5, y = 14.5, yend = 13,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.5)
  
  # IMG/M Output format
  p <- p +
    annotate("rect", xmin = 1, xmax = 8, ymin = 10.5, ymax = 12.5,
             fill = COLORS$output, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 4.5, y = 11.9,
             label = "Census Output Format", 
             size = 4, fontface = "bold", color = "white") +
    annotate("text", x = 4.5, y = 11.2, hjust = 0.5,
             label = "Taxon | OTU_count | Size_count\nPhylum/Family/Genus tables",
             size = 3.5, color = "white", lineheight = 1.1)
  
  # ============================================================================
  # RIGHT SIDE: NCBI ASSEMBLY DATA
  # ============================================================================
  
  # Main box - NCBI Assembly
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 18, ymax = 21,
             fill = COLORS$ncbi_primary, color = "white", linewidth = 2, alpha = 0.9) +
    annotate("text", x = 15.5, y = 20.2,
             label = "NCBI GenBank Assembly", 
             size = 5, fontface = "bold", color = "white") +
    annotate("text", x = 15.5, y = 19.5,
             label = "~3 Million Genomes", 
             size = 4, color = "white") +
    annotate("text", x = 15.5, y = 18.8,
             label = "All Domains of Life", 
             size = 4, color = "white") +
    annotate("text", x = 15.5, y = 18.3,
             label = "assembly_summary_genbank.txt", 
             size = 3, color = "white", fontface = "italic")
  
  # NCBI Data characteristics
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 14.5, ymax = 17.5,
             fill = "#3366ff", color = "white", linewidth = 1.5, alpha = 0.8) +
    annotate("text", x = 15.5, y = 17,
             label = "Data Structure:", 
             size = 4, fontface = "bold", color = "white") +
    annotate("text", x = 15.5, y = 16.4, hjust = 0.5,
             label = "• Strain-level genome assemblies\n• Assembly accessions\n• Species taxids (not aggregated)\n• Isolate vs MAG classification\n• Multiple strains per species",
             size = 3.5, color = "white", lineheight = 1.1)
  
  # Arrow down to parsing
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 14.5, yend = 13,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.5)

  # ============================================================================
  # NCBI PARSING PIPELINE (Multi-step process)
  # ============================================================================

  # Step 1: Taxonomic Mapping
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 11.5, ymax = 12.5,
             fill = COLORS$process, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 15.5, y = 12,
             label = "STEP 1: Taxonomic Mapping (taxonkit)",
             size = 3.5, fontface = "bold", color = "white")

  p <- p +
    annotate("text", x = 15.5, y = 11, hjust = 0.5,
             label = "taxid → Phylum/Family/Genus names",
             size = 3, color = COLORS$text_dark)

  # Arrow down
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 10.8, yend = 10.3,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.2)

  # Step 2: Species Aggregation
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 9, ymax = 10,
             fill = COLORS$process, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 15.5, y = 9.5,
             label = "STEP 2: Species-Level Aggregation",
             size = 3.5, fontface = "bold", color = "white")

  p <- p +
    annotate("text", x = 15.5, y = 8.5, hjust = 0.5,
             label = "Group by species_taxid\n(prevents strain inflation)",
             size = 3, color = COLORS$text_dark, lineheight = 1.1)

  # Arrow down
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 8.2, yend = 7.7,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.2)

  # Step 3: Taxonomic Rank Counting
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 6.5, ymax = 7.5,
             fill = COLORS$process, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 15.5, y = 7,
             label = "STEP 3: Hierarchical Counting",
             size = 3.5, fontface = "bold", color = "white")

  p <- p +
    annotate("text", x = 15.5, y = 6, hjust = 0.5,
             label = "Count genomes + species per taxon\nPhylum/Family/Genus levels",
             size = 3, color = COLORS$text_dark, lineheight = 1.1)

  # Arrow down
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 5.7, yend = 5.2,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.2)

  # Step 4: Isolate Classification
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 4, ymax = 5,
             fill = COLORS$process, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 15.5, y = 4.5,
             label = "STEP 4: Isolate vs MAG Classification",
             size = 3.5, fontface = "bold", color = "white")

  p <- p +
    annotate("text", x = 15.5, y = 3.5, hjust = 0.5,
             label = "Track cultured organisms\nvs environmental samples",
             size = 3, color = COLORS$text_dark, lineheight = 1.1)

  # Arrow down
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 3.2, yend = 2.7,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.2)

  # NCBI Output format (matching census format)
  p <- p +
    annotate("rect", xmin = 12, xmax = 19, ymin = 0.5, ymax = 2.5,
             fill = COLORS$output, color = "white", linewidth = 1.5, alpha = 0.9) +
    annotate("text", x = 15.5, y = 2.1,
             label = "NCBI Parsed Output Format",
             size = 4, fontface = "bold", color = "white") +
    annotate("text", x = 15.5, y = 1.4, hjust = 0.5,
             label = "Taxon | Genome_count | Species_count\nPhylum/Family/Genus tables\n✅ NOW COMPARABLE TO CENSUS",
             size = 3.5, color = "white", lineheight = 1.1)

  # ============================================================================
  # CENTER: THE CHALLENGE
  # ============================================================================

  p <- p +
    annotate("rect", xmin = 8.5, xmax = 11.5, ymin = 10, ymax = 12,
             fill = "#d24390", color = "white", linewidth = 2, alpha = 0.9) +
    annotate("text", x = 10, y = 11.5,
             label = "THE CHALLENGE",
             size = 4.5, fontface = "bold", color = "white") +
    annotate("text", x = 10, y = 10.5, hjust = 0.5,
             label = "Different data\nstructures require\ncareful parsing",
             size = 3.5, color = "white", lineheight = 1.1)

  # Arrows pointing to challenge
  p <- p +
    annotate("segment", x = 8, xend = 8.5, y = 11, yend = 11,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1, linetype = "dashed") +
    annotate("segment", x = 12, xend = 11.5, y = 11, yend = 11,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1, linetype = "dashed")

  # ============================================================================
  # BOTTOM: READY FOR MERGER (Next step preview)
  # ============================================================================

  # Arrow from IMG/M output to bottom
  p <- p +
    annotate("segment", x = 4.5, xend = 4.5, y = 10.5, yend = 3,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.5, linetype = "dashed")

  # Arrow from NCBI output to bottom
  p <- p +
    annotate("segment", x = 15.5, xend = 15.5, y = 0.5, yend = -0.5,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = COLORS$arrow, linewidth = 1.5, linetype = "dashed") +
    annotate("segment", x = 15.5, xend = 10, y = -0.5, yend = -0.5,
             color = COLORS$arrow, linewidth = 1.5, linetype = "dashed")

  # Arrow from IMG/M to merger point
  p <- p +
    annotate("segment", x = 4.5, xend = 10, y = 3, yend = -0.5,
             color = COLORS$arrow, linewidth = 1.5, linetype = "dashed")

  # Next step box
  p <- p +
    annotate("rect", xmin = 7, xmax = 13, ymin = -2.5, ymax = -0.8,
             fill = "#475093", color = "white", linewidth = 2, alpha = 0.9) +
    annotate("text", x = 10, y = -1.2,
             label = "READY FOR PART 2:",
             size = 4, fontface = "bold", color = "white") +
    annotate("text", x = 10, y = -1.9,
             label = "Triple-Anchor Merger Pipeline",
             size = 3.5, color = "white", fontface = "italic")

  # Key statistics boxes
  p <- p +
    annotate("rect", xmin = 0.5, xmax = 4, ymin = 0.5, ymax = 2.5,
             fill = "#7de0d0", color = "white", linewidth = 1.5, alpha = 0.8) +
    annotate("text", x = 2.25, y = 2.1,
             label = "Census Stats",
             size = 3.5, fontface = "bold", color = COLORS$text_dark) +
    annotate("text", x = 2.25, y = 1.4, hjust = 0.5,
             label = "287K OTUs (16S)\n4,578 genera\n95 phyla",
             size = 3, color = COLORS$text_dark, lineheight = 1.1)

  p <- p +
    annotate("rect", xmin = 16, xmax = 19.5, ymin = -2.5, ymax = -0.5,
             fill = "#3366ff", color = "white", linewidth = 1.5, alpha = 0.8) +
    annotate("text", x = 17.75, y = -0.9,
             label = "NCBI Stats",
             size = 3.5, fontface = "bold", color = "white") +
    annotate("text", x = 17.75, y = -1.6, hjust = 0.5,
             label = "~3M genomes\n142K species\n297 phyla",
             size = 3, color = "white", lineheight = 1.1)

  return(p)
}

# Generate the diagram
diagram <- create_part1_diagram()

# Save outputs
cat("Saving Part 1 diagram outputs...\n")

ggsave("part1_data_preparation.png", diagram, 
       width = 16, height = 12, dpi = 300, bg = "white")
cat("✅ Saved: part1_data_preparation.png (300 DPI)\n")

ggsave("part1_data_preparation_hires.png", diagram,
       width = 16, height = 12, dpi = 600, bg = "white")
cat("✅ Saved: part1_data_preparation_hires.png (600 DPI)\n")

ggsave("part1_data_preparation.pdf", diagram,
       width = 16, height = 12, bg = "white")
cat("✅ Saved: part1_data_preparation.pdf\n")

cat("\n🎉 Part 1 diagram generation complete!\n")

