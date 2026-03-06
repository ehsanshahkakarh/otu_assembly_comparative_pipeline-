#!/usr/bin/env Rscript
# Database Parsing and Structural Formation Pipeline Diagram
# Created: 2026-01-09
# Purpose: Publication-quality diagram showing data flow from raw databases to final analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

cat("Creating database parsing pipeline diagram...\n")

# Create the diagram using ggplot2
create_pipeline_diagram <- function() {
  
  # Define node positions (x, y coordinates)
  nodes <- data.frame(
    id = c(
      # Raw Data Sources (Top row)
      "16S_raw", "18S_raw", "NCBI_raw", "GTDB_raw", "EukProt_raw",
      
      # Parsing Layer (Second row)
      "16S_parse", "18S_parse", "NCBI_parse", "GTDB_parse", "EukProt_parse",
      
      # Intermediate Processing (Third row)
      "16S_counts", "18S_counts", "NCBI_counts", "GTDB_counts", "EukProt_counts",
      
      # Census Merge (Fourth row)
      "Census_merge",
      
      # Triple-Anchor Merger (Fifth row)
      "Triple_anchor",
      
      # Final Output (Bottom row)
      "Final_output"
    ),
    x = c(
      # Raw sources - evenly spaced
      1, 2, 3, 4, 5,
      # Parsing - same x positions
      1, 2, 3, 4, 5,
      # Counts - same x positions
      1, 2, 3, 4, 5,
      # Census merge - centered between 16S and 18S
      1.5,
      # Triple anchor - centered
      3,
      # Final output - centered
      3
    ),
    y = c(
      # Raw sources
      6, 6, 6, 6, 6,
      # Parsing
      5, 5, 5, 5, 5,
      # Counts
      4, 4, 4, 4, 4,
      # Census merge
      3,
      # Triple anchor
      2,
      # Final output
      1
    ),
    label = c(
      "16S Census\nOTU Clusters", "18S Census\nOTU Clusters", "NCBI\nAssembly", "GTDB\nMetadata", "EukProt\nProteomes",
      "16S Parser", "18S Parser", "NCBI Parser", "GTDB Parser", "EukProt Parser",
      "16S Counts\n(Phylum/Family/Genus)", "18S Counts\n(Division/Family/Genus)", 
      "NCBI Counts\n(Species-aggregated)", "GTDB Counts\n(Species-aggregated)", "EukProt Counts",
      "Census Merger\n(16S + 18S)",
      "Triple-Anchor Merger\n(Name + Accession + Lineage)",
      "Final Comparative Tables\n(Environmental vs Genomic)"
    ),
    color = c(
      # Raw sources - using eukaryotic colors
      "#416b7d", "#55d0ba", "#003ce1", "#c73de4", "#65417a",
      # Parsing - lighter versions
      "#5a8ba0", "#7de0d0", "#3366ff", "#d966f0", "#8a5a9a",
      # Counts - medium tones
      "#416b7d", "#55d0ba", "#003ce1", "#c73de4", "#65417a",
      # Census merge - special color
      "#cf8ac6",
      # Triple anchor - special color
      "#d24390",
      # Final output - success color
      "#663be6"
    ),
    stringsAsFactors = FALSE
  )
  
  # Define edges (connections between nodes)
  edges <- data.frame(
    from = c(
      # Raw to Parsing
      "16S_raw", "18S_raw", "NCBI_raw", "GTDB_raw", "EukProt_raw",
      # Parsing to Counts
      "16S_parse", "18S_parse", "NCBI_parse", "GTDB_parse", "EukProt_parse",
      # Counts to Census Merge
      "16S_counts", "18S_counts",
      # Counts and Census to Triple Anchor
      "Census_merge", "NCBI_counts", "GTDB_counts", "EukProt_counts",
      # Triple Anchor to Final
      "Triple_anchor"
    ),
    to = c(
      # Raw to Parsing
      "16S_parse", "18S_parse", "NCBI_parse", "GTDB_parse", "EukProt_parse",
      # Parsing to Counts
      "16S_counts", "18S_counts", "NCBI_counts", "GTDB_counts", "EukProt_counts",
      # Counts to Census Merge
      "Census_merge", "Census_merge",
      # Counts and Census to Triple Anchor
      "Triple_anchor", "Triple_anchor", "Triple_anchor", "Triple_anchor",
      # Triple Anchor to Final
      "Final_output"
    ),
    stringsAsFactors = FALSE
  )
  
  # Merge edge coordinates
  edges <- merge(edges, nodes[, c("id", "x", "y")], by.x = "from", by.y = "id")
  names(edges)[3:4] <- c("x_from", "y_from")
  edges <- merge(edges, nodes[, c("id", "x", "y")], by.x = "to", by.y = "id")
  names(edges)[5:6] <- c("x_to", "y_to")
  
  # Create the plot
  p <- ggplot() +
    # Draw edges (arrows)
    geom_segment(data = edges, 
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 color = "#333333", linewidth = 1.2, alpha = 0.7) +
    
    # Draw nodes (boxes)
    geom_tile(data = nodes, aes(x = x, y = y, fill = color),
              width = 0.8, height = 0.6, color = "white", linewidth = 2) +
    
    # Add labels
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              color = "white", size = 3.5, fontface = "bold", lineheight = 0.9) +
    
    # Styling
    scale_fill_identity() +
    theme_void() +
    coord_fixed(ratio = 1) +
    theme(plot.margin = margin(20, 20, 20, 20))
  
  return(p)
}

# Generate the diagram
diagram <- create_pipeline_diagram()

# Save outputs
cat("Saving diagram outputs...\n")

# PNG (standard resolution)
ggsave("database_pipeline_diagram.png", diagram, 
       width = 14, height = 10, dpi = 300, bg = "white")
cat("✅ Saved: database_pipeline_diagram.png (300 DPI)\n")

# PNG (high resolution for publication)
ggsave("database_pipeline_diagram_hires.png", diagram,
       width = 14, height = 10, dpi = 600, bg = "white")
cat("✅ Saved: database_pipeline_diagram_hires.png (600 DPI)\n")

# PDF (vector format)
ggsave("database_pipeline_diagram.pdf", diagram,
       width = 14, height = 10, bg = "white")
cat("✅ Saved: database_pipeline_diagram.pdf\n")

cat("\n🎉 Database pipeline diagram generation complete!\n")
cat("Files saved in: visuals/\n")

