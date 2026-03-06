#!/usr/bin/env Rscript
# Hierarchical Matching Algorithm Flowchart
# Created: 2026-01-05
# Purpose: Publication-quality diagram showing the matching algorithm hierarchy
# Matches aesthetic of alluvial and scatter plots

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(dplyr)
})

# Configuration
CONFIG <- list(
  # Output settings
  width = 16,
  height = 20,
  dpi = 300,
  
  # Colors matching your existing plots
  colors = list(
    tier1 = "#32CD32",        # Lime green (like Chloroflexota) - PRIMARY
    tier2 = "#ffb44c",        # Light orange (like Acidobacteriota) - FALLBACK
    tier3 = "#ff7200",        # Orange (like Planctomycetota) - VALIDATION
    tier4 = "#46bda3",        # Aqua (like Verrucomicrobiota) - AGGREGATION
    success = "#4c9b34",      # Dark green - matched
    fail = "#d19386",         # Light brown - census-only
    validation = "#c73de4",   # Magenta - GTDB-NCBI validation
    metrics = "#f51b7f",      # Bright pink - COINED METRICS
    start = "#416b7d",        # Dark blue-gray - input
    arrow = "#333333",        # Dark gray for arrows
    text = "#000000",         # Black for text
    bg = "#FFFFFF"            # White background
  ),
  
  # Text sizes
  text = list(
    title = 16,
    tier = 14,
    label = 12,
    detail = 10,
    small = 9
  )
)

# Create the diagram using ggplot2
create_matching_diagram <- function() {
  
  # Define node positions (x, y coordinates)
  nodes <- data.frame(
    id = c("start", "t1", "t2", "t3", "agg", "s1", "s2", "s3", "fail", "metrics", "validate"),
    x = c(5, 5, 5, 5, 5, 8.5, 8.5, 8.5, 1.5, 5, 9.5),
    y = c(19, 16.5, 13.5, 10.5, 7.5, 16.5, 13.5, 10.5, 10.5, 4, 14),
    width = c(3.5, 3.5, 3.5, 3.5, 3.5, 2.5, 2.5, 2.5, 2.5, 3.5, 2.8),
    height = c(1.2, 1.8, 1.8, 1.8, 1.5, 1.2, 1.2, 1.2, 1.5, 1.8, 2.5),
    color = c(CONFIG$colors$start, CONFIG$colors$tier1, CONFIG$colors$tier2, 
              CONFIG$colors$tier3, CONFIG$colors$tier4, CONFIG$colors$success,
              CONFIG$colors$success, CONFIG$colors$success, CONFIG$colors$fail,
              CONFIG$colors$metrics, CONFIG$colors$validation),
    stringsAsFactors = FALSE
  )
  
  # Define edges (connections)
  edges <- data.frame(
    from_id = c("start", "t1", "t1", "t2", "t2", "t3", "t3", "s1", "s2", "s3", 
                "agg", "validate", "validate", "validate"),
    to_id = c("t1", "s1", "t2", "s2", "t3", "s3", "fail", "agg", "agg", "agg",
              "metrics", "t1", "t2", "t3"),
    label = c("", "Match", "No Match", "Match", "No Match", "Match", "No Match",
              "", "", "", "", "Perfected", "Perfected", "Perfected"),
    type = c("main", "success", "fail", "success", "fail", "success", "fail",
             "main", "main", "main", "main", "validate", "validate", "validate"),
    stringsAsFactors = FALSE
  )
  
  # Merge node positions with edges
  edges <- edges %>%
    left_join(nodes %>% select(id, from_x = x, from_y = y), by = c("from_id" = "id")) %>%
    left_join(nodes %>% select(id, to_x = x, to_y = y), by = c("to_id" = "id"))
  
  # Create base plot
  p <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = CONFIG$colors$bg, color = NA),
      panel.background = element_rect(fill = CONFIG$colors$bg, color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_fixed(ratio = 1, xlim = c(0, 11), ylim = c(0, 20))
  
  # Add edges (arrows)
  for (i in 1:nrow(edges)) {
    edge <- edges[i, ]
    
    # Determine arrow color and style
    arrow_color <- if (edge$type == "validate") CONFIG$colors$validation else CONFIG$colors$arrow
    arrow_type <- if (edge$type == "validate") "dashed" else "solid"
    arrow_size <- if (edge$type == "validate") 0.8 else 1.0
    
    # Calculate arrow endpoints (adjust for node size)
    from_y_adj <- edge$from_y - 0.6
    to_y_adj <- edge$to_y + 0.9
    
    # Adjust x for side connections
    from_x <- edge$from_x
    to_x <- edge$to_x
    
    if (edge$to_id %in% c("s1", "s2", "s3")) {
      from_x <- edge$from_x + 1.75
      to_x <- edge$to_x - 1.25
    } else if (edge$to_id == "fail") {
      from_x <- edge$from_x - 1.75
      to_x <- edge$to_x + 1.25
    }
    
    p <- p + geom_segment(
      aes(x = from_x, y = from_y_adj, xend = to_x, yend = to_y_adj),
      arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
      color = arrow_color,
      size = arrow_size,
      linetype = arrow_type
    )
    
    # Add edge labels
    if (edge$label != "") {
      label_x <- (from_x + to_x) / 2
      label_y <- (from_y_adj + to_y_adj) / 2
      
      p <- p + annotate("text", x = label_x, y = label_y + 0.3,
                       label = edge$label, size = CONFIG$text$small,
                       fontface = "italic", color = arrow_color)
    }
  }
  
  # Add nodes (rectangles with rounded corners)
  for (i in 1:nrow(nodes)) {
    node <- nodes[i, ]
    
    p <- p + annotate("rect",
                     xmin = node$x - node$width/2, xmax = node$x + node$width/2,
                     ymin = node$y - node$height/2, ymax = node$y + node$height/2,
                     fill = node$color, color = "#333333", size = 1.2, alpha = 0.85)
  }
  
  # Add node labels
  p <- p +
    # START
    annotate("text", x = 5, y = 19.3, label = "Environmental Census Data",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 18.7, label = "(16S/18S OTUs with taxids)",
             size = CONFIG$text$small, color = "white") +
    
    # TIER 1
    annotate("text", x = 5, y = 17.1, label = "🥇 TIER 1: Taxid Matching",
             size = CONFIG$text$tier, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 16.6, label = "PRIMARY METHOD",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 16.1, label = "85-90% match rate",
             size = CONFIG$text$detail, color = "white") +
    annotate("text", x = 5, y = 15.7, label = "census_taxid == ncbi_taxid",
             size = CONFIG$text$small, fontface = "italic", color = "white") +

    # TIER 2
    annotate("text", x = 5, y = 14.1, label = "🥈 TIER 2: Lineage Matching",
             size = CONFIG$text$tier, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 13.6, label = "FALLBACK METHOD",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 13.1, label = "10-15% match rate",
             size = CONFIG$text$detail, color = "white") +
    annotate("text", x = 5, y = 12.7, label = "Handles misnaming",
             size = CONFIG$text$small, fontface = "italic", color = "white") +

    # TIER 3
    annotate("text", x = 5, y = 11.1, label = "🥉 TIER 3: Direct Name",
             size = CONFIG$text$tier, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 10.6, label = "VALIDATION METHOD",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 10.1, label = "<5% match rate",
             size = CONFIG$text$detail, color = "white") +
    annotate("text", x = 5, y = 9.7, label = "Edge cases",
             size = CONFIG$text$small, fontface = "italic", color = "white") +

    # TIER 4 - AGGREGATION
    annotate("text", x = 5, y = 8.1, label = "TIER 4: Hierarchical Aggregation",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 7.5, label = "Sum all descendant genomes",
             size = CONFIG$text$detail, color = "white") +

    # SUCCESS NODES
    annotate("text", x = 8.5, y = 16.8, label = "✅ Matched",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 8.5, y = 16.2, label = "Highest Confidence",
             size = CONFIG$text$small, color = "white") +

    annotate("text", x = 8.5, y = 13.8, label = "✅ Matched",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 8.5, y = 13.2, label = "High Confidence",
             size = CONFIG$text$small, color = "white") +

    annotate("text", x = 8.5, y = 10.8, label = "✅ Matched",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 8.5, y = 10.2, label = "Moderate Confidence",
             size = CONFIG$text$small, color = "white") +

    # FAIL NODE
    annotate("text", x = 1.5, y = 11, label = "❌ Census-Only",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 1.5, y = 10.5, label = "No genomic data",
             size = CONFIG$text$small, color = "white") +
    annotate("text", x = 1.5, y = 10, label = "novelty_factor = ∞",
             size = CONFIG$text$small, fontface = "italic", color = "white") +

    # METRICS NODE
    annotate("text", x = 5, y = 4.8, label = "Calculate COINED METRICS",
             size = CONFIG$text$tier, fontface = "bold", color = "white") +
    annotate("text", x = 5, y = 4.3, label = "🔴 novelty_factor",
             size = CONFIG$text$label, color = "white") +
    annotate("text", x = 5, y = 3.8, label = "🔵 overrepresentation_factor",
             size = CONFIG$text$label, color = "white") +
    annotate("text", x = 5, y = 3.3, label = "(Original Contributions)",
             size = CONFIG$text$small, fontface = "italic", color = "white") +

    # VALIDATION NODE
    annotate("text", x = 9.5, y = 15, label = "GTDB-NCBI",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 9.5, y = 14.5, label = "Validation Pipeline",
             size = CONFIG$text$label, fontface = "bold", color = "white") +
    annotate("text", x = 9.5, y = 13.9, label = "Proves complete",
             size = CONFIG$text$small, color = "white") +
    annotate("text", x = 9.5, y = 13.5, label = "overlap",
             size = CONFIG$text$small, color = "white") +
    annotate("text", x = 9.5, y = 13.1, label = "Validates algorithm",
             size = CONFIG$text$small, color = "white") +
    annotate("text", x = 9.5, y = 12.7, label = "effectiveness",
             size = CONFIG$text$small, color = "white") +

    # Title
    annotate("text", x = 5.5, y = 19.8, label = "Hierarchical Matching Algorithm (Perfected via GTDB-NCBI Validation)",
             size = CONFIG$text$title, fontface = "bold", color = CONFIG$colors$text, hjust = 0.5)

  return(p)
}

# Main execution
cat("Creating hierarchical matching algorithm diagram...\n")

# Create the diagram
diagram <- create_matching_diagram()

# Save as PNG
png_file <- "hierarchical_matching_algorithm.png"
ggsave(png_file, diagram, width = CONFIG$width, height = CONFIG$height, dpi = CONFIG$dpi, bg = "white")
cat(sprintf("✅ Saved PNG: %s\n", png_file))

# Save as PDF
pdf_file <- "hierarchical_matching_algorithm.pdf"
ggsave(pdf_file, diagram, width = CONFIG$width, height = CONFIG$height, device = cairo_pdf, bg = "white")
cat(sprintf("✅ Saved PDF: %s\n", pdf_file))

# Save high-res version for publication
png_file_hires <- "hierarchical_matching_algorithm_hires.png"
ggsave(png_file_hires, diagram, width = CONFIG$width, height = CONFIG$height, dpi = 600, bg = "white")
cat(sprintf("✅ Saved high-res PNG (600 DPI): %s\n", png_file_hires))

cat("\n🎉 Diagram generation complete!\n")
cat(sprintf("   - Standard PNG: %s (%d DPI)\n", png_file, CONFIG$dpi))
cat(sprintf("   - High-res PNG: %s (600 DPI)\n", png_file_hires))
cat(sprintf("   - PDF: %s\n", pdf_file))
cat("\nAll files saved in: visuals/\n")

