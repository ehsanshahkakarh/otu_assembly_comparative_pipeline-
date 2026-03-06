#!/usr/bin/env Rscript
# Hierarchical Matching Algorithm Flowchart (Base R Graphics)
# Created: 2026-01-05
# Purpose: Publication-quality diagram using base R graphics (no dependencies)

# Configuration
CONFIG <- list(
  # Output settings
  width = 16,
  height = 20,
  
  # Colors matching your existing plots
  colors = list(
    tier1 = "#32CD32",        # Lime green - PRIMARY
    tier2 = "#ffb44c",        # Light orange - FALLBACK
    tier3 = "#ff7200",        # Orange - VALIDATION
    tier4 = "#46bda3",        # Aqua - AGGREGATION
    success = "#4c9b34",      # Dark green - matched
    fail = "#d19386",         # Light brown - census-only
    validation = "#c73de4",   # Magenta - GTDB-NCBI validation
    metrics = "#f51b7f",      # Bright pink - COINED METRICS
    start = "#416b7d",        # Dark blue-gray - input
    arrow = "#333333",        # Dark gray for arrows
    text = "#000000"          # Black for text
  )
)

# Function to draw a rounded rectangle
draw_box <- function(x, y, width, height, col, border = "#333333", lwd = 2) {
  rect(x - width/2, y - height/2, x + width/2, y + height/2,
       col = col, border = border, lwd = lwd)
}

# Function to draw an arrow
draw_arrow <- function(x0, y0, x1, y1, col = "#333333", lwd = 2, lty = 1) {
  arrows(x0, y0, x1, y1, col = col, lwd = lwd, lty = lty, length = 0.15, angle = 20)
}

# Function to add multi-line text
add_text <- function(x, y, labels, cex = 1, col = "white", font = 1) {
  n <- length(labels)
  y_positions <- seq(y + (n-1) * 0.15, y - (n-1) * 0.15, length.out = n)
  for (i in 1:n) {
    text(x, y_positions[i], labels[i], cex = cex, col = col, font = font)
  }
}

# Create the diagram
create_diagram <- function() {
  
  # Node positions (x, y)
  nodes <- list(
    start = list(x = 5, y = 18, w = 3.5, h = 1),
    t1 = list(x = 5, y = 15.5, w = 3.5, h = 1.5),
    t2 = list(x = 5, y = 12.5, w = 3.5, h = 1.5),
    t3 = list(x = 5, y = 9.5, w = 3.5, h = 1.5),
    agg = list(x = 5, y = 6.5, w = 3.5, h = 1.2),
    s1 = list(x = 8.5, y = 15.5, w = 2.2, h = 1),
    s2 = list(x = 8.5, y = 12.5, w = 2.2, h = 1),
    s3 = list(x = 8.5, y = 9.5, w = 2.2, h = 1),
    fail = list(x = 1.5, y = 9.5, w = 2.2, h = 1.2),
    metrics = list(x = 5, y = 3.5, w = 3.5, h = 1.5),
    validate = list(x = 9.5, y = 13, w = 2.5, h = 2)
  )
  
  # Draw arrows first (so they appear behind boxes)
  # Main flow
  draw_arrow(nodes$start$x, nodes$start$y - 0.5, nodes$t1$x, nodes$t1$y + 0.75)
  draw_arrow(nodes$t1$x, nodes$t1$y - 0.75, nodes$t2$x, nodes$t2$y + 0.75)
  draw_arrow(nodes$t2$x, nodes$t2$y - 0.75, nodes$t3$x, nodes$t3$y + 0.75)
  draw_arrow(nodes$t3$x, nodes$t3$y - 0.75, nodes$fail$x + 1.1, nodes$fail$y + 0.6)
  
  # Success paths
  draw_arrow(nodes$t1$x + 1.75, nodes$t1$y, nodes$s1$x - 1.1, nodes$s1$y, col = CONFIG$colors$success)
  draw_arrow(nodes$t2$x + 1.75, nodes$t2$y, nodes$s2$x - 1.1, nodes$s2$y, col = CONFIG$colors$success)
  draw_arrow(nodes$t3$x + 1.75, nodes$t3$y, nodes$s3$x - 1.1, nodes$s3$y, col = CONFIG$colors$success)
  
  # Aggregation paths
  draw_arrow(nodes$s1$x, nodes$s1$y - 0.5, nodes$agg$x + 1.5, nodes$agg$y + 0.6)
  draw_arrow(nodes$s2$x, nodes$s2$y - 0.5, nodes$agg$x + 1, nodes$agg$y + 0.6)
  draw_arrow(nodes$s3$x, nodes$s3$y - 0.5, nodes$agg$x + 0.5, nodes$agg$y + 0.6)
  
  # Metrics path
  draw_arrow(nodes$agg$x, nodes$agg$y - 0.6, nodes$metrics$x, nodes$metrics$y + 0.75)
  
  # Validation arrows (dashed)
  draw_arrow(nodes$validate$x - 1.25, nodes$validate$y + 0.7, nodes$t1$x + 1.75, nodes$t1$y + 0.5,
             col = CONFIG$colors$validation, lty = 2, lwd = 1.5)
  draw_arrow(nodes$validate$x - 1.25, nodes$validate$y, nodes$t2$x + 1.75, nodes$t2$y + 0.3,
             col = CONFIG$colors$validation, lty = 2, lwd = 1.5)
  draw_arrow(nodes$validate$x - 1.25, nodes$validate$y - 0.7, nodes$t3$x + 1.75, nodes$t3$y + 0.1,
             col = CONFIG$colors$validation, lty = 2, lwd = 1.5)
  
  # Draw boxes
  draw_box(nodes$start$x, nodes$start$y, nodes$start$w, nodes$start$h, CONFIG$colors$start)
  draw_box(nodes$t1$x, nodes$t1$y, nodes$t1$w, nodes$t1$h, CONFIG$colors$tier1)
  draw_box(nodes$t2$x, nodes$t2$y, nodes$t2$w, nodes$t2$h, CONFIG$colors$tier2)
  draw_box(nodes$t3$x, nodes$t3$y, nodes$t3$w, nodes$t3$h, CONFIG$colors$tier3)
  draw_box(nodes$agg$x, nodes$agg$y, nodes$agg$w, nodes$agg$h, CONFIG$colors$tier4)
  draw_box(nodes$s1$x, nodes$s1$y, nodes$s1$w, nodes$s1$h, CONFIG$colors$success)
  draw_box(nodes$s2$x, nodes$s2$y, nodes$s2$w, nodes$s2$h, CONFIG$colors$success)
  draw_box(nodes$s3$x, nodes$s3$y, nodes$s3$w, nodes$s3$h, CONFIG$colors$success)
  draw_box(nodes$fail$x, nodes$fail$y, nodes$fail$w, nodes$fail$h, CONFIG$colors$fail)
  draw_box(nodes$metrics$x, nodes$metrics$y, nodes$metrics$w, nodes$metrics$h, CONFIG$colors$metrics)
  draw_box(nodes$validate$x, nodes$validate$y, nodes$validate$w, nodes$validate$h, CONFIG$colors$validation)
  
  # Add text labels
  # START
  add_text(nodes$start$x, nodes$start$y, 
           c("Environmental Census Data", "(16S/18S OTUs with taxids)"),
           cex = 1.1, font = 2)
  
  # TIER 1
  add_text(nodes$t1$x, nodes$t1$y,
           c("TIER 1: Taxid Matching", "PRIMARY METHOD", "85-90% match rate", "census_taxid == ncbi_taxid"),
           cex = c(1.2, 1.0, 0.9, 0.8), font = c(2, 2, 1, 3))
  
  # TIER 2
  add_text(nodes$t2$x, nodes$t2$y,
           c("TIER 2: Lineage Matching", "FALLBACK METHOD", "10-15% match rate", "Handles misnaming"),
           cex = c(1.2, 1.0, 0.9, 0.8), font = c(2, 2, 1, 3))
  
  # TIER 3
  add_text(nodes$t3$x, nodes$t3$y,
           c("TIER 3: Direct Name", "VALIDATION METHOD", "<5% match rate", "Edge cases"),
           cex = c(1.2, 1.0, 0.9, 0.8), font = c(2, 2, 1, 3))
  
  # AGGREGATION
  add_text(nodes$agg$x, nodes$agg$y,
           c("TIER 4: Hierarchical Aggregation", "Sum all descendant genomes"),
           cex = c(1.0, 0.9), font = c(2, 1))
  
  # SUCCESS NODES
  add_text(nodes$s1$x, nodes$s1$y, c("Matched", "Highest Confidence"), cex = c(1.0, 0.8), font = c(2, 1))
  add_text(nodes$s2$x, nodes$s2$y, c("Matched", "High Confidence"), cex = c(1.0, 0.8), font = c(2, 1))
  add_text(nodes$s3$x, nodes$s3$y, c("Matched", "Moderate Confidence"), cex = c(1.0, 0.8), font = c(2, 1))
  
  # FAIL NODE
  add_text(nodes$fail$x, nodes$fail$y,
           c("Census-Only", "No genomic data", "novelty_factor = inf"),
           cex = c(1.0, 0.8, 0.8), font = c(2, 1, 3))

