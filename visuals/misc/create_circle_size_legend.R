#!/usr/bin/env Rscript
# Circle Size Legend Generator
# Created: 2025-01-26
# Purpose: Create a standalone circle size legend that matches both 16S and 18S mega visuals
# This legend can be manually added to figures in external applications

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(cowplot)
})

# Set working directory to script location
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "create_circle_size_legend.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
setwd(script_dir)

# Configuration matching the mega visual scripts
config <- list(
  output_dir = file.path("final_visualizations"),
  size_range = c(10, 22),  # Same as both 16S and 18S scripts
  dpi = 300,
  legend_width = 4,   # Wider to accommodate labels moved far right
  legend_height = 2.5, # Compact height for 3 circles
  background_color = "white"
)

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Function to create representative circle size data (simplified to 3 examples)
create_representative_data <- function() {
  # Create 3 representative Genome/Isolate ratios that span the typical range
  # Based on analysis of actual data ranges from both 16S and 18S datasets

  # Three key ratios: low, medium, high
  genome_isolate_ratios <- c(1.0, 10.0, 100.0)

  # Apply the same transformation as both mega visual scripts
  circle_size_raw <- sqrt(genome_isolate_ratios)

  # Normalize using the same approach (3 + 9 * normalized)
  min_raw <- min(circle_size_raw)
  max_raw <- max(circle_size_raw)
  circle_sizes <- 3 + 9 * (circle_size_raw - min_raw) / (max_raw - min_raw)

  # Create data frame
  data.frame(
    Genome_Isolate_Ratio = genome_isolate_ratios,
    Circle_Size = circle_sizes,
    Label = paste0(genome_isolate_ratios, "×")
  )
}

# Create the circle size legend
create_circle_size_legend <- function() {
  # Get representative data
  legend_data <- create_representative_data()
  
  cat("Circle size legend data:\n")
  print(legend_data)
  
  # Create compact legend with circles vertically stacked and labels far to the right
  legend_plot <- ggplot(legend_data, aes(x = 0, y = rev(seq_along(Genome_Isolate_Ratio)))) +
    geom_point(aes(size = Circle_Size), color = "black", alpha = 0.9, stroke = 0.2) +
    geom_text(aes(label = Label), hjust = -0.5, size = 10, fontface = "bold", color = "black") +
    scale_size_continuous(range = config$size_range, guide = "none") +
    scale_x_continuous(limits = c(-0.5, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, length(legend_data$Genome_Isolate_Ratio) + 0.5),
                      expand = c(0, 0)) +
    labs(title = "Genome/Isolate Ratio") +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0,
                               margin = margin(b = 3)),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  return(legend_plot)
}

# Create reference comparison plot for fidelity checking
create_reference_plot <- function() {
  # Create a simple scatter plot with the same 3 ratios for consistency
  ref_data <- data.frame(
    x = c(10, 100, 1000),
    y = c(10, 100, 1000),
    ratio = c(1, 10, 100),
    size_raw = sqrt(c(1, 10, 100))
  )
  
  # Apply same transformation
  min_raw <- min(ref_data$size_raw)
  max_raw <- max(ref_data$size_raw)
  ref_data$circle_size <- 3 + 9 * (ref_data$size_raw - min_raw) / (max_raw - min_raw)
  
  ref_plot <- ggplot(ref_data, aes(x = x, y = y)) +
    geom_point(aes(size = circle_size), color = "black", alpha = 0.8) +
    geom_text(aes(label = paste0(ratio, "× ratio")),
              hjust = -0.2, vjust = -0.5, size = 10, fontface = "bold") +
    scale_size_continuous(range = config$size_range, guide = "none") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Reference Plot for Size Fidelity Check",
         subtitle = "Compare circle sizes with your main plots",
         x = "X axis (log scale)", y = "Y axis (log scale)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.background = element_rect(fill = config$background_color, color = NA),
      panel.background = element_rect(fill = config$background_color, color = NA)
    )
  
  return(ref_plot)
}

# Main function
main <- function() {
  cat("Creating Circle Size Legend\n")
  cat("===========================\n")
  
  # Create the legend
  legend_plot <- create_circle_size_legend()
  
  # Create reference plot
  reference_plot <- create_reference_plot()
  
  # Save the legend
  legend_file_png <- file.path(config$output_dir, "circle_size_legend.png")
  legend_file_pdf <- file.path(config$output_dir, "circle_size_legend.pdf")
  
  cat(paste("Saving circle size legend to:", legend_file_png, "\n"))
  ggsave(legend_file_png, legend_plot,
         width = config$legend_width, height = config$legend_height,
         dpi = config$dpi, bg = config$background_color, limitsize = FALSE,
         units = "in")
  
  ggsave(legend_file_pdf, legend_plot,
         width = config$legend_width, height = config$legend_height,
         dpi = config$dpi, bg = config$background_color, limitsize = FALSE,
         units = "in")
  
  # Save the reference plot
  ref_file_png <- file.path(config$output_dir, "circle_size_reference.png")
  ref_file_pdf <- file.path(config$output_dir, "circle_size_reference.pdf")
  
  cat(paste("Saving reference plot to:", ref_file_png, "\n"))
  ggsave(ref_file_png, reference_plot,
         width = 8, height = 6,
         dpi = config$dpi, bg = config$background_color, limitsize = FALSE,
         units = "in")
  
  ggsave(ref_file_pdf, reference_plot,
         width = 8, height = 6,
         dpi = config$dpi, bg = config$background_color, limitsize = FALSE,
         units = "in")
  
  cat("✅ Circle size legend creation complete!\n")
  cat(paste("   Legend files:", legend_file_png, "and", legend_file_pdf, "\n"))
  cat(paste("   Reference files:", ref_file_png, "and", ref_file_pdf, "\n"))
  cat(paste("   Legend dimensions:", config$legend_width, "x", config$legend_height, "inches\n"))
  cat(paste("   DPI:", config$dpi, "\n"))
  
  # Print scaling information
  cat("\n📏 SCALING FIDELITY GUIDE:\n")
  cat("========================\n")
  cat("1. Open your main plot and the reference plot side by side\n")
  cat("2. Compare circle sizes for the same ratios (1×, 10×, 100×)\n")
  cat("3. If circles don't match, adjust legend scale proportionally\n")
  cat("4. The legend uses the same size_range (10-22) and transformation\n")
  cat("5. At 300 DPI, circles should match exactly when imported at 100% scale\n")
}

# Run the legend creation
if (!interactive()) main()
