#!/usr/bin/env Rscript

# Color Palette Viewer for Mega Comprehensive Stacked Visual
# This script reads JSON color files and creates a visual preview

library(jsonlite)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Function to read JSON color file
read_color_json <- function(json_file) {
  if (!file.exists(json_file)) {
    stop(paste("JSON file not found:", json_file))
  }
  return(fromJSON(json_file))
}

# Function to create color palette visualization
create_color_plot <- function(colors, title, labels = NULL) {
  n_colors <- length(colors)
  
  # Create data frame for plotting
  df <- data.frame(
    x = rep(1, n_colors),
    y = 1:n_colors,
    color = colors,
    label = if(is.null(labels)) paste("Color", 1:n_colors) else labels,
    stringsAsFactors = FALSE
  )
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y, fill = color)) +
    geom_tile(width = 0.8, height = 0.8, color = "white", size = 1) +
    geom_text(aes(label = paste(label, "\n", color)), 
              color = "white", fontface = "bold", size = 3) +
    scale_fill_identity() +
    scale_y_continuous(breaks = 1:n_colors, labels = df$label, 
                       trans = "reverse") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    ) +
    coord_fixed(ratio = 0.3)
  
  return(p)
}

# Main function to create all visualizations
main <- function() {
  cat("Creating color palette visualizations...\n")
  
  # Read JSON files
  bacteria_data <- read_color_json("bacteria_colors.json")
  archaea_data <- read_color_json("archaea_colors.json")
  eukaryota_data <- read_color_json("eukaryota_colors.json")
  
  # Extract color vectors
  bacteria_colors <- unlist(bacteria_data$bacteria_colors)
  archaea_colors <- unlist(archaea_data$archaea_colors)
  eukaryota_colors <- unlist(eukaryota_data$eukaryota_colors)
  
  # Create individual plots
  bacteria_plot <- create_color_plot(
    bacteria_colors, 
    "Bacteria Colors (16 colors)",
    names(bacteria_colors)
  )
  
  archaea_plot <- create_color_plot(
    archaea_colors, 
    "Archaea Colors (5 colors)",
    names(archaea_colors)
  )
  
  eukaryota_plot <- create_color_plot(
    eukaryota_colors, 
    "Eukaryota Colors (13 colors)",
    names(eukaryota_colors)
  )
  
  # Arrange plots side by side
  combined_plot <- grid.arrange(
    bacteria_plot, archaea_plot, eukaryota_plot,
    ncol = 3,
    top = "Color Palettes for Mega Comprehensive Stacked Visual"
  )
  
  # Save the visualization
  output_file <- "color_palette_preview.png"
  ggsave(output_file, combined_plot, 
         width = 18, height = 10, dpi = 300, bg = "white")
  
  cat(paste("Color palette preview saved to:", output_file, "\n"))
  
  # Print summary
  cat("\nColor Summary:\n")
  cat("=============\n")
  cat(paste("Bacteria:", length(bacteria_colors), "colors\n"))
  cat(paste("Archaea:", length(archaea_colors), "colors\n"))
  cat(paste("Eukaryota:", length(eukaryota_colors), "colors\n"))
  cat(paste("Total:", length(bacteria_colors) + length(archaea_colors) + length(eukaryota_colors), "colors\n"))
  
  return(combined_plot)
}

# Run the main function if script is executed directly
if (!interactive()) {
  main()
}
