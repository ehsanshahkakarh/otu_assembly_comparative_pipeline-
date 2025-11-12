#!/usr/bin/env Rscript

# Simple Color Palette Viewer using base R graphics
# This script reads JSON color files and creates a visual preview

library(jsonlite)

# Function to read JSON color file
read_color_json <- function(json_file) {
  if (!file.exists(json_file)) {
    stop(paste("JSON file not found:", json_file))
  }
  return(fromJSON(json_file))
}

# Function to create color palette visualization using base R
create_color_plot <- function() {
  cat("Creating color palette visualizations...\n")
  
  # Read JSON files
  bacteria_data <- read_color_json("bacteria_colors.json")
  archaea_data <- read_color_json("archaea_colors.json")
  eukaryota_data <- read_color_json("eukaryota_colors.json")
  
  # Extract color vectors
  bacteria_colors <- unlist(bacteria_data$bacteria_colors)
  archaea_colors <- unlist(archaea_data$archaea_colors)
  eukaryota_colors <- unlist(eukaryota_data$eukaryota_colors)
  
  # Create PNG file
  png("color_palette_preview.png", width = 1800, height = 1000, res = 150)
  
  # Set up the layout
  par(mfrow = c(1, 3), mar = c(4, 2, 3, 2), oma = c(0, 0, 3, 0))
  
  # Plot bacteria colors
  plot(1, type = "n", xlim = c(0, 2), ylim = c(0, length(bacteria_colors) + 1),
       xlab = "", ylab = "", main = "Bacteria Colors (16)", axes = FALSE)
  for (i in 1:length(bacteria_colors)) {
    rect(0.2, i - 0.4, 1.8, i + 0.4, col = bacteria_colors[i], border = "white", lwd = 2)
    text(1, i, paste(names(bacteria_colors)[i], "\n", bacteria_colors[i]), 
         col = "white", cex = 0.7, font = 2)
  }
  
  # Plot archaea colors
  plot(1, type = "n", xlim = c(0, 2), ylim = c(0, length(archaea_colors) + 1),
       xlab = "", ylab = "", main = "Archaea Colors (5)", axes = FALSE)
  for (i in 1:length(archaea_colors)) {
    rect(0.2, i - 0.4, 1.8, i + 0.4, col = archaea_colors[i], border = "white", lwd = 2)
    text(1, i, paste(names(archaea_colors)[i], "\n", archaea_colors[i]), 
         col = "white", cex = 0.8, font = 2)
  }
  
  # Plot eukaryota colors
  plot(1, type = "n", xlim = c(0, 2), ylim = c(0, length(eukaryota_colors) + 1),
       xlab = "", ylab = "", main = "Eukaryota Colors (13)", axes = FALSE)
  for (i in 1:length(eukaryota_colors)) {
    rect(0.2, i - 0.4, 1.8, i + 0.4, col = eukaryota_colors[i], border = "white", lwd = 2)
    text(1, i, paste(names(eukaryota_colors)[i], "\n", eukaryota_colors[i]), 
         col = "white", cex = 0.7, font = 2)
  }
  
  # Add main title
  mtext("Color Palettes for Mega Comprehensive Stacked Visual", 
        outer = TRUE, cex = 1.5, font = 2)
  
  dev.off()
  
  cat("Color palette preview saved to: color_palette_preview.png\n")
  
  # Print summary
  cat("\nColor Summary:\n")
  cat("=============\n")
  cat(paste("Bacteria:", length(bacteria_colors), "colors\n"))
  cat(paste("Archaea:", length(archaea_colors), "colors\n"))
  cat(paste("Eukaryota:", length(eukaryota_colors), "colors\n"))
  cat(paste("Total:", length(bacteria_colors) + length(archaea_colors) + length(eukaryota_colors), "colors\n"))
}

# Run the main function
create_color_plot()
