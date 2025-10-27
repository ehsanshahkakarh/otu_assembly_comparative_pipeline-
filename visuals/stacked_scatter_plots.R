#!/usr/bin/env Rscript

# Stacked Scatter Plot Visualization Script
# Created: 2025-10-26
# Purpose: Generate stacked scatter plots for Bacteria, Archaea, and Eukaryota
# showing top 10 novelty_factor and overrepresentation_factor taxa

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(gridExtra)
  library(grid)
})

# Set up paths
path_18s <- "../Eukcensus_merge/18s_merged/csv_results"
path_16s <- "../Eukcensus_merge/16s_merged/csv_results"

# Function to read and process data
read_and_process_data <- function(file_path, domain_filter = NULL) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Filter by domain if specified
  if (!is.null(domain_filter) && "domain" %in% colnames(data)) {
    data <- data[data$domain == domain_filter, ]
  }
  
  # Filter for matched entries only
  if ("match_status" %in% colnames(data)) {
    data <- data[data$match_status == "matched", ]
  }
  
  # Remove rows with missing novelty_factor or overrepresentation_factor
  data <- data[!is.na(data$novelty_factor) & !is.na(data$overrepresentation_factor), ]
  
  # Filter for factors > 1.0
  data <- data[data$novelty_factor > 1.0 & data$overrepresentation_factor > 1.0, ]
  
  return(data)
}

# Function to get top taxa for annotation
get_top_taxa <- function(data, taxon_col, n = 10) {
  if (nrow(data) == 0) return(data.frame())
  
  # Get top 10 by novelty_factor
  top_novelty <- data %>%
    arrange(desc(novelty_factor)) %>%
    head(n) %>%
    mutate(annotation_type = "novelty")
  
  # Get top 10 by overrepresentation_factor
  top_overrep <- data %>%
    arrange(desc(overrepresentation_factor)) %>%
    head(n) %>%
    mutate(annotation_type = "overrepresentation")
  
  # Combine and remove duplicates
  top_taxa <- rbind(top_novelty, top_overrep) %>%
    distinct(!!sym(taxon_col), .keep_all = TRUE)
  
  return(top_taxa)
}

# Function to get domain-specific colors
get_domain_colors <- function(domain_name) {
  colors <- list(
    "Bacteria" = "#2E8B57",      # Sea green
    "Archaea" = "#C0342B",       # Red
    "Eukaryota" = "#4682B4"      # Steel blue
  )
  return(colors[[domain_name]])
}

# Function to create scatter plot
create_scatter_plot <- function(data, taxon_col, title, domain_name) {
  if (is.null(data) || nrow(data) == 0) {
    # Return empty plot
    return(ggplot() +
           ggtitle(paste(title, "(No data)")) +
           theme_minimal() +
           theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
  }

  # Get top taxa for annotation
  top_taxa <- get_top_taxa(data, taxon_col)

  # Determine circle size column
  size_col <- if ("isolate_percentage" %in% colnames(data)) "isolate_percentage" else "size_percentage"

  # Get domain-specific color
  domain_color <- get_domain_colors(domain_name)

  # Create base plot
  p <- ggplot(data, aes(x = overrepresentation_factor, y = novelty_factor)) +
    geom_point(aes(size = !!sym(size_col)), alpha = 0.6, color = domain_color, stroke = 0.5) +
    scale_size_continuous(name = "Size %", range = c(3, 15),
                         guide = guide_legend(override.aes = list(alpha = 1, color = domain_color))) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    labs(
      title = title,
      x = "Overrepresentation Factor",
      y = "Novelty Factor"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )

  # Add annotations for top taxa
  if (nrow(top_taxa) > 0) {
    # Create annotation labels with factor values (showing the higher of the two factors)
    top_taxa$max_factor <- pmax(top_taxa$novelty_factor, top_taxa$overrepresentation_factor)
    top_taxa$label <- paste0(top_taxa[[taxon_col]],
                            " (", format(round(top_taxa$max_factor, 1), nsmall = 1), "×)")

    p <- p +
      geom_point(data = top_taxa, aes(x = overrepresentation_factor, y = novelty_factor, size = !!sym(size_col)),
                 color = "black", fill = "yellow", shape = 21, stroke = 1.5, alpha = 0.9) +
      geom_text_repel(data = top_taxa,
                      aes(x = overrepresentation_factor, y = novelty_factor, label = label),
                      size = 3.5, color = "black", fontface = "bold",
                      segment.color = "black", segment.size = 0.4,
                      box.padding = 0.4, point.padding = 0.3,
                      max.overlaps = 25, force = 3, min.segment.length = 0)
  }

  return(p)
}

# Read data files
cat("Reading data files...\n")

# 16S data (Bacteria and Archaea)
bacteria_phylum <- read_and_process_data(file.path(path_16s, "16s_ncbi_merged_clean_phylum.csv"), "Bacteria")
bacteria_family <- read_and_process_data(file.path(path_16s, "16s_ncbi_merged_clean_family.csv"), "Bacteria")
archaea_phylum <- read_and_process_data(file.path(path_16s, "16s_ncbi_merged_clean_phylum.csv"), "Archaea")
archaea_family <- read_and_process_data(file.path(path_16s, "16s_ncbi_merged_clean_family.csv"), "Archaea")

# 18S data (Eukaryota)
eukaryota_phylum <- read_and_process_data(file.path(path_18s, "18s_ncbi_merged_clean_phylum.csv"), "Eukaryota")
eukaryota_family <- read_and_process_data(file.path(path_18s, "18s_ncbi_merged_clean_family.csv"), "Eukaryota")

# Create plots
cat("Creating scatter plots...\n")

# Bacteria plots
bacteria_phylum_plot <- create_scatter_plot(bacteria_phylum, "phylum", "Bacteria - Phylum Level", "Bacteria")
bacteria_family_plot <- create_scatter_plot(bacteria_family, "family", "Bacteria - Family Level", "Bacteria")

# Archaea plots
archaea_phylum_plot <- create_scatter_plot(archaea_phylum, "phylum", "Archaea - Phylum Level", "Archaea")
archaea_family_plot <- create_scatter_plot(archaea_family, "family", "Archaea - Family Level", "Archaea")

# Eukaryota plots
eukaryota_phylum_plot <- create_scatter_plot(eukaryota_phylum, "phylum", "Eukaryota - Phylum Level", "Eukaryota")
eukaryota_family_plot <- create_scatter_plot(eukaryota_family, "family", "Eukaryota - Family Level", "Eukaryota")

# Create stacked layout
cat("Arranging plots...\n")

# Create domain labels for left side
domain_labels <- list(
  textGrob("Bacteria", rot = 90, gp = gpar(fontsize = 14, fontface = "bold", col = "#2E8B57")),
  textGrob("Archaea", rot = 90, gp = gpar(fontsize = 14, fontface = "bold", col = "#C0342B")),
  textGrob("Eukaryota", rot = 90, gp = gpar(fontsize = 14, fontface = "bold", col = "#4682B4"))
)

# Arrange plots in a 3x2 grid (3 domains x 2 taxonomic levels)
stacked_plot <- grid.arrange(
  bacteria_phylum_plot, bacteria_family_plot,
  archaea_phylum_plot, archaea_family_plot,
  eukaryota_phylum_plot, eukaryota_family_plot,
  ncol = 2, nrow = 3,
  top = textGrob("Taxonomic Novelty and Overrepresentation Analysis\nTop 10 Taxa by Novelty and Overrepresentation Factors",
                 gp = gpar(fontsize = 18, fontface = "bold")),
  bottom = textGrob("Phylum Level (left column) | Family Level (right column)\nHighlighted taxa show top 10 by each factor",
                    gp = gpar(fontsize = 12, fontface = "italic"))
)

# Create output directory if it doesn't exist
output_dir <- "scatter_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the plot
output_file <- file.path(output_dir, "stacked_scatter_plots.pdf")
cat(paste("Saving plot to:", output_file, "\n"))

pdf(output_file, width = 18, height = 22)
grid.draw(stacked_plot)
dev.off()

# Also save as PNG
png_file <- file.path(output_dir, "stacked_scatter_plots.png")
png(png_file, width = 1800, height = 2200, res = 100)
grid.draw(stacked_plot)
dev.off()

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Bacteria Phylum entries:", ifelse(is.null(bacteria_phylum), 0, nrow(bacteria_phylum)), "\n")
cat("Bacteria Family entries:", ifelse(is.null(bacteria_family), 0, nrow(bacteria_family)), "\n")
cat("Archaea Phylum entries:", ifelse(is.null(archaea_phylum), 0, nrow(archaea_phylum)), "\n")
cat("Archaea Family entries:", ifelse(is.null(archaea_family), 0, nrow(archaea_family)), "\n")
cat("Eukaryota Phylum entries:", ifelse(is.null(eukaryota_phylum), 0, nrow(eukaryota_phylum)), "\n")
cat("Eukaryota Family entries:", ifelse(is.null(eukaryota_family), 0, nrow(eukaryota_family)), "\n")

cat("\nStacked scatter plots saved successfully!\n")
cat("Files created:\n")
cat(paste("- ", output_file, "\n"))
cat(paste("- ", png_file, "\n"))
cat("\nPlots show taxa with novelty_factor > 1.0 AND overrepresentation_factor > 1.0\n")
cat("Yellow highlighted points represent top 10 taxa by each factor\n")
