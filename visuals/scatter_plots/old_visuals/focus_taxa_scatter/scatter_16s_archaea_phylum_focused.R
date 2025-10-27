#!/usr/bin/env Rscript
# 16S Archaea Phylum Focused Scatter Plot
# Created: 2025-08-18
# Purpose: Clean, focused scatter plot for 16S Archaea at Phylum level only

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "scatter_16s_archaea_phylum_focused.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Clean configuration for 16S Archaea Phylum analysis
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "merged_output", "16s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 38,
  plot_height = 20,  # Expanded height for phyla legend
  dpi = 300,
  top_n = 10,
  text_size = 13,
  size_range = c(8, 40),  # Large circles to fill white space
  color_palette = "professional",
  # Options: "professional" (default), "Set2", "Set3", "Dark2", "rainbow", "jama", "cosmic_hallmarks_light", "tol"
  use_weighted_ratio = TRUE  # TOGGLE: Set to FALSE to use simple ratios, TRUE for weighted ratios
)

# Verify critical paths exist
if (!dir.exists(config$data_dir)) {
  stop(paste("Data directory not found:", config$data_dir))
}
if (!dir.exists(config$ncbi_data_dir)) {
  stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
}
cat(paste("Data directory verified:", config$data_dir, "\n"))
cat(paste("NCBI data directory verified:", config$ncbi_data_dir, "\n"))

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Generate color palette based on config
get_color_palette <- function(n_colors) {
  palette_type <- config$color_palette
  
  if (palette_type == "professional") {
    # Professional, discrete color palette - maximally distinguishable muted tones
    professional_colors <- c(
      "#2E5984",  # Deep navy blue
      "#8B4513",  # Saddle brown  
      "#FF8C00",  # Dark orange
      "#F0A0C0",  # Light pink
      "#CD853F",  # Peru (orange-brown)
      "#2F4F4F",  # Dark slate gray
      "#DC143C",  # Crimson red
      "#B8860B",  # Dark goldenrod
      "#8FBC8F",  # Dark sea green
      "#A0522D"   # Sienna (reddish brown)
    )
    
    # If we need more than 10 colors, extend with contrasting variants
    if (n_colors > 10) {
      # Create highly contrasting additional colors
      additional_colors <- c(
        "#4682B4", "#D2691E", "#9ACD32", "#DA70D6", "#20B2AA",
        "#FF6347", "#4169E1", "#32CD32", "#FF1493", "#00CED1"
      )
      professional_colors <- c(professional_colors, additional_colors)
    }
    
    # Return the requested number of colors
    return(professional_colors[1:min(n_colors, length(professional_colors))])
    
  } else if (palette_type == "tol") {
    # Requires palr package: install.packages("palr")
    if (requireNamespace("palr", quietly = TRUE)) {
      colors <- palr::tol_colors(n_colors)  # Tol color palette
    } else {
      cat("Warning: palr package not found, falling back to professional\n")
      return(get_color_palette(n_colors))  # Recursive call with professional palette
    }
  } else {
    # For any other palette type, fall back to professional
    cat(paste("Warning: Unknown palette type", palette_type, ", using professional palette\n"))
    return(get_color_palette(n_colors))  # Recursive call with professional palette
  }
  
  return(colors[1:n_colors])
}

# Load and process 16S archaea phylum data
load_archaea_phylum_data <- function() {
  filepath <- file.path(config$data_dir, "16s_ncbi_merged_clean_phylum.csv")
  if (!file.exists(filepath)) {
    cat(paste("File not found:", filepath, "\n"))
    cat("Available files in data directory:\n")
    print(list.files(config$data_dir, pattern = "*.csv"))
    stop(paste("16S NCBI merged file not found:", filepath))
  }

  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")

  # Rename columns to match expected format - keep original CSV column names
  colnames(data)[1] <- "Taxon"
  # Keep original column names from CSV for consistency
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count

  # Filter for ARCHAEA only
  archaea_data <- data %>% filter(domain == "Archaea")
  cat(paste("Found", nrow(archaea_data), "archaeal phyla\n"))
  cat(paste("Total phyla in dataset:", nrow(data), "\n"))

  # Filter for phyla with meaningful data (Census OR NCBI data > 0 for visualization)
  # Include both matched and ncbi_only entries for archaea
  data_with_data <- archaea_data %>% filter(Census_OTU_Count > 0 | NCBI_Species_Count > 0)
  cat(paste("Found", nrow(data_with_data), "archaeal phyla with Census or NCBI data\n"))

  # For ratio calculations, we need both Census AND NCBI data > 0
  data_with_both <- data_with_data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  cat(paste("Found", nrow(data_with_both), "archaeal phyla with both Census and NCBI data for ratios\n"))

  # Use all data for visualization, but only calculate ratios for those with both
  data <- data_with_data

  # Since these are already phyla, use the Taxon name as the Phylum
  data$Phylum <- data$Taxon

  # Report final phylum distribution
  phylum_counts <- table(data$Phylum)
  cat("Phylum distribution:\n")
  print(phylum_counts)

  # Calculate key metrics directly from raw counts
  # Handle cases where denominators might be 0
  # Novelty Ratio: How many Census OTUs per NCBI Species (Census richness vs NCBI richness)
  data$Novelty_Ratio <- ifelse(data$NCBI_Species_Count > 0,
                               data$Census_OTU_Count / data$NCBI_Species_Count,
                               NA)
  # Coverage Factor: How many NCBI Species per Census OTU (NCBI coverage of Census diversity)
  data$Coverage_Factor <- ifelse(data$Census_OTU_Count > 0,
                                data$NCBI_Species_Count / data$Census_OTU_Count,
                                NA)

  # EXPERIMENTAL: Weighted ratios - easy to toggle on/off
  if (config$use_weighted_ratio) {
    cat("Using weighted ratios: ratio * log10(denominator + 1)\n")
    data$Weighted_Novelty <- data$Novelty_Ratio * log10(data$NCBI_Species_Count + 1)
    data$Weighted_Coverage <- data$Coverage_Factor * log10(data$Census_OTU_Count + 1)
    # Use weighted ratios for ranking
    ranking_novelty <- data$Weighted_Novelty
    ranking_coverage <- data$Weighted_Coverage
  } else {
    cat("Using simple ratios\n")
    # Use simple ratios for ranking
    ranking_novelty <- data$Novelty_Ratio
    ranking_coverage <- data$Coverage_Factor
  }

  # Circle sizing based on genome to isolate ratio (how many genomes per isolate)
  # Calculate directly from raw counts - don't use pre-calculated percentages
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))

  # Use square root transformation for balanced size differences
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)

  # Normalize to a larger range (3-20) to fill more white space
  min_size <- min(data$Circle_Size_Raw, na.rm = TRUE)
  max_size <- max(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 3 + 17 * (data$Circle_Size_Raw - min_size) / (max_size - min_size)

  # Filter for meaningful ratios (> 1.0) before ranking, excluding NA values
  high_novelty_data <- data[!is.na(ranking_novelty) & ranking_novelty > 1.0, ]
  high_coverage_data <- data[!is.na(ranking_coverage) & ranking_coverage > 1.0, ]

  cat(paste("Phyla with novelty ratio > 1.0:", nrow(high_novelty_data), "\n"))
  cat(paste("Phyla with coverage ratio > 1.0:", nrow(high_coverage_data), "\n"))

  # Identify top taxa using the selected ranking method, but only from those > 1.0
  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE

  # Mark top novelty (only from those with ratio > 1.0)
  if (nrow(high_novelty_data) > 0) {
    # Get novelty values for the filtered data and rank them
    high_novelty_indices <- which(data$Taxon %in% high_novelty_data$Taxon)
    high_novelty_values <- ranking_novelty[high_novelty_indices]
    novelty_ranks <- rank(-high_novelty_values)

    # Get top N taxa from the high novelty group
    top_n_count <- min(config$top_n, length(high_novelty_values))
    top_novelty_mask <- novelty_ranks <= top_n_count
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_mask]

    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
    cat(paste("Selected top", length(top_novelty_taxa), "novelty phyla:", paste(top_novelty_taxa, collapse = ", "), "\n"))

    # Debug: Show all novelty ratios > 1.0 sorted
    cat("All phyla with novelty > 1.0 (sorted by novelty ratio):\n")
    novelty_debug <- high_novelty_data[order(-high_novelty_data$Novelty_Ratio), c("Taxon", "Novelty_Ratio")]
    print(head(novelty_debug, 15))
  }

  # Mark top coverage (only from those with ratio > 1.0)
  if (nrow(high_coverage_data) > 0) {
    # Get coverage values for the filtered data and rank them
    high_coverage_indices <- which(data$Taxon %in% high_coverage_data$Taxon)
    high_coverage_values <- ranking_coverage[high_coverage_indices]
    coverage_ranks <- rank(-high_coverage_values)

    # Get top N taxa from the high coverage group
    top_n_count <- min(config$top_n, length(high_coverage_values))
    top_coverage_mask <- coverage_ranks <= top_n_count
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_mask]

    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
    cat(paste("Selected top", length(top_coverage_taxa), "coverage phyla:", paste(top_coverage_taxa, collapse = ", "), "\n"))
  }

  # Debug: Print top coverage and novelty taxa
  cat("\n=== TOP COVERAGE PHYLA ===\n")
  top_coverage <- data %>% filter(Is_Top_Coverage) %>% arrange(desc(Coverage_Factor)) %>%
    select(Taxon, Coverage_Factor, Census_OTU_Count, NCBI_Species_Count)
  print(top_coverage)

  cat("\n=== TOP NOVELTY PHYLA ===\n")
  top_novelty <- data %>% filter(Is_Top_Novelty) %>% arrange(desc(Novelty_Ratio)) %>%
    select(Taxon, Novelty_Ratio, Census_OTU_Count, NCBI_Species_Count)
  print(top_novelty)

  return(data)
}

# Generate summary tables
generate_summary <- function(data, type = "coverage") {
  col <- if (type == "coverage") "Coverage_Factor" else "Novelty_Ratio"
  data %>% arrange(desc(.data[[col]])) %>% head(config$top_n) %>%
    select(Taxon, Coverage_Factor, Novelty_Ratio, Census_OTU_Count, NCBI_Species_Count, Phylum) %>%
    mutate(across(c(Coverage_Factor, Novelty_Ratio), ~ round(.x, 2)))
}

# Create professional scatter plot
create_archaea_phylum_scatter <- function(data) {
  # Base plot with professional styling - wider limits for phylum level
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 100000),  # Range for phylum-level data (wider)
                  expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 10000),  # Range for phylum-level data (wider)
                  expand = expansion(mult = c(0.05, 0.05)))

  # First add background points (light grey) for all non-top taxa - these go behind colored points
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "lightgrey", alpha = 0.6, stroke = 0.3)
  }

  # Then add highlighted points with colors by phyla (on top of background)
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]

  # Show all top archaea phyla without additional filtering
  meaningful_data <- top_data

  cat(paste("Showing all top archaea phyla:", nrow(meaningful_data), "out of", nrow(top_data), "\n"))
  top_data <- meaningful_data

  if (nrow(top_data) > 0) {
    # Since we're working with phyla directly, each taxon IS a phylum
    # Show only the meaningful phyla in the legend
    top_phyla <- unique(top_data$Taxon)  # Use Taxon since each row is a phylum
    top_phyla <- sort(top_phyla)
    n_phyla <- length(top_phyla)

    cat(paste("Showing", n_phyla, "meaningful archaea phyla with individual colors:", paste(top_phyla, collapse = ", "), "\n"))

    # Generate colors using configurable palette
    colors <- get_color_palette(n_phyla)

    # Add black outline layer first
    p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.5)

    # Add colored fill layer on top - color by taxon (which is the phylum name)
    p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Taxon, levels = top_phyla)),
                       alpha = 0.9) +
      scale_color_manual(values = setNames(colors, top_phyla), name = "Archaeal Phyla",
                        guide = guide_legend(override.aes = list(size = 12, alpha = 1, shape = 15),
                                           keywidth = unit(3.0, "cm"), keyheight = unit(2.5, "cm")))
  }

  # Professional styling
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "darkblue", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "right",
          # Professional styling with larger, darker axes
          plot.title = element_text(size = 64, hjust = 0.5, face = "bold", color = "grey50"),
          axis.title = element_text(size = 58, face = "bold", color = "grey50"),
          axis.text = element_text(size = 42, face = "bold", color = "grey50"),
          # Legend styling - optimized for taller figure
          legend.text = element_text(size = 48),
          legend.title = element_text(size = 56, face = "bold", color = "grey50"),
          legend.key.size = unit(4.0, "cm"),
          legend.key.width = unit(4.5, "cm"),
          legend.margin = margin(l = 50, t = 20, b = 20),  # Added top/bottom margins
          panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)) +
    labs(title = "16S EukCensus vs NCBI: Archaea Phylum",
         x = "16S EukCensus OTU Count", y = "NCBI Species Count",
         caption = paste(nrow(data), "archaeal phyla"))

  # Add circle size legend in top-left corner
  # Use the actual plot limits instead of data range for consistent positioning
  x_range <- c(1, 100000)  # Match the scale limits for phylum
  y_range <- c(1, 10000)   # Match the scale limits for phylum

  legend_x <- 10^(log10(x_range[1]) + 0.05 * diff(log10(x_range)))  # Left side of plot
  legend_y_top <- 10^(log10(y_range[2]) - 0.05 * diff(log10(y_range)))  # Close to top edge

  # Add circle size legend title on the left
  p <- p + annotate("text", x = legend_x, y = legend_y_top,
                   label = "Genomes/Isolates Ratio", hjust = 0, vjust = 0,
                   size = 10, fontface = "bold", color = "grey30")

  # Add example circles with realistic genome/isolate ratios from the data
  # Use realistic values that match actual data: 10x, 50x, 100x
  circle_ratios <- c(10, 50, 100)
  circle_labels <- paste0(circle_ratios, "x")

  for (i in 1:3) {
    y_pos <- 10^(log10(legend_y_top) - (i * 0.12 * diff(log10(y_range))))

    # Use simple, intuitive sizing based on the ratio values
    # Small circle for 10x, medium for 50x, large for 100x
    if (i == 1) {
      scaled_size <- 8   # Small circle for 10x
    } else if (i == 2) {
      scaled_size <- 20  # Medium circle for 50x
    } else {
      scaled_size <- 35  # Large circle for 100x
    }

    p <- p + annotate("point", x = legend_x, y = y_pos,
                     size = scaled_size, color = "black", alpha = 0.8, stroke = 1)

    # Place labels to the right of circles (since legend is on left side)
    label_x <- 10^(log10(legend_x) + 0.08 * diff(log10(x_range)))
    p <- p + annotate("text", x = label_x, y = y_pos,
                     label = circle_labels[i], hjust = 0, vjust = 0.5,
                     size = 6, color = "black", fontface = "bold")
  }

  # Add text labels with ggrepel for phylum-level data
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_data) > 0) {
    # Create label text for all top data
    top_data$label_text <- paste0(top_data$Taxon, " (",
                                 ifelse(top_data$Is_Top_Novelty, paste0(round(top_data$Novelty_Ratio, 1), "x nov"), ""),
                                 ifelse(top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty, paste0(round(top_data$Coverage_Factor, 1), "x cov"), ""),
                                 ifelse(top_data$Is_Top_Novelty & top_data$Is_Top_Coverage, paste0(round(top_data$Novelty_Ratio, 1), "x"), ""),
                                 ")")

    # Single ggrepel call for ALL annotations with better spacing
    p <- p + ggrepel::geom_text_repel(
      data = top_data,
      aes(label = label_text),
      color = "black", size = config$text_size, fontface = "bold",
      # Better settings to reduce overlaps while keeping labels close to points
      max.overlaps = Inf,
      max.iter = 5000,
      max.time = 3,
      force = 5,           # Reduced force to keep labels closer to points
      force_pull = 0.5,    # Higher pull to keep labels near points
      box.padding = 0.8,   # Smaller padding for tighter spacing
      point.padding = 0.5, # Smaller point padding
      min.segment.length = 0,
      segment.color = "grey40",
      segment.size = 0.6,
      segment.alpha = 0.9,
      xlim = c(NA, NA),
      ylim = c(NA, NA),
      seed = 123  # Different seed for potentially better positioning
    )
  }

  return(p)
}

# Main function
main <- function() {
  cat("16S Archaea Phylum Focused Analysis\n")
  cat("==================================\n")

  # Load and process data
  data <- load_archaea_phylum_data()

  # Create plot
  plot <- create_archaea_phylum_scatter(data)

  # Save plot and summaries
  ggsave(file.path(config$output_dir, "16s_archaea_phylum_focused.png"), plot,
         width = config$plot_width, height = config$plot_height, dpi = config$dpi, bg = "white")

  write.csv(generate_summary(data, "coverage"),
            file.path(config$output_dir, "16s_archaea_phylum_top_coverage.csv"), row.names = FALSE)
  write.csv(generate_summary(data, "novelty"),
            file.path(config$output_dir, "16s_archaea_phylum_top_novelty.csv"), row.names = FALSE)

  cat(paste("✅ Analysis complete! Processed", nrow(data), "archaeal phyla\n"))
  cat(paste("   Top novelty taxa:", sum(data$Is_Top_Novelty), "\n"))
  cat(paste("   Top coverage taxa:", sum(data$Is_Top_Coverage), "\n"))
  cat("\nGenerated files:\n")
  cat("  - 16s_archaea_phylum_focused.png\n")
  cat("  - 16s_archaea_phylum_top_coverage.csv\n")
  cat("  - 16s_archaea_phylum_top_novelty.csv\n")
}

# Run analysis
if (!interactive()) main()
