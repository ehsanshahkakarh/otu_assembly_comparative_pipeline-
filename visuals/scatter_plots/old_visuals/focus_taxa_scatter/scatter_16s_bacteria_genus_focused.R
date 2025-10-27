#!/usr/bin/env Rscript
# 16S Bacteria Genus Focused Scatter Plot
# Created: 2025-08-18
# Purpose: Clean, focused scatter plot for 16S Bacteria at Genus level only

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "scatter_16s_bacteria_genus_focused.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Clean configuration for 16S Bacteria Genus analysis
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "merged_output", "16s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 38,
  plot_height = 15,
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

# Simple archaea detection function
is_archaea <- function(taxon_name) {
  archaea_patterns <- c("Methanobrevibacter", "Methanococcus", "Methanospirillum", 
                       "Thermococcus", "Pyrococcus", "Sulfolobus", "Thermoplasma",
                       "Halobacterium", "Methanosarcina", "Archaeoglobus")
  any(sapply(archaea_patterns, function(pattern) grepl(pattern, taxon_name, ignore.case = TRUE)))
}

# Get phylum mapping for genus level
get_phylum_mapping <- function() {
  file_path <- file.path(config$ncbi_data_dir, "ncbi_genus_counts.csv")
  if (!file.exists(file_path)) {
    cat(paste("Warning: NCBI genus file not found:", file_path, "\n"))
    return(NULL)
  }

  data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat(paste("Loaded NCBI genus file with", nrow(data), "rows\n"))
  cat(paste("Columns:", paste(colnames(data), collapse = ", "), "\n"))
  
  # Check if we have the required columns
  if (!"lineage" %in% colnames(data)) {
    cat("Error: 'lineage' column not found in NCBI genus file\n")
    return(NULL)
  }
  if (!"lineage_ranks" %in% colnames(data)) {
    cat("Error: 'lineage_ranks' column not found in NCBI genus file\n")
    return(NULL)
  }
  
  # Check what the genus name column is called
  genus_col <- NULL
  for (col in c("genus_name", "genus", "taxon", "Taxon")) {
    if (col %in% colnames(data)) {
      genus_col <- col
      break
    }
  }
  
  if (is.null(genus_col)) {
    cat("Error: Could not find genus name column. Available columns:", paste(colnames(data), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat(paste("Using", genus_col, "as genus name column\n"))
  
  # Extract phylum from lineage using rank information
  phylum_map <- setNames(character(nrow(data)), data[[genus_col]])
  successful_mappings <- 0
  
  for (i in 1:nrow(data)) {
    lineage <- data$lineage[i]
    lineage_ranks <- data$lineage_ranks[i]
    
    if (!is.na(lineage) && lineage != "" && !is.na(lineage_ranks) && lineage_ranks != "") {
      lineage_parts <- strsplit(lineage, ";")[[1]]
      rank_parts <- strsplit(lineage_ranks, ";")[[1]]
      
      # Find phylum using rank information
      phylum_idx <- which(rank_parts == "phylum")[1]
      if (!is.na(phylum_idx) && phylum_idx <= length(lineage_parts)) {
        phylum_name <- trimws(lineage_parts[phylum_idx])
        phylum_map[data[[genus_col]][i]] <- phylum_name
        successful_mappings <- successful_mappings + 1
        
        # Show first few examples
        if (successful_mappings <= 3) {
          cat(paste("Example mapping:", data[[genus_col]][i], "->", phylum_name, "\n"))
        }
      } else {
        phylum_map[data[[genus_col]][i]] <- "Other"
      }
    } else {
      phylum_map[data[[genus_col]][i]] <- "Other"
    }
  }
  
  # Remove empty mappings
  phylum_map <- phylum_map[phylum_map != ""]
  cat(paste("Successfully mapped", successful_mappings, "genera to phyla\n"))
  cat(paste("Unique phyla found:", length(unique(phylum_map[phylum_map != "Other"])), "\n"))
  
  return(phylum_map)
}

# Generate color palette based on config
get_color_palette <- function(n_colors) {
  palette_type <- config$color_palette
  
  if (palette_type == "professional") {
    # Professional, discrete color palette - maximally distinguishable muted tones
    professional_colors <- c(
      "#2E5984",  # Deep navy blue
      "#8B4513",  # Saddle brown  
      "#FF8C00",  # Dark orange (for Chlorobiota - more orangey)
      "#F0A0C0",  # Light pink (for Gemmatimonadota)
      "#CD853F",  # Peru (orange-brown)
      "#2F4F4F",  # Dark slate gray
      "#DC143C",  # Crimson red (for Pseudomonadota)
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

# Load and process 16S bacteria genus data
load_bacteria_genus_data <- function() {
  filepath <- file.path(config$data_dir, "16s_ncbi_merged_clean_genus.csv")
  if (!file.exists(filepath)) {
    cat(paste("File not found:", filepath, "\n"))
    cat("Available files in data directory:\n")
    print(list.files(config$data_dir, pattern = "*.csv"))
    stop(paste("16S NCBI merged file not found:", filepath))
  }

  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")

  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  names(data)[names(data) == "census_otu_count"] <- "Census_OTU_Count"
  names(data)[names(data) == "ncbi_species_count"] <- "NCBI_Species_Count"
  names(data)[names(data) == "ncbi_genome_count"] <- "NCBI_Genome_Count"
  names(data)[names(data) == "isolate_count"] <- "Isolate_Count"

  # Filter for bacteria only (exclude archaea) - use domain column if available
  if ("domain" %in% colnames(data)) {
    data <- data %>% filter(domain == "Bacteria")
    cat(paste("Found", nrow(data), "bacteria genera after filtering by domain column\n"))
  } else {
    # Fallback to pattern matching if domain column not available
    archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
    data <- data %>% filter(!Taxon %in% archaea_taxa)
    cat(paste("Found", nrow(data), "bacteria genera after filtering archaea by pattern matching\n"))
  }

  # Filter for genera with both Census and NCBI data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  cat(paste("Found", nrow(data), "bacteria genera with both Census and NCBI data\n"))

  # Add phylum mapping
  phylum_map <- get_phylum_mapping()
  if (!is.null(phylum_map)) {
    data$Phylum <- ifelse(data$Taxon %in% names(phylum_map), phylum_map[data$Taxon], "Other")
  } else {
    data$Phylum <- "Other"
  }

  # Calculate key metrics
  data$Novelty_Ratio <- data$Census_OTU_Count / data$NCBI_Species_Count
  data$Coverage_Factor <- data$NCBI_Species_Count / data$Census_OTU_Count

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

  # Large circle sizing to fill white space effectively
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))

  # Use square root transformation for balanced size differences
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  
  # Normalize to a larger range (3-20) to fill more white space
  min_size <- min(data$Circle_Size_Raw, na.rm = TRUE)
  max_size <- max(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 3 + 17 * (data$Circle_Size_Raw - min_size) / (max_size - min_size)

  # Identify top taxa using the selected ranking method
  data$Is_Top_Novelty <- rank(-ranking_novelty) <= config$top_n
  data$Is_Top_Coverage <- rank(-ranking_coverage) <= config$top_n

  # Debug: Print top coverage and novelty taxa
  cat("\n=== TOP COVERAGE GENERA ===\n")
  top_coverage <- data %>% filter(Is_Top_Coverage) %>% arrange(desc(Coverage_Factor)) %>%
    select(Taxon, Coverage_Factor, Census_OTU_Count, NCBI_Species_Count)
  print(top_coverage)

  cat("\n=== TOP NOVELTY GENERA ===\n")
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
create_bacteria_genus_scatter <- function(data) {
  # Base plot with professional styling
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 1000),
                  expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 1000),
                  expand = expansion(mult = c(0.05, 0.05)))

  # Remove background points - only show highlighted taxa
  # Highlighted points with colors by phyla
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  if (nrow(top_data) > 0) {
    # Show only phyla of the annotated (top) taxa in the legend
    top_phyla <- unique(top_data$Phylum[top_data$Phylum != "Other"])
    top_phyla <- sort(top_phyla)
    n_phyla <- length(top_phyla)

    # Generate colors using configurable palette
    colors <- get_color_palette(n_phyla)

    # Add black outline layer first
    p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.5)

    # Add colored fill layer on top
    p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Phylum, levels = top_phyla)),
                       alpha = 0.9) +
      scale_color_manual(values = setNames(colors, top_phyla), name = "Phyla",
                        guide = guide_legend(override.aes = list(size = 12, alpha = 1, shape = 15),
                                           keywidth = unit(3.0, "cm"), keyheight = unit(2.5, "cm")))
  }

  # Background points (grey) - large to fill white space
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "white", alpha = 0.0)
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
          # Legend styling
          legend.text = element_text(size = 48),
          legend.title = element_text(size = 56, face = "bold", color = "grey50"),
          legend.key.size = unit(4.0, "cm"),
          legend.key.width = unit(4.5, "cm"),
          legend.margin = margin(l = 50),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)) +
    labs(title = "16S EukCensus vs NCBI: Bacteria Genus",
         x = "16S EukCensus OTU Count", y = "NCBI Species Count",
         caption = paste(nrow(data), "bacteria genera"))

  # Add circle size legend at the top right corner
  # Use the actual plot limits instead of data range
  x_range <- c(1, 10000)  # Match the scale limits
  y_range <- c(1, 1000)   # Match the scale limits

  legend_x <- 10^(log10(x_range[2]) - 0.25 * diff(log10(x_range)))  # More towards center from right edge
  legend_y_top <- 10^(log10(y_range[2]) - 0.06 * diff(log10(y_range)))  # A little higher from top edge

  # Add circle size legend title on the right
  p <- p + annotate("text", x = legend_x, y = legend_y_top,
                   label = "Genomes/Isolates Ratio", hjust = 1, vjust = 0,
                   size = 10, fontface = "bold", color = "grey30")

  # Add example circles with realistic genome/isolate ratios from the data
  # Use realistic values that match actual data: 10x, 50x, 100x
  circle_ratios <- c(10, 50, 100)
  circle_labels <- paste0(circle_ratios, "x")

  for (i in 1:3) {
    y_pos <- 10^(log10(legend_y_top) - (i * 0.12 * diff(log10(y_range))))

    # Use simple, intuitive sizing based on the ratio values
    # Small circle for 1x, medium for 50x, large for 150x
    if (i == 1) {
      scaled_size <- 8   # Small circle for 1x
    } else if (i == 2) {
      scaled_size <- 20  # Medium circle for 50x
    } else {
      scaled_size <- 35  # Large circle for 150x
    }

    p <- p + annotate("point", x = legend_x, y = y_pos,
                     size = scaled_size, color = "black", alpha = 0.8, stroke = 1)

    # Place labels to the left of circles (since legend is on right side)
    label_x <- 10^(log10(legend_x) - 0.08 * diff(log10(x_range)))
    p <- p + annotate("text", x = label_x, y = y_pos,
                     label = circle_labels[i], hjust = 1, vjust = 0.5,
                     size = 6, color = "black", fontface = "bold")
  }

  # Add ALL text labels - use single ggrepel call to ensure all 20 annotations appear
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_data) > 0) {
    # Create label text for all top data
    top_data$label_text <- paste0(top_data$Taxon, " (",
                                 ifelse(top_data$Is_Top_Novelty, paste0(round(top_data$Novelty_Ratio, 1), "x nov"), ""),
                                 ifelse(top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty, paste0(round(top_data$Coverage_Factor, 1), "x cov"), ""),
                                 ifelse(top_data$Is_Top_Novelty & top_data$Is_Top_Coverage, paste0(round(top_data$Novelty_Ratio, 1), "x"), ""),
                                 ")")

    # Single ggrepel call for ALL annotations to guarantee they all appear
    p <- p + ggrepel::geom_text_repel(
      data = top_data,
      aes(label = label_text),
      color = "black", size = config$text_size, fontface = "bold",
      # Aggressive settings to ensure ALL labels appear
      max.overlaps = Inf,
      max.iter = 10000,  # More iterations to find good positions
      max.time = 5,      # More time to solve positioning
      force = 20,        # High force to push labels apart
      force_pull = 0.1,  # Low pull to allow spreading into white space
      box.padding = 1.5,
      point.padding = 1.0,
      min.segment.length = 0,  # Always show segments
      segment.color = "grey40",
      segment.linewidth = 0.6,
      segment.alpha = 0.9,
      # Allow labels to spread into the middle white space
      xlim = c(NA, NA),  # No x limits
      ylim = c(NA, NA),  # No y limits
      # Use different seeds and directions for different groups
      seed = 42
    )
  }

  return(p)
}

# Main function
main <- function() {
  cat("16S Bacteria Genus Focused Analysis\n")
  cat("==================================\n")

  # Load and process data
  data <- load_bacteria_genus_data()

  # Create plot
  plot <- create_bacteria_genus_scatter(data)

  # Save plot and summaries
  ggsave(file.path(config$output_dir, "16s_bacteria_genus_focused.png"), plot,
         width = config$plot_width, height = config$plot_height, dpi = config$dpi, bg = "white")

  write.csv(generate_summary(data, "coverage"),
            file.path(config$output_dir, "16s_bacteria_genus_top_coverage.csv"), row.names = FALSE)
  write.csv(generate_summary(data, "novelty"),
            file.path(config$output_dir, "16s_bacteria_genus_top_novelty.csv"), row.names = FALSE)

  cat(paste("✅ Analysis complete! Processed", nrow(data), "bacteria genera\n"))
  cat(paste("   Top novelty taxa:", sum(data$Is_Top_Novelty), "\n"))
  cat(paste("   Top coverage taxa:", sum(data$Is_Top_Coverage), "\n"))
  cat("\nGenerated files:\n")
  cat("  - 16s_bacteria_genus_focused.png\n")
  cat("  - 16s_bacteria_genus_top_coverage.csv\n")
  cat("  - 16s_bacteria_genus_top_novelty.csv\n")
}

# Run analysis
if (!interactive()) main()
