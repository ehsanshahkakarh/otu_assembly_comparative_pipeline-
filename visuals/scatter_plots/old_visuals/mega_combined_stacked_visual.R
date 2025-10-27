#!/usr/bin/env Rscript
# Combined Mega Stacked Visual - 16S + 18S Scatter Plots
# Created: 2025-08-18
# Purpose: Create a unified mega-grid combining both 16S and 18S scatter plots

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "mega_combined_stacked_visual.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Configuration for mega visual
config <- list(
  data_dir_16s = file.path("..", "Eukcensus_merge", "merged_output", "16s_merged", "results"),
  data_dir_18s = file.path("..", "Eukcensus_merge", "merged_output", "18s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 120,  # Very wide for 3 columns to match individual plot scale
  plot_height = 40,  # Proportional height for 2 rows
  dpi = 300,
  top_n = 10,
  text_size = 6,     # Very small text for compact layout
  size_range = c(2, 10),  # Small circles for compact layout
  color_palette = "professional",
  use_weighted_ratio = TRUE
)

# Verify critical paths exist
if (!dir.exists(config$data_dir_16s)) {
  stop(paste("16S data directory not found:", config$data_dir_16s))
}
if (!dir.exists(config$data_dir_18s)) {
  stop(paste("18S data directory not found:", config$data_dir_18s))
}
if (!dir.exists(config$ncbi_data_dir)) {
  stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
}

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Master color palette for all taxa (both prokaryotes and eukaryotes)
get_master_color_palette <- function() {
  # Extended professional palette to cover all possible taxa
  master_colors <- c(
    "#2E5984",  "#8B4513",  "#FF8C00",  "#F0A0C0",  "#CD853F",
    "#2F4F4F",  "#DC143C",  "#B8860B",  "#8FBC8F",  "#A0522D",
    "#4682B4",  "#D2691E",  "#9ACD32",  "#DA70D6",  "#20B2AA",
    "#FF6347",  "#4169E1",  "#32CD32",  "#FF1493",  "#00CED1",
    "#FFD700",  "#FF69B4",  "#00FF7F",  "#8A2BE2",  "#FF4500",
    "#7B68EE",  "#00FA9A",  "#FF6B6B",  "#4ECDC4",  "#45B7D1",
    "#96CEB4",  "#FFEAA7",  "#DDA0DD",  "#98D8C8",  "#F7DC6F",
    "#BB8FCE",  "#85C1E9",  "#F8C471",  "#82E0AA",  "#F1948A",
    "#AED6F1",  "#A9DFBF",  "#F9E79F",  "#D7BDE2",  "#A3E4D7",
    "#FAD7A0",  "#D5A6BD",  "#AEB6BF",  "#D6EAF8",  "#EBDEF0"
  )
  
  return(master_colors)
}

# Simple archaea detection function
is_archaea <- function(taxon_name) {
  archaea_patterns <- c("Methanobrevibacter", "Methanococcus", "Methanospirillum",
                       "Thermococcus", "Pyrococcus", "Sulfolobus", "Thermoplasma",
                       "Halobacterium", "Methanosarcina", "Archaeoglobus")
  any(sapply(archaea_patterns, function(pattern) grepl(pattern, taxon_name, ignore.case = TRUE)))
}

# Standardize phylum names (consolidate subdivisions into main phyla)
standardize_phylum_names <- function(phylum_names) {
  # Map subdivision names to main phyla
  phylum_mapping <- c(
    "Methanobacteriota" = "Euryarchaeota",
    "Halobacteriota" = "Euryarchaeota",
    "Thermoplasmatota" = "Euryarchaeota"
  )

  # Apply mapping
  standardized <- phylum_names
  for (old_name in names(phylum_mapping)) {
    standardized[standardized == old_name] <- phylum_mapping[old_name]
  }

  return(standardized)
}

# Load 16S data function (simplified from mega_16s_stacked_visual.R)
load_16s_data <- function(level, domain) {
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)
  
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")
  
  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  
  # Filter by domain
  if ("domain" %in% colnames(data)) {
    data <- data %>% filter(domain == domain)
  } else {
    # Fallback to pattern matching for archaea
    if (domain == "Archaea") {
      archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
      data <- data %>% filter(Taxon %in% archaea_taxa)
    } else {
      archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
      data <- data %>% filter(!Taxon %in% archaea_taxa)
    }
  }
  
  # Filter for meaningful data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  
  # Add phylum information
  if (level == "phylum") {
    data$Phylum <- data$Taxon
  } else {
    data$Phylum <- "Unknown"  # Simplified for combined visual
  }

  # Standardize phylum names (consolidate Methanobacteriota -> Euryarchaeota, etc.)
  data$Phylum <- standardize_phylum_names(data$Phylum)
  
  # Calculate metrics
  data$Novelty_Ratio <- data$Census_OTU_Count / data$NCBI_Species_Count
  data$Coverage_Factor <- data$NCBI_Species_Count / data$Census_OTU_Count
  
  # Weighted ratios if enabled
  if (config$use_weighted_ratio) {
    data$Weighted_Novelty <- data$Novelty_Ratio * log10(data$NCBI_Species_Count + 1)
    data$Weighted_Coverage <- data$Coverage_Factor * log10(data$Census_OTU_Count + 1)
    ranking_novelty <- data$Weighted_Novelty
    ranking_coverage <- data$Weighted_Coverage
  } else {
    ranking_novelty <- data$Novelty_Ratio
    ranking_coverage <- data$Coverage_Factor
  }
  
  # Circle sizing
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))
  
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  min_size <- min(data$Circle_Size_Raw, na.rm = TRUE)
  max_size <- max(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 1 + 4 * (data$Circle_Size_Raw - min_size) / (max_size - min_size)
  
  # Identify top taxa
  high_novelty_data <- data[ranking_novelty > 1.0, ]
  high_coverage_data <- data[ranking_coverage > 1.0, ]
  
  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE
  
  # Mark top novelty
  if (nrow(high_novelty_data) > 0) {
    high_novelty_indices <- which(data$Taxon %in% high_novelty_data$Taxon)
    high_novelty_values <- ranking_novelty[high_novelty_indices]
    novelty_ranks <- rank(-high_novelty_values)
    top_n_count <- min(config$top_n, length(high_novelty_values))
    top_novelty_mask <- novelty_ranks <= top_n_count
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_mask]
    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
  }
  
  # Mark top coverage
  if (nrow(high_coverage_data) > 0) {
    high_coverage_indices <- which(data$Taxon %in% high_coverage_data$Taxon)
    high_coverage_values <- ranking_coverage[high_coverage_indices]
    coverage_ranks <- rank(-high_coverage_values)
    top_n_count <- min(config$top_n, length(high_coverage_values))
    top_coverage_mask <- coverage_ranks <= top_n_count
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_mask]
    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
  }
  
  return(data)
}

# Load 18S data function (simplified)
load_18s_data <- function(level) {
  filename <- paste0("18s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_18s, filename)
  
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")
  
  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  
  # Filter for meaningful data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  
  # Add division information (simplified)
  data$Division <- "Eukaryota"  # Simplified for combined visual
  
  # Calculate metrics (same as 16S)
  data$Novelty_Ratio <- data$Census_OTU_Count / data$NCBI_Species_Count
  data$Coverage_Factor <- data$NCBI_Species_Count / data$Census_OTU_Count
  
  # Weighted ratios if enabled
  if (config$use_weighted_ratio) {
    data$Weighted_Novelty <- data$Novelty_Ratio * log10(data$NCBI_Species_Count + 1)
    data$Weighted_Coverage <- data$Coverage_Factor * log10(data$Census_OTU_Count + 1)
    ranking_novelty <- data$Weighted_Novelty
    ranking_coverage <- data$Weighted_Coverage
  } else {
    ranking_novelty <- data$Novelty_Ratio
    ranking_coverage <- data$Coverage_Factor
  }
  
  # Circle sizing
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))
  
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  min_size <- min(data$Circle_Size_Raw, na.rm = TRUE)
  max_size <- max(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 1 + 4 * (data$Circle_Size_Raw - min_size) / (max_size - min_size)
  
  # Identify top taxa (same logic as 16S)
  high_novelty_data <- data[ranking_novelty > 1.0, ]
  high_coverage_data <- data[ranking_coverage > 1.0, ]
  
  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE
  
  # Mark top novelty
  if (nrow(high_novelty_data) > 0) {
    high_novelty_indices <- which(data$Taxon %in% high_novelty_data$Taxon)
    high_novelty_values <- ranking_novelty[high_novelty_indices]
    novelty_ranks <- rank(-high_novelty_values)
    top_n_count <- min(config$top_n, length(high_novelty_values))
    top_novelty_mask <- novelty_ranks <= top_n_count
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_mask]
    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
  }
  
  # Mark top coverage
  if (nrow(high_coverage_data) > 0) {
    high_coverage_indices <- which(data$Taxon %in% high_coverage_data$Taxon)
    high_coverage_values <- ranking_coverage[high_coverage_indices]
    coverage_ranks <- rank(-high_coverage_values)
    top_n_count <- min(config$top_n, length(high_coverage_values))
    top_coverage_mask <- coverage_ranks <= top_n_count
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_mask]
    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
  }
  
  return(data)
}

# Create individual scatter plot (minimal version for combined visual)
create_minimal_scatter <- function(data, dataset_type, level, domain = NULL) {
  # Base plot
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.05, 0.05)))

  # Get top data for highlighting
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]

  # Filter meaningful data
  if (config$use_weighted_ratio) {
    meaningful_data <- top_data[top_data$Weighted_Novelty > 1.0 | top_data$Weighted_Coverage > 1.0, ]
  } else {
    meaningful_data <- top_data[top_data$Novelty_Ratio > 1.0 | top_data$Coverage_Factor > 1.0, ]
  }

  # Apply exclusions for cleaner plots
  if (dataset_type == "16S") {
    excluded_taxa <- c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
    if ("Phylum" %in% colnames(meaningful_data)) {
      meaningful_data <- meaningful_data[!meaningful_data$Phylum %in% excluded_taxa, ]
    }
  }

  if (nrow(meaningful_data) > 0) {
    # Add black outline layer
    p <- p + geom_point(data = meaningful_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.2)

    # Add colored fill layer (simple gray for combined visual)
    p <- p + geom_point(data = meaningful_data,
                       aes(size = Circle_Size * 0.8),
                       color = "#4682B4", alpha = 0.7)
  }

  # Background points (minimal)
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "white", alpha = 0.0)
  }

  # Minimal styling for stacking
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "darkblue", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(
      # Remove titles and axis labels for stacking
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 10, face = "bold", color = "grey50"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
      plot.margin = margin(2, 2, 2, 2)
    )

  # Add minimal text annotations for top taxa (very selective)
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(meaningful_data) > 0) {
    # Only annotate top 3 taxa to avoid overcrowding
    if (nrow(meaningful_data) > 3) {
      if (config$use_weighted_ratio) {
        top_indices <- order(-meaningful_data$Weighted_Novelty)[1:3]
      } else {
        top_indices <- order(-meaningful_data$Novelty_Ratio)[1:3]
      }
      meaningful_data <- meaningful_data[top_indices, ]
    }

    # Create simplified labels
    meaningful_data$label_text <- paste0(meaningful_data$Taxon)

    p <- p + ggrepel::geom_text_repel(
      data = meaningful_data,
      aes(label = label_text),
      color = "black", size = config$text_size, fontface = "bold",
      max.overlaps = Inf,
      max.iter = 3000,
      force = 5,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.linewidth = 0.2,
      segment.alpha = 0.6,
      seed = 42
    )
  }

  return(p)
}

# Classify phyla by domain (for 16S data)
classify_phyla_by_domain <- function(phyla) {
  # Known archaeal phyla
  archaeal_phyla <- c("Euryarchaeota", "Nitrososphaerota", "Thermoproteota",
                     "Nanoarchaeota", "Promethearchaeota", "Halobacteriota",
                     "Thermoplasmatota", "Methanobacteriota")

  # Classify each phylum
  domains <- ifelse(phyla %in% archaeal_phyla, "Archaea", "Bacteria")
  return(domains)
}

# Main function to create combined mega visual
main <- function() {
  cat("Combined Mega Stacked Visual Creation\n")
  cat("====================================\n")

  # Define the grid structure (removed genus)
  levels <- c("phylum", "family")
  domains_16s <- c("Bacteria", "Archaea")

  # Create all plots
  plots <- list()

  cat("Creating 16S scatter plots...\n")
  for (level in levels) {
    for (domain in domains_16s) {
      cat(paste("Loading", domain, level, "data...\n"))
      data <- load_16s_data(level, domain)
      plot_key <- paste("16S", domain, level, sep = "_")
      plots[[plot_key]] <- create_minimal_scatter(data, "16S", level, domain)
    }
  }

  cat("Creating 18S scatter plots...\n")
  for (level in levels) {
    cat(paste("Loading eukaryote", level, "data...\n"))
    data <- load_18s_data(level)
    plot_key <- paste("18S", "Eukaryotes", level, sep = "_")
    plots[[plot_key]] <- create_minimal_scatter(data, "18S", level)
  }

  # Arrange plots in mega grid
  cat("Arranging plots in mega grid...\n")

  # Create headers
  dataset_headers <- list(
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "16S rRNA",
                       size = 20, fontface = "bold", color = "grey50") + theme_void(),
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "18S rRNA",
                       size = 20, fontface = "bold", color = "grey50") + theme_void()
  )

  domain_headers <- list(
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Bacteria",
                       size = 16, fontface = "bold", color = "grey50") + theme_void(),
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Archaea",
                       size = 16, fontface = "bold", color = "grey50") + theme_void(),
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Eukaryotes",
                       size = 16, fontface = "bold", color = "grey50") + theme_void()
  )

  level_headers <- list(
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Phyla",
                       size = 12, fontface = "bold", color = "grey50", angle = 90) + theme_void(),
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Family",
                       size = 12, fontface = "bold", color = "grey50", angle = 90) + theme_void()
  )

  # Create the mega grid structure
  # Top header row
  top_header <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    dataset_headers[[1]],     # 16S header (spans 2 columns)
    dataset_headers[[2]],     # 18S header (spans 1 column)
    ncol = 3, rel_widths = c(0.1, 2, 1)
  )

  # Domain header row
  domain_header <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    domain_headers[[1]],      # Bacteria
    domain_headers[[2]],      # Archaea
    domain_headers[[3]],      # Eukaryotes
    ncol = 4, rel_widths = c(0.1, 1, 1, 1)
  )

  # Data rows
  data_rows <- list()
  for (i in 1:length(levels)) {
    level <- levels[i]
    bacteria_plot <- plots[[paste("16S", "Bacteria", level, sep = "_")]]
    archaea_plot <- plots[[paste("16S", "Archaea", level, sep = "_")]]
    eukaryote_plot <- plots[[paste("18S", "Eukaryotes", level, sep = "_")]]

    data_rows[[i]] <- plot_grid(
      level_headers[[i]],
      bacteria_plot,
      archaea_plot,
      eukaryote_plot,
      ncol = 4, rel_widths = c(0.1, 1, 1, 1)
    )
  }

  # Combine all rows
  final_grid <- plot_grid(
    top_header,
    domain_header,
    data_rows[[1]],  # Phylum row
    data_rows[[2]],  # Family row
    ncol = 1, rel_heights = c(0.08, 0.06, 1, 1)
  )

  # Add overall title and axis labels
  title <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Comprehensive rRNA Database Coverage: EukCensus vs NCBI",
             size = 24, fontface = "bold", color = "grey50") +
    theme_void()

  x_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Environmental Count",
             size = 18, fontface = "bold", color = "grey50") +
    theme_void()

  y_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "GenBank Species Count",
             size = 18, fontface = "bold", color = "grey50", angle = 90) +
    theme_void()

  # Final assembly with titles
  complete_plot <- plot_grid(
    title,
    plot_grid(
      y_label,
      final_grid,
      ncol = 2, rel_widths = c(0.04, 1)
    ),
    x_label,
    ncol = 1, rel_heights = c(0.06, 1, 0.04)
  )

  # Save the mega visual
  output_file <- file.path(config$output_dir, "mega_combined_stacked_visual.png")
  cat(paste("Saving combined mega visual to:", output_file, "\n"))

  ggsave(output_file, complete_plot,
         width = config$plot_width, height = config$plot_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE)

  cat("✅ Combined mega visual creation complete!\n")
  cat(paste("   Generated:", output_file, "\n"))
  cat(paste("   Dimensions:", config$plot_width, "x", config$plot_height, "inches\n"))
  cat("   Layout: 16S (Bacteria + Archaea) + 18S (Eukaryotes) × Phylum/Family/Genus\n")
}

# Run the mega visual creation
if (!interactive()) main()
