#!/usr/bin/env Rscript
# Comprehensive Mega Stacked Visual - 16S Bacteria + 16S Archaea + 18S Eukaryota
# Created: 2025-10-26
# Purpose: Create a unified 3-column mega-grid with all three domains side-by-side

# ============================================================================
# COLOR EXCLUSION CONFIGURATION
# ============================================================================
# Specify colors to exclude from each domain's palette (without editing YAML)
# Use exact hex codes from the YAML file to exclude specific colors
# Leave empty vectors c() to use all colors from the palette

EXCLUDED_COLORS <- list(
  bacteria = c(),     # e.g., c("#4c9b34", "#72c859") to exclude first two greens
  archaea = c(),      # e.g., c("#f51b7f") to exclude bright pink
  eukaryota = c()     # e.g., c("#416b7d", "#69c1d4") to exclude first two blues
)

# Alternative: Exclude by position/index instead of hex codes
# Set to TRUE to use position-based exclusion instead of hex-based
USE_POSITION_EXCLUSION <- FALSE

EXCLUDED_POSITIONS <- list(
  bacteria = c(),     # e.g., c(1, 2) to exclude first two colors
  archaea = c(),      # e.g., c(1) to exclude first color
  eukaryota = c()     # e.g., c(1, 2) to exclude first two colors
)

# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(grid)
  library(cowplot)
  library(ggrepel)
  library(yaml)
})

# Shared helpers for deterministic, cross-figure taxon color assignment.
source("../../shared_config/color_mapping_functions.R")

# ============================================================================
# FIXED COLOR ASSIGNMENT - No repeats, uses all available colors
# ============================================================================

# Global map: taxon name → color (for consistency across plots)
TAXON_COLORS <- list()

# Simple function: assign color to a taxon
assign_taxon_color <- function(taxon, domain, color_config) {
  # Already seen this taxon? Return its color
  if (taxon %in% names(TAXON_COLORS)) {
    return(TAXON_COLORS[[taxon]])
  }

  # Check if taxon has hardcoded color in YAML
  color <- NULL
  if (domain == "Bacteria" && taxon %in% names(color_config$bacteria_colors)) {
    color <- color_config$bacteria_colors[[taxon]]
  } else if (domain == "Archaea" && taxon %in% names(color_config$archaea_colors)) {
    color <- color_config$archaea_colors[[taxon]]
  } else if (domain == "Eukaryota" && taxon %in% names(color_config$eukaryota_colors)) {
    color <- color_config$eukaryota_colors[[taxon]]
  }

  # Not hardcoded? Get from cross-domain pool (ALL colors from other domains)
  if (is.null(color)) {
    # Build pool: ALL colors from other domains (hardcoded + fallback)
    if (domain == "Bacteria") {
      # For Bacteria: use Archaea + Eukaryota colors
      pool <- c(
        unlist(color_config$archaea_colors),      # Archaea hardcoded
        unlist(color_config$eukaryota_colors),    # Eukaryota hardcoded
        unlist(color_config$fallback_colors$archaea),
        unlist(color_config$fallback_colors$eukaryota)
      )
    } else if (domain == "Eukaryota") {
      # For Eukaryota: use Bacteria + Archaea colors
      pool <- c(
        unlist(color_config$bacteria_colors),     # Bacteria hardcoded
        unlist(color_config$archaea_colors),      # Archaea hardcoded
        unlist(color_config$fallback_colors$bacteria),
        unlist(color_config$fallback_colors$archaea)
      )
    } else {  # Archaea
      # For Archaea: use Bacteria + Eukaryota colors
      pool <- c(
        unlist(color_config$bacteria_colors),     # Bacteria hardcoded
        unlist(color_config$eukaryota_colors),    # Eukaryota hardcoded
        unlist(color_config$fallback_colors$bacteria),
        unlist(color_config$fallback_colors$eukaryota)
      )
    }

    # Remove duplicates from pool
    pool <- unique(pool)

    # Get all colors already used
    used_colors <- unlist(TAXON_COLORS)

    # Find unused colors in the pool
    available_colors <- setdiff(pool, used_colors)

    # If we have available colors, use the first one
    # Otherwise, fall back to round-robin through the pool
    if (length(available_colors) > 0) {
      color <- available_colors[1]
    } else {
      # All colors used, start recycling
      index <- (length(TAXON_COLORS) %% length(pool)) + 1
      color <- pool[index]
    }
  }

  # Save mapping
  TAXON_COLORS[[taxon]] <<- color
  return(color)
}

# Initialize (just clear the list)
init_color_registry <- function() {
  TAXON_COLORS <<- list()
}

# ============================================================================
# Load shared taxonomic color mapping configuration
# This ensures consistency across all visualization types
load_shared_color_config <- function() {
  config_path <- file.path("..", "..", "shared_config", "taxonomic_color_mapping.yaml")
  if (!file.exists(config_path)) {
    stop("Shared color config not found at: ", config_path)
  }
  yaml::read_yaml(config_path)
}

# Configuration constants
PLOT_CONFIG <- list(
  # Plot dimensions and layout
  plot_width = 54,
  plot_height = 24,
  dpi = 300,

  # Data filtering
  top_n = 10,

  # Visual styling
  text_size = 11.5,  # Increased from 11 for better readability
  size_range = c(10, 22),

  # Circle aesthetics
  circle_shape = 21,
  circle_stroke = 0.6,
  circle_alpha = 0.9,
  bg_alpha = 0.3,

  # Grid layout
  col_spacing = 0.02,
  row_spacing = 0.05,
  legend_width_ratio = 0.15,
  main_plot_ratio = 0.85
)

# Layout configuration
LAYOUT_CONFIG <- list(
  col_widths = c(0.1, 1, 0.02, 1, 0.02, 1),  # Label, Bacteria, gap, Archaea, gap, Eukaryota
  row_heights = c(0.1, 1, 0.05, 1)  # Header, phylum, gap, family
)

# Function to calculate circle sizes based on isolate percentage
calculate_circle_size <- function(isolate_count, ncbi_genome_count) {
  isolate_percentage <- ifelse(ncbi_genome_count > 0, (isolate_count / ncbi_genome_count) * 100, 0)

  # Inverted scale: larger circles = poorly cultured (low isolate percentage)
  circle_size <- case_when(
    isolate_percentage == 0 ~ 25,
    isolate_percentage < 10 ~ 20,
    isolate_percentage < 50 ~ 15,
    TRUE ~ 10
  )

  return(circle_size)
}

# Function to create isolate percentage legend
create_isolate_legend <- function() {
  legend_data <- data.frame(
    isolate_percentage = c(0, 50, 100),
    circle_size = c(25, 15, 10),
    x = rep(1, 3),
    y = c(2, 1, 0.5),
    stringsAsFactors = FALSE
  )

  ggplot(legend_data, aes(x = x, y = y)) +
    geom_point(aes(size = circle_size), color = "black", fill = "gray",
               shape = PLOT_CONFIG$circle_shape, stroke = PLOT_CONFIG$circle_stroke, alpha = 0.8) +
    scale_size_identity() +
    theme_void() +
    xlim(0.5, 1.5) + ylim(0, 3)
}



# Helper function to extract phyla/divisions that appear in domain data (top taxa only)
extract_domain_phyla <- function(all_data, domain) {
  domain_phyla <- c()
  threshold <- 1.0  # Same threshold as used in plots
  top_n <- 10      # Same top_n as used in plots

  if (domain == "Eukaryota") {
    # For Eukaryota: only divisions that are in top 10 novelty or top 10 overrepresentation
    data_key <- paste("18S", domain, "phylum", sep = "_")
    if (data_key %in% names(all_data) && nrow(all_data[[data_key]]) > 0) {
      data <- all_data[[data_key]]

      # Get top novelty taxa (same logic as in plots)
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      # Get top overrepresentation taxa
      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      # Combine top taxa from both categories
      domain_phyla <- unique(c(top_novelty_taxa, top_coverage_taxa))
    }
  } else if (domain == "Bacteria") {
    # For Bacteria: only phyla that are in top 10 novelty or top 10 overrepresentation
    data_key <- paste("16S", domain, "phylum", sep = "_")
    if (data_key %in% names(all_data) && nrow(all_data[[data_key]]) > 0) {
      data <- all_data[[data_key]]

      # Get top novelty taxa
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      # Get top overrepresentation taxa
      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      # Combine top taxa from both categories
      domain_phyla <- unique(c(top_novelty_taxa, top_coverage_taxa))
    }
  } else if (domain == "Archaea") {
    # For Archaea: top phyla from phylum row + additional top phyla from family row

    # First get top phyla from phylum row
    phylum_data_key <- paste("16S", domain, "phylum", sep = "_")
    if (phylum_data_key %in% names(all_data) && nrow(all_data[[phylum_data_key]]) > 0) {
      data <- all_data[[phylum_data_key]]

      # Get top novelty and overrepresentation taxa from phylum level
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      domain_phyla <- c(domain_phyla, unique(c(top_novelty_taxa, top_coverage_taxa)))
    }

    # Then add any additional top phyla from family row that aren't already included
    family_data_key <- paste("16S", domain, "family", sep = "_")
    if (family_data_key %in% names(all_data) && nrow(all_data[[family_data_key]]) > 0) {
      data <- all_data[[family_data_key]]

      # Get top taxa from family level and extract their phyla
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      family_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        family_novelty_taxa <- head(novelty_candidates$Phylum, top_n)
      }

      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      family_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        family_coverage_taxa <- head(coverage_candidates$Phylum, top_n)
      }

      family_phyla <- unique(c(family_novelty_taxa, family_coverage_taxa))

      # Only add phyla that aren't already in the list and are actual phyla (not family names)
      color_config <- load_shared_color_config()
      additional_phyla <- family_phyla[!family_phyla %in% domain_phyla &
                                     family_phyla %in% names(color_config$archaea_colors)]
      domain_phyla <- c(domain_phyla, additional_phyla)
    }
  }

  # Remove unknowns and empty values
  domain_phyla <- unique(domain_phyla[!domain_phyla %in% c("Unknown", "", "Other", NA)])

  cat(paste("🔍", domain, "legend will show", length(domain_phyla), "top phyla:", paste(domain_phyla, collapse = ", "), "\n"))

  return(sort(domain_phyla))
}
# Function to create one long combined phyla legend for all domains (only phyla that appear in data)
create_combined_phyla_legend <- function(all_data) {
  color_config <- load_shared_color_config()

  # Extract phyla that actually appear in the data for each domain
  bacteria_phyla <- extract_domain_phyla(all_data, "Bacteria")
  archaea_phyla <- extract_domain_phyla(all_data, "Archaea")
  eukaryota_divisions <- extract_domain_phyla(all_data, "Eukaryota")

  # Combine all taxa with domain labels
  all_taxa <- c()
  all_colors <- c()
  all_domains <- c()

  # Add bacteria (only those that appear in data) - use dynamic assignment
  if (length(bacteria_phyla) > 0) {
    bacteria_colors <- sapply(bacteria_phyla, function(phylum) {
      assign_taxon_color(phylum, "Bacteria", color_config)
    })
    all_taxa <- c(all_taxa, bacteria_phyla)
    all_colors <- c(all_colors, unname(bacteria_colors))
    all_domains <- c(all_domains, rep("Bacteria", length(bacteria_phyla)))
  }

  # Add archaea (only those that appear in data) - use dynamic assignment
  if (length(archaea_phyla) > 0) {
    archaea_colors <- sapply(archaea_phyla, function(phylum) {
      assign_taxon_color(phylum, "Archaea", color_config)
    })
    all_taxa <- c(all_taxa, archaea_phyla)
    all_colors <- c(all_colors, unname(archaea_colors))
    all_domains <- c(all_domains, rep("Archaea", length(archaea_phyla)))
  }

  # Add eukaryota (only those that appear in data) - use dynamic assignment
  if (length(eukaryota_divisions) > 0) {
    euk_colors <- sapply(eukaryota_divisions, function(div) {
      assign_taxon_color(div, "Eukaryota", color_config)
    })
    all_taxa <- c(all_taxa, eukaryota_divisions)
    all_colors <- c(all_colors, unname(euk_colors))
    all_domains <- c(all_domains, rep("Eukaryota", length(eukaryota_divisions)))
  }

  # Create the combined legend
  combined_legend <- create_long_horizontal_legend(all_taxa, all_colors, all_domains)

  return(combined_legend)
}

# Function to create one long horizontal legend with domain separators
create_long_horizontal_legend <- function(all_taxa, all_colors, all_domains) {
  if (length(all_taxa) == 0) {
    return(ggplot() + theme_void())
  }

  # Create legend data
  legend_data <- data.frame(
    taxon = all_taxa,
    color = all_colors,
    domain = all_domains,
    x = seq_along(all_taxa),
    y = rep(1, length(all_taxa)),
    stringsAsFactors = FALSE
  )

  # Create domain boundary positions for separators
  domain_changes <- which(c(TRUE, all_domains[-1] != all_domains[-length(all_domains)]))
  domain_labels <- data.frame(
    domain = unique(all_domains),
    x_pos = domain_changes + (c(domain_changes[-1], length(all_taxa) + 1) - domain_changes) / 2,
    y_pos = rep(1.8, length(unique(all_domains))),
    stringsAsFactors = FALSE
  )

  cat(paste("🎨 Creating combined legend with", nrow(legend_data), "total taxa\n"))
  cat(paste("   Bacteria:", sum(all_domains == "Bacteria"), "taxa\n"))
  cat(paste("   Archaea:", sum(all_domains == "Archaea"), "taxa\n"))
  cat(paste("   Eukaryota:", sum(all_domains == "Eukaryota"), "taxa\n"))

  # Create the plot with clean format matching the eukaryota legend
  p <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_tile(aes(fill = color), color = "black", width = 0.8, height = 0.6, size = 0.5) +
    geom_text(aes(label = taxon), y = 0.3, angle = 45, hjust = 1, vjust = 1, size = 3) +
    scale_fill_identity() +
    theme_void() +
    ggtitle("Combined Phyla Legend: Bacteria | Archaea | Eukaryota") +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    xlim(0.5, length(all_taxa) + 0.5) + ylim(0, 1.5)

  # Add subtle domain separators (minimal visual impact)
  if (length(domain_changes) > 1) {
    separator_x <- domain_changes[-1] - 0.5
    for (sep_x in separator_x) {
      p <- p + geom_vline(xintercept = sep_x, color = "gray80", size = 0.5, linetype = "dotted")
    }
  }

  return(p)
}

# Helper function to create domain-specific horizontal legend with squares (kept for compatibility)
create_domain_specific_legend <- function(taxa_list, color_config, domain, title) {
  if (length(taxa_list) == 0) {
    return(ggplot() + theme_void() + ggtitle(title) +
           theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)))
  }

  # Get colors directly from the appropriate domain config
  if (domain == "Bacteria") {
    domain_colors <- color_config$bacteria_colors
  } else if (domain == "Archaea") {
    domain_colors <- color_config$archaea_colors
  } else if (domain == "Eukaryota") {
    domain_colors <- color_config$eukaryota_colors
  } else {
    domain_colors <- list()
  }

  # Create color vector - use assigned colors or fallback
  colors <- sapply(taxa_list, function(taxon) {
    if (taxon %in% names(domain_colors)) {
      return(domain_colors[[taxon]])
    } else {
      # Use shared color pool for unmapped taxa
      shared_pool <- unique(c(
        unlist(color_config$fallback_colors$bacteria),
        unlist(color_config$fallback_colors$archaea),
        unlist(color_config$fallback_colors$eukaryota)
      ))
      pool_index <- ((match(taxon, taxa_list) - 1) %% length(shared_pool)) + 1
      return(shared_pool[pool_index])
    }
  })

  # Create horizontal layout data
  legend_data <- data.frame(
    taxon = taxa_list,
    color = unname(colors),
    x = seq_along(taxa_list),
    y = rep(1, length(taxa_list)),
    stringsAsFactors = FALSE
  )

  cat(paste("🎨 Creating", domain, "legend with", nrow(legend_data), "taxa\n"))

  # Create horizontal legend with squares
  ggplot(legend_data, aes(x = x, y = y)) +
    geom_tile(aes(fill = color), color = "black", width = 0.8, height = 0.6, size = 0.5) +
    geom_text(aes(label = taxon), y = 0.3, angle = 45, hjust = 1, vjust = 1, size = 3) +
    scale_fill_identity() +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    xlim(0.5, length(taxa_list) + 0.5) + ylim(0, 1.5)
}

# Simple path configuration
# Output everything to the current directory (mega_scrip)
config <- list(
  data_dir_16s = "../../../Eukcensus_merge/16s_merged/csv_results",
  data_dir_18s = "../../../Eukcensus_merge/18s_merged/csv_results",
  output_dir = ".",  # Output to current directory (mega_scrip)
  source_data_dir = "source_data_v2",  # V2 version - separate directory
  ncbi_data_dir = "../../../ncbi_parse/csv_ncbi"
)

# Create output directories in the mega_scrip folder
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

# Function to save source data with organized sorting
save_source_data <- function(data, level, domain) {
  source_filename <- paste0(domain, "_", level, "_source_data.csv")
  source_filepath <- file.path(config$source_data_dir, source_filename)

  data_export <- data
  # Remove redundant domain and level columns - info is clear from filename

  # Organize data: novelty taxa first (sorted by novelty descending),
  # then overrepresentation taxa (sorted by coverage descending),
  # then remaining taxa (sorted by taxon name)
  novelty_data <- data_export[data_export$Is_Top_Novelty == TRUE, ]
  overrep_data <- data_export[data_export$Is_Top_Coverage == TRUE, ]
  other_data <- data_export[data_export$Is_Top_Novelty == FALSE & data_export$Is_Top_Coverage == FALSE, ]

  # Sort each group
  if (nrow(novelty_data) > 0) {
    novelty_data <- novelty_data[order(-novelty_data$novelty_factor), ]
  }
  if (nrow(overrep_data) > 0) {
    overrep_data <- overrep_data[order(-overrep_data$overrepresentation_factor), ]
  }
  if (nrow(other_data) > 0) {
    other_data <- other_data[order(other_data$Taxon), ]
  }

  # Combine in the desired order
  data_export <- rbind(novelty_data, overrep_data, other_data)

  write.csv(data_export, source_filepath, row.names = FALSE)
  cat(paste("Source data saved (organized by novelty/overrepresentation):", source_filepath, "\n"))
}

# Color palettes using shared taxonomic color mapping
get_color_palettes <- function() {
  # Load shared color configuration
  color_config <- load_shared_color_config()

  # Extract colors from shared config
  bacteria_colors_full <- unname(unlist(color_config$bacteria_colors))
  archaea_colors_full <- unname(unlist(color_config$archaea_colors))
  eukaryota_colors_full <- unname(unlist(color_config$eukaryota_colors))

  # Extract eukaryota division names (keys from the config)
  eukaryota_divisions <- names(color_config$eukaryota_colors)

  cat("🎨 LOADED SHARED COLOR CONFIG:\n")
  cat(paste("   Bacteria colors:", length(bacteria_colors_full), "\n"))
  cat(paste("   Archaea colors:", length(archaea_colors_full), "\n"))
  cat(paste("   Eukaryota colors:", length(eukaryota_colors_full), "\n"))
  cat(paste("   Eukaryota divisions:", paste(eukaryota_divisions, collapse = ", "), "\n"))

  # Apply exclusions
  if (USE_POSITION_EXCLUSION) {
    # Exclude by position/index
    bacteria_colors <- bacteria_colors_full[-EXCLUDED_POSITIONS$bacteria]
    archaea_colors <- archaea_colors_full[-EXCLUDED_POSITIONS$archaea]
    eukaryota_colors <- eukaryota_colors_full[-EXCLUDED_POSITIONS$eukaryota]
    eukaryota_divisions <- eukaryota_divisions[-EXCLUDED_POSITIONS$eukaryota]
  } else {
    # Exclude by hex code
    bacteria_colors <- bacteria_colors_full[!bacteria_colors_full %in% EXCLUDED_COLORS$bacteria]
    archaea_colors <- archaea_colors_full[!archaea_colors_full %in% EXCLUDED_COLORS$archaea]
    eukaryota_excluded_indices <- which(eukaryota_colors_full %in% EXCLUDED_COLORS$eukaryota)
    eukaryota_colors <- eukaryota_colors_full[!eukaryota_colors_full %in% EXCLUDED_COLORS$eukaryota]
    eukaryota_divisions <- eukaryota_divisions[-eukaryota_excluded_indices]
  }

  # Print exclusion summary
  if (length(EXCLUDED_COLORS$bacteria) > 0 || length(EXCLUDED_COLORS$archaea) > 0 ||
      length(EXCLUDED_COLORS$eukaryota) > 0 || length(EXCLUDED_POSITIONS$bacteria) > 0 ||
      length(EXCLUDED_POSITIONS$archaea) > 0 || length(EXCLUDED_POSITIONS$eukaryota) > 0) {

    cat("🎨 COLOR EXCLUSIONS APPLIED:\n")
    if (USE_POSITION_EXCLUSION) {
      if (length(EXCLUDED_POSITIONS$bacteria) > 0) {
        cat(paste("   Bacteria: Excluded positions", paste(EXCLUDED_POSITIONS$bacteria, collapse = ", "), "\n"))
      }
      if (length(EXCLUDED_POSITIONS$archaea) > 0) {
        cat(paste("   Archaea: Excluded positions", paste(EXCLUDED_POSITIONS$archaea, collapse = ", "), "\n"))
      }
      if (length(EXCLUDED_POSITIONS$eukaryota) > 0) {
        cat(paste("   Eukaryota: Excluded positions", paste(EXCLUDED_POSITIONS$eukaryota, collapse = ", "), "\n"))
      }
    } else {
      if (length(EXCLUDED_COLORS$bacteria) > 0) {
        cat(paste("   Bacteria: Excluded", paste(EXCLUDED_COLORS$bacteria, collapse = ", "), "\n"))
      }
      if (length(EXCLUDED_COLORS$archaea) > 0) {
        cat(paste("   Archaea: Excluded", paste(EXCLUDED_COLORS$archaea, collapse = ", "), "\n"))
      }
      if (length(EXCLUDED_COLORS$eukaryota) > 0) {
        cat(paste("   Eukaryota: Excluded", paste(EXCLUDED_COLORS$eukaryota, collapse = ", "), "\n"))
      }
    }
    cat(paste("   Final palette sizes: Bacteria =", length(bacteria_colors),
              ", Archaea =", length(archaea_colors),
              ", Eukaryota =", length(eukaryota_colors), "\n"))
  }

  return(list(
    bacteria = bacteria_colors,
    archaea = archaea_colors,
    eukaryota = setNames(eukaryota_colors, eukaryota_divisions)
  ))
}

# Helper function to filter divisions
filter_divisions <- function(data) {
  data %>% filter(!Division %in% c("Unknown", "", "Other", NA))
}

# Load 16S data function
load_16s_data <- function(level, domain) {
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)

  if (!file.exists(filepath)) stop(paste("File not found:", filepath))

  data <- read.csv(filepath, stringsAsFactors = FALSE) %>%
    filter(domain == !!domain, census_otu_count > 0, ncbi_species_count > 0)

  # Return empty data frame with required columns if no data after filtering
  if (nrow(data) == 0) {
    empty_df <- data.frame(
      Taxon = character(0),
      Circle_Size = numeric(0),
      Is_Top_Novelty = logical(0),
      Is_Top_Coverage = logical(0),
      stringsAsFactors = FALSE
    )
    return(empty_df)
  }

  # Standardize column names
  colnames(data)[1] <- "Taxon"
  data$Circle_Size <- calculate_circle_size(data$isolate_count, data$ncbi_genome_count)

  # Add phylum information
  data <- add_phylum_info(data, level, domain)

  # Identify top taxa (threshold > 1.0, up to top 10 that meet criteria)
  threshold <- 1.0

  # For novelty: only taxa above threshold, limited to top 10
  novelty_candidates <- data[data$novelty_factor > threshold, ]
  if (nrow(novelty_candidates) > 0) {
    novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
    top_novelty_taxa <- head(novelty_candidates$Taxon, PLOT_CONFIG$top_n)
    data$Is_Top_Novelty <- data$Taxon %in% top_novelty_taxa
  } else {
    data$Is_Top_Novelty <- FALSE
  }

  # For coverage: only taxa above threshold, limited to top 10
  coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
  if (nrow(coverage_candidates) > 0) {
    coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
    top_coverage_taxa <- head(coverage_candidates$Taxon, PLOT_CONFIG$top_n)
    data$Is_Top_Coverage <- data$Taxon %in% top_coverage_taxa
  } else {
    data$Is_Top_Coverage <- FALSE
  }

  return(data)
}

# Add phylum information
add_phylum_info <- function(data, level, domain) {
  if (level == "phylum") {
    data$Phylum <- data$Taxon
  } else {
    # For family level, try to get phylum mapping but be very conservative about filtering
    data$Phylum <- get_phylum_for_taxa(data$Taxon, level)

    # Always use family names as phylum for coloring to preserve all data
    # This ensures important families like Prochlorococcaceae are never lost
    cat(paste("ℹ️  Using", level, "names directly for coloring to preserve all taxa.\n"))
    data$Phylum <- data$Taxon  # Use family names directly - no filtering!

    # Optional: If you had valid phylum mappings, you could use them here
    # But for now, we prioritize data preservation over phylum accuracy
  }
  data$Phylum <- standardize_phylum_names(data$Phylum)
  return(data)
}

# Get phylum mapping for family level
get_phylum_for_taxa <- function(taxa, level) {
  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) return(rep("Unknown", length(taxa)))

  ncbi_data <- read.csv(file_path, stringsAsFactors = FALSE)
  taxon_col <- intersect(c("family_name", "genus_name", "family", "genus", "taxon", "Taxon"),
                        colnames(ncbi_data))[1]

  if (is.na(taxon_col) || !all(c("lineage", "lineage_ranks") %in% colnames(ncbi_data))) {
    return(rep("Unknown", length(taxa)))
  }

  # Extract phylum from lineage
  phylum_map <- setNames(rep("Unknown", nrow(ncbi_data)), ncbi_data[[taxon_col]])

  for (i in 1:nrow(ncbi_data)) {
    if (!is.na(ncbi_data$lineage[i]) && !is.na(ncbi_data$lineage_ranks[i])) {
      lineage <- strsplit(ncbi_data$lineage[i], ";")[[1]]
      ranks <- strsplit(ncbi_data$lineage_ranks[i], ";")[[1]]
      phylum_idx <- which(ranks == "phylum")
      if (length(phylum_idx) > 0 && phylum_idx[1] <= length(lineage)) {
        phylum_map[ncbi_data[[taxon_col]][i]] <- lineage[phylum_idx[1]]
      }
    }
  }

  result <- phylum_map[taxa]
  result[is.na(result)] <- "Unknown"
  return(unname(result))
}

# Standardize phylum names
standardize_phylum_names <- function(phylum_names) {
  phylum_mapping <- c("Methanobacteriota" = "Euryarchaeota")

  standardized <- phylum_names
  for (old_name in names(phylum_mapping)) {
    standardized[standardized == old_name] <- phylum_mapping[old_name]
  }
  return(standardized)
}

# Load 18S data function
load_18s_data <- function(level) {
  filename <- paste0("18s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_18s, filename)

  if (!file.exists(filepath)) stop(paste("File not found:", filepath))

  data <- read.csv(filepath, stringsAsFactors = FALSE) %>%
    filter(census_otu_count > 0, ncbi_species_count > 0)

  # Return empty data frame with required columns if no data after filtering
  if (nrow(data) == 0) {
    empty_df <- data.frame(
      Taxon = character(0),
      Circle_Size = numeric(0),
      Is_Top_Novelty = logical(0),
      Is_Top_Coverage = logical(0),
      stringsAsFactors = FALSE
    )
    return(empty_df)
  }

  # Standardize column names
  colnames(data)[1] <- "Taxon"
  data$Circle_Size <- calculate_circle_size(data$isolate_count, data$ncbi_genome_count)

  # Add division information
  if (level == "phylum") {
    data$Division <- data$Taxon
  } else {
    data$Division <- "Other"
    manual_overrides <- get_18s_manual_overrides()
    for (taxon in names(manual_overrides)) {
      if (taxon %in% data$Taxon) {
        data$Division[data$Taxon == taxon] <- manual_overrides[[taxon]]
      }
    }
  }

  data <- filter_divisions(data)
  # Return empty data frame with required columns if no data after division filtering
  if (nrow(data) == 0) {
    empty_df <- data.frame(
      Taxon = character(0),
      Circle_Size = numeric(0),
      Division = character(0),
      Is_Top_Novelty = logical(0),
      Is_Top_Coverage = logical(0),
      stringsAsFactors = FALSE
    )
    return(empty_df)
  }

  # Identify top taxa (threshold > 1.0, up to top 10 that meet criteria)
  threshold <- 1.0

  # For novelty: only taxa above threshold, limited to top 10
  novelty_candidates <- data[data$novelty_factor > threshold, ]
  if (nrow(novelty_candidates) > 0) {
    novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
    top_novelty_taxa <- head(novelty_candidates$Taxon, PLOT_CONFIG$top_n)
    data$Is_Top_Novelty <- data$Taxon %in% top_novelty_taxa
  } else {
    data$Is_Top_Novelty <- FALSE
  }

  # For coverage: only taxa above threshold, limited to top 10
  coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
  if (nrow(coverage_candidates) > 0) {
    coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
    top_coverage_taxa <- head(coverage_candidates$Taxon, PLOT_CONFIG$top_n)
    data$Is_Top_Coverage <- data$Taxon %in% top_coverage_taxa
  } else {
    data$Is_Top_Coverage <- FALSE
  }

  return(data)
}

# Manual overrides for 18S specific taxa (from 18S script)
get_18s_manual_overrides <- function() {
  list(
    # OPISTHOKONTA - Animals and Fungi
    "Insecta" = "Opisthokonta", "Mammalia" = "Opisthokonta", "Teleostei" = "Opisthokonta",
    "Amphibia" = "Opisthokonta", "Lepidosauria" = "Opisthokonta", "Arachnida" = "Opisthokonta",
    "Anthozoa" = "Opisthokonta", "Branchiopoda" = "Opisthokonta", "Caenogastropoda" = "Opisthokonta",
    "Digenea" = "Opisthokonta", "Rhabdocoela" = "Opisthokonta", "Evaginogenida" = "Opisthokonta",
    "Sordariomycetes" = "Opisthokonta", "Eurotiomycetes" = "Opisthokonta", "Saccharomycetales" = "Opisthokonta",
    "Dothideomycetes" = "Opisthokonta", "Agaricomycetes" = "Opisthokonta", "Tremellomycetes" = "Opisthokonta",
    "Lecanoromycetes" = "Opisthokonta", "Leotiomycetes" = "Opisthokonta", "Kickxellales" = "Opisthokonta",
    "Ustilaginomycetes" = "Opisthokonta", "Schizosaccharomycetes" = "Opisthokonta", "Dacrymycetes" = "Opisthokonta",
    "Lipomycetaceae" = "Opisthokonta", "Harpellales" = "Opisthokonta", "Spizellomycetaceae" = "Opisthokonta",
    "Rhizophydiaceae" = "Opisthokonta", "Haliphthorales" = "Opisthokonta",

    # ALVEOLATA - Ciliates, Dinoflagellates, Apicomplexans
    "Oxytrichidae" = "Alveolata", "Parameciidae" = "Alveolata", "Gregarinidae" = "Alveolata",
    "Actinocephalidae" = "Alveolata", "Theileriidae" = "Alveolata", "Sarcocystidae" = "Alveolata",
    "Eimeriidae" = "Alveolata", "Cryptosporidiidae" = "Alveolata", "Babesiidae" = "Alveolata",
    "Stylocephalidae" = "Alveolata", "Pseudocolliniidae" = "Alveolata", "Dysteriidae" = "Alveolata",
    "Adeleidae" = "Alveolata", "Xcellidae" = "Alveolata", "Corallicolidae" = "Alveolata",

    # STRAMENOPILES
    "Labyrinthulaceae" = "Stramenopiles", "Peronosporales" = "Stramenopiles", "Bicoecaceae" = "Stramenopiles",
    "Thaumatomonadidae" = "Stramenopiles", "Saprolegniales" = "Stramenopiles", "Thalassiosiraceae" = "Stramenopiles",
    "Thalassiosira" = "Stramenopiles", "Pythium" = "Stramenopiles",

    # DISCOBA
    "Trypanosomatidae" = "Discoba", "Neobodonidae" = "Discoba", "Bodonidae" = "Discoba",
    "Neovahlkampfiidae" = "Discoba", "Cercomonadidae" = "Discoba", "Psalteriomonadidae" = "Discoba",
    "Distigmidae" = "Discoba", "Spironemidae" = "Discoba", "Trachelomonas" = "Discoba",
    "Paracercomonadidae" = "Discoba", "Rhynchomonadidae" = "Discoba",

    # METAMONADA
    "Hexamitidae" = "Metamonada",

    # EVOSEA
    "Mastigamoebidae" = "Evosea", "Vannellidae" = "Evosea", "Paramoebidae" = "Evosea",
    "Vermamoebidae" = "Evosea", "Vampyrellidae" = "Evosea", "Balamuthiidae" = "Evosea",
    "Leptomyxidae" = "Evosea", "Hartmannulidae" = "Evosea", "Echinamoebidae" = "Evosea",
    "Plasmodiidae" = "Evosea", "Schizoplasmodiidae" = "Evosea", "Mastigamoeba" = "Evosea",

    # CHLOROPHYTA & STREPTOPHYTA
    "Coccomyxaceae" = "Chlorophyta", "Ceratiaceae" = "Chlorophyta", "Zea" = "Streptophyta"
  )
}

# Create individual scatter plot function
create_individual_scatter <- function(data, level, domain, master_colors) {
  save_source_data(data, level, domain)

  # Get top data for highlighting
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]

  # Base plot
  p <- ggplot(data, aes(x = census_otu_count, y = ncbi_species_count)) +
    scale_x_log10(labels = comma_format(), limits = c(1, 10000)) +
    scale_y_log10(labels = comma_format(), limits = c(1, 10000))

  # Background points
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "lightgray", fill = "lightgray",
                       shape = PLOT_CONFIG$circle_shape, alpha = PLOT_CONFIG$bg_alpha,
                       stroke = PLOT_CONFIG$circle_stroke)
  }

  # Top data points with proper color mapping using shared config
  if (nrow(top_data) > 0) {
    # Load shared color configuration
    color_config <- load_shared_color_config()

    # Determine color column and assign colors dynamically
    if (domain == "Eukaryota") {
      color_col <- "Division"
      plot_groups <- sort(unique(top_data$Division[!top_data$Division %in% c("Unknown", "", "Other", NA)]))
    } else {
      color_col <- "Phylum"
      plot_groups <- sort(unique(top_data$Phylum[!top_data$Phylum %in% c("Unknown", "", "Other", NA)]))
    }

    # Use shared mapping logic (same as alluvial) for exact cross-figure consistency.
    group_colors <- get_domain_colors(plot_groups, domain, color_config)

    # Print color mapping for verification
    cat(paste("🎨", domain, level, "color mapping:\n"))
    print_color_summary(names(group_colors), group_colors, domain)

    # Add colored points with phylum/division-based colors
    p <- p + geom_point(data = top_data,
                       aes_string(size = "Circle_Size",
                                 fill = paste0("factor(", color_col, ", levels = plot_groups)")),
                       color = "black",
                       shape = PLOT_CONFIG$circle_shape, alpha = PLOT_CONFIG$circle_alpha,
                       stroke = PLOT_CONFIG$circle_stroke) +
      scale_fill_manual(values = group_colors, guide = "none")
  }

  # Styling and theme
  p <- p +
    scale_size_continuous(range = PLOT_CONFIG$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, color = "grey50"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )

  # Organized annotations for top taxa with improved repelling
  if (nrow(top_data) > 0) {
    # Separate novelty and overrepresentation taxa
    novelty_data <- top_data[top_data$Is_Top_Novelty == TRUE, ]
    overrep_data <- top_data[top_data$Is_Top_Coverage == TRUE, ]

    # Sort novelty data by novelty factor (descending)
    if (nrow(novelty_data) > 0) {
      novelty_data <- novelty_data[order(-novelty_data$novelty_factor), ]

      # Different label styles for different levels
      if (level == "phylum") {
        # Phylum level: only show factor values
        novelty_data$label <- paste0(sprintf("%.1f", novelty_data$novelty_factor), "×")
      } else {
        # Family level: show taxon name and factor
        novelty_data$label <- paste0(novelty_data$Taxon, " (", sprintf("%.1f", novelty_data$novelty_factor), "×)")
      }

      # Add novelty annotations (repel DOWNWARD away from diagonal midline)
      p <- p + ggrepel::geom_text_repel(
        data = novelty_data,
        aes(label = label),
        size = PLOT_CONFIG$text_size * 0.9,  # Slightly smaller for better fit
        fontface = "bold",
        color = "black",
        max.overlaps = Inf,
        force = 4,                    # Strong repelling force
        force_pull = 0.1,            # Gentle pull toward points
        box.padding = 1.0,           # More padding around text boxes
        point.padding = 0.5,         # Padding around points
        segment.color = "gray40",    # Visible connector lines
        segment.size = 0.3,
        segment.alpha = 0.8,
        min.segment.length = 0,      # Always show connector lines
        direction = "both",          # Allow movement in both x and y
        nudge_y = -0.5,             # Strong DOWNWARD push away from midline
        nudge_x = 0,                # No horizontal bias
        xlim = c(NA, NA),           # No x limits
        ylim = c(NA, NA),           # No y limits
        seed = 42
      )
    }

    # Sort overrepresentation data by coverage factor (descending)
    if (nrow(overrep_data) > 0) {
      overrep_data <- overrep_data[order(-overrep_data$overrepresentation_factor), ]

      # Different label styles for different levels
      if (level == "phylum") {
        # Phylum level: only show factor values
        overrep_data$label <- paste0(sprintf("%.1f", overrep_data$overrepresentation_factor), "×")
      } else {
        # Family level: show taxon name and factor
        overrep_data$label <- paste0(overrep_data$Taxon, " (", sprintf("%.1f", overrep_data$overrepresentation_factor), "×)")
      }

      # Add overrepresentation annotations (repel UPWARD away from diagonal midline)
      p <- p + ggrepel::geom_text_repel(
        data = overrep_data,
        aes(label = label),
        size = PLOT_CONFIG$text_size * 0.9,  # Slightly smaller for better fit
        fontface = "bold",
        color = "black",
        max.overlaps = Inf,
        force = 4,                    # Strong repelling force
        force_pull = 0.1,            # Gentle pull toward points
        box.padding = 1.0,           # More padding around text boxes
        point.padding = 0.5,         # Padding around points
        segment.color = "gray40",    # Visible connector lines
        segment.size = 0.3,
        segment.alpha = 0.8,
        min.segment.length = 0,      # Always show connector lines
        direction = "both",          # Allow movement in both x and y
        nudge_y = 0.5,              # Strong UPWARD push away from midline
        nudge_x = 0,                # No horizontal bias
        xlim = c(NA, NA),           # No x limits
        ylim = c(NA, NA),           # No y limits
        seed = 123
      )
    }
  }


  return(p)
}




# Main function to create comprehensive mega visual
main <- function() {
  cat("Comprehensive Mega Stacked Visual Creation (V2 - SIMPLIFIED)\n")
  cat("=============================================================\n")

  # Initialize color mapping
  init_color_registry()
  cat("✅ Initialized color assignment (V2: FIXED - No repeats, cross-domain pool)\n")

  # Define the grid structure
  levels <- c("phylum", "family")
  domains_16s <- c("Bacteria", "Archaea")
  domain_18s <- "Eukaryota"

  # Get master color palettes
  master_colors <- get_color_palettes()

  # Collect all data
  all_data <- list()

  cat("Loading 16S data...\n")
  for (level in levels) {
    for (domain in domains_16s) {
      cat(paste("Loading", domain, level, "data...\n"))
      data <- load_16s_data(level, domain)
      all_data[[paste("16S", domain, level, sep = "_")]] <- data
    }
  }

  cat("Loading 18S data...\n")
  for (level in levels) {
    cat(paste("Loading eukaryote", level, "data...\n"))
    data <- load_18s_data(level)
    all_data[[paste("18S", domain_18s, level, sep = "_")]] <- data
  }

  # Create individual plots
  plots <- list()
  cat("Creating individual scatter plots...\n")

  for (level in levels) {
    # 16S plots
    for (domain in domains_16s) {
      data_key <- paste("16S", domain, level, sep = "_")
      data <- all_data[[data_key]]
      if (!is.null(data) && nrow(data) > 0) {
        cat(paste("Creating plot for", domain, level, "with", nrow(data), "taxa\n"))
        # Show first few taxa for verification
        cat(paste("  Sample taxa:", paste(head(data$Taxon, 3), collapse = ", "), "\n"))
        plots[[data_key]] <- create_individual_scatter(data, level, domain, master_colors)
      }
    }

    # 18S plot
    data_key <- paste("18S", domain_18s, level, sep = "_")
    data <- all_data[[data_key]]
    if (!is.null(data) && nrow(data) > 0) {
      cat(paste("Creating plot for", domain_18s, level, "\n"))
      plots[[data_key]] <- create_individual_scatter(data, level, domain_18s, master_colors)
    }
  }

  # Arrange plots in 2x3 grid (2 rows for levels, 3 columns for domains)
  cat("Arranging plots in comprehensive grid...\n")

  # Create row labels (left side)
  row_labels <- list()
  level_names <- c("Phyla", "Family")
  for (i in 1:length(levels)) {
    row_labels[[i]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = level_names[i],
               size = 12, fontface = "bold", color = "grey30", angle = 90) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  # Create column headers
  col_headers <- list()
  domain_names <- c("Bacteria", "Archaea", "Eukaryota")
  for (i in 1:length(domain_names)) {
    col_headers[[i]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = domain_names[i],
               size = 12, fontface = "bold", color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  # Create top row with column headers and spacing
  top_row <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    col_headers[[1]],         # Bacteria header
    ggplot() + theme_void(),  # Spacer between Bacteria and Archaea
    col_headers[[2]],         # Archaea header
    ggplot() + theme_void(),  # Spacer between Archaea and Eukaryota
    col_headers[[3]],         # Eukaryota header
    ncol = 6, rel_widths = c(0.1, 1, 0.05, 1, 0.05, 1)
  )

  # Create data rows with individual black frames around each plot
  data_rows <- list()
  for (i in 1:length(levels)) {
    level <- levels[i]

    # Get plots for this level
    bacteria_plot <- plots[[paste("16S", "Bacteria", level, sep = "_")]]
    archaea_plot <- plots[[paste("16S", "Archaea", level, sep = "_")]]
    eukaryota_plot <- plots[[paste("18S", "Eukaryota", level, sep = "_")]]

    # Add black frames around individual plots
    framed_bacteria_plot <- ggdraw(bacteria_plot) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))
    framed_archaea_plot <- ggdraw(archaea_plot) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))
    framed_eukaryota_plot <- ggdraw(eukaryota_plot) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))

    # Combine row label with framed plots and spacing
    data_rows[[i]] <- plot_grid(
      row_labels[[i]],
      framed_bacteria_plot,
      ggplot() + theme_void(),  # Spacer between Bacteria and Archaea
      framed_archaea_plot,
      ggplot() + theme_void(),  # Spacer between Archaea and Eukaryota
      framed_eukaryota_plot,
      ncol = 6, rel_widths = c(0.1, 1, 0.05, 1, 0.05, 1)
    )
  }

  # Create isolate percentage legend (keep on right side)
  isolate_legend <- create_isolate_legend()

  # Create one long combined phyla legend
  combined_phyla_legend <- create_combined_phyla_legend(all_data)

  # Combine main plots with isolate legend on the right
  main_plot_with_isolate <- plot_grid(
    plot_grid(
      top_row,                  # Column header row
      data_rows[[1]],           # Phylum row
      ggplot() + theme_void(),  # Small spacer
      data_rows[[2]],           # Family row
      ncol = 1, rel_heights = c(0.1, 1, 0.05, 1)
    ),
    isolate_legend,
    ncol = 2, rel_widths = c(0.85, 0.15)
  )

  # Combine main plot with long horizontal legend beneath
  complete_plot <- plot_grid(
    main_plot_with_isolate,
    ggplot() + theme_void(),  # Small spacer
    combined_phyla_legend,
    ncol = 1, rel_heights = c(0.75, 0.05, 0.2)
  )

  # Save the comprehensive mega visual (V2 - SIMPLIFIED COLOR ASSIGNMENT)
  output_file_png <- file.path(config$output_dir, "comprehensive_mega_stacked_visual_v2.png")
  output_file_pdf <- file.path(config$output_dir, "comprehensive_mega_stacked_visual_v2.pdf")

  cat(paste("Saving comprehensive mega visual (V2) to:", output_file_png, "\n"))
  cat(paste("Saving PDF version (V2) to:", output_file_pdf, "\n"))

  # Save PNG version
  ggsave(output_file_png, complete_plot,
         width = PLOT_CONFIG$plot_width, height = PLOT_CONFIG$plot_height,
         dpi = PLOT_CONFIG$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  # Save PDF version for Illustrator
  ggsave(output_file_pdf, complete_plot,
         width = PLOT_CONFIG$plot_width, height = PLOT_CONFIG$plot_height,
         bg = "white", limitsize = FALSE,
         units = "in")

  # Save the combined phyla legend as a separate file
  save_combined_legend(combined_phyla_legend)

  # Create master source data index
  create_source_data_index()

  cat("✅ Comprehensive mega visual creation complete (V2)!\n")
  cat(paste("   Main plot (V2):", output_file_png, "\n"))
  cat(paste("   Dimensions:", PLOT_CONFIG$plot_width, "x", PLOT_CONFIG$plot_height, "inches\n"))
  cat(paste("   Layout: 2 rows (phylum, family) × 3 columns (Bacteria, Archaea, Eukaryota)\n"))
  cat(paste("📊 Source data (V2) saved to:", config$source_data_dir, "\n"))
}

# Function to save the combined horizontal phyla legend
save_combined_legend <- function(combined_legend) {
  legend_dir <- file.path(config$output_dir, "phyla_legends")
  if (!dir.exists(legend_dir)) dir.create(legend_dir, recursive = TRUE)

  # Save the long combined legend (perfect for external manipulation) - V2
  combined_file <- file.path(legend_dir, "combined_phyla_legend_v2.png")
  ggsave(combined_file, combined_legend, width = 20, height = 4, dpi = 300, bg = "white")
  cat(paste("📊 Combined phyla legend (V2) saved:", combined_file, "\n"))
  cat(paste("    Format: Bacteria | Archaea | Eukaryota in one long horizontal strip\n"))
  cat(paste("    Perfect for external software manipulation!\n"))
}

# Function to create master index of all source data files
create_source_data_index <- function() {
  # Get all CSV files in source_data directory (excluding the index file itself)
  csv_files <- list.files(config$source_data_dir, pattern = "\\.csv$", full.names = FALSE)
  csv_files <- csv_files[!grepl("README", csv_files)]  # Exclude the index file itself

  if (length(csv_files) > 0) {
    # Separate source data files from legend files
    source_files <- csv_files[grepl("source_data", csv_files)]
    legend_files <- csv_files[grepl("legend", csv_files)]

    # Create descriptions
    source_descriptions <- sapply(source_files, function(f) {
      paste("Plot data for", gsub("_source_data.csv", "", f))
    })

    legend_descriptions <- sapply(legend_files, function(f) {
      paste("Color legend for", gsub("_legend.csv", "", f))
    })

    # Create index data frame
    index_data <- data.frame(
      File_Type = c(rep("Source_Data", length(source_files)), rep("Legend", length(legend_files))),
      Filename = c(source_files, legend_files),
      Description = c(source_descriptions, legend_descriptions),
      Created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )

    # Save index file
    index_filepath <- file.path(config$source_data_dir, "README_source_data_index.csv")
    write.csv(index_data, index_filepath, row.names = FALSE)
    cat(paste("📋 Source data index created:", index_filepath, "\n"))
    cat(paste("   Total files indexed:", nrow(index_data), "\n"))
  }
}



# Run the main visualization
if (!interactive()) {
  main()
}
