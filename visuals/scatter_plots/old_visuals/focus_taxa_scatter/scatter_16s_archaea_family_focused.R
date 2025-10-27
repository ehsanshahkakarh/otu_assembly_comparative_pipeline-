#!/usr/bin/env Rscript
# 16S Archaea Family Focused Scatter Plot
# Created: 2025-08-18
# Purpose: Clean, focused scatter plot for 16S Archaea at Family level only

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "scatter_16s_archaea_family_focused.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Clean configuration for 16S Archaea Family analysis
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "merged_output", "16s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 38,
  plot_height = 15,
  dpi = 300,
  top_n = 10,  # Top 10 coverage + top 10 novelty
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

# Get phylum mapping for family level
get_phylum_mapping <- function() {
  file_path <- file.path(config$ncbi_data_dir, "ncbi_family_counts.csv")
  if (!file.exists(file_path)) {
    cat(paste("Warning: NCBI family file not found:", file_path, "\n"))
    return(NULL)
  }

  data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat(paste("Loaded NCBI family file with", nrow(data), "rows\n"))
  cat(paste("Columns:", paste(colnames(data), collapse = ", "), "\n"))

  # Check if we have the required columns
  if (!"lineage" %in% colnames(data)) {
    cat("Error: 'lineage' column not found in NCBI family file\n")
    return(NULL)
  }
  if (!"lineage_ranks" %in% colnames(data)) {
    cat("Error: 'lineage_ranks' column not found in NCBI family file\n")
    return(NULL)
  }

  # Check what the family name column is called
  family_col <- NULL
  for (col in c("family_name", "family", "taxon", "Taxon")) {
    if (col %in% colnames(data)) {
      family_col <- col
      break
    }
  }

  if (is.null(family_col)) {
    cat("Error: Could not find family name column. Available columns:", paste(colnames(data), collapse = ", "), "\n")
    return(NULL)
  }

  cat(paste("Using", family_col, "as family name column\n"))

  # Extract phylum from lineage using rank information
  phylum_map <- setNames(character(nrow(data)), data[[family_col]])
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
        phylum_map[data[[family_col]][i]] <- phylum_name
        successful_mappings <- successful_mappings + 1

        # Show first few examples
        if (successful_mappings <= 3) {
          cat(paste("Example mapping:", data[[family_col]][i], "->", phylum_name, "\n"))
        }
      } else {
        phylum_map[data[[family_col]][i]] <- "Other"
      }
    } else {
      phylum_map[data[[family_col]][i]] <- "Other"
    }
  }

  cat(paste("Successfully mapped", successful_mappings, "families to phyla\n"))
  unique_phyla <- unique(phylum_map[phylum_map != "Other"])
  cat(paste("Unique phyla found:", length(unique_phyla), "\n"))
  cat(paste("Phyla:", paste(unique_phyla, collapse = ", "), "\n"))

  return(phylum_map)
}

# Manual overrides for specific taxa
get_manual_overrides <- function() {
  manual_overrides <- list(
    # Rename Methanobacteriota to Euryarchaeota for consistency with Census naming
    "Methanobacteriota" = "Euryarchaeota"
  )

  cat(paste("Manual overrides defined for", length(manual_overrides), "taxa\n"))
  for (taxon in names(manual_overrides)) {
    cat(paste("  Override:", taxon, "->", manual_overrides[[taxon]], "\n"))
  }

  return(manual_overrides)
}

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

# Load and process 16S archaea family data
load_archaea_family_data <- function() {
  filepath <- file.path(config$data_dir, "16s_ncbi_merged_clean_family.csv")
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
  cat(paste("Found", nrow(archaea_data), "archaeal families\n"))
  cat(paste("Total families in dataset:", nrow(data), "\n"))

  # Filter for families with meaningful data (Census OR NCBI data > 0 for visualization)
  # Include both matched and ncbi_only entries for archaea
  data_with_data <- archaea_data %>% filter(Census_OTU_Count > 0 | NCBI_Species_Count > 0)
  cat(paste("Found", nrow(data_with_data), "archaeal families with Census or NCBI data\n"))

  # For ratio calculations, we need both Census AND NCBI data > 0
  data_with_both <- data_with_data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  cat(paste("Found", nrow(data_with_both), "archaeal families with both Census and NCBI data for ratios\n"))

  # Use all data for visualization, but only calculate ratios for those with both
  data <- data_with_data

  # Add phylum mapping using NCBI data first, then apply manual overrides
  phylum_map <- get_phylum_mapping()
  if (!is.null(phylum_map)) {
    data$Phylum <- ifelse(data$Taxon %in% names(phylum_map), phylum_map[data$Taxon], "Other")
  } else {
    data$Phylum <- "Other"
  }

  # Apply manual overrides for specific taxa
  manual_overrides <- get_manual_overrides()
  if (!is.null(manual_overrides)) {
    for (taxon in names(manual_overrides)) {
      if (taxon %in% data$Taxon) {
        data$Phylum[data$Taxon == taxon] <- manual_overrides[[taxon]]
        cat(paste("Applied override:", taxon, "->", manual_overrides[[taxon]], "\n"))
      }
    }
  }

  # Handle NA and "Other" phylum assignments for archaea families
  # For archaea, if phylum is NA or "Other", try to infer from common archaea phyla
  na_or_other_mask <- is.na(data$Phylum) | data$Phylum == "Other" | data$Phylum == ""
  if (sum(na_or_other_mask) > 0) {
    cat(paste("Found", sum(na_or_other_mask), "families with NA/Other phylum assignments\n"))
    # For archaea families without clear phylum assignment, default to Euryarchaeota
    # since most archaea families in environmental samples belong to this phylum
    data$Phylum[na_or_other_mask] <- "Euryarchaeota"
    cat(paste("Assigned", sum(na_or_other_mask), "families to Euryarchaeota as default\n"))
  }

  # Global rename: Convert all "Methanobacteriota" to "Euryarchaeota" for consistency with Census naming
  methanobacteriota_mask <- data$Phylum == "Methanobacteriota"
  if (sum(methanobacteriota_mask) > 0) {
    cat(paste("Found", sum(methanobacteriota_mask), "families labeled as Methanobacteriota\n"))
    data$Phylum[methanobacteriota_mask] <- "Euryarchaeota"
    cat(paste("Renamed", sum(methanobacteriota_mask), "families from Methanobacteriota to Euryarchaeota\n"))
  }

  # Report final phylum distribution
  phylum_counts <- table(data$Phylum)
  cat("Final phylum distribution:\n")
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

  # STRICT FILTERING: Only consider families with BOTH novelty > 1.0 AND coverage > 1.0
  # This ensures we never show families with coverage < 1.0
  meaningful_families <- data[(!is.na(data$Novelty_Ratio) & data$Novelty_Ratio > 1.0) &
                             (!is.na(data$Coverage_Factor) & data$Coverage_Factor > 1.0), ]

  cat(paste("Families with BOTH novelty > 1.0 AND coverage > 1.0:", nrow(meaningful_families), "\n"))

  if (nrow(meaningful_families) == 0) {
    cat("No families meet the strict criteria (both novelty > 1.0 AND coverage > 1.0)\n")
    cat("Falling back to families with novelty > 1.0 OR coverage > 1.0\n")
    meaningful_families <- data[(!is.na(data$Novelty_Ratio) & data$Novelty_Ratio > 1.0) |
                               (!is.na(data$Coverage_Factor) & data$Coverage_Factor > 1.0), ]
    cat(paste("Families with novelty > 1.0 OR coverage > 1.0:", nrow(meaningful_families), "\n"))
  }

  # Filter for meaningful ratios (> 1.0) before ranking, excluding NA values
  high_novelty_data <- meaningful_families[!is.na(meaningful_families$Novelty_Ratio) & meaningful_families$Novelty_Ratio > 1.0, ]
  high_coverage_data <- meaningful_families[!is.na(meaningful_families$Coverage_Factor) & meaningful_families$Coverage_Factor > 1.0, ]

  cat(paste("Meaningful families with novelty ratio > 1.0:", nrow(high_novelty_data), "\n"))
  cat(paste("Meaningful families with coverage ratio > 1.0:", nrow(high_coverage_data), "\n"))

  # Identify top taxa using the selected ranking method, but only from meaningful families
  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE

  # Mark top novelty (only from meaningful families with ratio > 1.0)
  if (nrow(high_novelty_data) > 0) {
    # Get novelty values for the filtered data and rank them
    if (config$use_weighted_ratio) {
      ranking_values <- high_novelty_data$Weighted_Novelty
    } else {
      ranking_values <- high_novelty_data$Novelty_Ratio
    }
    novelty_ranks <- rank(-ranking_values)

    # Get top N taxa from the high novelty group
    top_n_count <- min(config$top_n, length(ranking_values))
    top_novelty_mask <- novelty_ranks <= top_n_count
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_mask]

    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
    cat(paste("Selected top", length(top_novelty_taxa), "novelty families:", paste(top_novelty_taxa, collapse = ", "), "\n"))
  }

  # Mark top coverage (only from meaningful families with ratio > 1.0)
  if (nrow(high_coverage_data) > 0) {
    # Get coverage values for the filtered data and rank them
    if (config$use_weighted_ratio) {
      ranking_values <- high_coverage_data$Weighted_Coverage
    } else {
      ranking_values <- high_coverage_data$Coverage_Factor
    }
    coverage_ranks <- rank(-ranking_values)

    # Get top N taxa from the high coverage group
    top_n_count <- min(config$top_n, length(ranking_values))
    top_coverage_mask <- coverage_ranks <= top_n_count
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_mask]

    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
    cat(paste("Selected top", length(top_coverage_taxa), "coverage families:", paste(top_coverage_taxa, collapse = ", "), "\n"))
  }

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
create_archaea_family_scatter <- function(data) {
  # Base plot with professional styling - wider limits for family level
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 1000),  # Range for family-level data
                  expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 1000),  # Range for family-level data
                  expand = expansion(mult = c(0.05, 0.05)))

  # Remove background points - only show highlighted taxa (already filtered for factors > 1.0)
  # Highlighted points with colors by families
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]

  cat(paste("Showing top archaea families (pre-filtered for factors > 1.0):", nrow(top_data), "\n"))

  if (nrow(top_data) > 0) {
    # Show only phyla of the meaningful annotated taxa in the legend
    top_phyla <- unique(top_data$Phylum[top_data$Phylum != "Other"])
    top_phyla <- sort(top_phyla)
    n_phyla <- length(top_phyla)

    # Generate colors using configurable palette
    colors <- get_color_palette(n_phyla)

    # Add black outline layer first
    p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.5)

    # Add colored fill layer on top - color by phylum, not individual family
    p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Phylum, levels = top_phyla)),
                       alpha = 0.9) +
      scale_color_manual(values = setNames(colors, top_phyla), name = "Archaeal Phyla",
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
          # Legend styling - optimized for taller figure
          legend.text = element_text(size = 48),
          legend.title = element_text(size = 56, face = "bold", color = "grey50"),
          legend.key.size = unit(4.0, "cm"),
          legend.key.width = unit(4.5, "cm"),
          legend.margin = margin(l = 50, t = 20, b = 20),  # Added top/bottom margins
          panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)) +
    labs(title = "16S EukCensus vs NCBI: Archaea Family",
         x = "16S EukCensus OTU Count", y = "NCBI Species Count",
         caption = paste(nrow(data), "archaeal families"))

  # Add circle size legend in top-left corner
  # Use the actual plot limits instead of data range for consistent positioning
  x_range <- c(1, 1000)  # Match the scale limits
  y_range <- c(1, 1000)   # Match the scale limits

  legend_x <- 10^(log10(x_range[1]) + 0.05 * diff(log10(x_range)))  # Left side of plot
  legend_y_top <- 10^(log10(y_range[2]) - 0.05 * diff(log10(y_range)))  # Move down - closer to plot area

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

  # Add text labels with ggrepel - split by diagonal position for better placement
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_data) > 0) {
    # Split data by diagonal position for better label placement
    top_data$log_x <- log10(top_data$Census_OTU_Count)
    top_data$log_y <- log10(top_data$NCBI_Species_Count)

    # Add specific nudging for different types of families
    top_data$nudge_x_val <- 0
    top_data$nudge_y_val <- 0

    # First, apply general nudging for coverage families (those with "cov" labels)
    # Coverage families: nudge up and to the left
    coverage_families <- top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty
    top_data$nudge_x_val[coverage_families] <- -0.2  # Pull coverage families left
    top_data$nudge_y_val[coverage_families] <- 0.2   # Pull coverage families up

    # Then apply specific manual positioning adjustments for problematic families
    # Special handling for Thermococcaceae - different positioning based on label type
    thermococcaceae_mask <- grepl("Thermococcaceae", top_data$Taxon)
    thermococcaceae_coverage <- thermococcaceae_mask & top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty
    thermococcaceae_other <- thermococcaceae_mask & !thermococcaceae_coverage

    # Thermococcaceae with coverage annotation - pull left and higher up
    top_data$nudge_x_val[thermococcaceae_coverage] <- -0.3  # Pull Thermococcaceae coverage left
    top_data$nudge_y_val[thermococcaceae_coverage] <- 0.4   # Pull Thermococcaceae coverage higher up

    # Other Thermococcaceae (novelty or both) - pull right and up
    top_data$nudge_x_val[thermococcaceae_other] <- 0.2   # Pull other Thermococcaceae right
    top_data$nudge_y_val[thermococcaceae_other] <- 0.15  # Pull other Thermococcaceae up

    top_data$nudge_x_val[grepl("Methanomassiliicoccaceae", top_data$Taxon)] <- -0.3  # Move Methanomassiliicoccaceae slightly to the right (from -0.5)
    top_data$nudge_y_val[grepl("Methanomassiliicoccaceae", top_data$Taxon)] <- -0.02 # Move Methanomassiliicoccaceae slightly upward (from -0.05)

    top_data$nudge_x_val[grepl("Methanocellaceae", top_data$Taxon)] <- -0.3  # Pull Methanocellaceae left
    top_data$nudge_y_val[grepl("Methanocellaceae", top_data$Taxon)] <- 0     # Keep Methanocellaceae at same height

    top_data$nudge_x_val[grepl("Haloferacaceae", top_data$Taxon)] <- 0       # Keep Haloferacaceae at same horizontal position
    top_data$nudge_y_val[grepl("Haloferacaceae", top_data$Taxon)] <- 0.3     # Pull Haloferacaceae up

    top_data$nudge_x_val[grepl("Natrialbaceae", top_data$Taxon)] <- 0        # Keep Natrialbaceae at same horizontal position
    top_data$nudge_y_val[grepl("Natrialbaceae", top_data$Taxon)] <- 0.25     # Pull Natrialbaceae up

    top_data$nudge_x_val[grepl("Methanobacteriaceae", top_data$Taxon)] <- 0  # Keep Methanobacteriaceae at same horizontal position
    top_data$nudge_y_val[grepl("Methanobacteriaceae", top_data$Taxon)] <- 0.3 # Pull Methanobacteriaceae up

    top_data$nudge_x_val[grepl("Thermoplasmataceae", top_data$Taxon)] <- 0.2  # Move Thermoplasmataceae to the right
    top_data$nudge_y_val[grepl("Thermoplasmataceae", top_data$Taxon)] <- -0.2 # Move Thermoplasmataceae down

    top_data$nudge_x_val[grepl("Methanospirillaceae", top_data$Taxon)] <- 0   # Keep Methanospirillaceae at same horizontal position
    top_data$nudge_y_val[grepl("Methanospirillaceae", top_data$Taxon)] <- -0.15 # Move Methanospirillaceae down for better legibility

    # Points above diagonal (more NCBI species than Census OTUs)
    above_diagonal <- top_data[top_data$log_y > top_data$log_x, ]
    # Points below diagonal (more Census OTUs than NCBI species)
    below_diagonal <- top_data[top_data$log_y <= top_data$log_x, ]

    cat(paste("Above diagonal:", nrow(above_diagonal), "families\n"))
    cat(paste("Below diagonal:", nrow(below_diagonal), "families\n"))

    # Labels for points above diagonal - use specific nudging
    if (nrow(above_diagonal) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = above_diagonal,
        aes(label = paste0(Taxon, " (",
                          ifelse(Is_Top_Novelty, paste0(round(Novelty_Ratio, 1), "x nov"), ""),
                          ifelse(Is_Top_Coverage & !Is_Top_Novelty, paste0(round(Coverage_Factor, 1), "x cov"), ""),
                          ifelse(Is_Top_Novelty & Is_Top_Coverage, paste0(round(Novelty_Ratio, 1), "x"), ""),
                          ")")),
        color = "black", size = config$text_size, fontface = "bold",
        box.padding = 2.0, point.padding = 1.2, force = 5,
        max.overlaps = Inf, min.segment.length = 0.2,
        segment.color = "grey40", segment.linewidth = 0.6, segment.alpha = 0.9,
        # Use specific nudging for each family
        nudge_x = above_diagonal$nudge_x_val,
        nudge_y = above_diagonal$nudge_y_val,
        direction = "both", force_pull = 2, seed = 42
      )
    }

    # Labels for points below diagonal - use specific nudging
    if (nrow(below_diagonal) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = below_diagonal,
        aes(label = paste0(Taxon, " (",
                          ifelse(Is_Top_Novelty, paste0(round(Novelty_Ratio, 1), "x nov"), ""),
                          ifelse(Is_Top_Coverage & !Is_Top_Novelty, paste0(round(Coverage_Factor, 1), "x cov"), ""),
                          ifelse(Is_Top_Novelty & Is_Top_Coverage, paste0(round(Novelty_Ratio, 1), "x"), ""),
                          ")")),
        color = "black", size = config$text_size, fontface = "bold",
        box.padding = 2.0, point.padding = 1.2, force = 5,
        max.overlaps = Inf, min.segment.length = 0.2,
        segment.color = "grey40", segment.linewidth = 0.6, segment.alpha = 0.9,
        # Use specific nudging for each family
        nudge_x = below_diagonal$nudge_x_val,
        nudge_y = below_diagonal$nudge_y_val,
        direction = "both", force_pull = 2, seed = 123
      )
    }
  }

  return(p)
}

# Main function
main <- function() {
  cat("16S Archaea Family Focused Analysis\n")
  cat("==================================\n")

  # Load and process data
  data <- load_archaea_family_data()

  # Create plot
  plot <- create_archaea_family_scatter(data)

  # Save plot and summaries
  ggsave(file.path(config$output_dir, "16s_archaea_family_focused.png"), plot,
         width = config$plot_width, height = config$plot_height, dpi = config$dpi, bg = "white")

  write.csv(generate_summary(data, "coverage"),
            file.path(config$output_dir, "16s_archaea_family_top_coverage.csv"), row.names = FALSE)
  write.csv(generate_summary(data, "novelty"),
            file.path(config$output_dir, "16s_archaea_family_top_novelty.csv"), row.names = FALSE)

  cat(paste("✅ Analysis complete! Processed", nrow(data), "archaeal families\n"))
  cat(paste("   Top novelty taxa:", sum(data$Is_Top_Novelty), "\n"))
  cat(paste("   Top coverage taxa:", sum(data$Is_Top_Coverage), "\n"))
  cat("\nGenerated files:\n")
  cat("  - 16s_archaea_family_focused.png\n")
  cat("  - 16s_archaea_family_top_coverage.csv\n")
  cat("  - 16s_archaea_family_top_novelty.csv\n")
}

# Run analysis
if (!interactive()) main()
