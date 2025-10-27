#!/usr/bin/env Rscript
# 16S Mega Stacked Visual - Combined Bacteria and Archaea Scatter Plots
# Created: 2025-08-18
# Purpose: Create a unified 2x3 grid of 16S scatter plots with shared color scheme

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(ggrepel)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "mega_16s_stacked_visual.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Configuration for mega visual
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "16s_merged", "csv_results"),
  output_dir = file.path("final_visualizations"),
  source_data_dir = file.path("source_data"),  # New: for source data output
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 32,   # Increased width to use space freed by removing title/labels
  plot_height = 20,  # Increased height to use space freed by removing title/labels
  legend_width = 8,  # Width for separate legend files
  legend_height = 16, # Height for separate legend files (taller for more phyla)
  dpi = 300,
  top_n = 10,
  text_size = 11,    # Clean text size for optimal readability
  size_range = c(10, 22),  # Match 18S script circle sizes
  color_palette = "professional",
  use_weighted_ratio = TRUE
)

# Verify critical paths exist
if (!dir.exists(config$data_dir)) {
  stop(paste("Data directory not found:", config$data_dir))
}
if (!dir.exists(config$ncbi_data_dir)) {
  stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
}

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Custom archaea color palette - prioritizing phyla with environmental census data
get_archaea_colors <- function() {
  # Prioritize archaea phyla that actually appear in environmental census data
  # Replaced Halobacteriota (0 census count) with Nanoarchaeota (highest novelty factor: 117.7×)
  archaea_color_map <- c(
    "Euryarchaeota" = "#55ee79",        # Most common archaea phylum - Color 1 (bright green)
    "Nitrososphaerota" = "#ef9e17",     # Important ammonia-oxidizing archaea - Color 2
    "Thermoproteota" = "#f61e5d",       # Hyperthermophilic archaea - Color 3
    "Nanoarchaeota" = "#4ad9d5"         # Highest novelty factor (117.7×) - Color 4
  )
  return(archaea_color_map)
}

# Master color palette for bacteria phyla - assigned to phyla that actually appear in plots
get_master_color_palette <- function() {
  # Colors assigned to the 16 most important bacteria phyla that appear in visualizations
  # Based on highest novelty factors and significance, excluding common/filtered phyla
  bacteria_color_map <- c(
    # Top novel phyla (highest factors)
    "Abditibacteriota" = "#87bd4b",        # 36.71× novel - highest factor
    "Calditrichota" = "#e4a98f",           # 35.00× novel
    "Planctomycetota" = "#923052",         # 26.34× novel - unique cell wall-less
    "Thermomicrobiota" = "#6cd490",        # 25.55× novel
    "Acidobacteriota" = "#985739",         # 24.66× novel - important soil bacteria
    "Gemmatimonadota" = "#9d9719",         # 19.97× novel
    "Chloroflexota" = "#b2f3fd",           # 18.45× novel - green non-sulfur bacteria
    "Verrucomicrobiota" = "#8d0571",       # 16.85× novel - environmental bacteria
    "Armatimonadota" = "#e99953",          # 15.34× novel
    "Myxococcota" = "#fae692",             # 14.31× novel
    "Bdellovibrionota" = "#164c3f",        # 10.71× novel
    "Thermodesulfobacteriota" = "#01414d", # 6.88× novel
    "Spirochaetota" = "#dcb94b",           # 6.82× novel - spiral-shaped bacteria
    "Bacteroidota" = "#de724d",            # 4.25× novel - important gut bacteria
    # Top overrepresented phyla
    "Chlorobiota" = "#458f9c",             # 18.67× overrepresented
    "Chlamydiota" = "#5b1f06",             # 9.50× overrepresented - dark brown for obligate parasites
    "Ignavibacteriota" = "#f16127"         # 1063.64× overrepresented - blood orange for distinction
  )

  return(bacteria_color_map)
}

# Simple archaea detection function
is_archaea <- function(taxon_name) {
  archaea_patterns <- c("Methanobrevibacter", "Methanococcus", "Methanospirillum",
                       "Thermococcus", "Pyrococcus", "Sulfolobus", "Thermoplasma",
                       "Halobacterium", "Methanosarcina", "Archaeoglobus")
  any(grepl(paste(archaea_patterns, collapse = "|"), taxon_name, ignore.case = TRUE))
}

# Function to extract and save source data for plotted taxa
extract_plotted_source_data <- function(data, level, domain) {
  # Apply same filtering logic as plotting functions
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Overrepresented, ]
  meaningful_data <- top_data[top_data$Novelty_Ratio > 1.0 | top_data$Overrepresented_Factor > 1.0, ]

  # Apply same exclusions as plots
  excluded_phyla <- c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
  meaningful_data <- meaningful_data[!meaningful_data$Phylum %in% excluded_phyla, ]

  if (nrow(meaningful_data) == 0) {
    return(NULL)
  }

  # All levels use the "Taxon" column for the taxonomic name
  level_col_name <- "Taxon"

  # Debug: print column names and dimensions
  cat(paste("Debug - Level:", level, "Column name:", level_col_name, "\n"))
  cat(paste("Debug - Data dimensions:", nrow(meaningful_data), "x", ncol(meaningful_data), "\n"))
  cat(paste("Debug - Available columns:", paste(colnames(meaningful_data), collapse = ", "), "\n"))

  # Check if the column exists
  if (!level_col_name %in% colnames(meaningful_data)) {
    cat(paste("Warning: Column", level_col_name, "not found. Available columns:", paste(colnames(meaningful_data), collapse = ", "), "\n"))
    return(NULL)
  }

  # Create source data with all requested columns
  source_data <- data.frame(
    Taxonomic_Level = level,
    Domain = domain,
    Taxon = meaningful_data[[level_col_name]],  # Use proper column name
    Phylum = meaningful_data$Phylum,
    Novelty_Factor = round(meaningful_data$Novelty_Ratio, 3),
    Overrepresentation_Factor = round(meaningful_data$Overrepresented_Factor, 3),
    Census_OTU_Count = meaningful_data$Census_OTU_Count,
    NCBI_Species_Count = meaningful_data$NCBI_Species_Count,
    NCBI_Genome_Count = meaningful_data$NCBI_Genome_Count,
    Isolate_Count = meaningful_data$Isolate_Count,
    Is_Top_Novelty = meaningful_data$Is_Top_Novelty,
    Is_Top_Overrepresented = meaningful_data$Is_Top_Overrepresented,
    stringsAsFactors = FALSE
  )

  return(source_data)
}

# Standardize phylum names (consolidate subdivisions into main phyla)
standardize_phylum_names <- function(phylum_names) {
  # Map subdivision names to main phyla - only Methanobacteriota is actually part of Euryarchaeota
  phylum_mapping <- c(
    "Methanobacteriota" = "Euryarchaeota"
    # Removed incorrect mappings: Halobacteriota and Thermoplasmatota are separate phyla
  )

  # Apply mapping using vectorized approach
  standardized <- phylum_names
  for (old_name in names(phylum_mapping)) {
    standardized[standardized == old_name] <- phylum_mapping[old_name]
  }

  return(standardized)
}

# Load data function (adapted from individual scripts)
load_16s_data <- function(level, domain) {
  # Determine file name based on level and domain
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir, filename)
  
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")
  
  # Rename columns to match expected format with safety checks
  colnames(data)[1] <- "Taxon"

  # Check for required columns and create if missing
  required_cols <- c("census_otu_count", "ncbi_species_count", "ncbi_genome_count", "isolate_count")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns:", paste(missing_cols, collapse = ", "), "- setting to 0"))
    data[missing_cols] <- 0
  }

  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  
  # Filter by domain
  if ("domain" %in% colnames(data)) {
    data <- data %>% filter(domain == !!domain)
  } else {
    # Fallback to pattern matching for archaea
    archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
    if (domain == "Archaea") {
      data <- data %>% filter(Taxon %in% archaea_taxa)
    } else {
      data <- data %>% filter(!Taxon %in% archaea_taxa)
    }
  }
  
  # Filter for meaningful data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  
  # Add phylum information
  if (level == "phylum") {
    data$Phylum <- data$Taxon
  } else {
    # For family and genus, we need to get phylum mapping
    data$Phylum <- get_phylum_for_taxa(data$Taxon, level)
    # Filter out any taxa that couldn't be mapped to known phyla
    data <- data %>% filter(Phylum != "Unknown" & Phylum != "")
  }

  # Standardize phylum names (consolidate Methanobacteriota -> Euryarchaeota, etc.)
  data$Phylum <- standardize_phylum_names(data$Phylum)
  
  # Use pre-calculated factors from merger scripts
  data$Novelty_Ratio <- data$novelty_factor
  data$Overrepresented_Factor <- data$overrepresentation_factor

  # Weighted ratios if enabled
  if (config$use_weighted_ratio) {
    data$Weighted_Novelty <- data$Novelty_Ratio * log10(data$NCBI_Species_Count + 1)
    data$Weighted_Overrepresented <- data$Overrepresented_Factor * log10(data$NCBI_Species_Count + 1)
    ranking_novelty <- data$Weighted_Novelty
    ranking_overrepresented <- data$Weighted_Overrepresented
  } else {
    ranking_novelty <- data$Novelty_Ratio
    ranking_overrepresented <- data$Overrepresented_Factor
  }
  
  # Circle sizing - simplified calculation
  data$Genome_Isolate_Ratio <- pmax(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), 1)

  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  size_range <- range(data$Circle_Size_Raw, na.rm = TRUE)
  # Ensure minimum circle size of 3 and maximum of 12 for better visibility
  data$Circle_Size <- 3 + 9 * (data$Circle_Size_Raw - size_range[1]) / diff(size_range)
  
  # Identify top taxa - simplified approach
  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Overrepresented <- FALSE

  # Mark top novelty taxa
  high_novelty_mask <- ranking_novelty > 1.0
  if (any(high_novelty_mask)) {
    novelty_ranks <- rank(-ranking_novelty[high_novelty_mask])
    top_novelty_indices <- which(high_novelty_mask)[novelty_ranks <= config$top_n]
    data$Is_Top_Novelty[top_novelty_indices] <- TRUE
  }

  # Mark top overrepresented taxa
  high_overrepresented_mask <- ranking_overrepresented > 1.0
  if (any(high_overrepresented_mask)) {
    overrepresented_ranks <- rank(-ranking_overrepresented[high_overrepresented_mask])
    top_overrepresented_indices <- which(high_overrepresented_mask)[overrepresented_ranks <= config$top_n]
    data$Is_Top_Overrepresented[top_overrepresented_indices] <- TRUE
  }
  
  return(data)
}

# Get phylum mapping for family and genus levels
get_phylum_for_taxa <- function(taxa, level) {
  # Load NCBI data to get phylum mapping
  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) {
    return(rep("Unknown", length(taxa)))
  }

  ncbi_data <- read.csv(file_path, stringsAsFactors = FALSE)

  # Find the taxon name column
  possible_cols <- c("family_name", "genus_name", "family", "genus", "taxon", "Taxon")
  taxon_col <- intersect(possible_cols, colnames(ncbi_data))[1]

  if (is.na(taxon_col) || !all(c("lineage", "lineage_ranks") %in% colnames(ncbi_data))) {
    return(rep("Unknown", length(taxa)))
  }

  # Extract phylum from lineage - fixed indexing
  phylum_map <- setNames(rep("", nrow(ncbi_data)), ncbi_data[[taxon_col]])

  valid_rows <- !is.na(ncbi_data$lineage) & !is.na(ncbi_data$lineage_ranks)
  if (any(valid_rows)) {
    valid_indices <- which(valid_rows)
    lineage_list <- strsplit(ncbi_data$lineage[valid_rows], ";")
    ranks_list <- strsplit(ncbi_data$lineage_ranks[valid_rows], ";")

    for (j in seq_along(valid_indices)) {
      i <- valid_indices[j]
      phylum_idx <- which(ranks_list[[j]] == "phylum")
      if (length(phylum_idx) > 0 && phylum_idx[1] <= length(lineage_list[[j]])) {
        phylum_map[ncbi_data[[taxon_col]][i]] <- lineage_list[[j]][phylum_idx[1]]
      }
    }
  }

  # Map taxa to phyla - vectorized approach
  result <- phylum_map[taxa]
  result[is.na(result) | result == ""] <- "Unknown"
  names(result) <- NULL

  return(result)
}

# ggrepel parameters for family-level annotations (original spacing)
get_repel_params <- function(color = "black", nudge_x = 0, nudge_y = 0) {
  list(
    color = color,
    size = config$text_size,
    fontface = "bold",
    max.overlaps = Inf,
    max.iter = 5000,
    force = 10,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.linewidth = 0.3,
    segment.alpha = 0.7,
    direction = "both",
    nudge_x = nudge_x,
    nudge_y = nudge_y
  )
}

# ggrepel parameters for phylum-level annotations (much more spacing from circles)
get_repel_params_phylum <- function(color = "black", nudge_x = 0, nudge_y = 0) {
  list(
    color = color,
    size = config$text_size,
    fontface = "bold",
    max.overlaps = Inf,
    max.iter = 8000,               # More iterations for better positioning
    force = 25,                    # Much stronger force to push labels further
    box.padding = 1.2,             # Much more padding between labels
    point.padding = 1.5,           # Much more padding around points to move labels away
    min.segment.length = 0,
    segment.color = "grey40",
    segment.linewidth = 0.3,
    segment.alpha = 0.7,
    direction = "both",
    nudge_x = nudge_x,
    nudge_y = nudge_y
  )
}

# Create individual scatter plot (without titles and axis labels for stacking)
create_individual_scatter <- function(data, level, domain, all_phyla, master_colors) {
  # Base plot with same limits as 18S script
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15))) +  # Increased expansion for label space
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15)))

  # Get top data for annotations (restore top N filtering for cleaner plots)
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Overrepresented, ]

  # Filter meaningful data for annotations only - use factor > 1.0 threshold
  meaningful_data <- top_data[top_data$Novelty_Ratio > 1.0 | top_data$Overrepresented_Factor > 1.0, ]

  # Apply exclusions for cleaner plots
  excluded_phyla <- c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
  meaningful_data <- meaningful_data[!meaningful_data$Phylum %in% excluded_phyla, ]

  # Background points first (draw behind meaningful data)
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Overrepresented, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "lightgray", alpha = 0.3)
  }

  if (nrow(meaningful_data) > 0) {
    # Get phyla present in this plot
    plot_phyla <- unique(meaningful_data$Phylum[meaningful_data$Phylum != "Unknown"])
    plot_phyla <- sort(plot_phyla)

    # Map to color scheme - allow all qualifying archaea phyla
    if (domain == "Archaea") {
      # Use custom archaea color map with fallback for unmapped phyla
      archaea_color_map <- get_archaea_colors()
      fallback_archaea_colors <- c("#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222")

      phyla_colors <- character(length(plot_phyla))
      for (i in seq_along(plot_phyla)) {
        phylum <- plot_phyla[i]
        if (phylum %in% names(archaea_color_map)) {
          phyla_colors[i] <- archaea_color_map[phylum]
        } else {
          # Use fallback colors for unmapped archaea phyla
          fallback_index <- ((i - 1) %% length(fallback_archaea_colors)) + 1
          phyla_colors[i] <- fallback_archaea_colors[fallback_index]
        }
      }
    } else {
      # Use bacteria color assignments with extended fallback for remaining phyla
      # Additional colors for bacteria phyla not in main palette
      extended_fallback_colors <- c(
        "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
        "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
        "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F"  # Fixed duplicate #85C1E9 → #C39BD3
      )

      phyla_colors <- ifelse(plot_phyla %in% names(master_colors),
                            master_colors[plot_phyla],
                            extended_fallback_colors[((seq_along(plot_phyla) - 1) %% length(extended_fallback_colors)) + 1])
    }

    # Add black outline layer
    p <- p + geom_point(data = meaningful_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.3)

    # Add colored fill layer
    p <- p + geom_point(data = meaningful_data,
                       aes(size = Circle_Size * 0.8, color = factor(Phylum, levels = plot_phyla)),
                       alpha = 0.9) +
      scale_color_manual(values = setNames(phyla_colors, plot_phyla),
                        name = "Phyla", guide = "none")  # Remove individual legends
  }

  # Minimal styling for stacking
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "darkblue", linetype = "dashed", alpha = 0.9, size = 2) +
    theme_minimal() +
    theme(
      # Remove titles and axis labels for stacking
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 32, face = "bold", color = "grey50"),  # X-axis text
      axis.text.y = element_text(size = 32, face = "bold", color = "grey50"), # Y-axis text
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line.x = element_blank(),  # Remove x-axis line
      axis.line.y = element_blank(),  # Remove y-axis line
      plot.margin = margin(1, 1, 1, 1),
      # Hide individual plot legends for master layout
      legend.position = "none"
    )

  # Add minimal text annotations for top taxa
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(meaningful_data) > 0) {
    # Create labels with factor numbers and nov/ovr designation for top taxa only
    # Use the original logic based on top novelty/overrepresented flags
    meaningful_data$designation <- ifelse(meaningful_data$Is_Top_Novelty, "Nov",
                                         ifelse(meaningful_data$Is_Top_Overrepresented & meaningful_data$Overrepresented_Factor > 1.0, "Ovr", "Nov"))

    # Always display the original factor values, regardless of weighting used for ranking
    factor_value <- ifelse(meaningful_data$Is_Top_Novelty,
                          meaningful_data$Novelty_Ratio,
                          meaningful_data$Overrepresented_Factor)

    # Final filter: only show taxa where the displayed factor value is > 1.0
    meaningful_data <- meaningful_data[factor_value > 1.0, ]
    factor_value <- factor_value[factor_value > 1.0]

    # Add asterisks for non-family taxa in family-level plots
    if (level == "family") {
      family_suffixes <- c("aceae", "idae", "ales", "ineae", "inae", "eae")
      suffix_pattern <- paste0("(", paste(family_suffixes, collapse = "|"), ")$")
      meaningful_data$is_family <- grepl(suffix_pattern, meaningful_data$Taxon, ignore.case = TRUE)
      meaningful_data$taxon_display <- ifelse(meaningful_data$is_family,
                                            meaningful_data$Taxon,
                                            paste0(meaningful_data$Taxon, "*"))
    } else {
      meaningful_data$taxon_display <- meaningful_data$Taxon
    }

    # Create labels - show full annotations for family level, only factors for phylum level
    if (level == "family") {
      meaningful_data$label_text <- paste0(meaningful_data$taxon_display, " (",
                                          round(factor_value, 1), "×)")
    } else {
      # For phylum level, show only factor values
      meaningful_data$label_text <- paste0("(", round(factor_value, 1), "×)")
    }

    # Simple consistent annotation settings with special positioning for Archaea family
    # All annotations in black color as discussed

    # Special positioning for Archaea family level to move away from central diagonal
    if (domain == "Archaea" && level == "family") {
      # Add all taxa in black - use different nudging for overrepresented vs novel
      overrepresented_data <- meaningful_data[meaningful_data$designation == "Ovr", ]
      novel_data <- meaningful_data[meaningful_data$designation == "Nov", ]

      # Add overrepresented taxa in black - push further up and left away from diagonal
      if (nrow(overrepresented_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -0.8, nudge_y = 0.8)
        repel_params$force <- 15  # Stronger repulsion force
        repel_params$box.padding <- 0.8  # More space between labels
        repel_params$point.padding <- 0.6  # More space around points
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = overrepresented_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Add novel taxa in black - push further down and right away from diagonal
      if (nrow(novel_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = 0.8, nudge_y = -0.8)
        repel_params$force <- 15  # Stronger repulsion force
        repel_params$box.padding <- 0.8  # More space between labels
        repel_params$point.padding <- 0.6  # More space around points
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = novel_data, aes(label = label_text)),
          repel_params
        ))
      }
    } else if (domain == "Bacteria" && level == "family") {
      # Enhanced spacing for Bacteria family level to prevent overlapping
      # Special handling for specific families that need custom positioning
      thermosynechococcaceae_data <- meaningful_data[grepl("Thermosynechococcaceae", meaningful_data$Taxon), ]
      chlorobiaceae_data <- meaningful_data[grepl("Chlorobiaceae", meaningful_data$Taxon), ]
      gloeomargarita_data <- meaningful_data[grepl("Gloeomargarita", meaningful_data$Taxon), ]
      physcisphaeraceae_data <- meaningful_data[grepl("Physcisphaeraceae", meaningful_data$Taxon), ]
      streptomycetaceae_data <- meaningful_data[grepl("Streptomycetaceae", meaningful_data$Taxon), ]
      cthoniobacteraceae_data <- meaningful_data[grepl("Cthoniobacteraceae", meaningful_data$Taxon), ]
      rhodothermaceae_data <- meaningful_data[grepl("Rhodothermaceae", meaningful_data$Taxon), ]
      merosscillaceae_data <- meaningful_data[grepl("Merosscillaceae", meaningful_data$Taxon), ]
      gemmatimonadaceae_data <- meaningful_data[grepl("Gemmatimonadaceae", meaningful_data$Taxon), ]
      other_data <- meaningful_data[!grepl("Thermosynechococcaceae|Chlorobiaceae|Gloeomargarita|Physcisphaeraceae|Streptomycetaceae|Cthoniobacteraceae|Rhodothermaceae|Merosscillaceae|Gemmatimonadaceae", meaningful_data$Taxon), ]

      # Handle Thermosynechococcaceae separately - move way up
      if (nrow(thermosynechococcaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = 0, nudge_y = 1.2)
        repel_params$force <- 15
        repel_params$max.iter <- 10000
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = thermosynechococcaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Chlorobiaceae separately - move up
      if (nrow(chlorobiaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = 0, nudge_y = 0.8)
        repel_params$force <- 15
        repel_params$max.iter <- 10000
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = chlorobiaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Gloeomargarita separately - move left to avoid overlap
      if (nrow(gloeomargarita_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -0.6, nudge_y = 0.2)
        repel_params$force <- 20
        repel_params$max.iter <- 12000
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = gloeomargarita_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Physcisphaeraceae separately - move up and right to spread out annotations
      if (nrow(physcisphaeraceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = 0.9, nudge_y = 0.7)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = physcisphaeraceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Streptomycetaceae separately - move far to the left
      if (nrow(streptomycetaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -1.2, nudge_y = 0)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = streptomycetaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Cthoniobacteraceae separately - move up and left to avoid Physcisphaeraceae
      if (nrow(cthoniobacteraceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -0.4, nudge_y = 0.8)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = cthoniobacteraceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Rhodothermaceae separately - move up and right
      if (nrow(rhodothermaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = 0.6, nudge_y = 0.9)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = rhodothermaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Merosscillaceae separately - move down and left
      if (nrow(merosscillaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -0.6, nudge_y = -0.8)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = merosscillaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle Gemmatimonadaceae separately - move up and slightly left
      if (nrow(gemmatimonadaceae_data) > 0) {
        repel_params <- get_repel_params(color = "black", nudge_x = -0.2, nudge_y = 1.0)
        repel_params$force <- 25
        repel_params$max.iter <- 15000
        repel_params$box.padding <- 1.0
        repel_params$point.padding <- 0.8
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = gemmatimonadaceae_data, aes(label = label_text)),
          repel_params
        ))
      }

      # Handle remaining taxa with enhanced spacing - all in black
      if (nrow(other_data) > 0) {
        repel_params <- get_repel_params(color = "black")
        repel_params$force <- 25  # Much stronger repulsion force
        repel_params$max.iter <- 15000  # More iterations for better positioning
        repel_params$box.padding <- 1.0  # Larger padding between labels
        repel_params$point.padding <- 0.8  # Larger padding around points
        p <- p + do.call(ggrepel::geom_text_repel, c(
          list(data = other_data, aes(label = label_text)),
          repel_params
        ))
      }
    } else {
      # Standard positioning for all other cases (Bacteria phylum, Archaea phylum)
      # Use different spacing based on taxonomic level
      if (level == "phylum") {
        repel_params <- get_repel_params_phylum(color = "black")
      } else {
        repel_params <- get_repel_params(color = "black")
      }
      p <- p + do.call(ggrepel::geom_text_repel, c(
        list(data = meaningful_data, aes(label = label_text)),
        repel_params
      ))
    }
  }



  return(p)
}

# Classify phyla by domain
classify_phyla_by_domain <- function(phyla) {
  # Known archaeal phyla
  archaeal_phyla <- c("Euryarchaeota", "Nitrososphaerota", "Thermoproteota",
                     "Nanoarchaeota", "Promethearchaeota", "Halobacteriota",
                     "Thermoplasmatota", "Methanobacteriota")

  # Classify each phylum
  domains <- ifelse(phyla %in% archaeal_phyla, "Archaea", "Bacteria")
  return(domains)
}

# Sort phyla by color gradient for visually appealing legend
sort_by_color_gradient <- function(phyla, colors) {
  if (length(phyla) == 0) return(phyla)

  # Convert hex colors to HSV for better sorting
  rgb_colors <- col2rgb(colors)
  hsv_colors <- rgb2hsv(rgb_colors)

  # Sort by hue first, then by saturation, then by value (brightness)
  # This creates a rainbow-like gradient arrangement
  color_order <- order(hsv_colors[1, ], hsv_colors[2, ], hsv_colors[3, ])

  return(phyla[color_order])
}

# Create master legend for all phyla with domain classification
create_master_legend <- function(all_phyla, master_colors) {
  # Classify phyla by domain
  domains <- classify_phyla_by_domain(all_phyla)

  # Separate phyla by domain
  archaea_phyla_unsorted <- all_phyla[domains == "Archaea"]
  bacteria_phyla_unsorted <- all_phyla[domains == "Bacteria"]

  # Create ordered list with domain headers - simplified approach
  legend_items <- c()
  legend_colors <- c()
  legend_shapes <- c()

  # Add Archaea section - include all qualifying archaea phyla sorted by color gradient
  if (length(archaea_phyla_unsorted) > 0) {
    archaea_color_map <- get_archaea_colors()
    fallback_archaea_colors <- c("#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222")

    # Get colors for all archaea phyla first
    archaea_colors_unsorted <- character(length(archaea_phyla_unsorted))
    for (i in seq_along(archaea_phyla_unsorted)) {
      phylum <- archaea_phyla_unsorted[i]
      if (phylum %in% names(archaea_color_map)) {
        archaea_colors_unsorted[i] <- archaea_color_map[phylum]
      } else {
        # Use fallback colors for unmapped archaea phyla
        fallback_index <- ((i - 1) %% length(fallback_archaea_colors)) + 1
        archaea_colors_unsorted[i] <- fallback_archaea_colors[fallback_index]
      }
    }

    # Sort by color gradient
    archaea_phyla <- sort_by_color_gradient(archaea_phyla_unsorted, archaea_colors_unsorted)
    archaea_colors <- archaea_colors_unsorted[match(archaea_phyla, archaea_phyla_unsorted)]

    legend_items <- c(legend_items, "ARCHAEA", archaea_phyla)
    legend_colors <- c(legend_colors, "white", archaea_colors)
    legend_shapes <- c(legend_shapes, NA, rep(15, length(archaea_phyla)))
  }

  # Add Bacteria section - sorted by color gradient
  if (length(bacteria_phyla_unsorted) > 0) {
    # Extended fallback colors for bacteria phyla not in main palette
    extended_fallback_colors <- c(
      "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
      "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
      "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F"  # Fixed duplicate #85C1E9 → #C39BD3
    )

    # Get colors for all bacteria phyla first
    bacteria_colors_unsorted <- ifelse(bacteria_phyla_unsorted %in% names(master_colors),
                                      master_colors[bacteria_phyla_unsorted],
                                      extended_fallback_colors[((seq_along(bacteria_phyla_unsorted) - 1) %% length(extended_fallback_colors)) + 1])

    # Sort by color gradient
    bacteria_phyla <- sort_by_color_gradient(bacteria_phyla_unsorted, bacteria_colors_unsorted)
    bacteria_colors <- bacteria_colors_unsorted[match(bacteria_phyla, bacteria_phyla_unsorted)]

    legend_items <- c(legend_items, "BACTERIA", bacteria_phyla)
    legend_colors <- c(legend_colors, "white", bacteria_colors)
    legend_shapes <- c(legend_shapes, NA, rep(15, length(bacteria_phyla)))
  }

  # Create separate legends for Archaea and Bacteria
  legends <- list()

  # Create Archaea legend if archaea phyla exist
  if (length(archaea_phyla_unsorted) > 0) {
    archaea_legend_items <- c("ARCHAEA", archaea_phyla)
    archaea_legend_colors <- c("white", archaea_colors)
    archaea_legend_shapes <- c(NA, rep(15, length(archaea_phyla)))

    archaea_dummy_data <- data.frame(
      x = rep(1, length(archaea_legend_items)),
      y = rep(1, length(archaea_legend_items)),
      Item = factor(archaea_legend_items, levels = archaea_legend_items)
    )

    archaea_legend_plot <- ggplot(archaea_dummy_data, aes(x = x, y = y, color = Item)) +
      geom_point(size = 8, alpha = 0.9) +
      scale_color_manual(values = setNames(archaea_legend_colors, archaea_legend_items),
                        name = NULL,
                        guide = guide_legend(
                          override.aes = list(size = c(NA, rep(8, length(archaea_legend_items)-1)),
                                            alpha = c(0, rep(1, length(archaea_legend_items)-1)),
                                            shape = archaea_legend_shapes),
                          keywidth = unit(1.5, "cm"),
                          keyheight = unit(1.0, "cm"),
                          ncol = 1
                        )) +
      theme_void() +
      theme(
        legend.text = element_text(size = 24, face = c("bold", rep("plain", length(archaea_legend_items)-1))),
        legend.title = element_blank(),
        legend.key.size = unit(2.0, "cm"),
        legend.spacing.y = unit(0.8, "cm"),
        plot.margin = margin(0, 0, 0, 0)
      )

    legends$archaea <- get_legend(archaea_legend_plot)
  }

  # Create Bacteria legend if bacteria phyla exist
  if (length(bacteria_phyla_unsorted) > 0) {
    bacteria_legend_items <- c("BACTERIA", bacteria_phyla)
    bacteria_legend_colors <- c("white", bacteria_colors)
    bacteria_legend_shapes <- c(NA, rep(15, length(bacteria_phyla)))

    bacteria_dummy_data <- data.frame(
      x = rep(1, length(bacteria_legend_items)),
      y = rep(1, length(bacteria_legend_items)),
      Item = factor(bacteria_legend_items, levels = bacteria_legend_items)
    )

    bacteria_legend_plot <- ggplot(bacteria_dummy_data, aes(x = x, y = y, color = Item)) +
      geom_point(size = 8, alpha = 0.9) +
      scale_color_manual(values = setNames(bacteria_legend_colors, bacteria_legend_items),
                        name = NULL,
                        guide = guide_legend(
                          override.aes = list(size = c(NA, rep(8, length(bacteria_legend_items)-1)),
                                            alpha = c(0, rep(1, length(bacteria_legend_items)-1)),
                                            shape = bacteria_legend_shapes),
                          keywidth = unit(1.5, "cm"),
                          keyheight = unit(1.0, "cm"),
                          ncol = 1
                        )) +
      theme_void() +
      theme(
        legend.text = element_text(size = 24, face = c("bold", rep("plain", length(bacteria_legend_items)-1))),
        legend.title = element_blank(),
        legend.key.size = unit(2.0, "cm"),
        legend.spacing.y = unit(0.8, "cm"),
        plot.margin = margin(0, 0, 0, 0)
      )

    legends$bacteria <- get_legend(bacteria_legend_plot)
  }

  return(legends)
}

# Main function to create mega visual
main <- function() {
  cat("16S Mega Stacked Visual Creation\n")
  cat("===============================\n")

  # Define the grid structure (removed genus)
  levels <- c("phylum", "family")
  domains <- c("Bacteria", "Archaea")

  # Collect all data to determine master phyla list
  all_data <- list()
  all_phyla_set <- character(0)

  cat("Loading all data to determine master phyla list...\n")
  for (level in levels) {
    for (domain in domains) {
      cat(paste("Loading", domain, level, "data...\n"))
      data <- load_16s_data(level, domain)
      all_data[[paste(domain, level, sep = "_")]] <- data

      # Get top data for annotations (same filtering as individual plots)
      top_data <- data[data$Is_Top_Novelty | data$Is_Top_Overrepresented, ]

      # Filter meaningful data for annotations only - use factor > 1.0 threshold
      meaningful_data <- top_data[top_data$Novelty_Ratio > 1.0 | top_data$Overrepresented_Factor > 1.0, ]

      # Apply exclusions for cleaner plots (same as individual scripts)
      excluded_phyla <- c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
      meaningful_data <- meaningful_data[!meaningful_data$Phylum %in% excluded_phyla, ]

      # Only collect phyla from taxa that actually get plotted with colors and annotations
      meaningful_phyla <- unique(meaningful_data$Phylum[meaningful_data$Phylum != "Unknown" & meaningful_data$Phylum != ""])
      all_phyla_set <- unique(c(all_phyla_set, meaningful_phyla))
    }
  }

  # Sort phyla alphabetically for consistent coloring
  all_phyla <- sort(all_phyla_set)
  master_colors <- get_master_color_palette()

  cat(paste("Found", length(all_phyla), "unique phyla across all datasets\n"))
  cat(paste("Phyla:", paste(all_phyla, collapse = ", "), "\n"))

  # Extract and save source data organized by rank and factor type
  cat("Extracting source data for plotted taxa...\n")
  if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

  # Process phylum and family levels only
  for (level in c("phylum", "family")) {
    # Combine data from both domains for this level
    combined_level_data <- list()

    for (domain in c("Bacteria", "Archaea")) {
      data_key <- paste(domain, level, sep = "_")
      if (data_key %in% names(all_data)) {
        source_data <- extract_plotted_source_data(all_data[[data_key]], level, domain)
        if (!is.null(source_data)) {
          combined_level_data[[domain]] <- source_data
          cat(paste("  ", domain, level, ":", nrow(source_data), "plotted taxa\n"))
        }
      }
    }

    # Combine domains and split by factor type
    if (length(combined_level_data) > 0) {
      all_level_data <- do.call(rbind, combined_level_data)

      # Split into novelty and overrepresented
      novelty_data <- all_level_data[all_level_data$Is_Top_Novelty == TRUE, ]
      overrepresented_data <- all_level_data[all_level_data$Is_Top_Overrepresented == TRUE, ]

      # Save separate files
      if (nrow(novelty_data) > 0) {
        novelty_file <- file.path(config$source_data_dir, paste0("16s_", level, "_novelty_source_data.csv"))
        write.csv(novelty_data, novelty_file, row.names = FALSE)
        cat(paste("  Novelty data saved:", novelty_file, "(", nrow(novelty_data), "taxa)\n"))
      }

      if (nrow(overrepresented_data) > 0) {
        overrepresented_file <- file.path(config$source_data_dir, paste0("16s_", level, "_overrepresented_source_data.csv"))
        write.csv(overrepresented_data, overrepresented_file, row.names = FALSE)
        cat(paste("  Overrepresented data saved:", overrepresented_file, "(", nrow(overrepresented_data), "taxa)\n"))
      }
    }
  }

  # Create individual plots
  plots <- list()
  plot_names <- character(0)

  cat("Creating individual scatter plots...\n")
  for (i in 1:length(levels)) {
    level <- levels[i]
    for (j in 1:length(domains)) {
      domain <- domains[j]
      data_key <- paste(domain, level, sep = "_")
      data <- all_data[[data_key]]

      cat(paste("Creating plot for", domain, level, "...\n"))
      plot <- create_individual_scatter(data, level, domain, all_phyla, master_colors)

      plot_key <- paste(domain, level, sep = "_")
      plots[[plot_key]] <- plot
      plot_names <- c(plot_names, plot_key)
    }
  }

  # Create and save separate legends
  cat("Creating and saving separate legends...\n")
  legends <- create_master_legend(all_phyla, master_colors)

  # Save Archaea legend if it exists
  if (!is.null(legends$archaea)) {
    archaea_pdf_file <- file.path(config$output_dir, "16s_archaea_legend.pdf")
    archaea_png_file <- file.path(config$output_dir, "16s_archaea_legend.png")

    cat(paste("Saving Archaea legend to:", archaea_png_file, "\n"))
    ggsave(archaea_png_file, legends$archaea,
           width = config$legend_width * 0.6, height = config$legend_height * 0.4,
           dpi = config$dpi, bg = "white", limitsize = FALSE,
           units = "in")

    ggsave(archaea_pdf_file, legends$archaea,
           width = config$legend_width * 0.6, height = config$legend_height * 0.4,
           dpi = config$dpi, bg = "white", limitsize = FALSE,
           units = "in")
  }

  # Save Bacteria legend if it exists
  if (!is.null(legends$bacteria)) {
    bacteria_pdf_file <- file.path(config$output_dir, "16s_bacteria_legend.pdf")
    bacteria_png_file <- file.path(config$output_dir, "16s_bacteria_legend.png")

    cat(paste("Saving Bacteria legend to:", bacteria_png_file, "\n"))
    ggsave(bacteria_png_file, legends$bacteria,
           width = config$legend_width, height = config$legend_height,
           dpi = config$dpi, bg = "white", limitsize = FALSE,
           units = "in")

    ggsave(bacteria_pdf_file, legends$bacteria,
           width = config$legend_width, height = config$legend_height,
           dpi = config$dpi, bg = "white", limitsize = FALSE,
           units = "in")
  }

  # Arrange plots in 2x2 grid with individual frames (2 rows for levels, 2 columns for domains) - NO LEGEND
  cat("Arranging plots in grid with individual frames...\n")

  # Create row labels (left side)
  row_labels <- list()
  for (i in 1:length(levels)) {
    level_label <- switch(levels[i],
                         "phylum" = "Phyla",
                         "family" = "Family")

    row_labels[[i]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = level_label,
               size = 12, fontface = "bold", color = "grey30", angle = 90) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  # Create column labels (top)
  col_labels <- list()
  for (j in 1:length(domains)) {
    col_labels[[j]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = domains[j],
               size = 12, fontface = "bold", color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  # Create the grid layout with spacing between plots
  # Top row: empty corner + column headers with spacer
  top_row <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    col_labels[[1]],          # Bacteria header
    ggplot() + theme_void(),  # Vertical spacer
    col_labels[[2]],          # Archaea header
    ncol = 4, rel_widths = c(0.1, 1, 0.03, 1)  # Reduced horizontal spacing to match vertical visually
  )

  # Data rows: row label + 2 framed plots with spacer between them
  data_rows <- list()
  for (i in 1:length(levels)) {
    level <- levels[i]
    bacteria_plot <- plots[[paste("Bacteria", level, sep = "_")]]
    archaea_plot <- plots[[paste("Archaea", level, sep = "_")]]

    # Add black frames around individual plots
    framed_bacteria_plot <- ggdraw(bacteria_plot) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))
    framed_archaea_plot <- ggdraw(archaea_plot) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))

    data_rows[[i]] <- plot_grid(
      row_labels[[i]],
      framed_bacteria_plot,
      ggplot() + theme_void(),  # Vertical spacer between plots
      framed_archaea_plot,
      ncol = 4, rel_widths = c(0.1, 1, 0.03, 1)  # Reduced horizontal spacing to match vertical visually
    )
  }

  # Combine all rows with spacing between plots - NO LEGEND
  main_grid <- plot_grid(
    top_row,
    data_rows[[1]],           # Phylum row
    ggplot() + theme_void(),  # Horizontal spacer
    data_rows[[2]],           # Family row
    ncol = 1, rel_heights = c(0.1, 1, 0.05, 1)  # Same spacing as horizontal (0.05)
  )

  # Final plot without title or axis labels (to be added elsewhere)
  complete_plot <- main_grid

  # Save the mega visual
  output_file <- file.path(config$output_dir, "16s_mega_stacked_visual.png")
  cat(paste("Saving mega visual to:", output_file, "\n"))

  ggsave(output_file, complete_plot,
         width = config$plot_width, height = config$plot_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  cat("✅ Mega visual creation complete!\n")
  cat(paste("   Generated main plot:", output_file, "\n"))
  if (!is.null(legends$archaea)) {
    cat(paste("   Generated Archaea legend PNG:", file.path(config$output_dir, "16s_archaea_legend.png"), "\n"))
    cat(paste("   Generated Archaea legend PDF:", file.path(config$output_dir, "16s_archaea_legend.pdf"), "\n"))
  }
  if (!is.null(legends$bacteria)) {
    cat(paste("   Generated Bacteria legend PNG:", file.path(config$output_dir, "16s_bacteria_legend.png"), "\n"))
    cat(paste("   Generated Bacteria legend PDF:", file.path(config$output_dir, "16s_bacteria_legend.pdf"), "\n"))
  }
  cat(paste("   Main plot dimensions:", config$plot_width, "x", config$plot_height, "inches\n"))
  cat(paste("   Total phyla in legends:", length(all_phyla), "\n"))
}

# Run the mega visual creation
if (!interactive()) main()
