#!/usr/bin/env Rscript
# Bacteria-Specific Stacked Visual
# Created: 2025-10-26
# Purpose: Create stacked bacterial plots (phylum + family) with custom phyla color assignments

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(grid)
  library(cowplot)
  library(ggrepel)
  library(ggtext)
})

# ============================================================================
# BACTERIA PHYLA COLOR ASSIGNMENT SECTION
# ============================================================================
# Assign specific colors to bacterial phyla from the BAC color scheme
# BAC color scheme: ["#548877","#82dd83","#ff7200","#ccdd75","#a58c26","#e5b726",
#                   "#b1caae","#c37231","#68a200","#4bca8f","#dbb579","#ffb06b",
#                   "#386529","#74c4d2","#7cce54","#b8e08f","#45ffea","#994417"]

# BAC color scheme: ["#548877","#82dd83","#ff7200","#ccdd75","#a58c26","#e5b726",
#                   "#b1caae","#c37231","#68a200","#4bca8f","#dbb579","#ffb06b",
#                   "#386529","#74c4d2","#7cce54","#b8e08f","#45ffea","#994417"]

# Master color palette for bacteria (using your BAC colors)
get_master_bacteria_color_palette <- function() {
  bacteria_colors <- c("#548877", "#82dd83", "#ff7200", "#ccdd75", "#a58c26",
                      "#e5b726", "#b1caae", "#c37231", "#68a200", "#4bca8f",
                      "#dbb579", "#ffb06b", "#386529", "#74c4d2", "#7cce54",
                      "#b8e08f", "#45ffea", "#994417")

  return(bacteria_colors)
}

# ============================================================================
# AESTHETIC CONFIGURATION SECTION (MATCHING COMPREHENSIVE SCRIPT)
# ============================================================================

# Circle aesthetics configuration (from comprehensive script)
CIRCLE_AESTHETICS <- list(
  shape = 21,           # Fillable circles
  stroke = 0.6,         # Thick black outline
  alpha = 0.9,          # Transparency
  size_multiplier = 0.8,
  bg_color = "lightgray",
  bg_fill = "lightgray",
  bg_alpha = 0.3,
  outline_color = "black",
  outline_fill = "black"
)

# Plot theme configuration (from comprehensive script)
PLOT_THEME <- list(
  axis_text_size = 32,
  axis_text_color = "grey50",
  grid_major_color = "grey90",
  grid_major_size = 0.3,
  diagonal_color = "black",
  diagonal_type = "dashed",
  diagonal_alpha = 0.9,
  diagonal_size = 3,
  plot_margin = margin(1, 1, 1, 1)
)

# Configuration for stacked bacteria visual
config <- list(
  data_dir_16s = file.path("..", "..", "Eukcensus_merge", "16s_merged", "csv_results"),
  output_dir = file.path("bacteria_stacked_plots"),
  source_data_dir = file.path("source_data"),
  ncbi_data_dir = file.path("..", "..", "00ncbi_parse", "csv_ncbi"),
  plot_width = 18,   # Width for stacked layout
  plot_height = 20,  # Height for 2-row stack (phylum + family)
  dpi = 300,
  top_n = 10,
  text_size = 11,
  size_range = c(10, 22),
  domain = "Bacteria"
)

# Layout configuration for stacked visual
LAYOUT_CONFIG <- list(
  row_spacing = 0.05,     # Space between phylum and family rows
  row_heights = c(0.1, 1, 0.05, 1),  # Header, phylum, gap, family
  legend_width_ratio = 0.15,  # 15% for legend
  main_plot_ratio = 0.85      # 85% for main plots
)

# Create output directories
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

# Function to calculate inverted circle sizes (from comprehensive script)
calculate_isolate_circle_size <- function(isolate_count, ncbi_genome_count) {
  isolate_percentage <- ifelse(ncbi_genome_count > 0,
                              (isolate_count / ncbi_genome_count) * 100,
                              0)
  
  circle_size <- case_when(
    isolate_percentage == 0 ~ 30,
    isolate_percentage < 1 ~ 28,
    isolate_percentage < 5 ~ 26,
    isolate_percentage < 10 ~ 24,
    isolate_percentage < 25 ~ 22,
    isolate_percentage < 50 ~ 20,
    isolate_percentage < 75 ~ 15,
    TRUE ~ 10
  )
  
  return(circle_size)
}

# Load bacteria data function
load_bacteria_data <- function(level) {
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)
  
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")
  
  # Filter for Bacteria only
  cat(paste("Before filtering:", nrow(data), "total taxa\n"))
  data <- data %>% filter(domain == "Bacteria")
  cat(paste("After filtering for Bacteria:", nrow(data), "taxa\n"))
  
  if (nrow(data) == 0) {
    cat("No bacterial data found\n")
    return(data.frame())
  }
  
  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  
  # Filter for meaningful data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  
  # Use pre-calculated factors directly (no weighted ratios)
  data$Novelty_Ratio <- data$novelty_factor
  data$Coverage_Factor <- data$overrepresentation_factor

  # Use original factors for ranking (no weighting)
  ranking_novelty <- data$Novelty_Ratio
  ranking_coverage <- data$Coverage_Factor

  # Circle sizing
  data$Circle_Size <- calculate_isolate_circle_size(data$Isolate_Count, data$NCBI_Genome_Count)

  # Add phylum information
  data <- add_phylum_info(data, level)

  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE
  
  # Filter for meaningful data (same as focus_taxa_scatter scripts)
  high_novelty_data <- data[data$Novelty_Ratio > 1.0, ]
  high_coverage_data <- data[data$Coverage_Factor > 1.0, ]

  # Mark top novelty taxa (using focus_taxa_scatter logic)
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
    cat(paste("Selected top", length(top_novelty_taxa), "novelty taxa:", paste(top_novelty_taxa, collapse = ", "), "\n"))
  }

  # Mark top coverage taxa (using focus_taxa_scatter logic)
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
    cat(paste("Selected top", length(top_coverage_taxa), "coverage taxa:", paste(top_coverage_taxa, collapse = ", "), "\n"))
  }

  cat(paste("Final bacterial", level, "data:", nrow(data), "taxa,",
            sum(data$Is_Top_Novelty), "top novelty,",
            sum(data$Is_Top_Coverage), "top coverage\n"))
  
  return(data)
}

# Add phylum information function
add_phylum_info <- function(data, level) {
  if (level == "phylum") {
    data$Phylum <- data$Taxon
  } else {
    data$Phylum <- get_phylum_for_taxa(data$Taxon, level)
    data <- data %>% filter(Phylum != "Unknown" & Phylum != "")
  }
  
  data$Phylum <- standardize_phylum_names(data$Phylum)
  return(data)
}

# Get phylum mapping for family level
get_phylum_for_taxa <- function(taxa, level) {
  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) {
    return(rep("Unknown", length(taxa)))
  }
  
  ncbi_data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  possible_cols <- c("family_name", "genus_name", "family", "genus", "taxon", "Taxon")
  taxon_col <- intersect(possible_cols, colnames(ncbi_data))[1]
  
  if (is.na(taxon_col) || !all(c("lineage", "lineage_ranks") %in% colnames(ncbi_data))) {
    return(rep("Unknown", length(taxa)))
  }
  
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
  
  result <- phylum_map[taxa]
  result[is.na(result) | result == ""] <- "Unknown"
  names(result) <- NULL
  
  return(result)
}

# Standardize phylum names (for bacteria, we shouldn't map to archaea!)
standardize_phylum_names <- function(phylum_names) {
  # Only map bacterial subdivision names to main bacterial phyla
  # Remove the archaea mapping since this is bacteria-only script
  phylum_mapping <- c(
    # Add bacterial-specific mappings if needed
    # "Subdivision_name" = "Main_phylum_name"
  )

  standardized <- phylum_names
  for (old_name in names(phylum_mapping)) {
    standardized[standardized == old_name] <- phylum_mapping[old_name]
  }

  return(standardized)
}

# Function to create isolate percentage legend (from comprehensive script)
create_isolate_legend <- function() {
  legend_data <- data.frame(
    isolate_percentage = c(0, 50, 100),
    circle_size = c(30, 20, 10),
    x = rep(1, 3),
    y = c(2, 1, 0.5),
    stringsAsFactors = FALSE
  )

  legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_point(aes(size = circle_size),
               color = CIRCLE_AESTHETICS$outline_color,
               fill = CIRCLE_AESTHETICS$outline_fill,
               shape = CIRCLE_AESTHETICS$shape,
               stroke = CIRCLE_AESTHETICS$stroke,
               alpha = 0.8) +
    scale_size_identity() +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5)) +
    xlim(0.5, 1.5) +
    ylim(0, 3) +
    labs(title = NULL)

  return(legend_plot)
}

# Create individual bacteria scatter plot (for stacking - no legends, minimal styling)
create_bacteria_scatter <- function(data, level) {

  # Save comprehensive source data with annotation details
  source_filename <- paste0("bacteria_", level, "_source_data.csv")
  source_filepath <- file.path(config$source_data_dir, source_filename)

  data_export <- data
  data_export$Domain <- "Bacteria"
  data_export$Level <- level
  data_export$Plot_Type <- ifelse(data_export$Is_Top_Novelty, "Novel",
                                 ifelse(data_export$Is_Top_Coverage, "Overrepresented", "Background"))

  # Add annotation details for figure reference
  data_export$Will_Be_Annotated <- data_export$Is_Top_Novelty | data_export$Is_Top_Coverage
  data_export$Annotation_Reason <- case_when(
    data_export$Is_Top_Novelty & data_export$Is_Top_Coverage ~ "Both Novel and Overrepresented",
    data_export$Is_Top_Novelty ~ "High Novelty Factor",
    data_export$Is_Top_Coverage ~ "High Coverage Factor",
    TRUE ~ "Background (not annotated)"
  )

  # Add ranking information
  data_export$Novelty_Rank <- rank(-data_export$Novelty_Ratio)
  data_export$Coverage_Rank <- rank(-data_export$Coverage_Factor)
  data_export$NCBI_Count_Rank <- rank(-data_export$NCBI_Species_Count)
  data_export$Census_Count_Rank <- rank(-data_export$Census_OTU_Count)

  write.csv(data_export, source_filepath, row.names = FALSE)
  cat(paste("📊 Source data saved:", source_filepath, "\n"))

  # Save separate annotation-only file for easy reference
  annotation_data <- data_export[data_export$Will_Be_Annotated, ]
  if (nrow(annotation_data) > 0) {
    annotation_filename <- paste0("bacteria_", level, "_annotations_only.csv")
    annotation_filepath <- file.path(config$source_data_dir, annotation_filename)
    write.csv(annotation_data, annotation_filepath, row.names = FALSE)
    cat(paste("📝 Annotation data saved:", annotation_filepath, "\n"))
  }

  # Base plot with expanded limits
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15))) +
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15)))

  # Get top data for highlighting
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  cat(paste("Taxa to be annotated (", nrow(top_data), "):", paste(top_data$Taxon, collapse = ", "), "\n"))

  # Background points first
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = CIRCLE_AESTHETICS$bg_color,
                       fill = CIRCLE_AESTHETICS$bg_fill,
                       shape = CIRCLE_AESTHETICS$shape,
                       alpha = CIRCLE_AESTHETICS$bg_alpha,
                       stroke = CIRCLE_AESTHETICS$stroke)
  }

  if (nrow(top_data) > 0) {
    # Get unique phyla and assign colors (same method as comprehensive script)
    plot_phyla <- sort(unique(top_data$Phylum[!top_data$Phylum %in% c("Unknown", "", "Other", NA)]))
    cat(paste("Plot phyla found:", paste(plot_phyla, collapse = ", "), "\n"))

    # Get master color palette and assign sequentially (like comprehensive script)
    master_colors <- get_master_bacteria_color_palette()

    # Assign colors sequentially to phyla (same as comprehensive script method)
    # Handle case where we have more phyla than colors by cycling through colors
    color_indices <- ((1:length(plot_phyla) - 1) %% length(master_colors)) + 1
    phyla_colors <- master_colors[color_indices]
    names(phyla_colors) <- plot_phyla

    cat("Final color mapping:\n")
    for (i in 1:length(phyla_colors)) {
      cat(paste("  ", names(phyla_colors)[i], "->", phyla_colors[i], "\n"))
    }

    # Add colored points with black outline
    p <- p + geom_point(data = top_data,
                       aes(size = Circle_Size,
                           fill = factor(Phylum, levels = plot_phyla)),
                       color = "black",
                       shape = CIRCLE_AESTHETICS$shape,
                       alpha = CIRCLE_AESTHETICS$alpha,
                       stroke = CIRCLE_AESTHETICS$stroke) +
      scale_fill_manual(values = phyla_colors,
                        name = "Bacterial Phyla",
                        guide = "none")  # Remove legend for stacking
  }

  # Minimal styling for stacking (matching comprehensive script)
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1,
                color = PLOT_THEME$diagonal_color,
                linetype = PLOT_THEME$diagonal_type,
                alpha = PLOT_THEME$diagonal_alpha,
                size = PLOT_THEME$diagonal_size) +
    theme_minimal() +
    theme(
      # Remove titles and axis labels for stacking
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = PLOT_THEME$axis_text_size, face = "bold", color = PLOT_THEME$axis_text_color),
      axis.text.y = element_text(size = PLOT_THEME$axis_text_size, face = "bold", color = PLOT_THEME$axis_text_color),
      panel.grid.major = element_line(color = PLOT_THEME$grid_major_color, linewidth = PLOT_THEME$grid_major_size),
      panel.grid.minor = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = PLOT_THEME$plot_margin,
      # Hide individual plot legends for master layout
      legend.position = "none"
    )

  # Add directional annotations with aggressive ggrepel settings
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_data) > 0) {
    # Get factor values for display
    factor_value <- ifelse(top_data$Is_Top_Novelty,
                          top_data$Novelty_Ratio,
                          top_data$Coverage_Factor)

    # Create labels based on taxonomic level
    if (level == "phylum") {
      # For phylum: only show factor values (no taxon names)
      top_data$label_text <- paste0("(", sprintf("%.1f", factor_value), "×)")
    } else {
      # For family: show taxon names with factor values
      top_data$label_text <- paste0(top_data$Taxon, " (",
                                   sprintf("%.1f", factor_value), "×)")
    }

    # Split data by PRIMARY category for directional annotation (handle overlaps)
    # Priority: If both novelty AND coverage, treat as novelty (down-right)
    novel_data <- top_data[top_data$Is_Top_Novelty, ]  # All novelty taxa (including dual)
    overrep_only_data <- top_data[top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty, ]  # Coverage-only taxa

    # Novel taxa (including dual): DOWN-RIGHT from x=y line
    if (nrow(novel_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = novel_data,
        aes(label = label_text),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        max.overlaps = Inf,
        max.iter = 10000,
        max.time = 5,
        force = 20,
        force_pull = 0.1,
        box.padding = 1.5,
        point.padding = 1.0,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.linewidth = 0.5,
        segment.alpha = 0.8,
        direction = "both",
        nudge_x = 0.3,    # Right
        nudge_y = -0.6    # Down from x=y line
      )
    }

    # Coverage-only taxa: UP-LEFT from x=y line
    if (nrow(overrep_only_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = overrep_only_data,
        aes(label = label_text),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        max.overlaps = Inf,
        max.iter = 10000,
        max.time = 5,
        force = 20,
        force_pull = 0.1,
        box.padding = 1.5,
        point.padding = 1.0,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.linewidth = 0.5,
        segment.alpha = 0.8,
        direction = "both",
        nudge_x = -0.3,   # Left
        nudge_y = 0.6     # Up from x=y line
      )
    }

    # Debug: Verify all top taxa are being annotated
    total_annotated <- nrow(novel_data) + nrow(overrep_only_data)
    cat(paste("Annotation check:", total_annotated, "of", nrow(top_data), "top taxa will be annotated\n"))
    cat(paste("  Novel (down-right):", nrow(novel_data), "taxa\n"))
    cat(paste("  Coverage-only (up-left):", nrow(overrep_only_data), "taxa\n"))
  }

  return(p)
}

# Main function to create bacteria stacked visual
main <- function() {
  cat("Bacteria Stacked Visual Creation\n")
  cat("===============================\n")

  # Set working directory to script location
  script_dir <- dirname(normalizePath(ifelse(interactive(),
                                            file.path(getwd(), "bacteria_specific_scatter_visual.R"),
                                            commandArgs(trailingOnly = FALSE)[4])))
  setwd(script_dir)
  cat(paste("Working directory set to:", getwd(), "\n"))

  # Verify paths exist
  if (!dir.exists(config$data_dir_16s)) {
    stop(paste("16S data directory not found:", config$data_dir_16s))
  }
  if (!dir.exists(config$ncbi_data_dir)) {
    stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
  }

  # Define levels to process
  levels <- c("phylum", "family")

  # Collect all data and create plots
  plots <- list()
  all_data <- list()

  cat("Loading bacterial data...\n")
  for (level in levels) {
    cat(paste("Loading bacterial", level, "data...\n"))
    data <- load_bacteria_data(level)

    if (nrow(data) > 0) {
      all_data[[level]] <- data
      cat(paste("Creating", level, "plot with", nrow(data), "taxa\n"))
      plots[[level]] <- create_bacteria_scatter(data, level)
    } else {
      cat(paste("No data available for", level, "level\n"))
    }
  }

  if (length(plots) == 0) {
    stop("No plots could be created - no data available")
  }

  # Create stacked layout (matching comprehensive script style)
  cat("Arranging plots in stacked layout...\n")

  # Create row labels (left side)
  row_labels <- list()
  level_names <- c("Phyla", "Family")
  for (i in 1:length(levels)) {
    if (levels[i] %in% names(plots)) {
      row_labels[[i]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = level_names[i],
                 size = 12, fontface = "bold", color = "grey30", angle = 90) +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
    }
  }

  # Create header
  header <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Bacteria",
             size = 12, fontface = "bold", color = "grey50") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  # Create top row with header
  top_row <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    header,                   # Bacteria header
    ncol = 2, rel_widths = c(0.1, 1)
  )

  # Create data rows with black frames
  data_rows <- list()
  for (i in 1:length(levels)) {
    level <- levels[i]
    if (level %in% names(plots)) {
      # Add black frame around plot
      framed_plot <- ggdraw(plots[[level]]) +
        theme(plot.background = element_rect(color = "black", fill = NA, size = 3))

      # Combine row label with framed plot
      data_rows[[i]] <- plot_grid(
        row_labels[[i]],
        framed_plot,
        ncol = 2, rel_widths = c(0.1, 1)
      )
    }
  }

  # Create isolate percentage legend
  isolate_legend <- create_isolate_legend()

  # Combine all rows
  main_plot <- plot_grid(
    top_row,                  # Header row
    data_rows[[1]],           # Phylum row
    ggplot() + theme_void(),  # Small spacer
    data_rows[[2]],           # Family row
    ncol = 1, rel_heights = LAYOUT_CONFIG$row_heights
  )

  # Add legend to the right margin
  complete_plot <- plot_grid(
    main_plot,
    isolate_legend,
    ncol = 2, rel_widths = c(LAYOUT_CONFIG$main_plot_ratio, LAYOUT_CONFIG$legend_width_ratio)
  )

  # Save the stacked visual
  output_file <- file.path(config$output_dir, "bacteria_stacked_visual.png")
  cat(paste("Saving bacteria stacked visual to:", output_file, "\n"))

  ggsave(output_file, complete_plot,
         width = config$plot_width, height = config$plot_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  # Also save as PDF
  pdf_file <- file.path(config$output_dir, "bacteria_stacked_visual.pdf")
  ggsave(pdf_file, complete_plot,
         width = config$plot_width, height = config$plot_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  # Create color legend reference
  create_color_legend_reference()

  cat("✅ Bacteria stacked visual creation complete!\n")
  cat(paste("   Main plot:", output_file, "\n"))
  cat(paste("   Dimensions:", config$plot_width, "x", config$plot_height, "inches\n"))
  cat(paste("   Layout: 2 rows (phylum, family) × 1 column (Bacteria only)\n"))
  cat(paste("📊 Source data saved to:", config$source_data_dir, "\n"))
}

# Function to create color legend reference
create_color_legend_reference <- function() {
  # Create a reference file showing the BAC color palette
  bacteria_colors <- get_master_bacteria_color_palette()

  legend_data <- data.frame(
    Color_Index = 1:length(bacteria_colors),
    Color = bacteria_colors,
    Note = "Colors assigned sequentially to phyla as they appear in data",
    stringsAsFactors = FALSE
  )

  legend_file <- file.path(config$source_data_dir, "bacteria_color_palette.csv")
  write.csv(legend_data, legend_file, row.names = FALSE)
  cat(paste("🎨 Color palette saved:", legend_file, "\n"))
}

# Execute main function if script is run directly
if (!interactive()) {
  main()
}
