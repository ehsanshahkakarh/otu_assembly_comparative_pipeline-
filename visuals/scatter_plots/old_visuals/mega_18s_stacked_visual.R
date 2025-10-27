#!/usr/bin/env Rscript
# 18S Mega Stacked Visual - Combined Eukaryote Scatter Plots
# Created: 2025-08-18
# Purpose: Create a unified 2x1 grid of 18S scatter plots with shared color scheme

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

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(), 
                                          file.path(getwd(), "mega_18s_stacked_visual.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Circle aesthetics configuration (matching comprehensive version)
CIRCLE_AESTHETICS <- list(
  shape = 21,           # Filled circles with outline
  stroke = 0.8,         # Outline thickness
  alpha = 0.9,          # Transparency

  # Size scaling
  size_multiplier = 0.8,  # For colored layer (slightly smaller than outline)

  # Background circle settings
  bg_color = "lightgray",
  bg_fill = "lightgray",
  bg_alpha = 0.3,

  # Outline layer settings
  outline_color = "black",
  outline_fill = "black"
)

# Plot theme configuration (matching comprehensive version)
PLOT_THEME <- list(
  # Text sizes
  axis_text_size = 32,
  axis_text_color = "grey50",

  # Grid settings
  grid_major_color = "grey90",
  grid_major_size = 0.3,

  # Diagonal line settings
  diagonal_color = "black",
  diagonal_type = "dashed",
  diagonal_alpha = 0.9,
  diagonal_size = 3,

  # Margins
  plot_margin = margin(1, 1, 1, 1)
)

# Configuration for mega visual (paths relative to script location)
config <- list(
  data_dir_18s = file.path("..", "..", "Eukcensus_merge", "18s_merged", "csv_results"),
  output_dir = file.path("final_visualizations"),
  source_data_dir = file.path("source_data"),
  ncbi_data_dir = file.path("..", "..", "00ncbi_parse", "csv_ncbi"),
  plot_width = 18,   # Updated to match comprehensive proportions
  plot_height = 20,  # Same height as comprehensive
  legend_width = 8,  # Width for separate legend files
  legend_height = 12, # Height for separate legend files
  dpi = 300,
  top_n = 10,
  text_size = 11,    # Clean text size for optimal readability
  size_range = c(10, 22),  # Large circles for better visibility (matching comprehensive)
  color_palette = "professional"
)

# Verify critical paths exist
if (!dir.exists(config$data_dir_18s)) {
  stop(paste("18S data directory not found:", config$data_dir_18s))
}
if (!dir.exists(config$ncbi_data_dir)) {
  stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
}

# Create output directories
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

# Function to save source data and legends (matching comprehensive version)
save_source_data_and_legends <- function(data, level, domain, master_colors) {
  # Save the source data used for plotting
  source_filename <- paste0(domain, "_", level, "_source_data.csv")
  source_filepath <- file.path(config$source_data_dir, source_filename)

  # Add metadata columns for clarity
  data_export <- data
  data_export$Domain <- domain
  data_export$Level <- level
  data_export$Plot_Type <- ifelse(data_export$Is_Top_Novelty, "Novel",
                                 ifelse(data_export$Is_Top_Coverage, "Overrepresented", "Background"))

  write.csv(data_export, source_filepath, row.names = FALSE)
  cat(paste("  📊 Source data saved:", source_filepath, "\n"))

  # Save legend information for Eukaryota divisions
  color_col <- "Division"
  plot_groups <- sort(unique(data$Division[!data$Division %in% c("Unknown", "", "Other", NA)]))
  legend_type <- "Divisions"

  # Create legend data frame
  if (length(plot_groups) > 0) {
    group_colors <- sapply(plot_groups, function(grp) {
      if (grp %in% names(master_colors)) {
        master_colors[grp]
      } else {
        fallback_colors <- c("#E74C3C", "#1D8348", "#117864", "#D35400", "#922B21",
                            "#6C3483", "#1B4F72", "#758C6C", "#85C1E9", "#F1C40F", "#BB8FCE", "#48C9B0")
        fallback_colors[((match(grp, plot_groups) - 1) %% length(fallback_colors)) + 1]
      }
    })

    legend_data <- data.frame(
      Taxon = plot_groups,
      Color = as.character(group_colors),
      Domain = domain,
      Level = level,
      Type = legend_type,
      stringsAsFactors = FALSE
    )

    legend_filename <- paste0(domain, "_", level, "_legend.csv")
    legend_filepath <- file.path(config$source_data_dir, legend_filename)
    write.csv(legend_data, legend_filepath, row.names = FALSE)
    cat(paste("  🎨 Legend data saved:", legend_filepath, "\n"))
  }
}

# Master color palette for eukaryotic divisions
get_master_color_palette <- function() {
  provided_colors <- c("6934b4","fda2da","bf00bf","bee869","ffef44","ff7f00",
                      "4e7abc","859d9a","0000ff","6f2c2e","fb9a9a","ff0002")

  plotted_divisions <- c("Opisthokonta", "Alveolata", "Rhizaria", "Discoba",
                        "Stramenopiles", "Evosea", "Streptophyta", "Chlorophyta",
                        "Tubulinea", "Metamonada", "Discosea", "Rhodophyta")

  # Create color mapping
  eukaryota_color_map <- setNames(paste0("#", provided_colors[1:length(plotted_divisions)]),
                                 plotted_divisions)
  return(eukaryota_color_map)
}

# Allowed eukaryotic divisions constant
ALLOWED_DIVISIONS <- c("Alveolata", "Chlorophyta", "Discoba", "Discosea", "Evosea",
                      "Metamonada", "Opisthokonta", "Rhizaria", "Rhodophyta",
                      "Stramenopiles", "Streptophyta", "Tubulinea")

# Helper function to filter divisions
filter_divisions <- function(data) {
  data_before <- nrow(data)
  data <- data %>% filter(!Division %in% c("Unknown", "", "Other", NA))
  cat(paste("Filtered from", data_before, "to", nrow(data), "taxa after removing unmapped divisions\n"))
  return(data)
}

# Function to extract and save source data for plotted taxa
extract_18s_plotted_source_data <- function(data, level) {
  # Apply same filtering logic as plotting functions
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  meaningful_data <- top_data  # 18S shows ALL top taxa without additional filtering

  if (nrow(meaningful_data) == 0) {
    return(NULL)
  }

  # All levels use the "Taxon" column for the taxonomic name
  level_col_name <- "Taxon"

  # Create source data with all requested columns
  source_data <- data.frame(
    Taxonomic_Level = level,
    Domain = "Eukaryota",
    Taxon = meaningful_data[[level_col_name]],  # Use proper column name
    Division = meaningful_data$Division,  # 18S uses Division instead of Phylum
    Novelty_Factor = round(meaningful_data$Novelty_Ratio, 3),
    Overrepresentation_Factor = round(meaningful_data$Coverage_Factor, 3),
    Census_OTU_Count = meaningful_data$Census_OTU_Count,
    NCBI_Species_Count = meaningful_data$NCBI_Species_Count,
    NCBI_Genome_Count = meaningful_data$NCBI_Genome_Count,
    Isolate_Count = meaningful_data$Isolate_Count,
    Is_Top_Novelty = meaningful_data$Is_Top_Novelty,
    Is_Top_Coverage = meaningful_data$Is_Top_Coverage,
    stringsAsFactors = FALSE
  )

  return(source_data)
}

# Load data function
load_18s_data <- function(level) {
  # Determine file name based on level
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

  # Add division information (exactly like individual 18S scripts)
  if (level == "phylum") {
    # For phylum level, each taxon IS a division (like individual scripts do)
    data$Division <- data$Taxon
    cat(paste("Phylum level: assigned", nrow(data), "taxa as divisions\n"))
  } else {
    # For family level, assign divisions based on taxonomic knowledge
    # Start with "Other" and then apply manual overrides
    data$Division <- "Other"
    cat(paste("Family level: starting with", nrow(data), "families, will apply manual overrides\n"))
  }

  # Apply manual overrides for specific taxa
  manual_overrides <- get_manual_overrides()
  if (!is.null(manual_overrides)) {
    for (taxon in names(manual_overrides)) {
      if (taxon %in% data$Taxon) {
        data$Division[data$Taxon == taxon] <- manual_overrides[[taxon]]
      }
    }
    cat(paste("Applied", sum(data$Taxon %in% names(manual_overrides)), "manual overrides\n"))
  }

  # Filter divisions and check for remaining data
  data <- filter_divisions(data)
  if (nrow(data) == 0) {
    cat("WARNING: No data remaining after division filtering!\n")
    return(data.frame())
  }

  cat(paste("Final divisions:", paste(head(unique(data$Division), 8), collapse = ", "), "\n"))
  
  # Use pre-calculated factors from merger scripts
  data$Novelty_Ratio <- data$novelty_factor
  data$Coverage_Factor <- data$overrepresentation_factor

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

  # Circle sizing - matching 16S script approach exactly
  data$Genome_Isolate_Ratio <- pmax(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), 1)

  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  size_range <- range(data$Circle_Size_Raw, na.rm = TRUE)
  # Ensure minimum circle size of 3 and maximum of 12 for better visibility (matching 16S)
  data$Circle_Size <- 3 + 9 * (data$Circle_Size_Raw - size_range[1]) / diff(size_range)

  # Identify top taxa - use original factors for thresholding, not weighted
  high_novelty_data <- data[data$Novelty_Ratio > 1.0, ]
  high_coverage_data <- data[data$Coverage_Factor > 1.0, ]

  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE

  # Mark top novelty - use weighted for ranking but original threshold
  if (nrow(high_novelty_data) > 0) {
    high_novelty_indices <- which(data$Taxon %in% high_novelty_data$Taxon)
    high_novelty_values <- ranking_novelty[high_novelty_indices]
    novelty_ranks <- rank(-high_novelty_values)
    top_n_count <- min(config$top_n, length(high_novelty_values))
    top_novelty_mask <- novelty_ranks <= top_n_count
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_mask]
    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
  }

  # Mark top coverage - use weighted for ranking but original threshold
  if (nrow(high_coverage_data) > 0) {
    high_coverage_indices <- which(data$Taxon %in% high_coverage_data$Taxon)
    high_coverage_values <- ranking_coverage[high_coverage_indices]
    coverage_ranks <- rank(-high_coverage_values)
    top_n_count <- min(config$top_n, length(high_coverage_values))
    top_coverage_mask <- coverage_ranks <= top_n_count
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_mask]
    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
  }

  # Summary output
  cat(paste("Final data:", nrow(data), "taxa,",
            sum(data$Is_Top_Novelty), "top novelty,",
            sum(data$Is_Top_Coverage), "top coverage\n"))

  return(data)
}



# Manual overrides for specific taxa
get_manual_overrides <- function() {
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

# Create individual scatter plot (without titles and axis labels for stacking)
create_individual_18s_scatter <- function(data, level, all_phyla, master_colors) {

  # Save source data and legend information
  save_source_data_and_legends(data, level, "Eukaryota", master_colors)

  # ggrepel parameters for consistent annotation (matching comprehensive version)
  get_repel_params <- function(color = "black", nudge_x = 0, nudge_y = 0) {
    list(
      color = color,
      size = config$text_size,
      fontface = "bold",
      max.overlaps = Inf,
      max.iter = 8000,           # Higher iteration count like comprehensive
      force = 15,                # Higher force like comprehensive
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

  # Base plot with expanded limits to accommodate labels
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15))) +  # Increased expansion for label space
    scale_y_log10(labels = comma_format(),
                  limits = c(1, 10000),
                  expand = expansion(mult = c(0.15, 0.15)))

  # Get top data for highlighting
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]

  # Show ALL top taxa (both novelty and coverage) - no additional filtering
  meaningful_data <- top_data

  # Background points first (draw behind meaningful data)
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = CIRCLE_AESTHETICS$bg_color,
                       fill = CIRCLE_AESTHETICS$bg_fill,
                       shape = CIRCLE_AESTHETICS$shape,
                       alpha = CIRCLE_AESTHETICS$bg_alpha,
                       stroke = CIRCLE_AESTHETICS$stroke)
  }

  if (nrow(meaningful_data) > 0) {
    # Get divisions present in this plot
    plot_divisions <- sort(unique(meaningful_data$Division[!meaningful_data$Division %in% c("Unknown", "", "Other", NA)]))

    # Map to master color scheme with fallback colors
    fallback_colors <- c("#E74C3C", "#1D8348", "#117864", "#D35400", "#922B21",
                        "#6C3483", "#1B4F72", "#758C6C", "#85C1E9", "#F1C40F", "#BB8FCE", "#48C9B0")

    division_colors <- sapply(plot_divisions, function(div) {
      if (div %in% names(master_colors)) {
        master_colors[div]
      } else {
        fallback_colors[((match(div, plot_divisions) - 1) %% length(fallback_colors)) + 1]
      }
    })

    # Add single layer with colored fill and thin black outline (matching comprehensive)
    p <- p + geom_point(data = meaningful_data,
                       aes(size = Circle_Size, fill = factor(Division, levels = plot_divisions)),
                       color = "black",  # Thin black outline
                       shape = CIRCLE_AESTHETICS$shape,
                       alpha = CIRCLE_AESTHETICS$alpha,
                       stroke = CIRCLE_AESTHETICS$stroke) +
      scale_fill_manual(values = setNames(division_colors, plot_divisions),
                       name = "Divisions", guide = "none")  # Remove individual legends
  }



  # Minimal styling for stacking using configuration
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

  # Add text annotations for top taxa
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(meaningful_data) > 0) {
    # Get factor values for display
    factor_value <- ifelse(meaningful_data$Is_Top_Novelty,
                          meaningful_data$Novelty_Ratio,
                          meaningful_data$Coverage_Factor)

    # Create labels - show full annotations for family level, only factors for phylum level (matching 16S)
    if (level == "family") {
      # Add asterisks for non-family taxa in family-level plots
      family_suffixes <- c("aceae", "idae", "ales", "ineae", "inae", "eae")
      meaningful_data$taxon_display <- ifelse(sapply(meaningful_data$Taxon, function(taxon) {
        any(sapply(family_suffixes, function(suffix) grepl(paste0(suffix, "$"), taxon, ignore.case = TRUE)))
      }), meaningful_data$Taxon, paste0(meaningful_data$Taxon, "*"))

      meaningful_data$full_label <- paste0(meaningful_data$taxon_display, " (",
                                          sprintf("%.1f", factor_value), "×)")
    } else {
      # For phylum level, show only factor values (matching 16S approach)
      meaningful_data$full_label <- paste0("(", sprintf("%.1f", factor_value), "×)")
    }

    # DIRECTIONAL ANNOTATION STRATEGY (matching comprehensive version):
    # Novel taxa -> DOWN and RIGHT (away from center diagonal)
    # Overrepresented taxa -> UP and LEFT (away from center diagonal)
    # This creates clear visual separation and prevents overlaps

    # Split data by category type
    novel_data <- meaningful_data[meaningful_data$Is_Top_Novelty, ]
    overrep_data <- meaningful_data[meaningful_data$Is_Top_Coverage, ]

    # NOVEL TAXA: Repel DOWN and to the RIGHT (away from center line)
    if (nrow(novel_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = novel_data,
        aes(label = full_label),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        max.overlaps = Inf,        # Never drop annotations
        max.iter = 8000,           # High iteration count for good positioning
        force = 15,                # Strong repulsion force
        force_pull = 2,            # Additional pulling force
        box.padding = 0.4,         # Space around text boxes
        point.padding = 0.3,       # Space around points
        min.segment.length = 0,    # Always show connector lines
        segment.color = "grey50",
        segment.linewidth = 0.4,
        segment.alpha = 0.8,
        direction = "both",        # Allow movement in both directions
        nudge_x = 0.3,            # Push RIGHT (away from center line)
        nudge_y = -0.6,           # Push DOWN
        expand = expansion(mult = 0.15)  # Extra plot space
      )
    }

    # OVERREPRESENTED TAXA: Repel UP and to the LEFT (away from center line)
    if (nrow(overrep_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = overrep_data,
        aes(label = full_label),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        max.overlaps = Inf,        # Never drop annotations
        max.iter = 8000,           # High iteration count for good positioning
        force = 15,                # Strong repulsion force
        force_pull = 2,            # Additional pulling force
        box.padding = 0.4,         # Space around text boxes
        point.padding = 0.3,       # Space around points
        min.segment.length = 0,    # Always show connector lines
        segment.color = "grey50",
        segment.linewidth = 0.4,
        segment.alpha = 0.8,
        direction = "both",        # Allow movement in both directions
        nudge_x = -0.3,           # Push LEFT (away from center line)
        nudge_y = 0.6,            # Push UP
        expand = expansion(mult = 0.15)  # Extra plot space
      )
    }
  }



  return(p)
}

# Create master legend for all eukaryotic divisions
create_master_18s_legend <- function(all_divisions, master_colors) {
  # Generate colors with fallback
  fallback_colors <- c("#E74C3C", "#1D8348", "#117864", "#D35400", "#922B21",
                      "#6C3483", "#1B4F72", "#758C6C", "#85C1E9", "#F1C40F", "#BB8FCE", "#48C9B0")

  legend_colors <- sapply(all_divisions, function(div) {
    if (div %in% names(master_colors)) {
      master_colors[div]
    } else {
      fallback_colors[((match(div, all_divisions) - 1) %% length(fallback_colors)) + 1]
    }
  })

  # Create legend plot
  dummy_data <- data.frame(
    x = rep(1, length(all_divisions)),
    y = rep(1, length(all_divisions)),
    Division = factor(all_divisions, levels = all_divisions)
  )

  legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = Division)) +
    geom_point(size = 8, alpha = 0.9) +
    scale_color_manual(values = setNames(legend_colors, all_divisions),
                      name = "Divisions",
                      guide = guide_legend(
                        override.aes = list(size = 8, alpha = 1, shape = 15),
                        keywidth = unit(2.0, "cm"),
                        keyheight = unit(1.5, "cm"),
                        ncol = 1
                      )) +
    theme_void() +
    theme(
      legend.text = element_text(size = 34, face = "bold"),
      legend.title = element_text(size = 38, face = "bold", color = "grey50"),
      legend.key.size = unit(3.2, "cm"),
      legend.spacing.y = unit(1.5, "cm"),
      plot.margin = margin(0, 0, 0, 0)
    )

  return(get_legend(legend_plot))
}

# Create size legend for circle sizes
create_size_legend <- function(sample_data) {
  # Create dummy data with representative circle sizes
  size_values <- c(min(sample_data$Circle_Size, na.rm = TRUE),
                   median(sample_data$Circle_Size, na.rm = TRUE),
                   max(sample_data$Circle_Size, na.rm = TRUE))

  # Calculate corresponding genome/isolate ratios for labels
  min_ratio <- min(sample_data$Genome_Isolate_Ratio, na.rm = TRUE)
  max_ratio <- max(sample_data$Genome_Isolate_Ratio, na.rm = TRUE)
  median_ratio <- median(sample_data$Genome_Isolate_Ratio, na.rm = TRUE)

  dummy_data <- data.frame(
    x = rep(1, 3),
    y = rep(1, 3),
    Circle_Size = size_values,
    Ratio_Label = c(sprintf("%.1f", min_ratio),
                    sprintf("%.1f", median_ratio),
                    sprintf("%.1f", max_ratio))
  )

  size_legend_plot <- ggplot(dummy_data, aes(x = x, y = y, size = Circle_Size)) +
    geom_point(color = "black", alpha = 0.8) +
    scale_size_continuous(
      range = config$size_range,
      name = "Genome/Isolate\nRatio",
      breaks = size_values,
      labels = dummy_data$Ratio_Label,
      guide = guide_legend(
        override.aes = list(color = "black", alpha = 0.8),
        title.position = "top",
        title.hjust = 0.5,
        keywidth = unit(1.5, "cm"),
        keyheight = unit(1.5, "cm"),
        ncol = 1
      )
    ) +
    theme_void() +
    theme(
      legend.title = element_text(size = 14, face = "bold", color = "grey50"),
      legend.text = element_text(size = 12, color = "grey50"),
      legend.key = element_blank(),
      legend.spacing.y = unit(0.5, "cm"),
      plot.margin = margin(0, 0, 0, 0)
    )

  # Extract just the legend
  size_legend <- get_legend(size_legend_plot)
  return(size_legend)
}

# Main function to create 18S mega visual
main <- function() {
  cat("18S Mega Stacked Visual Creation\n")
  cat("===============================\n")

  # Define the grid structure (only eukaryotes, 2 levels: division = phylum level)
  levels <- c("phylum", "family")  # phylum data = divisions for 18S

  # Collect all data to determine master phyla list
  all_data <- list()
  all_phyla_set <- character(0)

  cat("Loading all data to determine master divisions list...\n")
  for (level in levels) {
    data <- load_18s_data(level)
    if (nrow(data) == 0) {
      cat(paste("WARNING: No data loaded for level", level, "- skipping\n"))
      next
    }

    all_data[[level]] <- data

    # Collect divisions from top taxa
    top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
    meaningful_divisions <- unique(top_data$Division[!top_data$Division %in% c("Unknown", "", "Other", NA)])
    all_phyla_set <- unique(c(all_phyla_set, meaningful_divisions))

    cat(paste("Level", level, "- found", length(meaningful_divisions), "meaningful divisions\n"))
  }

  # Sort divisions alphabetically for consistent coloring
  all_divisions <- sort(all_phyla_set)  # Keep variable name for compatibility
  master_colors <- get_master_color_palette()

  cat(paste("Found", length(all_divisions), "unique divisions across all datasets\n"))
  cat(paste("Divisions:", paste(all_divisions, collapse = ", "), "\n"))

  # Extract and save source data organized by rank and factor type
  cat("Extracting source data for plotted taxa...\n")
  if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

  # Process phylum and family levels only
  for (level in c("phylum", "family")) {
    if (level %in% names(all_data)) {
      source_data <- extract_18s_plotted_source_data(all_data[[level]], level)
      if (!is.null(source_data)) {
        cat(paste("  18S", level, ":", nrow(source_data), "plotted taxa\n"))

        # Split into novelty and overrepresented (coverage)
        novelty_data <- source_data[source_data$Is_Top_Novelty == TRUE, ]
        overrepresented_data <- source_data[source_data$Is_Top_Coverage == TRUE, ]

        # Save separate files
        if (nrow(novelty_data) > 0) {
          novelty_file <- file.path(config$source_data_dir, paste0("18s_", level, "_novelty_source_data.csv"))
          write.csv(novelty_data, novelty_file, row.names = FALSE)
          cat(paste("  Novelty data saved:", novelty_file, "(", nrow(novelty_data), "taxa)\n"))
        }

        if (nrow(overrepresented_data) > 0) {
          overrepresented_file <- file.path(config$source_data_dir, paste0("18s_", level, "_overrepresented_source_data.csv"))
          write.csv(overrepresented_data, overrepresented_file, row.names = FALSE)
          cat(paste("  Overrepresented data saved:", overrepresented_file, "(", nrow(overrepresented_data), "taxa)\n"))
        }
      }
    }
  }

  # Create individual plots
  plots <- list()
  for (level in levels) {
    if (!level %in% names(all_data) || nrow(all_data[[level]]) == 0) {
      cat(paste("Skipping plot creation for", level, "- no data available\n"))
      next
    }

    cat(paste("Creating plot for", level, "level...\n"))
    plots[[level]] <- create_individual_18s_scatter(all_data[[level]], level, all_divisions, master_colors)
  }

  # Create and save master legend separately
  cat("Creating and saving master legend...\n")
  master_legend <- create_master_18s_legend(all_divisions, master_colors)

  # Save legend as separate PDF and PNG files
  legend_pdf_file <- file.path(config$output_dir, "18s_divisions_legend.pdf")
  legend_png_file <- file.path(config$output_dir, "18s_divisions_legend.png")

  cat(paste("Saving legend to:", legend_pdf_file, "\n"))
  ggsave(legend_pdf_file, master_legend,
         width = config$legend_width, height = config$legend_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  cat(paste("Saving legend to:", legend_png_file, "\n"))
  ggsave(legend_png_file, master_legend,
         width = config$legend_width, height = config$legend_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  # Arrange plots in 2x1 grid with column header (matching 16S structure)
  cat("Arranging plots in grid with Eukaryota header...\n")

  # Create row labels (left side)
  row_labels <- list()
  level_names <- c("Divisions", "Family")  # Use "Divisions" for eukaryotic phyla
  for (i in 1:length(levels)) {
    row_labels[[i]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = level_names[i],
               size = 12, fontface = "bold", color = "grey30", angle = 90) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  # Create column header (matching 16S approach)
  col_header <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Eukaryota",
             size = 12, fontface = "bold", color = "grey50") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  # Create top row with column header
  top_row <- plot_grid(
    ggplot() + theme_void(),  # Empty corner
    col_header,               # Eukaryota header
    ncol = 2, rel_widths = c(0.1, 1)
  )

  # Create data rows with individual black frames around each plot
  data_rows <- list()
  for (i in 1:length(levels)) {
    level <- levels[i]

    # Add black frame around individual plot
    framed_plot <- ggdraw(plots[[level]]) +
      theme(plot.background = element_rect(color = "black", fill = NA, size = 3))

    # Combine row label with framed plot
    data_rows[[i]] <- plot_grid(
      row_labels[[i]],
      framed_plot,
      ncol = 2, rel_widths = c(0.1, 1)
    )
  }

  # Combine all rows with header and spacing (matching 16S structure)
  complete_plot <- plot_grid(
    top_row,                  # Column header row
    data_rows[[1]],           # Division row (phylum level)
    ggplot() + theme_void(),  # Small spacer
    data_rows[[2]],           # Family row
    ncol = 1, rel_heights = c(0.1, 1, 0.05, 1)
  )

  # Save the mega visual
  output_file <- file.path(config$output_dir, "18s_mega_stacked_visual.png")
  cat(paste("Saving mega visual to:", output_file, "\n"))

  ggsave(output_file, complete_plot,
         width = config$plot_width, height = config$plot_height,
         dpi = config$dpi, bg = "white", limitsize = FALSE,
         units = "in")

  # Create master source data index
  create_source_data_index()

  cat("✅ 18S Mega visual creation complete!\n")
  cat(paste("   Main plot:", output_file, "\n"))
  cat(paste("   Legend files:", legend_pdf_file, "and", legend_png_file, "\n"))
  cat(paste("   Dimensions:", config$plot_width, "x", config$plot_height, "inches\n"))
  cat(paste("   Total divisions:", length(all_divisions), "\n"))
  cat(paste("📊 Source data saved to:", config$source_data_dir, "\n"))
}

# Function to create master index of all source data files (matching comprehensive version)
create_source_data_index <- function() {
  # Get all CSV files in source_data directory
  csv_files <- list.files(config$source_data_dir, pattern = "\\.csv$", full.names = FALSE)

  if (length(csv_files) > 0) {
    # Separate source data files from legend files
    source_files <- csv_files[grepl("source_data", csv_files)]
    legend_files <- csv_files[grepl("legend", csv_files)]

    # Create index data frame
    index_data <- data.frame(
      File_Type = c(rep("Source_Data", length(source_files)), rep("Legend", length(legend_files))),
      Filename = c(source_files, legend_files),
      Description = c(
        paste("Plot data for", gsub("_source_data.csv", "", source_files)),
        paste("Color legend for", gsub("_legend.csv", "", legend_files))
      ),
      Created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )

    # Save index file
    index_filepath <- file.path(config$source_data_dir, "README_18s_source_data_index.csv")
    write.csv(index_data, index_filepath, row.names = FALSE)
    cat(paste("📋 18S source data index created:", index_filepath, "\n"))
  }
}

# Run the mega visual creation
if (!interactive()) main()