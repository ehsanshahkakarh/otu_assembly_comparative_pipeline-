#!/usr/bin/env Rscript
# Comprehensive Mega Stacked Visual - 16S Bacteria + 16S Archaea + 18S Eukaryota
# Created: 2025-10-26
# Purpose: Create a unified 3-column mega-grid with all three domains side-by-side

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(grid)
  library(cowplot)
  library(ggrepel)
})

# Configuration constants
PLOT_CONFIG <- list(
  # Plot dimensions and layout
  plot_width = 54,
  plot_height = 20,
  dpi = 300,

  # Data filtering
  top_n = 10,

  # Visual styling
  text_size = 11,
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

# Function to create simple legend
create_legend <- function() {
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

# Simple path configuration
config <- list(
  data_dir_16s = "../../Eukcensus_merge/16s_merged/csv_results",
  data_dir_18s = "../../Eukcensus_merge/18s_merged/csv_results",
  output_dir = "final_visualizations",
  source_data_dir = "source_data",
  ncbi_data_dir = "../../00ncbi_parse/csv_ncbi"
)

# Create output directories
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

# Function to save source data
save_source_data <- function(data, level, domain) {
  source_filename <- paste0(domain, "_", level, "_source_data.csv")
  source_filepath <- file.path(config$source_data_dir, source_filename)

  data_export <- data
  data_export$Domain <- domain
  data_export$Level <- level

  write.csv(data_export, source_filepath, row.names = FALSE)
  cat(paste("Source data saved:", source_filepath, "\n"))
}

# Color palettes
get_color_palettes <- function() {
  bacteria_colors <- c(
    "#4c9b34", "#72c859", "#9fe18b", "#cfe99d",   # spaced greens
    "#ff7200", "#d58a2f", "#ffb44c", "#f5d24f", "#a67c28",  # warm band
    "#548877", "#46bda3", "#55e3ff", "#7ac7da",   # teals/aquas
    "#bfb1d3", "#80456e",                          # moved from euks (lavender, plum)
    "#994417"                                      # brown
  )

  archaea_colors <- c("#f51b7f", "#ff3f4d", "#d19386", "#8c2a50", "#f5c7bd")

  eukaryota_colors <- c("#416b7d", "#69c1d4", "#55d0ba", "#003ce1", "#c73de4",
                       "#65417a", "#68536c", "#cf8ac6", "#d24390", "#475093",
                       "#663be6", "#7a9dcd", "#2E8B57")

  eukaryota_divisions <- c("Opisthokonta", "Alveolata", "Rhizaria", "Discoba",
                          "Stramenopiles", "Evosea", "Streptophyta", "Chlorophyta",
                          "Metamonada", "Discosea", "Rhodophyta", "Tubulinea")

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

  if (nrow(data) == 0) return(data.frame())

  # Standardize column names
  colnames(data)[1] <- "Taxon"
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  data$Novelty_Ratio <- data$novelty_factor
  data$Coverage_Factor <- data$overrepresentation_factor
  data$Circle_Size <- calculate_circle_size(data$Isolate_Count, data$NCBI_Genome_Count)

  # Add phylum information
  data <- add_phylum_info(data, level, domain)

  # Identify top taxa (threshold > 1.0)
  threshold <- 1.0
  data$Is_Top_Novelty <- data$Novelty_Ratio > threshold &
                        data$Taxon %in% head(data[order(-data$Novelty_Ratio), ]$Taxon, PLOT_CONFIG$top_n)
  data$Is_Top_Coverage <- data$Coverage_Factor > threshold &
                         data$Taxon %in% head(data[order(-data$Coverage_Factor), ]$Taxon, PLOT_CONFIG$top_n)

  return(data)
}

# Add phylum information
add_phylum_info <- function(data, level, domain) {
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

  if (nrow(data) == 0) return(data.frame())

  # Standardize column names
  colnames(data)[1] <- "Taxon"
  data$Census_OTU_Count <- data$census_otu_count
  data$NCBI_Species_Count <- data$ncbi_species_count
  data$NCBI_Genome_Count <- data$ncbi_genome_count
  data$Isolate_Count <- data$isolate_count
  data$Novelty_Ratio <- data$novelty_factor
  data$Coverage_Factor <- data$overrepresentation_factor
  data$Circle_Size <- calculate_circle_size(data$Isolate_Count, data$NCBI_Genome_Count)

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
  if (nrow(data) == 0) return(data.frame())

  # Identify top taxa (threshold > 1.0)
  threshold <- 1.0
  data$Is_Top_Novelty <- data$Novelty_Ratio > threshold &
                        data$Taxon %in% head(data[order(-data$Novelty_Ratio), ]$Taxon, PLOT_CONFIG$top_n)
  data$Is_Top_Coverage <- data$Coverage_Factor > threshold &
                         data$Taxon %in% head(data[order(-data$Coverage_Factor), ]$Taxon, PLOT_CONFIG$top_n)

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
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(), limits = c(1, 10000)) +
    scale_y_log10(labels = comma_format(), limits = c(1, 10000))

  # Background points
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "lightgray", fill = "lightgray",
                       shape = PLOT_CONFIG$circle_shape, alpha = PLOT_CONFIG$bg_alpha,
                       stroke = PLOT_CONFIG$circle_stroke)
  }

  # Top data points with proper color mapping
  if (nrow(top_data) > 0) {
    # Determine color column and get appropriate colors
    if (domain == "Eukaryota") {
      color_col <- "Division"
      plot_groups <- sort(unique(top_data$Division[!top_data$Division %in% c("Unknown", "", "Other", NA)]))
      # Use eukaryota colors from master palette
      group_colors <- master_colors$eukaryota[plot_groups]
      # Fill in missing colors with fallbacks
      missing_colors <- is.na(group_colors)
      if (any(missing_colors)) {
        fallback_colors <- c("#E74C3C", "#1D8348", "#117864", "#D35400", "#922B21", "#6C3483")
        group_colors[missing_colors] <- fallback_colors[1:sum(missing_colors)]
      }
      names(group_colors) <- plot_groups
    } else {
      color_col <- "Phylum"
      plot_groups <- sort(unique(top_data$Phylum[!top_data$Phylum %in% c("Unknown", "", "Other", NA)]))
      # Use bacteria or archaea colors from master palette
      if (domain == "Bacteria") {
        group_colors <- master_colors$bacteria[1:length(plot_groups)]
      } else {
        group_colors <- master_colors$archaea[1:length(plot_groups)]
      }
      names(group_colors) <- plot_groups
    }

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

  # Simple annotations for top taxa
  if (nrow(top_data) > 0) {
    # Create simple labels
    factor_value <- pmax(top_data$Novelty_Ratio, top_data$Coverage_Factor)
    top_data$label <- paste0(top_data$Taxon, " (", sprintf("%.1f", factor_value), "×)")

    p <- p + ggrepel::geom_text_repel(
      data = top_data,
      aes(label = label),
      size = PLOT_CONFIG$text_size,
      fontface = "bold",
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5
    )
  }


  return(p)
}




# Main function to create comprehensive mega visual
main <- function() {
  cat("Comprehensive Mega Stacked Visual Creation\n")
  cat("==========================================\n")

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

  # Create isolate percentage legend
  isolate_legend <- create_legend()

  # Combine all rows with header and spacing
  main_plot <- plot_grid(
    top_row,                  # Column header row
    data_rows[[1]],           # Phylum row
    ggplot() + theme_void(),  # Small spacer
    data_rows[[2]],           # Family row
    ncol = 1, rel_heights = c(0.1, 1, 0.05, 1)
  )

  # Add legend to the right margin
  complete_plot <- plot_grid(
    main_plot,
    isolate_legend,
    ncol = 2, rel_widths = c(PLOT_CONFIG$main_plot_ratio, PLOT_CONFIG$legend_width_ratio)
  )

  # Save the comprehensive mega visual
  output_file_png <- file.path(config$output_dir, "comprehensive_mega_stacked_visual.png")
  output_file_pdf <- file.path(config$output_dir, "comprehensive_mega_stacked_visual.pdf")

  cat(paste("Saving comprehensive mega visual to:", output_file_png, "\n"))
  cat(paste("Saving PDF version to:", output_file_pdf, "\n"))

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

  # Create master source data index
  create_source_data_index()

  cat("✅ Comprehensive mega visual creation complete!\n")
  cat(paste("   Main plot:", output_file, "\n"))
  cat(paste("   Dimensions:", PLOT_CONFIG$plot_width, "x", PLOT_CONFIG$plot_height, "inches\n"))
  cat(paste("   Layout: 2 rows (phylum, family) × 3 columns (Bacteria, Archaea, Eukaryota)\n"))
  cat(paste("📊 Source data saved to:", config$source_data_dir, "\n"))
}

# Function to create master index of all source data files
create_source_data_index <- function() {
  # Get all CSV files in source_data directory
  csv_files <- list.files(config$source_data_dir, pattern = "\\.csv$", full.names = FALSE)

  if (length(csv_files) > 0) {
    # Separate source data files from legend files
    source_files <- csv_files[grepl("source_data", csv_files)]
    legend_files <- csv_files[grepl("legend", csv_files)]

    # Create index data frame
    total_files <- length(source_files) + length(legend_files)
    index_data <- data.frame(
      File_Type = c(rep("Source_Data", length(source_files)), rep("Legend", length(legend_files))),
      Filename = c(source_files, legend_files),
      Description = c(
        paste("Plot data for", gsub("_source_data.csv", "", source_files)),
        paste("Color legend for", gsub("_legend.csv", "", legend_files))
      ),
      Created = rep(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), total_files),
      stringsAsFactors = FALSE
    )

    # Save index file
    index_filepath <- file.path(config$source_data_dir, "README_source_data_index.csv")
    write.csv(index_data, index_filepath, row.names = FALSE)
    cat(paste("📋 Source data index created:", index_filepath, "\n"))
  }
}



# Run the main visualization
if (!interactive()) {
  main()
}
