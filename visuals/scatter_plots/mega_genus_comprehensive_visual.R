#!/usr/bin/env Rscript
# Comprehensive Mega Genus Visual - 16S Bacteria + 16S Archaea + 18S Eukaryota
# Created: 2025-10-17
# Purpose: Create a unified 3-column genus-level visualization with all three domains side-by-side

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
                                          file.path(getwd(), "mega_genus_comprehensive_visual.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Configuration for comprehensive mega genus visual
config <- list(
  data_dir_16s = file.path("..", "Eukcensus_merge", "16s_merged", "csv_results"),
  data_dir_18s = file.path("..", "Eukcensus_merge", "18s_merged", "csv_results"),
  output_dir = file.path("final_visualizations"),
  source_data_dir = file.path("source_data"),
  ncbi_data_dir = file.path("..", "00ncbi_parse", "csv_ncbi"),
  plot_width = 48,   # Wide for 3 columns - matches stacked visual for square plots
  plot_height = 20,  # Taller height for more square individual plots
  legend_width = 10,  # Width for separate legend files
  legend_height = 15, # Height for separate legend files
  dpi = 300,
  top_n = 15,        # Top 15 instead of top 10
  text_size = 11,    # Clean text size for optimal readability
  size_range = c(10, 22),  # Large circles for better visibility
  color_palette = "professional"
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

# Create output directories
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$source_data_dir, recursive = TRUE, showWarnings = FALSE)

# Master color palette for bacteria (from 16S script)
get_master_16s_color_palette <- function() {
  bacteria_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                      "#c49c94", "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5")
  
  archaea_colors <- c("#C0342B", "#5247A6", "#1D8890", "#B78933", "#8B4513", 
                     "#2F4F4F", "#800080", "#008B8B", "#B22222", "#556B2F",
                     "#A0522D", "#BC8F8F", "#F5DEB3", "#FFE4B5", "#FFDAB9")
  
  return(list(bacteria = bacteria_colors, archaea = archaea_colors))
}

# Master color palette for eukaryotic divisions (from 18S script)
get_master_18s_color_palette <- function() {
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

# Helper function to filter divisions (for 18S)
filter_divisions <- function(data) {
  data_before <- nrow(data)

  # Debug: Check coverage factors before filtering
  cat(paste("DEBUG: Before division filtering - Max overrepresentation_factor:", max(data$overrepresentation_factor, na.rm = TRUE), "\n"))
  cat(paste("DEBUG: Before division filtering - Taxa with overrepresentation_factor > 1.0:", sum(data$overrepresentation_factor > 1.0, na.rm = TRUE), "\n"))
  cat(paste("DEBUG: Sample high coverage taxa before filtering:", paste(head(data$Taxon[data$overrepresentation_factor > 1.0], 3), collapse = ", "), "\n"))

  data <- data %>% filter(!Division %in% c("Unknown", "", "Other", NA))
  cat(paste("Filtered from", data_before, "to", nrow(data), "taxa after removing unmapped divisions\n"))

  # Debug: Check coverage factors after filtering
  if (nrow(data) > 0) {
    cat(paste("DEBUG: After division filtering - Max overrepresentation_factor:", max(data$overrepresentation_factor, na.rm = TRUE), "\n"))
    cat(paste("DEBUG: After division filtering - Taxa with overrepresentation_factor > 1.0:", sum(data$overrepresentation_factor > 1.0, na.rm = TRUE), "\n"))
  }

  return(data)
}

# Load 16S data function (16S files contain both domains, filter by domain column)
load_16s_data <- function(level, domain) {
  # 16S files are not separated by domain - load the combined file and filter
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)

  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }

  cat(paste("Loading", domain, level, "data from:", filepath, "\n"))
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  cat(paste("Loaded", nrow(data), "total rows\n"))

  # Filter by domain (Bacteria or Archaea) - add debug output
  cat(paste("Before filtering:", nrow(data), "total taxa\n"))
  target_domain <- domain  # Store the target domain to avoid variable name conflict
  data <- data %>% filter(domain == target_domain)
  cat(paste("After filtering for", target_domain, ":", nrow(data), "taxa\n"))

  # Show sample taxa for verification
  if (nrow(data) > 0) {
    sample_taxa <- head(data$Taxon, 3)
    cat(paste("Sample", target_domain, "taxa:", paste(sample_taxa, collapse = ", "), "\n"))
  }

  # Standardize column names - for genus level, the taxon column is "genus"
  if ("genus" %in% colnames(data)) data$Taxon <- data$genus
  if ("taxon" %in% colnames(data)) data$Taxon <- data$taxon
  if ("census_otu_count" %in% colnames(data)) data$Census_OTU_Count <- data$census_otu_count
  if ("ncbi_species_count" %in% colnames(data)) data$NCBI_Species_Count <- data$ncbi_species_count
  if ("ncbi_genome_count" %in% colnames(data)) data$NCBI_Genome_Count <- data$ncbi_genome_count
  if ("isolate_count" %in% colnames(data)) data$Isolate_Count <- data$isolate_count

  # Filter for meaningful data
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)

  # Use pre-calculated factors from merger scripts directly (no duplication needed)
  # Debug: Check factor column assignment
  cat(paste("DEBUG: After factor assignment - novelty_factor sample:", paste(head(data$novelty_factor, 3), collapse = ", "), "\n"))
  cat(paste("DEBUG: After factor assignment - overrepresentation_factor sample:", paste(head(data$overrepresentation_factor, 3), collapse = ", "), "\n"))

  # For ranking, use the pre-calculated factors directly (no weighting needed)
  ranking_novelty <- data$novelty_factor
  ranking_coverage <- data$overrepresentation_factor

  # Circle sizing - genome/isolate ratio
  data$Genome_Isolate_Ratio <- pmax(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), 1)
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  size_range <- range(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 3 + 9 * (data$Circle_Size_Raw - size_range[1]) / diff(size_range)

  # Add phylum information
  data <- add_phylum_info(data, level, domain)

  # Debug: Check factor ranges
  cat(paste("Novelty factor range:", round(min(data$novelty_factor, na.rm = TRUE), 2), "to", round(max(data$novelty_factor, na.rm = TRUE), 2), "\n"))
  cat(paste("Coverage factor range:", round(min(data$overrepresentation_factor, na.rm = TRUE), 2), "to", round(max(data$overrepresentation_factor, na.rm = TRUE), 2), "\n"))

  # Identify top taxa using pre-calculated factors with 1.0 threshold
  high_novelty_data <- data[data$novelty_factor > 1.0, ]
  high_coverage_data <- data[data$overrepresentation_factor > 1.0, ]

  cat(paste("Taxa with novelty > 1.0:", nrow(high_novelty_data), "\n"))
  cat(paste("Taxa with coverage > 1.0:", nrow(high_coverage_data), "\n"))

  # Debug: Check the Taxon column
  if (nrow(high_novelty_data) > 0) {
    cat(paste("DEBUG: high_novelty_data Taxon column class:", class(high_novelty_data$Taxon), "\n"))
    cat(paste("DEBUG: Sample high novelty taxa:", paste(head(high_novelty_data$Taxon, 3), collapse = ", "), "\n"))
  }

  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE

  if (nrow(high_novelty_data) > 0) {
    top_novelty_indices <- order(high_novelty_data$novelty_factor, decreasing = TRUE)[1:min(config$top_n, nrow(high_novelty_data))]
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_indices]
    cat(paste("DEBUG: top_novelty_indices length:", length(top_novelty_indices), "\n"))
    cat(paste("DEBUG: top_novelty_taxa length:", length(top_novelty_taxa), "\n"))
    cat(paste("DEBUG: Sample top novelty taxa:", paste(head(top_novelty_taxa, 3), collapse = ", "), "\n"))
    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
    cat(paste("Selected top", length(top_novelty_taxa), "novelty taxa\n"))
  }

  if (nrow(high_coverage_data) > 0) {
    top_coverage_indices <- order(high_coverage_data$overrepresentation_factor, decreasing = TRUE)[1:min(config$top_n, nrow(high_coverage_data))]
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_indices]
    cat(paste("DEBUG: top_coverage_indices length:", length(top_coverage_indices), "\n"))
    cat(paste("DEBUG: top_coverage_taxa length:", length(top_coverage_taxa), "\n"))
    cat(paste("DEBUG: Sample top coverage taxa:", paste(head(top_coverage_taxa, 3), collapse = ", "), "\n"))
    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
    cat(paste("Selected top", length(top_coverage_taxa), "coverage taxa\n"))
  }

  # Summary output
  cat(paste("Final", domain, level, "data:", nrow(data), "taxa,",
            sum(data$Is_Top_Novelty), "top novelty,",
            sum(data$Is_Top_Coverage), "top coverage\n"))

  return(data)
}

# Add phylum information for 16S data
add_phylum_info <- function(data, level, domain) {
  if (level == "genus") {
    # For genus level, assign a simple color grouping instead of complex phylum mapping
    # This avoids the complexity of mapping thousands of genera to phyla
    data$Phylum <- "Mixed"  # Simple placeholder for genus-level coloring
    cat(paste("Genus level: assigned", nrow(data), "genera to 'Mixed' group for coloring\n"))
  }

  return(data)
}

# Get phylum mapping for genus level (from mega_16s_stacked_visual.R)
get_phylum_for_taxa <- function(taxa, level) {
  # Load NCBI data to get phylum mapping
  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) {
    return(rep("Unknown", length(taxa)))
  }

  ncbi_data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat(paste("Loaded NCBI", level, "file with", nrow(ncbi_data), "rows\n"))
  cat(paste("Columns:", paste(head(colnames(ncbi_data), 10), collapse = ", "), "\n"))

  # Find the taxon name column - for genus level, the column is "genus"
  if (level == "genus") {
    possible_cols <- c("genus", "genus_name", "taxon", "Taxon")
  } else {
    possible_cols <- c("family_name", "genus_name", "family", "genus", "taxon", "Taxon")
  }
  taxon_col <- intersect(possible_cols, colnames(ncbi_data))[1]
  cat(paste("Using taxon column:", taxon_col, "\n"))

  if (is.na(taxon_col) || !all(c("lineage", "lineage_ranks") %in% colnames(ncbi_data))) {
    cat("ERROR: Required columns not found\n")
    return(rep("Unknown", length(taxa)))
  }

  # Extract phylum from lineage
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

  # Debug output
  mapped_count <- sum(result != "Unknown")
  cat(paste("Successfully mapped", mapped_count, "out of", length(taxa), "taxa to phyla\n"))
  if (mapped_count > 0) {
    unique_phyla <- unique(result[result != "Unknown"])
    cat(paste("Found phyla:", paste(head(unique_phyla, 5), collapse = ", "), "\n"))
  }

  return(result)
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

# Load 18S data function (adapted from updated 18S script)
load_18s_data <- function(level) {
  # Determine file name based on level
  filename <- paste0("18s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_18s, filename)

  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }

  cat(paste("Loading 18S", level, "data from:", filepath, "\n"))
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  cat(paste("DEBUG: Loaded", nrow(data), "total rows from CSV\n"))

  # Standardize column names - for genus level, the taxon column is "genus"
  if ("genus" %in% colnames(data)) data$Taxon <- data$genus
  if ("taxon" %in% colnames(data)) data$Taxon <- data$taxon
  if ("census_otu_count" %in% colnames(data)) data$Census_OTU_Count <- data$census_otu_count
  if ("ncbi_species_count" %in% colnames(data)) data$NCBI_Species_Count <- data$ncbi_species_count
  if ("ncbi_genome_count" %in% colnames(data)) data$NCBI_Genome_Count <- data$ncbi_genome_count
  if ("isolate_count" %in% colnames(data)) data$Isolate_Count <- data$isolate_count

  # Debug: Check if high overrepresentation genera are present before filtering
  high_overrep_genera <- c("Oryza", "Colletotrichum", "Coemansia", "Cryptococcus", "Aspergillus", "Saccharomyces")
  present_before_filter <- high_overrep_genera[high_overrep_genera %in% data$Taxon]
  cat(paste("DEBUG: High overrepresentation genera present before filtering:", paste(present_before_filter, collapse = ", "), "\n"))

  # Filter for meaningful data
  data_before_filter <- nrow(data)
  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
  cat(paste("DEBUG: After Census_OTU_Count > 0 & NCBI_Species_Count > 0 filter:", nrow(data), "rows (was", data_before_filter, ")\n"))

  # Debug: Check if high overrepresentation genera survived the filtering
  present_after_filter <- high_overrep_genera[high_overrep_genera %in% data$Taxon]
  cat(paste("DEBUG: High overrepresentation genera present after filtering:", paste(present_after_filter, collapse = ", "), "\n"))

  # Debug: Show sample of the 135 genera that survived
  cat(paste("DEBUG: Sample of 135 surviving genera:", paste(head(data$Taxon, 10), collapse = ", "), "\n"))

  # Add division information (exactly like individual 18S scripts)
  if (level == "genus") {
    # For genus level, assign divisions based on taxonomic knowledge
    # Start with "Other" and then apply manual overrides
    data$Division <- "Other"
    cat(paste("Genus level: starting with", nrow(data), "genera, will apply manual overrides\n"))

  # Debug: Check if high overrepresentation genera are present before manual overrides
  high_overrep_genera <- c("Oryza", "Colletotrichum", "Coemansia", "Cryptococcus", "Aspergillus", "Saccharomyces")
  present_high_overrep <- high_overrep_genera[high_overrep_genera %in% data$Taxon]
  cat(paste("DEBUG: High overrepresentation genera present before manual overrides:", paste(present_high_overrep, collapse = ", "), "\n"))
  if (length(present_high_overrep) > 0) {
    sample_factors <- data$overrepresentation_factor[data$Taxon %in% present_high_overrep[1:min(3, length(present_high_overrep))]]
    cat(paste("DEBUG: Sample overrepresentation factors:", paste(sample_factors, collapse = ", "), "\n"))
  }
  }

  # Apply manual overrides for specific taxa
  manual_overrides <- get_18s_manual_overrides()
  if (!is.null(manual_overrides)) {
    applied_count <- 0
    for (taxon in names(manual_overrides)) {
      if (taxon %in% data$Taxon) {
        data$Division[data$Taxon == taxon] <- manual_overrides[[taxon]]
        applied_count <- applied_count + 1
        cat(paste("  Applied override:", taxon, "->", manual_overrides[[taxon]], "\n"))
      }
    }
    cat(paste("Applied", applied_count, "manual overrides\n"))
  }

  # Filter divisions and check for remaining data
  data <- filter_divisions(data)
  if (nrow(data) == 0) {
    cat("WARNING: No data remaining after division filtering!\n")
    return(data.frame())
  }

  cat(paste("Final divisions:", paste(head(unique(data$Division), 8), collapse = ", "), "\n"))

  # Use pre-calculated factors from merger scripts directly (no duplication needed)
  # Debug: Check factor column assignment for 18S
  cat(paste("DEBUG 18S: After factor assignment - novelty_factor sample:", paste(head(data$novelty_factor, 3), collapse = ", "), "\n"))
  cat(paste("DEBUG 18S: After factor assignment - overrepresentation_factor sample:", paste(head(data$overrepresentation_factor, 3), collapse = ", "), "\n"))
  cat(paste("DEBUG 18S: Max overrepresentation_factor before filtering:", max(data$overrepresentation_factor, na.rm = TRUE), "\n"))
  cat(paste("DEBUG 18S: Taxa with overrepresentation_factor > 1.0 before filtering:", sum(data$overrepresentation_factor > 1.0, na.rm = TRUE), "\n"))

  # For ranking, use the pre-calculated factors directly (no weighting needed)
  ranking_novelty <- data$novelty_factor
  ranking_coverage <- data$overrepresentation_factor

  # Circle sizing - matching 16S script approach exactly
  data$Genome_Isolate_Ratio <- pmax(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), 1)
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)
  size_range <- range(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 3 + 9 * (data$Circle_Size_Raw - size_range[1]) / diff(size_range)

  # Debug: Check factor ranges
  cat(paste("18S Novelty factor range:", round(min(data$novelty_factor, na.rm = TRUE), 2), "to", round(max(data$novelty_factor, na.rm = TRUE), 2), "\n"))
  cat(paste("18S Coverage factor range:", round(min(data$overrepresentation_factor, na.rm = TRUE), 2), "to", round(max(data$overrepresentation_factor, na.rm = TRUE), 2), "\n"))

  # Identify top taxa using pre-calculated factors with 1.0 threshold
  high_novelty_data <- data[data$novelty_factor > 1.0, ]
  high_coverage_data <- data[data$overrepresentation_factor > 1.0, ]

  cat(paste("18S Taxa with novelty > 1.0:", nrow(high_novelty_data), "\n"))
  cat(paste("18S Taxa with coverage > 1.0:", nrow(high_coverage_data), "\n"))

  data$Is_Top_Novelty <- FALSE
  data$Is_Top_Coverage <- FALSE

  if (nrow(high_novelty_data) > 0) {
    top_novelty_indices <- order(high_novelty_data$novelty_factor, decreasing = TRUE)[1:min(config$top_n, nrow(high_novelty_data))]
    top_novelty_taxa <- high_novelty_data$Taxon[top_novelty_indices]
    data$Is_Top_Novelty[data$Taxon %in% top_novelty_taxa] <- TRUE
    cat(paste("18S Selected top", length(top_novelty_taxa), "novelty taxa\n"))
  }

  if (nrow(high_coverage_data) > 0) {
    top_coverage_indices <- order(high_coverage_data$overrepresentation_factor, decreasing = TRUE)[1:min(config$top_n, nrow(high_coverage_data))]
    top_coverage_taxa <- high_coverage_data$Taxon[top_coverage_indices]
    data$Is_Top_Coverage[data$Taxon %in% top_coverage_taxa] <- TRUE
    cat(paste("18S Selected top", length(top_coverage_taxa), "coverage taxa\n"))
  }

  # Summary output
  cat(paste("Final 18S", level, "data:", nrow(data), "taxa,",
            sum(data$Is_Top_Novelty), "top novelty,",
            sum(data$Is_Top_Coverage), "top coverage\n"))

  return(data)
}

# Manual overrides for 18S eukaryotic taxa (from 18S script)
get_18s_manual_overrides <- function() {
  manual_overrides <- list(
    "Kickxellales" = "Opisthokonta",
    "Lecanoromycetes" = "Opisthokonta",
    "Lipomycetaceae" = "Opisthokonta",
    "Amphibia" = "Opisthokonta",
    "Lepidosauria" = "Opisthokonta",
    "Neovahlkampfiidae" = "Discoba",
    "Labyrinthulaceae" = "Stramenopiles",
    # Add genus-level overrides
    "Neovahlkampfia" = "Discoba",
    "Labyrinthula" = "Stramenopiles",
    "Mastigamoeba" = "Evosea",
    "Stylonychia" = "Alveolata",
    "Trachelomonas" = "Discoba",
    "Gibberella" = "Opisthokonta",
    "Zea" = "Streptophyta",
    "Thalassiosira" = "Stramenopiles",
    "Pythium" = "Stramenopiles",
    # High overrepresentation genera
    "Oryza" = "Streptophyta",           # 23.0 overrepresentation
    "Colletotrichum" = "Opisthokonta",  # 10.4 overrepresentation
    "Coemansia" = "Opisthokonta",       # 12.167 overrepresentation
    "Cryptococcus" = "Opisthokonta",    # 5.0 overrepresentation
    "Aspergillus" = "Opisthokonta",     # 5.479 overrepresentation
    "Neurospora" = "Opisthokonta",      # 2.615 overrepresentation
    "Penicillium" = "Opisthokonta",     # 1.74 overrepresentation
    "Saccharomyces" = "Opisthokonta",   # 1.615 overrepresentation
    "Candida" = "Opisthokonta"          # 1.345 overrepresentation
  )

  cat(paste("Manual overrides defined for", length(manual_overrides), "taxa\n"))
  for (taxon in names(manual_overrides)) {
    cat(paste("  Override:", taxon, "->", manual_overrides[[taxon]], "\n"))
  }

  return(manual_overrides)
}

# Function to load family mapping from census data
load_family_mapping <- function() {
  # Load 18S taxonomic combinations for genus-family mapping
  combinations_file <- "../18S_censusparse/metadata/sanity_check/taxonomic_combinations_detailed.csv"

  if (!file.exists(combinations_file)) {
    cat("Warning: Taxonomic combinations file not found:", combinations_file, "\n")
    return(data.frame(genus = character(), family = character(), stringsAsFactors = FALSE))
  }

  # Read the file, skipping comment lines
  combinations <- read.csv(combinations_file, comment.char = "#", stringsAsFactors = FALSE)

  # Extract genus-family mapping
  family_map <- combinations[, c("genus", "family")]
  family_map <- family_map[!is.na(family_map$genus) & !is.na(family_map$family), ]

  # Remove duplicates (keep first occurrence)
  family_map <- family_map[!duplicated(family_map$genus), ]

  cat("Loaded family mapping for", nrow(family_map), "genera\n")
  return(family_map)
}

# Function to add family information to data
add_family_info <- function(data, family_mapping, gene_type) {
  if (nrow(data) == 0) return(data)

  # For 18S data, use the family mapping directly
  if (gene_type == "18s") {
    data$family <- family_mapping$family[match(data$genus, family_mapping$genus)]
  } else {
    # For 16S data (bacteria/archaea), family info might not be available in the same way
    # Set to NA for now, could be enhanced later with 16S family data
    data$family <- NA
  }

  # Fill missing families with "Unknown"
  data$family[is.na(data$family)] <- "Unknown"

  return(data)
}

# ggrepel parameters for genus-level annotations
get_repel_params <- function(color = "black", nudge_x = 0, nudge_y = 0) {
  list(
    color = color,
    size = config$text_size,
    fontface = "bold",
    max.overlaps = Inf,
    max.iter = 10000,        # More iterations for better positioning
    force = 15,              # Stronger repulsion force for genus-level data
    box.padding = 1.2,       # Much more space around labels
    point.padding = 1.0,     # More space around points
    min.segment.length = 0,
    segment.color = "grey40",
    segment.linewidth = 0.3,
    segment.alpha = 0.7,
    direction = "both",
    nudge_x = nudge_x,
    nudge_y = nudge_y
  )
}

# Create individual scatter plot for each domain
create_scatter_plot <- function(data, domain, level, master_colors) {
  cat(paste("Creating", domain, level, "plot with", nrow(data), "taxa\n"))

  if (nrow(data) == 0) {
    cat(paste("No data for", domain, level, "- creating empty plot\n"))
    return(ggplot() + theme_void())
  }

  # Debug: Check data ranges before setting limits
  cat(paste("DEBUG:", domain, "- Census_OTU_Count range:", round(min(data$Census_OTU_Count, na.rm = TRUE), 2), "to", round(max(data$Census_OTU_Count, na.rm = TRUE), 2), "\n"))
  cat(paste("DEBUG:", domain, "- NCBI_Species_Count range:", round(min(data$NCBI_Species_Count, na.rm = TRUE), 2), "to", round(max(data$NCBI_Species_Count, na.rm = TRUE), 2), "\n"))

  # Base plot with automatic limits based on actual data range
  x_range <- range(data$Census_OTU_Count, na.rm = TRUE)
  y_range <- range(data$NCBI_Species_Count, na.rm = TRUE)

  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  limits = c(max(0.1, x_range[1] * 0.1), x_range[2] * 10),  # Dynamic range based on data
                  expand = expansion(mult = c(0.25, 0.25))) +  # More expansion for annotation space
    scale_y_log10(labels = comma_format(),
                  limits = c(max(0.1, y_range[1] * 0.1), y_range[2] * 10),  # Dynamic range based on data
                  expand = expansion(mult = c(0.25, 0.25)))  # More expansion for annotation space

  # Get top data for highlighting
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  meaningful_data <- top_data

  # Debug: Check data structure
  cat(paste("DEBUG:", domain, "- Total data points:", nrow(data), "\n"))
  cat(paste("DEBUG:", domain, "- Background points:", nrow(data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]), "\n"))
  cat(paste("DEBUG:", domain, "- Top data points:", nrow(top_data), "\n"))
  cat(paste("DEBUG:", domain, "- Circle_Size range:", round(min(data$Circle_Size, na.rm = TRUE), 2), "to", round(max(data$Circle_Size, na.rm = TRUE), 2), "\n"))

  # Background points first (draw behind meaningful data)
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "lightgray", alpha = 0.3)
  }

  if (nrow(meaningful_data) > 0) {
    # For genus level, use simple color scheme based on novelty/coverage categories
    # Create color grouping based on top category
    meaningful_data$Color_Group <- ifelse(meaningful_data$Is_Top_Novelty, "Novel", "Overrepresented")

    # Simple color scheme for genus level
    genus_colors <- c("Novel" = "#2E86AB", "Overrepresented" = "#A23B72")
    plot_groups <- c("Novel", "Overrepresented")
    group_colors <- genus_colors

    # Add black outline layer
    p <- p + geom_point(data = meaningful_data, aes(size = Circle_Size),
                       color = "black", alpha = 0.9, stroke = 0.3)

    # Add colored fill layer
    p <- p + geom_point(data = meaningful_data,
                       aes(size = Circle_Size * 0.8, color = factor(Color_Group, levels = plot_groups)),
                       alpha = 0.9) +
      scale_color_manual(values = setNames(group_colors, plot_groups),
                        name = "Category",
                        guide = "none")  # Remove individual legends
  }

  # Minimal styling for stacking with solid background
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", alpha = 0.9, linewidth = 1.5) +
    theme_minimal() +
    theme(
      # Remove titles and axis labels for stacking
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = config$text_size, color = "black"),
      panel.background = element_rect(fill = "white", color = NA),  # Solid white background
      plot.background = element_rect(fill = "white", color = NA),   # Solid white plot background
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      plot.margin = margin(5, 5, 5, 5)
    )

  # Add annotations if there are meaningful data points
  if (nrow(meaningful_data) > 0) {
    # Determine which factor to display based on top category
    factor_value <- ifelse(meaningful_data$Is_Top_Novelty, meaningful_data$novelty_factor, meaningful_data$overrepresentation_factor)

    # Create labels - show full annotations for genus level
    meaningful_data$full_label <- paste0(meaningful_data$Taxon, " (",
                                        sprintf("%.1f", factor_value), "×)")

    # Updated annotation strategy:
    # Novel taxa (top 15) -> South (bottom) of plot
    # Overrepresented taxa (top 15) -> North (up and away from center)

    # Split data by category type
    novel_data <- meaningful_data[meaningful_data$Is_Top_Novelty, ]
    overrep_data <- meaningful_data[meaningful_data$Is_Top_Coverage, ]

    # Novel taxa: Always repel SOUTH (toward bottom of plot)
    if (nrow(novel_data) > 0) {
      novel_params <- get_repel_params(color = "black", nudge_y = -0.8)  # Strong downward push
      p <- p + do.call(ggrepel::geom_text_repel, c(
        list(data = novel_data, aes(label = full_label)),
        novel_params
      ))
    }

    # Overrepresented taxa: Always repel UP and away from center
    if (nrow(overrep_data) > 0) {
      overrep_params <- get_repel_params(color = "black", nudge_y = 0.8)  # Strong upward push
      p <- p + do.call(ggrepel::geom_text_repel, c(
        list(data = overrep_data, aes(label = full_label)),
        overrep_params
      ))
    }
  }

  return(p)
}

# Main function to create comprehensive mega genus visual
main <- function() {
  cat("Comprehensive Mega Genus Visual Creation\n")
  cat("========================================\n")

  # Define the grid structure - single row with 3 domains
  level <- "genus"
  domains_16s <- c("Bacteria", "Archaea")
  domain_18s <- "Eukaryota"

  # Get master color palettes
  master_colors_16s <- get_master_16s_color_palette()
  master_colors_18s <- get_master_18s_color_palette()

  # Load all data
  cat("\nLoading data for all domains...\n")
  bacteria_data <- load_16s_data(level, "Bacteria")
  archaea_data <- load_16s_data(level, "Archaea")
  eukaryota_data <- load_18s_data(level)

  # Create individual plots
  cat("\nCreating individual plots...\n")
  bacteria_plot <- create_scatter_plot(bacteria_data, "Bacteria", level, master_colors_16s)
  archaea_plot <- create_scatter_plot(archaea_data, "Archaea", level, master_colors_16s)
  eukaryota_plot <- create_scatter_plot(eukaryota_data, "Eukaryota", level, master_colors_18s)

  # Create column headers
  bacteria_header <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Bacteria", size = 16, fontface = "bold") +
    theme_void() + theme(plot.margin = margin(5, 5, 5, 5))

  archaea_header <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Archaea", size = 16, fontface = "bold") +
    theme_void() + theme(plot.margin = margin(5, 5, 5, 5))

  eukaryota_header <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Eukaryota", size = 16, fontface = "bold") +
    theme_void() + theme(plot.margin = margin(5, 5, 5, 5))

  # Create row label
  genus_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Genus", size = 16, fontface = "bold", angle = 90) +
    theme_void() + theme(plot.margin = margin(5, 5, 5, 5))

  # Combine into comprehensive layout
  cat("\nCombining plots into comprehensive layout...\n")

  # Header row with domain labels - adjusted for more square individual plots
  header_row <- plot_grid(
    NULL, bacteria_header, NULL, archaea_header, NULL, eukaryota_header,
    ncol = 6, rel_widths = c(0.08, 1.2, 0.05, 1.2, 0.05, 1.2)  # More square proportions for each plot
  )

  # Data row with plots - adjusted for more square individual plots
  data_row <- plot_grid(
    genus_label, bacteria_plot, NULL, archaea_plot, NULL, eukaryota_plot,
    ncol = 6, rel_widths = c(0.08, 1.2, 0.05, 1.2, 0.05, 1.2)  # More square proportions for each plot
  )

  # Combine header and data rows with solid background
  comprehensive_plot <- plot_grid(
    header_row, data_row,
    ncol = 1, rel_heights = c(0.1, 1)
  ) +
  theme(plot.background = element_rect(fill = "white", color = NA))

  # Save comprehensive plot
  output_file <- file.path(config$output_dir, "comprehensive_mega_genus_visual.png")
  cat(paste("Saving comprehensive plot to:", output_file, "\n"))

  ggsave(output_file, comprehensive_plot,
         width = config$plot_width, height = config$plot_height, dpi = config$dpi, units = "in")

  cat(paste("Comprehensive mega genus visual saved successfully!\n"))
  cat(paste("Dimensions:", config$plot_width, "×", config$plot_height, "inches\n"))
  cat(paste("Top taxa per category:", config$top_n, "\n"))

  # Export source data for all domains
  cat("\nExporting source data...\n")

  # Combine all top taxa data
  bacteria_top <- bacteria_data[bacteria_data$Is_Top_Novelty | bacteria_data$Is_Top_Coverage, ]
  archaea_top <- archaea_data[archaea_data$Is_Top_Novelty | archaea_data$Is_Top_Coverage, ]
  eukaryota_top <- eukaryota_data[eukaryota_data$Is_Top_Novelty | eukaryota_data$Is_Top_Coverage, ]

  # Load family information from taxonomic combinations
  cat("Loading family information from census data...\n")
  family_mapping <- load_family_mapping()

  # Add family information to each dataset
  bacteria_top <- add_family_info(bacteria_top, family_mapping, "16s")
  archaea_top <- add_family_info(archaea_top, family_mapping, "16s")
  eukaryota_top <- add_family_info(eukaryota_top, family_mapping, "18s")

  # Define essential columns to keep (remove unnecessary columns, include family)
  essential_cols <- c("genus", "family", "census_otu_count", "ncbi_genome_count", "ncbi_species_count",
                     "isolate_count", "isolate_percentage", "novelty_factor",
                     "overrepresentation_factor", "domain", "Division", "Is_Top_Novelty",
                     "Is_Top_Coverage", "Domain")

  if (nrow(bacteria_top) > 0) {
    bacteria_top$Domain <- "Bacteria"
    # Keep only essential columns
    bacteria_export <- bacteria_top[, intersect(essential_cols, colnames(bacteria_top))]
    bacteria_csv_file <- file.path(config$source_data_dir, "comprehensive_genus_bacteria_top_taxa.csv")
    bacteria_tsv_file <- file.path(config$source_data_dir, "comprehensive_genus_bacteria_top_taxa.tsv")
    write.csv(bacteria_export, bacteria_csv_file, row.names = FALSE)
    write.table(bacteria_export, bacteria_tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(paste("Bacteria top taxa exported to:", bacteria_csv_file, "\n"))
    cat(paste("Bacteria top taxa exported to:", bacteria_tsv_file, "\n"))
  }

  if (nrow(archaea_top) > 0) {
    archaea_top$Domain <- "Archaea"
    # Keep only essential columns
    archaea_export <- archaea_top[, intersect(essential_cols, colnames(archaea_top))]
    archaea_csv_file <- file.path(config$source_data_dir, "comprehensive_genus_archaea_top_taxa.csv")
    archaea_tsv_file <- file.path(config$source_data_dir, "comprehensive_genus_archaea_top_taxa.tsv")
    write.csv(archaea_export, archaea_csv_file, row.names = FALSE)
    write.table(archaea_export, archaea_tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(paste("Archaea top taxa exported to:", archaea_csv_file, "\n"))
    cat(paste("Archaea top taxa exported to:", archaea_tsv_file, "\n"))
  }

  if (nrow(eukaryota_top) > 0) {
    eukaryota_top$Domain <- "Eukaryota"
    # Keep only essential columns
    eukaryota_export <- eukaryota_top[, intersect(essential_cols, colnames(eukaryota_top))]
    eukaryota_csv_file <- file.path(config$source_data_dir, "comprehensive_genus_eukaryota_top_taxa.csv")
    eukaryota_tsv_file <- file.path(config$source_data_dir, "comprehensive_genus_eukaryota_top_taxa.tsv")
    write.csv(eukaryota_export, eukaryota_csv_file, row.names = FALSE)
    write.table(eukaryota_export, eukaryota_tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(paste("Eukaryota top taxa exported to:", eukaryota_csv_file, "\n"))
    cat(paste("Eukaryota top taxa exported to:", eukaryota_tsv_file, "\n"))
  }

  cat("\nComprehensive mega genus visual creation complete!\n")
}

# Execute main function
if (!interactive()) {
  main()
}
