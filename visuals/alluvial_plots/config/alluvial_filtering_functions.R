#!/usr/bin/env Rscript
# ============================================================================
# Alluvial Plot Filtering Functions
# ============================================================================
# Helper functions to load configuration and apply filtering logic
# for 18S and 16S alluvial plots
#
# Author: Alluvial Config Team  
# Date: 2025-01-17
# ============================================================================

# Load required libraries
if (!require(yaml, quietly = TRUE)) {
  stop("Package 'yaml' is required but not installed")
}
if (!require(dplyr, quietly = TRUE)) {
  stop("Package 'dplyr' is required but not installed")
}

# Load alluvial filtering configuration
load_alluvial_config <- function(config_path = "config/alluvial_filtering_config.yaml") {
  # Try multiple possible paths
  possible_paths <- c(
    config_path,
    file.path("config", "alluvial_filtering_config.yaml"),
    file.path("..", "config", "alluvial_filtering_config.yaml"),
    file.path("alluvial_plots", "config", "alluvial_filtering_config.yaml")
  )
  
  config_file <- NULL
  for (path in possible_paths) {
    if (file.exists(path)) {
      config_file <- path
      break
    }
  }
  
  if (is.null(config_file)) {
    stop(paste("Could not find alluvial filtering config. Tried paths:", 
               paste(possible_paths, collapse = ", ")))
  }
  
  cat(paste("Loading alluvial config from:", config_file, "\n"))
  return(yaml.load_file(config_file))
}

# Load and validate data based on configuration
load_alluvial_data <- function(domain_config, config_name) {
  cat(paste("\n=== Loading data for", config_name, "===\n"))
  
  # Load merged data
  merged_path <- domain_config$data_sources$merged_data
  if (!file.exists(merged_path)) {
    stop(paste("ERROR: Cannot find merged data file at:", merged_path))
  }
  
  merged_data <- read.csv(merged_path, stringsAsFactors = FALSE)
  cat(paste("Merged data loaded:", nrow(merged_data), "rows\n"))
  
  # Load census data if specified
  census_data <- NULL
  if (!is.null(domain_config$data_sources$census_data)) {
    census_path <- domain_config$data_sources$census_data
    if (file.exists(census_path)) {
      census_data <- read.csv(census_path, stringsAsFactors = FALSE)
      cat(paste("Census data loaded:", nrow(census_data), "rows\n"))
    } else {
      cat(paste("Warning: Census data file not found at:", census_path, "\n"))
    }
  }
  
  return(list(merged = merged_data, census = census_data))
}

# Apply basic data filtering
apply_basic_filtering <- function(data, domain_config, global_config) {
  cat("Applying basic data filtering...\n")
  
  # Filter for domain
  if (!is.null(domain_config$domain_filter)) {
    data <- data %>% filter(domain == domain_config$domain_filter)
    cat(paste("After domain filtering:", nrow(data), "rows\n"))
  }
  
  # Filter for matched status
  data <- data %>% filter(match_status == 'matched')
  cat(paste("After match status filtering:", nrow(data), "rows\n"))
  
  # Remove NA and empty values
  required_cols <- domain_config$required_columns
  for (col in required_cols) {
    if (col %in% colnames(data)) {
      data <- data %>% filter(!is.na(.data[[col]]) & .data[[col]] != "" & .data[[col]] != "N/A")
    }
  }
  cat(paste("After NA/empty filtering:", nrow(data), "rows\n"))
  
  # Apply minimum thresholds
  if ("census_otu_count" %in% colnames(data)) {
    data <- data %>% filter(census_otu_count >= global_config$min_otu_count)
  }
  if ("ncbi_genome_count" %in% colnames(data)) {
    data <- data %>% filter(ncbi_genome_count >= global_config$min_genome_count)
  }
  if ("ncbi_species_count" %in% colnames(data)) {
    data <- data %>% filter(ncbi_species_count >= global_config$min_species_count)
  }
  
  cat(paste("After threshold filtering:", nrow(data), "rows\n"))
  
  # Exclude specific taxa
  if (!is.null(domain_config$exclude_taxa)) {
    data <- data %>% filter(!phylum %in% domain_config$exclude_taxa)
    cat(paste("After taxa exclusion:", nrow(data), "rows\n"))
  }
  
  return(data)
}

# Extract .U. entries from census data
extract_u_entries <- function(census_data, domain_config, global_config) {
  if (is.null(census_data) || !domain_config$u_entries$enabled) {
    return(data.frame())
  }
  
  cat("Extracting .U. (unidentified) entries...\n")
  
  # Filter for .U. entries
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= domain_config$u_entries$min_otu_count)
  
  # Apply domain patterns
  for (pattern in domain_config$u_entries$domain_patterns) {
    if (startsWith(pattern, "!")) {
      # Exclusion pattern
      exclude_pattern <- substring(pattern, 2)
      u_entries <- u_entries %>% filter(!grepl(exclude_pattern, Name_to_use))
    } else {
      # Inclusion pattern
      u_entries <- u_entries %>% filter(grepl(pattern, Name_to_use))
    }
  }
  
  if (nrow(u_entries) == 0) {
    cat("No .U. entries found matching criteria\n")
    return(data.frame())
  }
  
  # Convert to standard format
  u_entries_formatted <- u_entries %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      match_status = "u_entry",
      domain = domain_config$domain_filter
    )
  
  cat(paste("Found", nrow(u_entries_formatted), ".U. entries\n"))
  return(u_entries_formatted)
}

# Apply filtering strategy
apply_filtering_strategy <- function(data, domain_config) {
  strategy <- domain_config$filtering$strategy
  cat(paste("Applying filtering strategy:", strategy, "\n"))
  
  if (strategy == "top_abundance") {
    return(apply_top_abundance_filter(data, domain_config))
  } else if (strategy == "threshold_based") {
    return(apply_threshold_filter(data, domain_config))
  } else if (strategy == "hybrid") {
    return(apply_hybrid_filter(data, domain_config))
  } else {
    stop(paste("Unknown filtering strategy:", strategy))
  }
}

# Top abundance filtering
apply_top_abundance_filter <- function(data, domain_config) {
  top_n <- domain_config$filtering$top_n %||% domain_config$filtering$hybrid$min_top_abundance %||% 8
  
  selected_data <- data %>%
    arrange(desc(census_size_count)) %>%
    head(top_n)
  
  cat(paste("Selected top", nrow(selected_data), "taxa by abundance\n"))
  return(selected_data)
}

# Threshold-based filtering  
apply_threshold_filter <- function(data, domain_config) {
  novelty_thresh <- domain_config$filtering$novelty_threshold %||% 1.0
  overrep_thresh <- domain_config$filtering$overrepresentation_threshold %||% 1.0
  
  selected_data <- data %>%
    filter(novelty_factor > novelty_thresh | overrepresentation_factor > overrep_thresh)
  
  cat(paste("Selected", nrow(selected_data), "taxa meeting threshold criteria\n"))
  return(selected_data)
}

# Hybrid filtering (combines abundance and threshold approaches)
apply_hybrid_filter <- function(data, domain_config) {
  hybrid_config <- domain_config$filtering$hybrid
  
  # Get top taxa by abundance
  min_top <- hybrid_config$min_top_abundance %||% 6
  top_abundance <- data %>%
    arrange(desc(census_size_count)) %>%
    head(min_top)
  
  # Get high-factor taxa
  novelty_thresh <- hybrid_config$novelty_threshold %||% 10.0
  overrep_thresh <- hybrid_config$overrepresentation_threshold %||% 1.5
  
  high_factor <- data %>%
    filter(novelty_factor > novelty_thresh | overrepresentation_factor > overrep_thresh) %>%
    filter(!phylum %in% top_abundance$phylum)  # Don't duplicate
  
  # Combine and limit
  combined <- bind_rows(top_abundance, high_factor)
  max_total <- hybrid_config$max_total_taxa %||% 12
  
  if (nrow(combined) > max_total) {
    # Prioritize by a combined score
    combined <- combined %>%
      mutate(combined_score = log10(census_size_count) + novelty_factor + overrepresentation_factor) %>%
      arrange(desc(combined_score)) %>%
      head(max_total) %>%
      select(-combined_score)
  }
  
  cat(paste("Hybrid selection:", nrow(top_abundance), "by abundance +", 
            nrow(high_factor), "by factors =", nrow(combined), "total\n"))
  
  return(combined)
}

# Create "Other" category from remaining data
create_other_category <- function(all_data, selected_data, global_config) {
  if (!global_config$include_other_category) {
    return(data.frame())
  }

  other_data <- all_data %>%
    filter(!phylum %in% selected_data$phylum)

  if (nrow(other_data) == 0) {
    return(data.frame())
  }

  # Aggregate other data
  other_summary <- data.frame(
    phylum = global_config$other_category_name,
    census_otu_count = sum(other_data$census_otu_count, na.rm = TRUE),
    census_size_count = sum(other_data$census_size_count, na.rm = TRUE),
    ncbi_genome_count = sum(other_data$ncbi_genome_count, na.rm = TRUE),
    ncbi_species_count = sum(other_data$ncbi_species_count, na.rm = TRUE),
    match_status = "aggregated",
    domain = unique(selected_data$domain)[1],
    stringsAsFactors = FALSE
  )

  cat(paste("Created 'Other' category with", nrow(other_data), "aggregated taxa\n"))
  return(other_summary)
}

# Calculate percentages for alluvial data
calculate_percentages <- function(data) {
  # Calculate total counts for percentage calculations
  total_otu <- sum(data$census_otu_count, na.rm = TRUE)
  total_size <- sum(data$census_size_count, na.rm = TRUE)
  total_genomes <- sum(data$ncbi_genome_count, na.rm = TRUE)
  total_species <- sum(data$ncbi_species_count, na.rm = TRUE)

  data <- data %>%
    mutate(
      otu_percentage = ifelse(total_otu > 0, (census_otu_count / total_otu) * 100, 0),
      size_percentage = ifelse(total_size > 0, (census_size_count / total_size) * 100, 0),
      genome_percentage = ifelse(total_genomes > 0, (ncbi_genome_count / total_genomes) * 100, 0),
      species_percentage = ifelse(total_species > 0, (ncbi_species_count / total_species) * 100, 0)
    )

  return(data)
}

# Main function to process alluvial data
process_alluvial_data <- function(domain_key, config_path = NULL) {
  cat(paste("\n", paste(rep("=", 60), collapse = ""), "\n"))
  cat(paste("PROCESSING ALLUVIAL DATA FOR:", toupper(domain_key), "\n"))
  cat(paste(paste(rep("=", 60), collapse = ""), "\n"))

  # Load configuration
  config <- load_alluvial_config(config_path)
  global_config <- config$global
  domain_config <- config[[domain_key]]

  if (is.null(domain_config)) {
    stop(paste("No configuration found for domain:", domain_key))
  }

  # Load data
  data_list <- load_alluvial_data(domain_config, domain_key)
  merged_data <- data_list$merged
  census_data <- data_list$census

  # Apply basic filtering
  filtered_data <- apply_basic_filtering(merged_data, domain_config, global_config)

  # Extract .U. entries
  u_entries <- extract_u_entries(census_data, domain_config, global_config)

  # Combine matched data with .U. entries
  combined_data <- bind_rows(filtered_data, u_entries)

  # Apply filtering strategy
  selected_data <- apply_filtering_strategy(combined_data, domain_config)

  # Create "Other" category
  other_category <- create_other_category(combined_data, selected_data, global_config)

  # Final dataset
  final_data <- bind_rows(selected_data, other_category)

  # Calculate percentages
  final_data <- calculate_percentages(final_data)

  cat(paste("\nFINAL DATASET SUMMARY:\n"))
  cat(paste("- Total taxa selected:", nrow(final_data), "\n"))
  cat(paste("- Includes 'Other' category:", nrow(other_category) > 0, "\n"))
  cat(paste("- Total OTUs:", sum(final_data$census_otu_count, na.rm = TRUE), "\n"))
  cat(paste("- Total sequences:", sum(final_data$census_size_count, na.rm = TRUE), "\n"))
  cat(paste("- Total genomes:", sum(final_data$ncbi_genome_count, na.rm = TRUE), "\n"))
  cat(paste("- Total species:", sum(final_data$ncbi_species_count, na.rm = TRUE), "\n"))

  return(list(
    data = final_data,
    config = list(global = global_config, domain = domain_config),
    summary = list(
      total_taxa = nrow(final_data),
      has_other = nrow(other_category) > 0,
      total_otus = sum(final_data$census_otu_count, na.rm = TRUE),
      total_sequences = sum(final_data$census_size_count, na.rm = TRUE),
      total_genomes = sum(final_data$ncbi_genome_count, na.rm = TRUE),
      total_species = sum(final_data$ncbi_species_count, na.rm = TRUE)
    )
  ))
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x
