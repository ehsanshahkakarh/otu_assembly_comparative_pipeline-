#!/usr/bin/env Rscript
# Source Data Generator - Minimal Dependencies Version
# ===================================================
# This script generates only the source data files without creating visualizations
# Avoids problematic libraries like ggplot2/scales that have GLIBCXX issues

cat("Source Data Generation (Minimal Dependencies)\n")
cat("==============================================\n")

# Simple path configuration
config <- list(
  data_dir_16s = "../../Eukcensus_merge/16s_merged/csv_results",
  data_dir_18s = "../../Eukcensus_merge/18s_merged/csv_results",
  source_data_dir = "source_data",
  ncbi_data_dir = "../../00ncbi_parse/csv_ncbi"
)

# Create output directory
if (!dir.exists(config$source_data_dir)) dir.create(config$source_data_dir, recursive = TRUE)

# Configuration constants
PLOT_CONFIG <- list(
  top_n = 10
)

# Function to calculate circle size based on isolate percentage
calculate_circle_size <- function(isolate_count, genome_count) {
  isolate_pct <- ifelse(genome_count > 0, (isolate_count / genome_count) * 100, 0)
  
  size <- ifelse(isolate_pct == 0, 25,
          ifelse(isolate_pct < 10, 20,
          ifelse(isolate_pct < 50, 15, 10)))
  
  return(size)
}

# Function to add phylum information
add_phylum_info <- function(data, level, domain) {
  if (level == "phylum") {
    data$Division <- data$Taxon
  } else {
    # For family level, try to get phylum from NCBI data
    ncbi_file <- file.path(config$ncbi_data_dir, paste0("ncbi_", tolower(domain), "_", level, "_counts.csv"))
    
    if (file.exists(ncbi_file)) {
      ncbi_data <- read.csv(ncbi_file, stringsAsFactors = FALSE)
      if ("phylum" %in% colnames(ncbi_data)) {
        # Match by taxon name
        match_idx <- match(data$Taxon, ncbi_data[[1]])
        data$Division <- ifelse(is.na(match_idx), "Unknown", ncbi_data$phylum[match_idx])
      } else {
        data$Division <- "Unknown"
      }
    } else {
      data$Division <- "Unknown"
    }
  }
  
  return(data)
}

# Function to save source data
save_source_data <- function(data, domain, level) {
  if (nrow(data) == 0) {
    cat(paste("  No data to save for", domain, level, "\n"))
    return()
  }
  
  filename <- paste0(domain, "_", level, "_source_data.csv")
  filepath <- file.path(config$source_data_dir, filename)
  
  write.csv(data, filepath, row.names = FALSE)
  
  # Report statistics
  novelty_count <- sum(data$Is_Top_Novelty == TRUE, na.rm = TRUE)
  coverage_count <- sum(data$Is_Top_Coverage == TRUE, na.rm = TRUE)
  
  cat(paste("  ✓ Saved:", filename, "\n"))
  cat(paste("    Rows:", nrow(data), "| Cols:", ncol(data), "\n"))
  cat(paste("    Novelty flags:", novelty_count, "| Coverage flags:", coverage_count, "\n"))
}

# Function to load and process 16S data
load_16s_data <- function(level, domain) {
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)
  
  cat(paste("Loading", domain, level, "data from:", filename, "\n"))
  
  if (!file.exists(filepath)) {
    cat(paste("  File not found:", filepath, "\n"))
    return(data.frame())
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data <- data[data$domain == domain & data$census_otu_count > 0 & data$ncbi_species_count > 0, ]
  
  # Remove any duplicate columns and unnecessary match detail columns
  duplicate_cols <- c("Census_OTU_Count", "NCBI_Species_Count", "NCBI_Genome_Count",
                     "Isolate_Count", "Novelty_Ratio", "Coverage_Factor",
                     "Domain", "Level", "coverage_percentage")
  match_detail_cols <- c("direct_matches", "taxid_matches", "lineage_matches")
  cols_to_remove <- c(duplicate_cols, match_detail_cols)[c(duplicate_cols, match_detail_cols) %in% colnames(data)]
  if (length(cols_to_remove) > 0) {
    data <- data[, !colnames(data) %in% cols_to_remove, drop = FALSE]
  }
  
  if (nrow(data) == 0) return(data.frame())
  
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
  
  # Filter to only include flagged taxa (those that will appear in figures)
  novelty_data <- data[data$Is_Top_Novelty == TRUE, ]
  coverage_data <- data[data$Is_Top_Coverage == TRUE & data$Is_Top_Novelty == FALSE, ]

  # Sort the filtered data
  if (nrow(novelty_data) > 0) novelty_data <- novelty_data[order(-novelty_data$novelty_factor), ]
  if (nrow(coverage_data) > 0) coverage_data <- coverage_data[order(-coverage_data$overrepresentation_factor), ]

  # Only include taxa that will appear in figures
  data <- rbind(novelty_data, coverage_data)

  # Move total_matches column to the end if it exists
  if ("total_matches" %in% colnames(data)) {
    total_matches_col <- data$total_matches
    data <- data[, !colnames(data) %in% "total_matches", drop = FALSE]
    data$total_matches <- total_matches_col
  }

  cat(paste("  Processed:", nrow(data), "taxa\n"))
  return(data)
}

# Function to load and process 18S data
load_18s_data <- function(level) {
  filename <- paste0("18s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_18s, filename)
  
  cat(paste("Loading Eukaryota", level, "data from:", filename, "\n"))
  
  if (!file.exists(filepath)) {
    cat(paste("  File not found:", filepath, "\n"))
    return(data.frame())
  }
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data <- data[data$census_otu_count > 0 & data$ncbi_species_count > 0, ]
  
  # Remove any duplicate columns and unnecessary match detail columns
  duplicate_cols <- c("Census_OTU_Count", "NCBI_Species_Count", "NCBI_Genome_Count",
                     "Isolate_Count", "Novelty_Ratio", "Coverage_Factor",
                     "Domain", "Level", "coverage_percentage")
  match_detail_cols <- c("direct_matches", "taxid_matches", "lineage_matches")
  cols_to_remove <- c(duplicate_cols, match_detail_cols)[c(duplicate_cols, match_detail_cols) %in% colnames(data)]
  if (length(cols_to_remove) > 0) {
    data <- data[, !colnames(data) %in% cols_to_remove, drop = FALSE]
  }
  
  if (nrow(data) == 0) return(data.frame())
  
  # Standardize column names
  colnames(data)[1] <- "Taxon"
  data$Circle_Size <- calculate_circle_size(data$isolate_count, data$ncbi_genome_count)
  
  # Add division information
  data <- add_phylum_info(data, level, "Eukaryota")
  
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
  
  # Filter to only include flagged taxa (those that will appear in figures)
  novelty_data <- data[data$Is_Top_Novelty == TRUE, ]
  coverage_data <- data[data$Is_Top_Coverage == TRUE & data$Is_Top_Novelty == FALSE, ]

  # Sort the filtered data
  if (nrow(novelty_data) > 0) novelty_data <- novelty_data[order(-novelty_data$novelty_factor), ]
  if (nrow(coverage_data) > 0) coverage_data <- coverage_data[order(-coverage_data$overrepresentation_factor), ]

  # Only include taxa that will appear in figures
  data <- rbind(novelty_data, coverage_data)

  # Move total_matches column to the end if it exists
  if ("total_matches" %in% colnames(data)) {
    total_matches_col <- data$total_matches
    data <- data[, !colnames(data) %in% "total_matches", drop = FALSE]
    data$total_matches <- total_matches_col
  }

  cat(paste("  Processed:", nrow(data), "taxa\n"))
  return(data)
}

# Main execution
main <- function() {
  cat("Generating source data files...\n")
  
  # Define the grid structure
  levels <- c("phylum", "family")
  domains_16s <- c("Bacteria", "Archaea")
  domain_18s <- "Eukaryota"
  
  # Process 16S data
  cat("\nProcessing 16S data...\n")
  for (level in levels) {
    for (domain in domains_16s) {
      cat(paste("\n--- Processing", domain, level, "---\n"))
      data <- load_16s_data(level, domain)
      save_source_data(data, domain, level)
    }
  }
  
  # Process 18S data
  cat("\nProcessing 18S data...\n")
  for (level in levels) {
    cat(paste("\n--- Processing Eukaryota", level, "---\n"))
    data <- load_18s_data(level)
    save_source_data(data, "Eukaryota", level)
  }
  
  cat("\n✅ Source data generation complete!\n")
  cat(paste("   Files saved to:", config$source_data_dir, "\n"))
}

# Run the main function
main()
