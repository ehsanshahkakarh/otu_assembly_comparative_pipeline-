#!/usr/bin/env Rscript
# Generate Source Data for Mega Visual Figures
# Created: 2025-01-26
# Purpose: Extract genome/isolate ratio data for all mega visual plots

library(dplyr)
library(readr)

# Configuration
config <- list(
  # Data directories
  data_16s = "Eukcensus_merge/16s_merged/csv_results",
  data_18s = "Eukcensus_merge/18s_merged/csv_results", 
  output_dir = "final_visuals",
  
  # Output file names
  output_16s_bacteria_phylum = "source_data_16s_bacteria_phylum.csv",
  output_16s_bacteria_family = "source_data_16s_bacteria_family.csv",
  output_16s_archaea_phylum = "source_data_16s_archaea_phylum.csv", 
  output_16s_archaea_family = "source_data_16s_archaea_family.csv",
  output_18s_divisions = "source_data_18s_divisions.csv",
  output_18s_family = "source_data_18s_family.csv"
)

# Create output directory
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
  cat("Created output directory:", config$output_dir, "\n")
}

# Function to calculate genome/isolate ratio
calculate_genome_isolate_ratio <- function(data) {
  data$Genome_Isolate_Ratio <- ifelse(data$isolate_count > 0,
                                      data$ncbi_genome_count / data$isolate_count,
                                      max(data$ncbi_genome_count / pmax(data$isolate_count, 1), na.rm = TRUE))
  return(data)
}

# Function to process 16S data
process_16s_data <- function(level, domain_filter) {
  file_path <- file.path(config$data_16s, paste0("16s_ncbi_merged_clean_", level, ".csv"))
  
  if (!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  cat("Processing 16S", level, "data for", domain_filter, "\n")
  
  # Load data
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # Filter by domain
  data <- data %>% filter(domain == domain_filter)
  
  # Calculate genome/isolate ratio
  data <- calculate_genome_isolate_ratio(data)
  
  # Select relevant columns for source data
  source_data <- data %>%
    select(
      Taxon = all_of(level),
      Domain = domain,
      Census_OTU_Count = census_otu_count,
      Census_Size_Count = census_size_count,
      NCBI_Genome_Count = ncbi_genome_count,
      NCBI_Species_Count = ncbi_species_count,
      Isolate_Count = isolate_count,
      Isolate_Percentage = isolate_percentage,
      Genome_Isolate_Ratio,
      Coverage_Percentage = coverage_percentage,
      Match_Status = match_status
    ) %>%
    arrange(desc(Genome_Isolate_Ratio))
  
  return(source_data)
}

# Function to process 18S data
process_18s_data <- function(level) {
  file_path <- file.path(config$data_18s, paste0("18s_ncbi_merged_clean_", level, ".csv"))
  
  if (!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  cat("Processing 18S", level, "data\n")
  
  # Load data
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # Calculate genome/isolate ratio
  data <- calculate_genome_isolate_ratio(data)
  
  # Select relevant columns for source data
  source_data <- data %>%
    select(
      Taxon = all_of(level),
      Domain = domain,
      Census_OTU_Count = census_otu_count,
      Census_Size_Count = census_size_count,
      NCBI_Genome_Count = ncbi_genome_count,
      NCBI_Species_Count = ncbi_species_count,
      Isolate_Count = isolate_count,
      Isolate_Percentage = isolate_percentage,
      Genome_Isolate_Ratio,
      Coverage_Percentage = coverage_percentage,
      Match_Status = match_status
    ) %>%
    arrange(desc(Genome_Isolate_Ratio))
  
  return(source_data)
}

# Generate all source data files
cat("=== Generating Source Data Files ===\n\n")

# 16S Bacteria Phylum
bacteria_phylum <- process_16s_data("phylum", "Bacteria")
if (!is.null(bacteria_phylum)) {
  output_path <- file.path(config$output_dir, config$output_16s_bacteria_phylum)
  write_csv(bacteria_phylum, output_path)
  cat("✓ Saved:", output_path, "(", nrow(bacteria_phylum), "taxa )\n")
}

# 16S Bacteria Family
bacteria_family <- process_16s_data("family", "Bacteria")
if (!is.null(bacteria_family)) {
  output_path <- file.path(config$output_dir, config$output_16s_bacteria_family)
  write_csv(bacteria_family, output_path)
  cat("✓ Saved:", output_path, "(", nrow(bacteria_family), "taxa )\n")
}

# 16S Archaea Phylum
archaea_phylum <- process_16s_data("phylum", "Archaea")
if (!is.null(archaea_phylum)) {
  output_path <- file.path(config$output_dir, config$output_16s_archaea_phylum)
  write_csv(archaea_phylum, output_path)
  cat("✓ Saved:", output_path, "(", nrow(archaea_phylum), "taxa )\n")
}

# 16S Archaea Family
archaea_family <- process_16s_data("family", "Archaea")
if (!is.null(archaea_family)) {
  output_path <- file.path(config$output_dir, config$output_16s_archaea_family)
  write_csv(archaea_family, output_path)
  cat("✓ Saved:", output_path, "(", nrow(archaea_family), "taxa )\n")
}

# 18S Divisions (Phylum level)
divisions_data <- process_18s_data("division")
if (!is.null(divisions_data)) {
  output_path <- file.path(config$output_dir, config$output_18s_divisions)
  write_csv(divisions_data, output_path)
  cat("✓ Saved:", output_path, "(", nrow(divisions_data), "taxa )\n")
}

# 18S Family
family_18s <- process_18s_data("family")
if (!is.null(family_18s)) {
  output_path <- file.path(config$output_dir, config$output_18s_family)
  write_csv(family_18s, output_path)
  cat("✓ Saved:", output_path, "(", nrow(family_18s), "taxa )\n")
}

cat("\n=== Source Data Generation Complete ===\n")
cat("All files saved to:", config$output_dir, "\n")
cat("Files contain genome/isolate ratios and all relevant metrics for each taxon.\n")
