#!/usr/bin/env Rscript

# Debug script to check data structure
library(dplyr)

# Load data
merged_data_path <- "../../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv"

data_18s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

cat("=== DATA STRUCTURE DEBUG ===\n")
cat("Merged data columns:", paste(colnames(data_18s), collapse = ", "), "\n")
cat("Census data columns:", paste(colnames(census_division_data), collapse = ", "), "\n")

# Process eukaryotic data
process_eukaryotic_data <- function() {
  matched_data <- data_18s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Eukaryota")
  
  return(matched_data)
}

matched_data <- process_eukaryotic_data()
cat("Matched data rows:", nrow(matched_data), "\n")
cat("Matched data columns:", paste(colnames(matched_data), collapse = ", "), "\n")

# Check first few rows
cat("First 3 rows of matched data:\n")
print(head(matched_data[, c("phylum", "ncbi_genome_count", "census_size_count", "census_otu_count", "ncbi_species_count")], 3))

# Create simple test data
test_data <- data.frame(
  alluvium = rep(1:3, each = 4),
  phylum = rep(c("Test1", "Test2", "Test3"), each = 4),
  x = rep(c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"), 3),
  stratum = rep(c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"), 3),
  percentage = c(10, 20, 30, 40, 15, 25, 35, 45, 5, 15, 25, 35),
  stringsAsFactors = FALSE
)

cat("Test data structure:\n")
print(test_data)

cat("Test data x levels:", paste(unique(test_data$x), collapse = ", "), "\n")
cat("Test data stratum levels:", paste(unique(test_data$stratum), collapse = ", "), "\n")
