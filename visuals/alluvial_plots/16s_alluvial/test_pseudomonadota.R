#!/usr/bin/env Rscript

# Simple test to check Pseudomonadota in R data processing
library(dplyr)

PROCESS_DOMAIN <- "Bacteria"

# Load data
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)

cat("=== R Data Loading Test ===\n")
cat("Total entries:", nrow(data_16s), "\n")

# Filter for bacteria
bacteria_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(domain == "Bacteria")

cat("Bacteria entries:", nrow(bacteria_data), "\n")

# Check for Pseudomonadota
pseudo <- bacteria_data %>% filter(phylum == "Pseudomonadota")
cat("Pseudomonadota entries:", nrow(pseudo), "\n")

if (nrow(pseudo) > 0) {
  cat("Pseudomonadota data:\n")
  cat("  Genomes:", pseudo$ncbi_genome_count, "\n")
  cat("  Species:", pseudo$ncbi_species_count, "\n")
  cat("  OTUs:", pseudo$census_otu_count, "\n")
  cat("  Sequences:", pseudo$census_size_count, "\n")
} else {
  cat("ERROR: Pseudomonadota not found!\n")
  cat("Available phyla:\n")
  phyla_counts <- bacteria_data %>% count(phylum, sort = TRUE)
  print(head(phyla_counts, 10))
}

cat("=== Test Complete ===\n")
