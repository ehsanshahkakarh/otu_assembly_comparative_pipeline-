#!/usr/bin/env Rscript

# Quick check of what phyla are in the plot data
library(dplyr)
library(yaml)

PROCESS_DOMAIN <- "Bacteria"

# Load data
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv"

data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

# Process domain-specific data
matched_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
  filter(domain == PROCESS_DOMAIN)

# Get .U. entries
u_entries <- census_division_data %>%
  filter(grepl("Bacteria\\.U\\.", Name_to_use)) %>%
  filter(otu_count >= 10) %>%
  select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
  mutate(
    ncbi_genome_count = 0,
    ncbi_species_count = 0,
    domain = PROCESS_DOMAIN,
    match_status = "census_only"
  )

# Combine
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)

# Select top 12
top_phyla <- combined_data %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    otu_pct = (census_otu_count / total_otu_count) * 100,
    size_pct = (census_size_count / total_size_count) * 100,
    total_representation = genome_pct + species_pct + otu_pct + size_pct
  ) %>%
  arrange(desc(total_representation)) %>%
  head(12)

cat("\n=== TOP 12 BACTERIA PHYLA ===\n")
for (i in 1:nrow(top_phyla)) {
  cat(sprintf("%2d. %-30s | Genome: %6.2f%% | Species: %6.2f%% | OTU: %6.2f%% | Size: %6.2f%% | Total: %7.2f\n",
              i,
              top_phyla$phylum[i],
              top_phyla$genome_pct[i],
              top_phyla$species_pct[i],
              top_phyla$otu_pct[i],
              top_phyla$size_pct[i],
              top_phyla$total_representation[i]))
}

cat("\n=== CHECKING PSEUDOMONADOTA ===\n")
pseudo <- top_phyla %>% filter(phylum == "Pseudomonadota")
if (nrow(pseudo) > 0) {
  cat("✅ Pseudomonadota IS in top 12\n")
  cat(sprintf("   Genome: %.2f%%, Species: %.2f%%, OTU: %.2f%%, Size: %.2f%%\n",
              pseudo$genome_pct, pseudo$species_pct, pseudo$otu_pct, pseudo$size_pct))
} else {
  cat("❌ Pseudomonadota NOT in top 12\n")
}

