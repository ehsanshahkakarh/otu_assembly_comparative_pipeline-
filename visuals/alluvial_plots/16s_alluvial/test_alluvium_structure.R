#!/usr/bin/env Rscript

# Test the alluvium structure to see if Pseudomonadota flows are correct

# Simulate the data structure
long_data <- data.frame(
  alluvium = c(1, 1, 1, 1,  # Pseudomonadota
               2, 2, 2, 2,  # Bacillota
               3, 3, 3, 3), # Campylobacterota
  phylum = c(rep("1. Pseudomonadota", 4),
             rep("2. Bacillota", 4),
             rep("3. Campylobacterota", 4)),
  x = rep(c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"), 3),
  stratum = rep(c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"), 3),
  percentage = c(
    57.83, 28.68, 22.93, 39.53,  # Pseudomonadota - DOMINANT
    23.65, 16.11, 21.31, 17.06,  # Bacillota
    5.84, 0.90, 1.01, 0.89       # Campylobacterota
  ),
  stringsAsFactors = FALSE
)

cat("=== ORIGINAL DATA ===\n")
print(long_data)

# Apply the preprocessing from the script
library(dplyr)

long_data_f <- long_data %>%
  group_by(phylum) %>%
  mutate(alluvium = cur_group_id()) %>%  # This is line 241 in the script
  ungroup()

cat("\n=== AFTER PREPROCESSING (line 241) ===\n")
print(long_data_f)

cat("\n=== CHECKING ALLUVIUM IDS ===\n")
cat("Pseudomonadota alluvium IDs:", unique(long_data_f$alluvium[long_data_f$phylum == "1. Pseudomonadota"]), "\n")
cat("Bacillota alluvium IDs:", unique(long_data_f$alluvium[long_data_f$phylum == "2. Bacillota"]), "\n")
cat("Campylobacterota alluvium IDs:", unique(long_data_f$alluvium[long_data_f$phylum == "3. Campylobacterota"]), "\n")

# Check if percentages are preserved
cat("\n=== PSEUDOMONADOTA PERCENTAGES ===\n")
pseudo_data <- long_data_f %>% filter(phylum == "1. Pseudomonadota")
print(pseudo_data[, c("x", "percentage", "alluvium")])

cat("\n=== SUMMARY ===\n")
cat("Total rows:", nrow(long_data_f), "\n")
cat("Unique phyla:", length(unique(long_data_f$phylum)), "\n")
cat("Unique alluvium IDs:", length(unique(long_data_f$alluvium)), "\n")
cat("Pseudomonadota rows:", sum(long_data_f$phylum == "1. Pseudomonadota"), "\n")
cat("Pseudomonadota total percentage:", sum(pseudo_data$percentage), "\n")

