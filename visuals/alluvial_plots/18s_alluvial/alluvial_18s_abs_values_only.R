#!/usr/bin/env Rscript
#
# 18S Eukaryotic Alluvial Plot Generator - CLEAN MERGED DATA APPROACH
# ===================================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Total Species
#
# Uses ONLY the merged 18S data for clean, consistent analysis
# Shows top divisions by total representation for Eukaryota domain
# Uses ABSOLUTE VALUES ONLY for both node heights and flow widths
# Advanced alluvial aesthetics with optimized flow guidance
#
# Author: Absolute Alluvial Team
# Date: 2025-01-20
#

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(RColorBrewer)
library(viridis)
library(scales)
library(tidyr)

cat("=== 18S rRNA Eukaryotic Alluvial Plot Generator (Clean Merged Data) ===\n\n")

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv"

# Load the 18S rRNA merged data
cat("Loading 18S rRNA merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 18S merged data file at: ", merged_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
data_18s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
cat("18S rRNA data loaded:", nrow(data_18s), "rows\n")

# Load 18S division data for .U. entries
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 18S census data file at: ", census_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

# Function to process eukaryotic data
process_eukaryotic_data <- function() {
  cat("\n=== Processing Eukaryota Domain ===\n")

  # Filter for matched phyla for Eukaryota domain
  matched_data <- data_18s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Eukaryota")

  cat("Matched eukaryotic divisions found:", nrow(matched_data), "\n")

  return(matched_data)
}

# Master eukaryotic color palette from mega 18S script
get_eukaryotic_colors <- function() {
  eukaryotic_color_map <- c(
    "Opisthokonta" = "#2E5984",      # Deep navy blue - most abundant
    "Streptophyta" = "#8B4513",      # Saddle brown - land plants
    "Stramenopiles" = "#FF8C00",     # Dark orange - diverse protists
    "Alveolata" = "#F0A0C0",         # Light pink - ciliates, dinoflagellates
    "Chlorophyta" = "#CD853F",       # Peru - green algae
    "Discoba" = "#2F4F4F",           # Dark slate gray - flagellates
    "Rhizaria" = "#DC143C",          # Crimson red - amoeboid protists
    "Rhodophyta" = "#B8860B",        # Dark goldenrod - red algae
    "Metamonada" = "#8FBC8F",        # Dark sea green - anaerobic flagellates
    "Cryptista" = "#A0522D",         # Sienna - cryptophytes
    "Amoebozoa" = "#4682B4",         # Steel blue - amoebas
    "Haptista" = "#D2691E"           # Chocolate - haptophytes
  )
  return(eukaryotic_color_map)
}

# Extended fallback colors for eukaryotic divisions
get_extended_eukaryotic_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

# Extract eukaryotic .U. entries from 18S census data
get_eukaryotic_u_entries <- function(census_data) {
  # Filter for eukaryotic .U. entries
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 10) %>%  # Lower threshold for eukaryotes
    filter(grepl("Eukaryota", Name_to_use) |
           (!grepl("Bacteria", Name_to_use) & !grepl("Archaea", Name_to_use))) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = "Eukaryota",
      match_status = "census_only"
    )

  return(u_entries)
}

u_entries <- get_eukaryotic_u_entries(census_division_data)
cat("Eukaryotic .U. entries found:", nrow(u_entries), "\n")

# Process the eukaryotic data
matched_data <- process_eukaryotic_data()
# Calculate totals from merged file + census totals for proper percentages
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)
total_size_count <- sum(census_division_data$size_count, na.rm = TRUE)  # Use full census for sequences
total_otu_count <- sum(census_division_data$otu_count, na.rm = TRUE)    # Use full census for OTUs

cat("Total counts (Eukaryota matched + full census for sequences/OTUs):\n")
cat("  Eukaryota Genomes:", scales::comma(total_genome_count), "\n")
cat("  18S Sequences:", scales::comma(total_size_count), "\n")
cat("  18S OTUs:", scales::comma(total_otu_count), "\n")
cat("  Eukaryota Species:", scales::comma(total_species_count), "\n")

# Calculate percentages for matched data
matched_data <- matched_data %>%
  mutate(
    otu_percentage = (census_otu_count / total_otu_count) * 100,
    size_percentage = (census_size_count / total_size_count) * 100,
    genome_pct_db = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100
  )

# Calculate percentages for .U. entries
if (nrow(u_entries) > 0) {
  u_entries <- u_entries %>%
    mutate(
      otu_percentage = (census_otu_count / total_otu_count) * 100,
      size_percentage = (census_size_count / total_size_count) * 100,
      genome_pct_db = 0,  # No genomes for .U. entries
      species_pct = 0     # No species for .U. entries
    )
}

# Combine matched eukaryotic data with eukaryotic .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Select top 8 entries by size_percentage (eukaryotic divisions + eukaryotic .U. entries)
n_top <- min(8, nrow(combined_data))
top_phyla <- combined_data %>%
  arrange(desc(size_percentage)) %>%
  head(n_top)

cat(paste("Top", n_top, "entries selected (eukaryotic divisions + eukaryotic .U.):", nrow(top_phyla), "\n"))

# Calculate "Other" category
other_data <- combined_data %>%
  filter(!phylum %in% top_phyla$phylum)

other_genome_count <- sum(other_data$ncbi_genome_count, na.rm = TRUE)
other_size_count <- sum(other_data$census_size_count, na.rm = TRUE)
other_otu_count <- sum(other_data$census_otu_count, na.rm = TRUE)
other_species_count <- sum(other_data$ncbi_species_count, na.rm = TRUE)

cat("Other category totals:\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Create long format data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Add top phyla data (including .U. entries)
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],        # Node 1: NCBI Total Genomes
      top_phyla$census_size_count[i],        # Node 2: 18S EukCensus Sequences
      top_phyla$census_otu_count[i],         # Node 3: 18S EukCensus OTUs
      top_phyla$ncbi_species_count[i]        # Node 4: NCBI Total Species
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
  stratum = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
  stringsAsFactors = FALSE
)

long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Fix aesthetic issues
cat("Applying advanced alluvial preprocessing...\n")

# 1. Ensure each phylum has a row at every axis (fill missing with 0)
long_data_f <- long_data %>%
  complete(x, phylum, fill = list(absolute_count = 0)) %>%
  group_by(phylum) %>%
  mutate(alluvium = cur_group_id()) %>%  # stable id per phylum
  ungroup()

# 2. Order phyla by size at the first axis (prettier strata stacking)
first_axis <- sort(unique(long_data_f$x))[1]
sizes_first <- long_data_f %>%
  filter(x == first_axis) %>%
  arrange(desc(absolute_count)) %>%
  select(phylum) %>% pull()

long_data_f <- long_data_f %>%
  mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

cat("Advanced preprocessing complete - data optimized for clean alluvial flows\n")

# Create professional color palette using eukaryotic colors
phyla_names <- unique(long_data_f$phylum)
n_colors <- length(phyla_names)

# Get eukaryotic color mappings
eukaryotic_colors <- get_eukaryotic_colors()
extended_colors <- get_extended_eukaryotic_colors()

# Assign colors to phyla
colors <- character(n_colors)
names(colors) <- phyla_names

for (i in seq_along(phyla_names)) {
  phylum_name <- phyla_names[i]
  clean_phylum <- gsub("^\\d+\\. ", "", phylum_name)

  if (phylum_name == "Other") {
    colors[i] <- "#808080"  # Gray for Other
  } else if (grepl("\\.U\\.", clean_phylum)) {
    colors[i] <- "#ffcc99"  # Light orange for eukaryotic .U. entries
  } else if (clean_phylum %in% names(eukaryotic_colors)) {
    colors[i] <- eukaryotic_colors[clean_phylum]
  } else {
    # Use extended fallback colors for unmapped phyla
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
  }
}

cat("Color mapping complete for", length(colors), "phyla\n")

# Create ADVANCED alluvial plot with optimized aesthetics
p_abs <- ggplot(
  long_data_f,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_count, fill = phylum)
) +
  # Emphasized flows with forward guidance to reduce crossings; earlier knot to tighten curves
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  # Minimal nodes (strata) - flows meet directly at axis lines
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = colors, name = "Division") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # No node titles - clean minimal appearance
    axis.text.x = element_blank(),
    # Y-axis ticks - larger and bold
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Legend styling - moved much further right
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1.8, "cm"),
    legend.margin = margin(l = 80),
    # Plot titles - larger
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 24, hjust = 0.5, color = "gray40"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    y = NULL,
    x = NULL
  )

# Save advanced alluvial plot with optimized dimensions
ggsave("alluvial_18s_abs_values_only.png", p_abs, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("alluvial_18s_abs_values_only.pdf", p_abs, width = 24, height = 10, dpi = 300, bg = "white")

cat("\n=== 18S Eukaryotic Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_18s_abs_values_only.png\n")
cat("  - alluvial_18s_abs_values_only.pdf\n")
cat("\n18S eukaryotic alluvial plot generated with:\n")
cat("  - Clean merged data approach\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional eukaryotic color scheme\n")
cat("  - Thin nodes and elegant flows\n")
cat("  - Legend positioned far right\n")
cat("  - Node titles moved down\n\n")
