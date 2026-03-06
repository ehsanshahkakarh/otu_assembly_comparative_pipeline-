#!/usr/bin/env Rscript

# ============================================================================
# 16S Prokaryotic STACKED Alluvial Plot - ABSOLUTE + PERCENTAGE
# ============================================================================
# This script creates stacked alluvial plots for 16S prokaryotic data
# showing both Bacteria and Archaea with ABSOLUTE (top) and PERCENTAGE (bottom)
# values, with a single shared legend.
#
# Data Flow: Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species
#
# Key Features:
# - Bacteria absolute values (top) + Bacteria percentage (bottom)
# - Archaea absolute values (top) + Archaea percentage (bottom)
# - Single shared legend for all plots
# - Uses patchwork for plot composition
# - Professional bacteria and archaea color schemes
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)
library(patchwork)  # For combining plots

# Load shared color mapping functions
source("../../../shared_config/color_mapping_functions.R")

cat("=== 16S Prokaryotic STACKED Alluvial Plot (Absolute + Percentage) ===\n")
cat("Creating stacked visualization: Bacteria & Archaea\n")
cat("Each domain: Absolute (top) + Percentage (bottom)\n")
cat("Single shared legend for all plots\n\n")

# Create output directories if they don't exist
if (!dir.exists("../figures")) {
  dir.create("../figures")
  cat("Created '../figures/' directory for plot outputs\n")
}

# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv"

# Load merged 16S data
cat("Loading 16S merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 16S merged data file at: ", merged_data_path)
}
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)

# Load census division data for .U. entries
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 16S census data file at: ", census_data_path)
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("16S merged data:", nrow(data_16s), "rows\n")
cat("16S census data:", nrow(census_division_data), "rows\n\n")

# Load shared color configuration
cat("Loading colors from shared configuration...\n")
color_config <- load_taxonomic_colors("../../../shared_config/taxonomic_color_mapping.yaml")

# Get bacterial colors from shared configuration
get_bacterial_colors <- function() {
  bacteria_list <- color_config$bacteria_colors
  bacteria_colors <- character(length(bacteria_list))
  names(bacteria_colors) <- names(bacteria_list)
  
  for (i in 1:length(bacteria_list)) {
    bacteria_colors[names(bacteria_list)[i]] <- bacteria_list[[i]]
  }
  
  return(bacteria_colors)
}

# Get archaeal colors from shared configuration
get_archaeal_colors <- function() {
  archaea_list <- color_config$archaea_colors
  archaea_colors <- character(length(archaea_list))
  names(archaea_colors) <- names(archaea_list)
  
  for (i in 1:length(archaea_list)) {
    archaea_colors[names(archaea_list)[i]] <- archaea_list[[i]]
  }
  
  return(archaea_colors)
}

bacterial_colors <- get_bacterial_colors()
archaeal_colors <- get_archaeal_colors()
cat("Loaded", length(bacterial_colors), "bacterial colors from shared configuration\n")
cat("Loaded", length(archaeal_colors), "archaeal colors from shared configuration\n\n")

# Extended color palette for unmapped phyla
extended_colors <- c(
  "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
  "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
  "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9"
)

# ============================================================================
# BACTERIA PROCESSING
# ============================================================================

cat("Processing BACTERIA data...\n")

# Process bacteria-specific data
process_bacteria_data <- function() {
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Bacteria")  # BACTERIA ONLY
  
  return(matched_data)
}

# Get bacteria .U. entries
get_bacteria_u_entries <- function(census_data) {
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 50) %>%
    filter(grepl("Bacteria", Name_to_use) |
           grepl("Proteobacteria", Name_to_use) |
           grepl("Firmicutes", Name_to_use) |
           grepl("Bacteroidetes", Name_to_use) |
           grepl("Actinobacteria", Name_to_use) |
           (!grepl("Eukaryota", Name_to_use) & !grepl("Archaea", Name_to_use))) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = "Bacteria",
      match_status = "census_only"
    )
  
  return(u_entries)
}

matched_bacteria <- process_bacteria_data()
u_bacteria <- get_bacteria_u_entries(census_division_data)

cat("Bacteria data processing results:\n")
cat("  - Matched bacteria phyla:", nrow(matched_bacteria), "\n")
cat("  - Bacteria .U. entries:", nrow(u_bacteria), "\n\n")

# Combine bacteria data
combined_bacteria <- bind_rows(matched_bacteria, u_bacteria)

# Calculate bacteria totals
total_bacteria_genome <- sum(matched_bacteria$ncbi_genome_count, na.rm = TRUE)
total_bacteria_species <- sum(matched_bacteria$ncbi_species_count, na.rm = TRUE)
total_bacteria_otu <- sum(combined_bacteria$census_otu_count, na.rm = TRUE)
total_bacteria_size <- sum(combined_bacteria$census_size_count, na.rm = TRUE)

cat("Bacteria totals:\n")
cat("  - Genomes:", scales::comma(total_bacteria_genome), "\n")
cat("  - Species:", scales::comma(total_bacteria_species), "\n")
cat("  - OTUs:", scales::comma(total_bacteria_otu), "\n")
cat("  - Sequences:", scales::comma(total_bacteria_size), "\n\n")

# Calculate bacteria percentages
combined_bacteria <- combined_bacteria %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_bacteria_genome) * 100,
    species_pct = (ncbi_species_count / total_bacteria_species) * 100,
    size_pct = (census_size_count / total_bacteria_size) * 100,
    otu_pct = (census_otu_count / total_bacteria_otu) * 100,
    total_representation = genome_pct + species_pct + size_pct + otu_pct
  )

# Select top 12 bacteria phyla
n_top_bacteria <- 12
top_bacteria <- combined_bacteria %>%
  arrange(desc(total_representation)) %>%
  head(n_top_bacteria)

cat("Top", n_top_bacteria, "bacterial phyla selected\n\n")

# ============================================================================
# ARCHAEA PROCESSING
# ============================================================================

cat("Processing ARCHAEA data...\n")

# Process archaea-specific data
process_archaea_data <- function() {
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Archaea")  # ARCHAEA ONLY

  return(matched_data)
}

# Get archaea .U. entries
get_archaea_u_entries <- function(census_data) {
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 10) %>%
    filter(grepl("Archaea", Name_to_use)) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = "Archaea",
      match_status = "census_only"
    )

  return(u_entries)
}

matched_archaea <- process_archaea_data()
u_archaea <- get_archaea_u_entries(census_division_data)

cat("Archaea data processing results:\n")
cat("  - Matched archaea phyla:", nrow(matched_archaea), "\n")
cat("  - Archaea .U. entries:", nrow(u_archaea), "\n\n")

# Combine archaea data
combined_archaea <- bind_rows(matched_archaea, u_archaea)

# Calculate archaea totals
total_archaea_genome <- sum(matched_archaea$ncbi_genome_count, na.rm = TRUE)
total_archaea_species <- sum(matched_archaea$ncbi_species_count, na.rm = TRUE)
total_archaea_otu <- sum(combined_archaea$census_otu_count, na.rm = TRUE)
total_archaea_size <- sum(combined_archaea$census_size_count, na.rm = TRUE)

cat("Archaea totals:\n")
cat("  - Genomes:", scales::comma(total_archaea_genome), "\n")
cat("  - Species:", scales::comma(total_archaea_species), "\n")
cat("  - OTUs:", scales::comma(total_archaea_otu), "\n")
cat("  - Sequences:", scales::comma(total_archaea_size), "\n\n")

# Calculate archaea percentages
combined_archaea <- combined_archaea %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_archaea_genome) * 100,
    species_pct = (ncbi_species_count / total_archaea_species) * 100,
    size_pct = (census_size_count / total_archaea_size) * 100,
    otu_pct = (census_otu_count / total_archaea_otu) * 100,
    total_representation = genome_pct + species_pct + size_pct + otu_pct
  )

# Select top 6-8 archaea phyla
n_top_archaea <- min(8, nrow(combined_archaea))
top_archaea <- combined_archaea %>%
  arrange(desc(total_representation)) %>%
  head(n_top_archaea)

cat("Top", n_top_archaea, "archaeal phyla selected\n\n")

# ============================================================================
# CREATE BACTERIA ALLUVIAL PLOT
# ============================================================================

cat("Creating bacteria alluvial plot...\n")

# Create bacteria long format data
bacteria_long <- data.frame()

for (i in 1:nrow(top_bacteria)) {
  bacteria_long <- rbind(bacteria_long, data.frame(
    phylum = top_bacteria$phylum[i],
    x = "Genbank_Genome_%",
    percentage = top_bacteria$genome_pct[i],
    alluvium = i,
    domain = "Bacteria",
    stringsAsFactors = FALSE
  ))
  bacteria_long <- rbind(bacteria_long, data.frame(
    phylum = top_bacteria$phylum[i],
    x = "IMG_Genome_%",
    percentage = top_bacteria$size_pct[i],
    alluvium = i,
    domain = "Bacteria",
    stringsAsFactors = FALSE
  ))
  bacteria_long <- rbind(bacteria_long, data.frame(
    phylum = top_bacteria$phylum[i],
    x = "16S_OTU_%",
    percentage = top_bacteria$otu_pct[i],
    alluvium = i,
    domain = "Bacteria",
    stringsAsFactors = FALSE
  ))
  bacteria_long <- rbind(bacteria_long, data.frame(
    phylum = top_bacteria$phylum[i],
    x = "Genbank_Species_%",
    percentage = top_bacteria$species_pct[i],
    alluvium = i,
    domain = "Bacteria",
    stringsAsFactors = FALSE
  ))
}

# Add "Other" for bacteria
other_bacteria_alluvium <- nrow(top_bacteria) + 1
other_bacteria_genome_pct <- max(0, 100 - sum(top_bacteria$genome_pct, na.rm = TRUE))
other_bacteria_size_pct <- max(0, 100 - sum(top_bacteria$size_pct, na.rm = TRUE))
other_bacteria_otu_pct <- max(0, 100 - sum(top_bacteria$otu_pct, na.rm = TRUE))
other_bacteria_species_pct <- max(0, 100 - sum(top_bacteria$species_pct, na.rm = TRUE))

bacteria_long <- rbind(bacteria_long, data.frame(
  phylum = "Other",
  x = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  percentage = c(other_bacteria_genome_pct, other_bacteria_size_pct, other_bacteria_otu_pct, other_bacteria_species_pct),
  alluvium = other_bacteria_alluvium,
  domain = "Bacteria",
  stringsAsFactors = FALSE
))

bacteria_long$x <- factor(bacteria_long$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))

# Assign colors to bacteria phyla
bacteria_phyla_names <- unique(bacteria_long$phylum)
bacteria_plot_colors <- character(length(bacteria_phyla_names))
names(bacteria_plot_colors) <- bacteria_phyla_names

for (i in 1:length(bacteria_phyla_names)) {
  phylum_name <- as.character(bacteria_phyla_names[i])

  if (phylum_name == "Other") {
    bacteria_plot_colors[i] <- "#CCCCCC"
  } else {
    clean_phylum <- gsub("\\.U\\..*$", "", phylum_name)
    if (clean_phylum %in% names(bacterial_colors)) {
      bacteria_plot_colors[i] <- bacterial_colors[clean_phylum]
    } else {
      fallback_index <- ((i - 1) %% length(extended_colors)) + 1
      bacteria_plot_colors[i] <- extended_colors[fallback_index]
    }
  }
}

# Prepare bacteria absolute values data
bacteria_long_abs <- bacteria_long %>%
  mutate(
    absolute_value = case_when(
      x == "Genbank_Genome_%" ~ percentage / 100 * total_bacteria_genome,
      x == "IMG_Genome_%" ~ percentage / 100 * total_bacteria_size,
      x == "16S_OTU_%" ~ percentage / 100 * total_bacteria_otu,
      x == "Genbank_Species_%" ~ percentage / 100 * total_bacteria_species,
      TRUE ~ 0
    )
  )

# Create bacteria ABSOLUTE plot
p_bacteria_abs <- ggplot(
  bacteria_long_abs,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_value, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = bacteria_plot_colors, name = "Phylum") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "16S_OTU_%" = "16S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(1.0, "cm"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  labs(
    title = "Absolute Values",
    y = "Count",
    x = NULL
  )

# Create bacteria PERCENTAGE plot
p_bacteria_pct <- ggplot(
  bacteria_long,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = bacteria_plot_colors, name = "Phylum") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "16S_OTU_%" = "16S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(1.0, "cm"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.margin = margin(5, 10, 10, 10)
  ) +
  labs(
    title = "Percentage Values",
    y = "Percentage (%)",
    x = NULL
  )

# ============================================================================
# CREATE ARCHAEA ALLUVIAL PLOT
# ============================================================================

cat("Creating archaea alluvial plot...\n")

# Create archaea long format data
archaea_long <- data.frame()

for (i in 1:nrow(top_archaea)) {
  archaea_long <- rbind(archaea_long, data.frame(
    phylum = top_archaea$phylum[i],
    x = "Genbank_Genome_%",
    percentage = top_archaea$genome_pct[i],
    alluvium = i,
    domain = "Archaea",
    stringsAsFactors = FALSE
  ))
  archaea_long <- rbind(archaea_long, data.frame(
    phylum = top_archaea$phylum[i],
    x = "IMG_Genome_%",
    percentage = top_archaea$size_pct[i],
    alluvium = i,
    domain = "Archaea",
    stringsAsFactors = FALSE
  ))
  archaea_long <- rbind(archaea_long, data.frame(
    phylum = top_archaea$phylum[i],
    x = "16S_OTU_%",
    percentage = top_archaea$otu_pct[i],
    alluvium = i,
    domain = "Archaea",
    stringsAsFactors = FALSE
  ))
  archaea_long <- rbind(archaea_long, data.frame(
    phylum = top_archaea$phylum[i],
    x = "Genbank_Species_%",
    percentage = top_archaea$species_pct[i],
    alluvium = i,
    domain = "Archaea",
    stringsAsFactors = FALSE
  ))
}

# Add "Other" for archaea
other_archaea_alluvium <- nrow(top_archaea) + 1
other_archaea_genome_pct <- max(0, 100 - sum(top_archaea$genome_pct, na.rm = TRUE))
other_archaea_size_pct <- max(0, 100 - sum(top_archaea$size_pct, na.rm = TRUE))
other_archaea_otu_pct <- max(0, 100 - sum(top_archaea$otu_pct, na.rm = TRUE))
other_archaea_species_pct <- max(0, 100 - sum(top_archaea$species_pct, na.rm = TRUE))

archaea_long <- rbind(archaea_long, data.frame(
  phylum = "Other",
  x = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  percentage = c(other_archaea_genome_pct, other_archaea_size_pct, other_archaea_otu_pct, other_archaea_species_pct),
  alluvium = other_archaea_alluvium,
  domain = "Archaea",
  stringsAsFactors = FALSE
))

archaea_long$x <- factor(archaea_long$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))

# Assign colors to archaea phyla
archaea_phyla_names <- unique(archaea_long$phylum)
archaea_plot_colors <- character(length(archaea_phyla_names))
names(archaea_plot_colors) <- archaea_phyla_names

for (i in 1:length(archaea_phyla_names)) {
  phylum_name <- as.character(archaea_phyla_names[i])

  if (phylum_name == "Other") {
    archaea_plot_colors[i] <- "#CCCCCC"
  } else {
    clean_phylum <- gsub("\\.U\\..*$", "", phylum_name)
    if (clean_phylum %in% names(archaeal_colors)) {
      archaea_plot_colors[i] <- archaeal_colors[clean_phylum]
    } else {
      fallback_index <- ((i - 1) %% length(extended_colors)) + 1
      archaea_plot_colors[i] <- extended_colors[fallback_index]
    }
  }
}

# Prepare archaea absolute values data
archaea_long_abs <- archaea_long %>%
  mutate(
    absolute_value = case_when(
      x == "Genbank_Genome_%" ~ percentage / 100 * total_archaea_genome,
      x == "IMG_Genome_%" ~ percentage / 100 * total_archaea_size,
      x == "16S_OTU_%" ~ percentage / 100 * total_archaea_otu,
      x == "Genbank_Species_%" ~ percentage / 100 * total_archaea_species,
      TRUE ~ 0
    )
  )

# Create archaea ABSOLUTE plot
p_archaea_abs <- ggplot(
  archaea_long_abs,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_value, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = archaea_plot_colors, name = "Phylum") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "16S_OTU_%" = "16S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(1.0, "cm"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  labs(
    title = "Absolute Values",
    y = "Count",
    x = NULL
  )

# Create archaea PERCENTAGE plot
p_archaea_pct <- ggplot(
  archaea_long,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = archaea_plot_colors, name = "Phylum") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "16S_OTU_%" = "16S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(1.0, "cm"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.margin = margin(5, 10, 10, 10)
  ) +
  labs(
    title = "Percentage Values",
    y = "Percentage (%)",
    x = NULL
  )

# ============================================================================
# STACK PLOTS WITH SHARED LEGEND
# ============================================================================

cat("\nStacking bacteria and archaea plots with shared legend...\n")

# Combine all phyla names and colors for shared legend
all_phyla <- c(bacteria_phyla_names, archaea_phyla_names)
all_colors <- c(bacteria_plot_colors, archaea_plot_colors)

# Remove duplicates (e.g., "Other" appears in both)
unique_phyla <- unique(all_phyla)
unique_colors <- character(length(unique_phyla))
names(unique_colors) <- unique_phyla

for (phylum in unique_phyla) {
  # Use bacteria color if available, otherwise archaea color
  if (phylum %in% names(bacteria_plot_colors)) {
    unique_colors[phylum] <- bacteria_plot_colors[phylum]
  } else {
    unique_colors[phylum] <- archaea_plot_colors[phylum]
  }
}

# Create SEPARATE stacked figures for Bacteria and Archaea
cat("\nCreating separate stacked figures for Bacteria and Archaea...\n")

# ============================================================================
# BACTERIA STACKED FIGURE (Absolute + Percentage)
# ============================================================================

# Update bacteria plots to use bacteria colors only
p_bacteria_abs_final <- p_bacteria_abs + scale_fill_manual(values = bacteria_plot_colors, name = "Phylum")
p_bacteria_pct_final <- p_bacteria_pct + scale_fill_manual(values = bacteria_plot_colors, name = "Phylum")

# Stack bacteria plots vertically with shared legend
p_bacteria_stacked <- p_bacteria_abs_final / p_bacteria_pct_final +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "16S Bacterial Diversity Flow: Absolute vs Percentage",
    theme = theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5))
  )

# Save bacteria stacked plot
cat("\nSaving bacteria stacked alluvial plot...\n")
ggsave("../figures/alluvial_16s_bacteria_stacked_abs_pct.png", p_bacteria_stacked, width = 24, height = 16, dpi = 300, bg = "white")
ggsave("../figures/alluvial_16s_bacteria_stacked_abs_pct.pdf", p_bacteria_stacked, width = 24, height = 16, dpi = 300, bg = "white")

# ============================================================================
# ARCHAEA STACKED FIGURE (Absolute + Percentage)
# ============================================================================

# Update archaea plots to use archaea colors only
p_archaea_abs_final <- p_archaea_abs + scale_fill_manual(values = archaea_plot_colors, name = "Phylum")
p_archaea_pct_final <- p_archaea_pct + scale_fill_manual(values = archaea_plot_colors, name = "Phylum")

# Stack archaea plots vertically with shared legend
p_archaea_stacked <- p_archaea_abs_final / p_archaea_pct_final +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "16S Archaeal Diversity Flow: Absolute vs Percentage",
    theme = theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5))
  )

# Save archaea stacked plot
cat("\nSaving archaea stacked alluvial plot...\n")
ggsave("../figures/alluvial_16s_archaea_stacked_abs_pct.png", p_archaea_stacked, width = 24, height = 16, dpi = 300, bg = "white")
ggsave("../figures/alluvial_16s_archaea_stacked_abs_pct.pdf", p_archaea_stacked, width = 24, height = 16, dpi = 300, bg = "white")

cat("\n=== 16S Prokaryotic STACKED Alluvial Plots Created Successfully ===\n")
cat("Files saved:\n")
cat("\nBACTERIA:\n")
cat("  - ../figures/alluvial_16s_bacteria_stacked_abs_pct.png\n")
cat("  - ../figures/alluvial_16s_bacteria_stacked_abs_pct.pdf\n")
cat("\nARCHAEA:\n")
cat("  - ../figures/alluvial_16s_archaea_stacked_abs_pct.png\n")
cat("  - ../figures/alluvial_16s_archaea_stacked_abs_pct.pdf\n")
cat("\nEach stacked figure contains:\n")
cat("  - Absolute values (top)\n")
cat("  - Percentage values (bottom)\n")
cat("  - Single shared legend\n")
cat("  - Professional color schemes\n")
cat("  - Advanced alluvial aesthetics\n\n")

