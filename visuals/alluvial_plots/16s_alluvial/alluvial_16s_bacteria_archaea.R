#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - BACTERIA & ARCHAEA
# ============================================================
#
# Creates alluvial plots showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#
# Generates separate plots for both Bacteria and Archaea domains
# Shows top phyla by total representation for each domain
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

cat("=== 16S rRNA Prokaryotic Alluvial Plot Generator (Bacteria & Archaea) ===\n\n")

# Load the 16S rRNA merged data
cat("Loading 16S rRNA merged data...\n")
data_16s <- read.csv("../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("16S rRNA data loaded:", nrow(data_16s), "rows\n")

# Load 16S division data for .U. entries
census_division_data <- read.csv("../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv", stringsAsFactors = FALSE)

# Master bacterial color palette from mega 16S script
get_bacterial_colors <- function() {
  bacteria_color_map <- c(
    "Abditibacteriota" = "#87bd4b",        # 36.71× novel - highest factor
    "Calditrichota" = "#e4a98f",           # 35.00× novel
    "Planctomycetota" = "#923052",         # 26.34× novel - unique cell wall-less
    "Thermomicrobiota" = "#6cd490",        # 25.55× novel
    "Acidobacteriota" = "#985739",         # 24.66× novel - soil bacteria
    "Gemmatimonadota" = "#9d9719",         # 23.33× novel
    "Chloroflexota" = "#b2f3fd",           # 18.45× novel - photosynthetic
    "Verrucomicrobiota" = "#8d0571",       # 16.86× novel - purple for distinction
    "Armatimonadota" = "#e99953",          # 15.00× novel
    "Myxococcota" = "#fae692",             # 14.29× novel - social bacteria
    "Bdellovibrionota" = "#164c3f",        # 13.33× novel - predatory bacteria
    "Thermodesulfobacteriota" = "#01414d", # 12.50× novel - thermophilic sulfate reducers
    "Spirochaetota" = "#dcb94b",           # 11.11× novel - spiral-shaped
    "Bacteroidota" = "#de724d",            # 4.25× novel - major gut bacteria
    "Chlorobiota" = "#458f9c",             # 18.67× overrepresented
    "Chlamydiota" = "#5b1f06",             # 9.50× overrepresented - dark brown for obligate parasites
    "Ignavibacteriota" = "#f16127"         # 1063.64× overrepresented - blood orange for distinction
  )
  return(bacteria_color_map)
}

# Master archaea color palette from mega 16S script
get_archaea_colors <- function() {
  archaea_color_map <- c(
    "Euryarchaeota" = "#55ee79",        # Most common archaea phylum - bright green
    "Nitrososphaerota" = "#ef9e17",     # Important ammonia-oxidizing archaea
    "Thermoproteota" = "#f61e5d",       # Hyperthermophilic archaea
    "Nanoarchaeota" = "#4ad9d5"         # Highest novelty factor (117.7×)
  )
  return(archaea_color_map)
}

# Extended fallback colors
get_extended_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

# Excluded bacterial phyla that use fallback colors
get_excluded_bacterial_phyla <- function() {
  c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
}

# Extract domain-specific .U. entries from 16S census data
get_domain_u_entries <- function(census_data, domain_name) {
  if (domain_name == "Bacteria") {
    # Filter for bacterial .U. entries
    u_entries <- census_data %>%
      filter(grepl("\\.U\\.", Name_to_use)) %>%
      filter(otu_count >= 50) %>%  # Minimum threshold for inclusion
      filter(grepl("Bacteria", Name_to_use) | 
             grepl("Proteobacteria", Name_to_use) |
             grepl("Firmicutes", Name_to_use) |
             grepl("Bacteroidetes", Name_to_use) |
             grepl("Actinobacteria", Name_to_use) |
             (!grepl("Eukaryota", Name_to_use) & !grepl("Archaea", Name_to_use)))
    
    if (nrow(u_entries) > 0) {
      u_data <- data.frame(
        phylum = "Bacteria.U.phylum",
        census_otu_count = sum(u_entries$otu_count, na.rm = TRUE),
        census_size_count = sum(u_entries$size_count, na.rm = TRUE),
        ncbi_genome_count = 0,
        ncbi_species_count = 0,
        domain = "Bacteria",
        stringsAsFactors = FALSE
      )
      return(u_data)
    }
  } else if (domain_name == "Archaea") {
    # Filter for archaea .U. entries
    u_entries <- census_data %>%
      filter(grepl("\\.U\\.", Name_to_use)) %>%
      filter(otu_count >= 10) %>%  # Lower threshold for archaea (less abundant)
      filter(grepl("Archaea", Name_to_use))
    
    if (nrow(u_entries) > 0) {
      u_data <- data.frame(
        phylum = "Archaea.U.phylum",
        census_otu_count = sum(u_entries$otu_count, na.rm = TRUE),
        census_size_count = sum(u_entries$size_count, na.rm = TRUE),
        ncbi_genome_count = 0,
        ncbi_species_count = 0,
        domain = "Archaea",
        stringsAsFactors = FALSE
      )
      return(u_data)
    }
  }
  
  return(data.frame())
}

# Function to process domain-specific data and create alluvial plot
process_domain_alluvial <- function(domain_name) {
  cat(paste("\n=== Processing", domain_name, "Domain ===\n"))
  
  # Filter for matched phyla for specific domain
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == domain_name)
  
  cat("Matched", tolower(domain_name), "phyla found:", nrow(matched_data), "\n")
  
  # Get .U. entries for this domain
  u_entries <- get_domain_u_entries(census_division_data, domain_name)
  if (nrow(u_entries) > 0) {
    cat(paste(domain_name, ".U. entries found:", nrow(u_entries), "\n"))
    matched_data <- rbind(matched_data, u_entries)
  }
  
  # Determine number of top phyla to show
  n_top <- if (domain_name == "Bacteria") 12 else min(8, nrow(matched_data))
  
  # Calculate total representation score and select top phyla
  matched_data$total_representation <- matched_data$ncbi_genome_count + 
                                      matched_data$census_otu_count + 
                                      matched_data$census_size_count + 
                                      matched_data$ncbi_species_count
  
  # Sort by total representation and take top N
  top_phyla <- matched_data %>%
    arrange(desc(total_representation)) %>%
    head(n_top)
  
  cat("Top", n_top, tolower(domain_name), "phyla selected\n")
  
  # Calculate "Other" category
  other_data <- matched_data %>%
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
  
  return(list(
    top_phyla = top_phyla,
    other_counts = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
    domain = domain_name
  ))
}

# Function to create alluvial plot for a domain
create_domain_alluvial <- function(domain_data) {
  domain_name <- domain_data$domain
  top_phyla <- domain_data$top_phyla
  other_counts <- domain_data$other_counts

  cat(paste("\n=== Creating", domain_name, "Alluvial Plot ===\n"))

  # Create long format data for alluvial plot
  long_data <- data.frame()

  # Add top phyla data
  for (i in 1:nrow(top_phyla)) {
    phylum_data <- data.frame(
      alluvium = rep(i, 4),
      phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
      x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
      stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
      absolute_count = c(
        top_phyla$ncbi_genome_count[i],
        top_phyla$census_size_count[i],
        top_phyla$census_otu_count[i],
        top_phyla$ncbi_species_count[i]
      ),
      stringsAsFactors = FALSE
    )
    long_data <- rbind(long_data, phylum_data)
  }

  # Add "Other" category
  other_data <- data.frame(
    alluvium = rep(nrow(top_phyla) + 1, 4),
    phylum = rep("Other", 4),
    x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    absolute_count = other_counts,
    stringsAsFactors = FALSE
  )

  long_data <- rbind(long_data, other_data)

  # Fix x-axis ordering
  long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))
  long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))

  # ADVANCED ALLUVIAL PREPROCESSING
  long_data_f <- long_data %>%
    complete(x, phylum, fill = list(absolute_count = 0)) %>%
    group_by(phylum) %>%
    mutate(alluvium = cur_group_id()) %>%
    ungroup()

  # Order phyla by size at first axis
  first_axis <- sort(unique(long_data_f$x))[1]
  sizes_first <- long_data_f %>%
    filter(x == first_axis) %>%
    arrange(desc(absolute_count)) %>%
    select(phylum) %>% pull()

  long_data_f <- long_data_f %>%
    mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

  # Get colors for this domain
  colors <- get_domain_colors(unique(long_data_f$phylum), domain_name)

  # Create plot
  p <- ggplot(
    long_data_f,
    aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_count, fill = phylum)
  ) +
    geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.15,
                  lode.guidance = "forward", knot.pos = 0.35) +
    geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
                 linewidth = 0.4, width = 0.12) +
    scale_fill_manual(values = colors, name = "Phylum") +
    scale_x_discrete(expand = c(0.02, 0.02)) +
    scale_y_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = expansion(mult = c(0, 0.02))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 24, face = "bold", angle = 0, hjust = 0.5, margin = margin(t = 20)),
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 24, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 16, face = "bold"),
      legend.key.size = unit(1.8, "cm"),
      legend.margin = margin(l = 80),
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

  return(list(plot = p, data = long_data_f, domain = domain_name))
}

# Function to get domain-specific colors
get_domain_colors <- function(phyla_names, domain_name) {
  n_colors <- length(phyla_names)
  colors <- character(n_colors)
  names(colors) <- phyla_names

  if (domain_name == "Bacteria") {
    bacterial_colors <- get_bacterial_colors()
    extended_colors <- get_extended_colors()
    excluded_phyla <- get_excluded_bacterial_phyla()

    for (i in seq_along(phyla_names)) {
      phylum_name <- phyla_names[i]
      clean_phylum <- gsub("^\\d+\\. ", "", phylum_name)

      if (phylum_name == "Other") {
        colors[i] <- "#808080"
      } else if (grepl("\\.U\\.", clean_phylum)) {
        colors[i] <- "#ff9999"
      } else if (clean_phylum %in% excluded_phyla) {
        fallback_index <- ((i - 1) %% length(extended_colors)) + 1
        colors[i] <- extended_colors[fallback_index]
      } else if (clean_phylum %in% names(bacterial_colors)) {
        colors[i] <- bacterial_colors[clean_phylum]
      } else {
        fallback_index <- ((i - 1) %% length(extended_colors)) + 1
        colors[i] <- extended_colors[fallback_index]
      }
    }
  } else if (domain_name == "Archaea") {
    archaea_colors <- get_archaea_colors()
    extended_archaea_colors <- c("#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222", "#FF4500", "#32CD32", "#9932CC")

    for (i in seq_along(phyla_names)) {
      phylum_name <- phyla_names[i]
      clean_phylum <- gsub("^\\d+\\. ", "", phylum_name)

      if (phylum_name == "Other") {
        colors[i] <- "#808080"
      } else if (grepl("\\.U\\.", clean_phylum)) {
        colors[i] <- "#ffcc99"
      } else if (clean_phylum %in% names(archaea_colors)) {
        colors[i] <- archaea_colors[clean_phylum]
      } else {
        fallback_index <- ((i - 1) %% length(extended_archaea_colors)) + 1
        colors[i] <- extended_archaea_colors[fallback_index]
      }
    }
  }

  return(colors)
}

# Process both domains
bacteria_data <- process_domain_alluvial("Bacteria")
archaea_data <- process_domain_alluvial("Archaea")

# Create alluvial plots
bacteria_plot <- create_domain_alluvial(bacteria_data)
archaea_plot <- create_domain_alluvial(archaea_data)

# Save plots
ggsave("alluvial_16s_bacteria_abs_values.png", bacteria_plot$plot, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("alluvial_16s_bacteria_abs_values.pdf", bacteria_plot$plot, width = 24, height = 10, dpi = 300, bg = "white")

ggsave("alluvial_16s_archaea_abs_values.png", archaea_plot$plot, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("alluvial_16s_archaea_abs_values.pdf", archaea_plot$plot, width = 24, height = 10, dpi = 300, bg = "white")

cat("\n=== Alluvial Plots Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_16s_bacteria_abs_values.png/pdf\n")
cat("  - alluvial_16s_archaea_abs_values.png/pdf\n")
cat("\nBoth bacteria and archaea alluvial plots generated with optimized aesthetics!\n")
