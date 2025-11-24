#!/usr/bin/env Rscript
#
# 16S Archaeal Alluvial Plot Generator - ABSOLUTE VALUES
# =====================================================
#
# Creates alluvial plots showing the flow from:
# Genbank Genome Count → IMG Genome Count → 16S OTU Count → Genbank Species Count
#
# Generates plots specifically for Archaea domain with proper archaea-specific:
# - Phyla organization and filtering
# - Percentage calculations based on archaea-only totals
# - Color schemes optimized for archaea phyla
# - Advanced alluvial aesthetics with optimized flow guidance
#
# Author: Archaea Alluvial Team
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
library(yaml)

cat("=== 16S rRNA Archaeal Alluvial Plot Generator (Absolute Values) ===\n\n")

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv"

# Load the 16S rRNA merged data
cat("Loading 16S rRNA merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 16S merged data file at: ", merged_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
cat("16S rRNA data loaded:", nrow(data_16s), "rows\n")

# Load 16S division data for .U. entries (archaeal domain)
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 16S census data file at: ", census_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)
cat("16S census data loaded:", nrow(census_division_data), "rows\n")

# Function to process archaea-specific data
process_archaea_data <- function() {
  cat("\n=== Processing Archaea Domain ===\n")

  # Filter for matched archaea phyla ONLY
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Archaea")  # ARCHAEA ONLY

  cat("Matched archaea phyla found:", nrow(matched_data), "\n")
  
  # Display archaea phyla found
  if (nrow(matched_data) > 0) {
    cat("Archaea phyla in dataset:\n")
    archaea_phyla <- unique(matched_data$phylum)
    for (i in 1:length(archaea_phyla)) {
      phylum_data <- matched_data[matched_data$phylum == archaea_phyla[i], ]
      total_genomes <- sum(phylum_data$ncbi_genome_count, na.rm = TRUE)
      total_species <- sum(phylum_data$ncbi_species_count, na.rm = TRUE)
      total_otus <- sum(phylum_data$census_otu_count, na.rm = TRUE)
      total_sequences <- sum(phylum_data$census_size_count, na.rm = TRUE)
      
      cat(sprintf("  %d. %-20s | Genomes: %6s | Species: %6s | OTUs: %8s | Sequences: %8s\n",
                  i, archaea_phyla[i],
                  format(total_genomes, big.mark = ","),
                  format(total_species, big.mark = ","),
                  format(total_otus, big.mark = ","),
                  format(total_sequences, big.mark = ",")))
    }
  }

  return(matched_data)
}

# Process archaea data
matched_data <- process_archaea_data()

# Extract archaea-specific .U. entries from 16S census data
get_archaea_u_entries <- function(census_data) {
  # Filter for archaea .U. entries
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 10) %>%  # Lower threshold for archaea (less abundant)
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

u_entries <- get_archaea_u_entries(census_division_data)
cat(paste("Archaea .U. entries found:", nrow(u_entries), "\n"))

# Combine matched archaea data with archaea .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals from ARCHAEA-ONLY data for proper percentages
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)

# For sequences and OTUs, use COMBINED archaea data (matched + .U. entries)
# This matches the approach used in 18S scripts
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)

cat(paste("Total counts (ARCHAEA-ONLY calculations):\n"))
cat(paste("  Archaea Genomes:", scales::comma(total_genome_count), "\n"))
cat("  Archaea 16S Sequences:", scales::comma(total_size_count), "\n")
cat("  Archaea 16S OTUs:", scales::comma(total_otu_count), "\n")
cat(paste("  Archaea Species:", scales::comma(total_species_count), "\n"))

# Calculate percentages for matched data (archaea-specific totals)
matched_data <- matched_data %>%
  mutate(
    otu_percentage = (census_otu_count / total_otu_count) * 100,
    size_percentage = (census_size_count / total_size_count) * 100,
    genome_pct_db = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100
  )

# Calculate percentages for .U. entries (archaea-specific totals)
if (nrow(u_entries) > 0) {
  u_entries <- u_entries %>%
    mutate(
      otu_percentage = (census_otu_count / total_otu_count) * 100,
      size_percentage = (census_size_count / total_size_count) * 100,
      genome_pct_db = 0,  # No genomes for .U. entries
      species_pct = 0     # No species for .U. entries
    )
}

# Select top archaea entries by total representation (archaea-appropriate number)
top_n <- 6  # Fewer archaea phyla than bacteria

top_phyla <- combined_data %>%
  mutate(
    total_representation = genome_pct_db + species_pct + otu_percentage + size_percentage
  ) %>%
  arrange(desc(total_representation)) %>%
  head(top_n)

cat("\n=== DEBUG: Top Archaea Phyla Selected ===\n")
for (i in 1:nrow(top_phyla)) {
  entry_type <- if (top_phyla$match_status[i] == "census_only") "[.U.]" else "[matched]"
  cat(sprintf("%2d. %-25s %6s | Genomes: %8s | Species: %8s | OTUs: %8s | Sequences: %8s\n",
              i,
              top_phyla$phylum[i],
              entry_type,
              format(top_phyla$ncbi_genome_count[i], big.mark = ","),
              format(top_phyla$ncbi_species_count[i], big.mark = ","),
              format(top_phyla$census_otu_count[i], big.mark = ","),
              format(top_phyla$census_size_count[i], big.mark = ",")))
}
cat("=== END DEBUG ===\n\n")

# Calculate "Other" category from remaining archaea data
other_data <- combined_data %>%
  filter(!phylum %in% top_phyla$phylum)

other_genome_count <- sum(other_data$ncbi_genome_count, na.rm = TRUE)
other_size_count <- sum(other_data$census_size_count, na.rm = TRUE)
other_otu_count <- sum(other_data$census_otu_count, na.rm = TRUE)
other_species_count <- sum(other_data$ncbi_species_count, na.rm = TRUE)

cat("\nOther category (remaining", nrow(other_data), "archaeal entries):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n\n")

# Create long format data for alluvial plot
long_data <- data.frame()

# Add top phyla data (handling .U. entries correctly)
for (i in 1:nrow(top_phyla)) {
  # Check if this is a .U. entry (census_only)
  is_u_entry <- top_phyla$match_status[i] == "census_only"

  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
    stratum = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
    absolute_count = c(
      if (is_u_entry) 0 else top_phyla$ncbi_genome_count[i],        # Node 1: 0 for .U. entries
      top_phyla$census_size_count[i],                               # Node 2: IMG Genome Count (16S sequences)
      top_phyla$census_otu_count[i],                                # Node 3: 16S OTU Count
      if (is_u_entry) 0 else top_phyla$ncbi_species_count[i]        # Node 4: 0 for .U. entries
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
  stratum = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"))
long_data$stratum <- factor(long_data$stratum, levels = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Fix aesthetic issues
cat("Applying advanced alluvial preprocessing...\n")

# Handle zero flows properly for ggalluvial
# Simple solution: Replace zero values with minimal visible values (1 count)
# This keeps the alluvium structure intact while making zeros barely visible

long_data_f <- long_data

# Replace zero values with minimal flow width (1 count)
long_data_f$absolute_count[long_data_f$absolute_count == 0] <- 1

cat("Advanced preprocessing complete - replaced zero values with 1-count minimal flows\n")

# Debug: Check the data structure and counts for Genbank nodes
cat("\n=== DEBUG: Genbank Node Data Check ===\n")
genbank_genome_data <- long_data[long_data$x == "Genbank_Genome_Count", ]
genbank_species_data <- long_data[long_data$x == "Genbank_Species_Count", ]

cat("Genbank Genome Node flows:\n")
for (i in 1:nrow(genbank_genome_data)) {
  if (genbank_genome_data$absolute_count[i] > 0) {
    cat(sprintf("  %-25s: %8s genomes\n",
                genbank_genome_data$phylum[i],
                format(genbank_genome_data$absolute_count[i], big.mark = ",")))
  }
}

cat("\nGenbank Species Node flows:\n")
for (i in 1:nrow(genbank_species_data)) {
  if (genbank_species_data$absolute_count[i] > 0) {
    cat(sprintf("  %-25s: %8s species\n",
                genbank_species_data$phylum[i],
                format(genbank_species_data$absolute_count[i], big.mark = ",")))
  }
}
cat("=== END DEBUG ===\n\n")

# Load shared color configuration
source("../../shared_config/color_mapping_functions.R")

# Load archaea colors from shared configuration
cat("Loading archaea colors from shared configuration...\n")
color_config <- load_taxonomic_colors("../../shared_config/taxonomic_color_mapping.yaml")

# Extract archaea colors and convert to named vector
archaea_colors_list <- color_config$archaea_colors
archaea_colors <- character(length(archaea_colors_list))
names(archaea_colors) <- names(archaea_colors_list)
for (name in names(archaea_colors_list)) {
  archaea_colors[name] <- as.character(archaea_colors_list[[name]])
}

# Extended fallback colors for archaea phyla not in main palette
get_extended_archaea_colors <- function() {
  c("#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222", "#FF4500", "#32CD32", "#9932CC")
}

# Generate colors for the plot using shared configuration
extended_colors <- get_extended_archaea_colors()

# Map colors to phyla in the plot
unique_phyla <- unique(long_data_f$phylum)
colors <- character(length(unique_phyla))
names(colors) <- unique_phyla

for (i in 1:length(unique_phyla)) {
  phylum_name <- unique_phyla[i]

  # Skip "Other" category for now
  if (phylum_name == "Other") {
    next
  }

  # Extract base phylum name (remove numbering)
  base_phylum <- gsub("^\\d+\\. ", "", phylum_name)

  # Check for .U. entries first and assign distinctive color
  if (grepl("\\.U\\.", base_phylum)) {
    colors[phylum_name] <- "#9ACD32"  # Yellow-green - distinctive for unclassified entries
  } else if (base_phylum %in% names(archaea_colors)) {
    # Try to match with archaea color palette
    colors[phylum_name] <- archaea_colors[base_phylum]
  } else {
    # Use extended colors for unmapped phyla
    color_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[phylum_name] <- extended_colors[color_index]
  }
}

# Set "Other" category to light gray
colors["Other"] <- "#CCCCCC"

cat("Color mapping for archaea phyla:\n")
for (i in 1:length(colors)) {
  cat(sprintf("  %-25s: %s\n", names(colors)[i], colors[i]))
}
cat("\n")

# Prepare node annotations data
max_count <- max(c(total_genome_count, total_size_count, total_otu_count, total_species_count))
node_labels <- data.frame(
  x = factor(c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
             levels = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count")),
  label = c(
    paste0("Genbank Genomes\n", scales::comma(total_genome_count)),
    paste0("IMG Sequences\n", scales::comma(total_size_count)),
    paste0("16S OTUs\n", scales::comma(total_otu_count)),
    paste0("Genbank Species\n", scales::comma(total_species_count))
  ),
  y = max_count * 1.05,  # Position above the plot
  stringsAsFactors = FALSE
)

# Create ADVANCED absolute alluvial plot with optimized aesthetics
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
  # Add node labels with counts
  geom_text(data = node_labels, aes(x = x, y = y, label = label, fill = NULL, stratum = NULL, alluvium = NULL),
            inherit.aes = FALSE, size = 6, fontface = "bold", hjust = 0.5, vjust = 0) +
  scale_fill_manual(values = colors, name = "Archaea Division") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    limits = c(0, max_count * 1.1),  # Increased to accommodate node labels
    expand = expansion(mult = c(0, 0))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # No node titles - clean minimal appearance
    axis.text.x = element_blank(),
    # Y-axis ticks - larger and bold
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
    axis.title.x = element_blank(),
    # Legend positioning and styling
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(l = 20),
    # Plot margins and background
    plot.margin = margin(20, 20, 20, 20),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Title styling
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15))
  ) +
  labs(
    title = "16S Archaeal Taxonomic Coverage Flow (Absolute Values)",
    subtitle = paste("Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species | Top", nrow(top_phyla), "Archaea Phyla"),
    y = "Absolute Count"
  )

# Save the plot
output_png <- "alluvial_16s_archaea_abs_values_only.png"
output_pdf <- "alluvial_16s_archaea_abs_values_only.pdf"

cat("Saving archaea alluvial plot...\n")
ggsave(output_png, p_abs, width = 24, height = 10, dpi = 300, bg = "white")
ggsave(output_pdf, p_abs, width = 24, height = 10, dpi = 300, bg = "white")

# Check which phyla have NCBI data for summary
ncbi_visible <- top_phyla %>%
  filter(ncbi_genome_count > 0 | ncbi_species_count > 0)

if (nrow(ncbi_visible) == 0) {
  cat("\n⚠️  WARNING: No archaea phyla have NCBI genome or species data!\n")
  cat("   This suggests archaea are severely underrepresented in NCBI databases.\n")
  cat("   The plot will show flows only through 16S census data (sequences/OTUs).\n\n")
} else {
  cat("\n✅ Archaea phyla with NCBI representation:\n")
  for (i in 1:nrow(ncbi_visible)) {
    cat(sprintf("   - %-20s: %s genomes, %s species\n",
                ncbi_visible$phylum[i],
                format(ncbi_visible$ncbi_genome_count[i], big.mark = ","),
                format(ncbi_visible$ncbi_species_count[i], big.mark = ",")))
  }
}

cat("\n📊 Summary (Archaea NCBI data):\n")
cat("   - Total archaea phyla:", nrow(matched_data), "\n")
cat("   - Top entries displayed:", nrow(top_phyla), "\n")
cat("   - Phyla with NCBI genomes:", nrow(ncbi_visible), "\n")
cat("   - Total Archaea Genomes:", scales::comma(total_genome_count), "\n")
cat("   - Total 16S Sequences:", scales::comma(total_size_count), "\n")
cat("   - Total 16S OTUs:", scales::comma(total_otu_count), "\n")
cat("   - Total Archaea Species:", scales::comma(total_species_count), "\n")

cat("\n⚠️  If you only see 1-2 phyla in the plot, this is expected for archaea:\n")
cat("   1. Archaea are less abundant than bacteria in most environments\n")
cat("   2. Archaea are underrepresented in NCBI genome databases\n")
cat("   3. Most archaea diversity is captured in 16S census data, not genomes\n\n")

# Create detailed flow annotations with absolute values and percentages
cat("Creating detailed flow annotations file...\n")

# Calculate totals for percentage calculations
total_genomes <- sum(combined_data$ncbi_genome_count, na.rm = TRUE)
total_sequences <- sum(combined_data$census_size_count, na.rm = TRUE)
total_otus <- sum(combined_data$census_otu_count, na.rm = TRUE)
total_species <- sum(combined_data$ncbi_species_count, na.rm = TRUE)

# Create detailed annotations for each taxon at each node
flow_annotations <- data.frame()

# Add annotations for top phyla
for (i in 1:nrow(top_phyla)) {
  taxon_name <- top_phyla$phylum[i]

  # Node 1: Genbank Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Genome_Count",
    Node_Order = 1,
    Absolute_Count = top_phyla$ncbi_genome_count[i],
    Percentage = round((top_phyla$ncbi_genome_count[i] / total_genomes) * 100, 2),
    Flow_Width = top_phyla$ncbi_genome_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 2: IMG Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "IMG_Genome_Count",
    Node_Order = 2,
    Absolute_Count = top_phyla$census_size_count[i],
    Percentage = round((top_phyla$census_size_count[i] / total_sequences) * 100, 2),
    Flow_Width = top_phyla$census_size_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 3: 16S OTUs
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "16S_OTU_Count",
    Node_Order = 3,
    Absolute_Count = top_phyla$census_otu_count[i],
    Percentage = round((top_phyla$census_otu_count[i] / total_otus) * 100, 2),
    Flow_Width = top_phyla$census_otu_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 4: Genbank Species
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Species_Count",
    Node_Order = 4,
    Absolute_Count = top_phyla$ncbi_species_count[i],
    Percentage = round((top_phyla$ncbi_species_count[i] / total_species) * 100, 2),
    Flow_Width = top_phyla$ncbi_species_count[i],
    stringsAsFactors = FALSE
  ))
}

# Add "Other" category annotations
flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "Genbank_Genome_Count",
  Node_Order = 1,
  Absolute_Count = other_genome_count,
  Percentage = round((other_genome_count / total_genomes) * 100, 2),
  Flow_Width = other_genome_count,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "IMG_Genome_Count",
  Node_Order = 2,
  Absolute_Count = other_size_count,
  Percentage = round((other_size_count / total_sequences) * 100, 2),
  Flow_Width = other_size_count,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "16S_OTU_Count",
  Node_Order = 3,
  Absolute_Count = other_otu_count,
  Percentage = round((other_otu_count / total_otus) * 100, 2),
  Flow_Width = other_otu_count,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "Genbank_Species_Count",
  Node_Order = 4,
  Absolute_Count = other_species_count,
  Percentage = round((other_species_count / total_species) * 100, 2),
  Flow_Width = other_species_count,
  stringsAsFactors = FALSE
))

# Sort by node order and then by flow width (descending)
flow_annotations <- flow_annotations %>%
  arrange(Node_Order, desc(Flow_Width))

write.table(flow_annotations, "alluvial_16s_archaea_abs_flow_annotations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_Count", "IMG_Genome_Count", "16S_OTU_Count", "Genbank_Species_Count"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    "Genbank Total Genomes (Archaea)",
    "IMG Genome Count (16S sequences)",
    "16S OTU Count",
    "Genbank Total Species (Archaea)"
  ),
  Total_Count = c(
    scales::comma(total_genomes),
    scales::comma(total_sequences),
    scales::comma(total_otus),
    scales::comma(total_species)
  ),
  Data_Type = c("Absolute Count", "Absolute Count", "Absolute Count", "Absolute Count"),
  stringsAsFactors = FALSE
)

write.table(node_descriptions, "alluvial_16s_archaea_abs_node_descriptions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary statistics file
summary_stats <- data.frame(
  Metric = c(
    "Total_Taxa_Shown",
    "Top_Taxa_Count",
    "Other_Category_Included",
    "Total_Archaea_Genomes",
    "Total_16S_Sequences",
    "Total_16S_OTUs",
    "Total_Archaea_Species",
    "Filtering_Method",
    "Color_System"
  ),
  Value = c(
    nrow(top_phyla) + 1,  # +1 for Other category
    nrow(top_phyla),
    "Yes",
    scales::comma(total_genomes),
    scales::comma(total_sequences),
    scales::comma(total_otus),
    scales::comma(total_species),
    paste("Top", nrow(top_phyla), "by total representation"),
    "Archaea color mapping"
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, "alluvial_16s_archaea_abs_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 16S Archaea Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_16s_archaea_abs_values_only.png\n")
cat("  - alluvial_16s_archaea_abs_values_only.pdf\n")
cat("  - alluvial_16s_archaea_abs_flow_annotations.tsv\n")
cat("  - alluvial_16s_archaea_abs_node_descriptions.tsv\n")
cat("  - alluvial_16s_archaea_abs_summary.tsv\n")

cat("\n16S Archaea alluvial plot generated with:\n")
cat("  - Archaea-specific data filtering and totals\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional archaea color scheme\n")
cat("  - Detailed flow annotations (TSV format)\n")
cat("  - Node descriptions and summary statistics\n")
cat("  - Archaea-only percentage calculations\n")
