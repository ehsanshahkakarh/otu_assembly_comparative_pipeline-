#!/usr/bin/env Rscript

# ============================================================================
# 16S Archaeal Alluvial Plot - PERCENTAGE VALUES (Archaea-Specific Approach)
# ============================================================================
# This script creates percentage-based alluvial plots for 16S archaeal data
# using ONLY archaea-specific data filtering and percentage calculations.
# 
# Data Flow: Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species
# 
# Key Features:
# - Uses only archaea domain data for percentage calculations
# - Advanced alluvial preprocessing for clean flows
# - Professional archaea color schemes
# - Optimized aesthetics (thin nodes, elegant flows)
# - Percentage normalization based on archaea-only totals
# - Archaea-specific phyla organization and filtering
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)

cat("=== 16S Archaea Alluvial Plot (Percentage Values) ===\n")
cat("Using archaea-specific data approach for reliable visualization\n\n")

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv"

# Load merged 16S data (ONLY source needed)
cat("Loading 16S merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 16S merged data file at: ", merged_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)

# Load census division data for .U. entries
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 16S census data file at: ", census_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("16S merged data:", nrow(data_16s), "rows\n")
cat("16S census data:", nrow(census_division_data), "rows\n\n")

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
  # Filter for archaea .U. entries with specific patterns
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

# Process the data
matched_data <- process_archaea_data()
u_entries <- get_archaea_u_entries(census_division_data)

cat("Archaea data processing results:\n")
cat("  - Matched archaea phyla:", nrow(matched_data), "\n")
cat("  - Archaea .U. entries:", nrow(u_entries), "\n")

# Display archaea phyla found
if (nrow(matched_data) > 0) {
  cat("\nArchaea phyla in dataset:\n")
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

# Combine matched data with .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals for percentage calculations (ARCHAEA-ONLY)
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)

# For sequences and OTUs, use COMBINED archaea data (matched + .U. entries)
# This matches the approach used in 18S scripts
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)

cat("\nTotal counts for percentage calculations (ARCHAEA-ONLY):\n")
cat("  - Total Archaea Genomes:", scales::comma(total_genome_count), "\n")
cat("  - Total Archaea Species:", scales::comma(total_species_count), "\n")
cat("  - Total 16S OTUs:", scales::comma(total_otu_count), "\n")
cat("  - Total 16S Sequences:", scales::comma(total_size_count), "\n\n")

# Select top entries by total representation (archaea-appropriate number)
top_n <- 6  # Fewer archaea phyla than bacteria

top_phyla <- combined_data %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    otu_pct = (census_otu_count / total_otu_count) * 100,
    size_pct = (census_size_count / total_size_count) * 100,
    total_representation = genome_pct + species_pct + otu_pct + size_pct
  ) %>%
  arrange(desc(total_representation)) %>%
  head(top_n)

cat("Top archaea phyla selected by total representation:\n")
for (i in 1:nrow(top_phyla)) {
  entry_type <- if (top_phyla$match_status[i] == "census_only") "[.U.]" else "[matched]"
  cat(sprintf("  %d. %-20s %6s | Total rep: %6.2f%%\n",
              i, top_phyla$phylum[i], entry_type, top_phyla$total_representation[i]))
}

# Calculate "Other" category percentages
other_data <- combined_data %>%
  filter(!phylum %in% top_phyla$phylum)

other_genome_count <- sum(other_data$ncbi_genome_count, na.rm = TRUE)
other_species_count <- sum(other_data$ncbi_species_count, na.rm = TRUE)
other_otu_count <- sum(other_data$census_otu_count, na.rm = TRUE)
other_size_count <- sum(other_data$census_size_count, na.rm = TRUE)

other_genome_pct <- (other_genome_count / total_genome_count) * 100
other_species_pct <- (other_species_count / total_species_count) * 100
other_otu_pct <- (other_otu_count / total_otu_count) * 100
other_size_pct <- (other_size_count / total_size_count) * 100

cat(sprintf("\nOther category: %.2f%% total representation\n", 
            other_genome_pct + other_species_pct + other_otu_pct + other_size_pct))

# Create long format data for alluvial plot
long_data <- data.frame()

# Add top phyla data (handling .U. entries correctly)
for (i in 1:nrow(top_phyla)) {
  # Check if this is a .U. entry (census_only)
  is_u_entry <- top_phyla$match_status[i] == "census_only"

  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
    stratum = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
    percentage = c(
      if (is_u_entry) 0 else top_phyla$genome_pct[i],    # Node 1: 0% for .U. entries
      top_phyla$size_pct[i],                             # Node 2: IMG Genome % (16S sequences)
      top_phyla$otu_pct[i],                              # Node 3: 16S OTU %
      if (is_u_entry) 0 else top_phyla$species_pct[i]    # Node 4: 0% for .U. entries
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  stratum = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)

long_data <- rbind(long_data, other_data)

# Fix x-axis ordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))
long_data$stratum <- factor(long_data$stratum, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Fix aesthetic issues
cat("Applying advanced alluvial preprocessing...\n")

# Handle zero flows properly for ggalluvial
# Simple solution: Replace zero values with minimal visible values (0.1)
# This keeps the alluvium structure intact while making zeros barely visible

long_data_f <- long_data

# Replace zero values with minimal flow width (0.1%)
long_data_f$percentage[long_data_f$percentage == 0] <- 0.1

cat("Advanced preprocessing complete - replaced zero values with 0.1% minimal flows\n")

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

# Extended color palette for unmapped phyla
get_extended_colors <- function() {
  extended_colors <- c(
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5"
  )
  return(extended_colors)
}

# Generate colors for the plot using shared configuration
extended_colors <- get_extended_colors()

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
node_labels <- data.frame(
  x = factor(c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
             levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%")),
  label = c(
    paste0("Genbank Genomes\n", scales::comma(total_genome_count)),
    paste0("IMG Sequences\n", scales::comma(total_size_count)),
    paste0("16S OTUs\n", scales::comma(total_otu_count)),
    paste0("Genbank Species\n", scales::comma(total_species_count))
  ),
  y = 102,  # Position above the plot
  stringsAsFactors = FALSE
)

# Create ADVANCED percentage alluvial plot with optimized aesthetics
p_pct <- ggplot(
  long_data_f,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)
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
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, 110),  # Increased to accommodate node labels
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
    title = "16S Archaeal Taxonomic Coverage Flow (Percentage Values)",
    subtitle = paste("Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species | Top", nrow(top_phyla), "Archaea Phyla"),
    y = "Percentage (%)"
  )

# Save the plot
output_png <- "alluvial_16s_archaea_pct_values_only.png"
output_pdf <- "alluvial_16s_archaea_pct_values_only.pdf"

cat("Saving archaea percentage alluvial plot...\n")
ggsave(output_png, p_pct, width = 24, height = 10, dpi = 300, bg = "white")
ggsave(output_pdf, p_pct, width = 24, height = 10, dpi = 300, bg = "white")

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
    cat(sprintf("   - %-20s: %.2f%% genomes, %.2f%% species\n",
                ncbi_visible$phylum[i],
                ncbi_visible$genome_pct[i],
                ncbi_visible$species_pct[i]))
  }
}

cat("\n📊 Summary (Archaea percentage data):\n")
cat("   - Total archaea phyla:", nrow(matched_data), "\n")
cat("   - Top entries displayed:", nrow(top_phyla), "\n")
cat("   - Phyla with NCBI genomes:", nrow(ncbi_visible), "\n")
cat("   - Percentage calculations based on archaea-only totals\n")

# Create detailed flow annotations with percentages and absolute values
cat("\nCreating detailed flow annotations file...\n")

# Create detailed annotations for each taxon at each node
flow_annotations <- data.frame()

# Add annotations for top phyla
for (i in 1:nrow(top_phyla)) {
  taxon_name <- top_phyla$phylum[i]

  # Node 1: Genbank Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Genome_%",
    Node_Order = 1,
    Absolute_Count = top_phyla$ncbi_genome_count[i],
    Percentage = round(top_phyla$genome_pct[i], 2),
    Flow_Width = top_phyla$genome_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 2: IMG Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "IMG_Genome_%",
    Node_Order = 2,
    Absolute_Count = top_phyla$census_size_count[i],
    Percentage = round(top_phyla$size_pct[i], 2),
    Flow_Width = top_phyla$size_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 3: 16S OTUs
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "16S_OTU_%",
    Node_Order = 3,
    Absolute_Count = top_phyla$census_otu_count[i],
    Percentage = round(top_phyla$otu_pct[i], 2),
    Flow_Width = top_phyla$otu_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 4: Genbank Species
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Species_%",
    Node_Order = 4,
    Absolute_Count = top_phyla$ncbi_species_count[i],
    Percentage = round(top_phyla$species_pct[i], 2),
    Flow_Width = top_phyla$species_pct[i],
    stringsAsFactors = FALSE
  ))
}

# Add "Other" category annotations
flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "Genbank_Genome_%",
  Node_Order = 1,
  Absolute_Count = other_genome_count,
  Percentage = round(other_genome_pct, 2),
  Flow_Width = other_genome_pct,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "IMG_Genome_%",
  Node_Order = 2,
  Absolute_Count = other_size_count,
  Percentage = round(other_size_pct, 2),
  Flow_Width = other_size_pct,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "16S_OTU_%",
  Node_Order = 3,
  Absolute_Count = other_otu_count,
  Percentage = round(other_otu_pct, 2),
  Flow_Width = other_otu_pct,
  stringsAsFactors = FALSE
))

flow_annotations <- rbind(flow_annotations, data.frame(
  Taxon = "Other",
  Node = "Genbank_Species_%",
  Node_Order = 4,
  Absolute_Count = other_species_count,
  Percentage = round(other_species_pct, 2),
  Flow_Width = other_species_pct,
  stringsAsFactors = FALSE
))

# Sort by node order and then by flow width (descending)
flow_annotations <- flow_annotations %>%
  arrange(Node_Order, desc(Flow_Width))

write.table(flow_annotations, "alluvial_16s_archaea_pct_flow_annotations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    "Genbank Total Genomes (Archaea)",
    "IMG Genome Count (16S sequences)",
    "16S OTU Count",
    "Genbank Total Species (Archaea)"
  ),
  Total_Count = c(
    scales::comma(total_genome_count),
    scales::comma(total_size_count),
    scales::comma(total_otu_count),
    scales::comma(total_species_count)
  ),
  Data_Type = c("Percentage", "Percentage", "Percentage", "Percentage"),
  stringsAsFactors = FALSE
)

write.table(node_descriptions, "alluvial_16s_archaea_pct_node_descriptions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

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
    "Color_System",
    "Percentage_Base"
  ),
  Value = c(
    nrow(top_phyla) + 1,  # +1 for Other category
    nrow(top_phyla),
    "Yes",
    scales::comma(total_genome_count),
    scales::comma(total_size_count),
    scales::comma(total_otu_count),
    scales::comma(total_species_count),
    paste("Top", nrow(top_phyla), "by total representation"),
    "Archaea color mapping",
    "Archaea-only totals"
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, "alluvial_16s_archaea_pct_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 16S Archaea Percentage Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_16s_archaea_pct_values_only.png\n")
cat("  - alluvial_16s_archaea_pct_values_only.pdf\n")
cat("  - alluvial_16s_archaea_pct_flow_annotations.tsv\n")
cat("  - alluvial_16s_archaea_pct_node_descriptions.tsv\n")
cat("  - alluvial_16s_archaea_pct_summary.tsv\n")

cat("\n16S Archaea percentage alluvial plot generated with:\n")
cat("  - Archaea-specific data filtering and percentage calculations\n")
cat("  - Advanced alluvial aesthetics with optimized flow guidance\n")
cat("  - Professional archaea color scheme\n")
cat("  - Detailed flow annotations (TSV format)\n")
cat("  - Node descriptions and summary statistics\n")
cat("  - Percentage normalization based on archaea-only totals\n")
cat("  - Clean minimal appearance with thin nodes\n")

cat("\n⚠️  Note: If flows appear thin, this is expected for archaea:\n")
cat("   1. Archaea are less abundant than bacteria in most environments\n")
cat("   2. Archaea are underrepresented in NCBI genome databases\n")
cat("   3. Most archaea diversity is captured in 16S census data, not genomes\n")
cat("   4. Percentage calculations are based on archaea-only totals for accuracy\n")
