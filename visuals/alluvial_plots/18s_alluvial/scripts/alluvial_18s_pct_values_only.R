#!/usr/bin/env Rscript

# ============================================================================
# 18S Eukaryotic Alluvial Plot - PERCENTAGE VALUES (Clean Merged Data Approach)
# ============================================================================
# This script creates percentage-based alluvial plots for 18S eukaryotic data
# using ONLY the clean merged data files for consistency and reliability.
#
# Data Flow: NCBI Eukaryota Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Eukaryota Species
#
# Key Features:
# - Uses only merged CSV files (no raw assembly dependencies)
# - Advanced alluvial preprocessing for clean flows
# - Professional eukaryotic color scheme
# - Optimized aesthetics (thin nodes, elegant flows)
# - Percentage normalization for better comparison
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)

# Load shared color mapping functions
source("../../shared_config/color_mapping_functions.R")

cat("=== 18S Eukaryotic Alluvial Plot (Percentage Values) ===\n")
cat("Using clean merged data approach for reliable visualization\n\n")

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../../misc/Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv"

# Load merged 18S data (ONLY source needed)
cat("Loading 18S merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 18S merged data file at: ", merged_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
data_18s <- read.csv(merged_data_path, stringsAsFactors = FALSE)

# Load census division data for .U. entries
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 18S census data file at: ", census_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("  - 18S merged entries:", nrow(data_18s), "\n")
cat("  - Census division entries:", nrow(census_division_data), "\n\n")

# Process eukaryotic data (filter for Eukaryota domain only)
process_eukaryotic_data <- function() {
  matched_data <- data_18s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Eukaryota")

  return(matched_data)
}

# Load shared color configuration
cat("Loading shared color configuration...\n")
color_config <- load_taxonomic_colors()

# Extract eukaryotic .U. entries from census data
get_eukaryotic_u_entries <- function(census_data) {
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 10) %>%
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

# Process data
matched_data <- process_eukaryotic_data()
u_entries <- get_eukaryotic_u_entries(census_division_data)

cat("Processed eukaryotic data:\n")
cat("  - Matched entries:", nrow(matched_data), "\n")
cat("  - .U. entries:", nrow(u_entries), "\n\n")

# Combine matched data with .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals for percentage calculations
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)

cat("Total counts for percentage calculations:\n")
cat("  - Total Eukaryota Genomes:", scales::comma(total_genome_count), "\n")
cat("  - Total Eukaryota Species:", scales::comma(total_species_count), "\n")
cat("  - Total 18S OTUs:", scales::comma(total_otu_count), "\n")
cat("  - Total 18S Sequences:", scales::comma(total_size_count), "\n\n")

# Select top 8 entries by total representation (genome + sequence coverage)
top_phyla <- combined_data %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    otu_pct = (census_otu_count / total_otu_count) * 100,
    size_pct = (census_size_count / total_size_count) * 100,
    total_representation = genome_pct + species_pct + otu_pct + size_pct
  ) %>%
  arrange(desc(total_representation)) %>%
  head(8)  # Top 8 for eukaryotic diversity

cat("Top 8 eukaryotic divisions selected for visualization\n\n")

# Calculate "Other" percentages
top_genome_pct <- sum(top_phyla$genome_pct, na.rm = TRUE)
top_species_pct <- sum(top_phyla$species_pct, na.rm = TRUE)
top_otu_pct <- sum(top_phyla$otu_pct, na.rm = TRUE)
top_size_pct <- sum(top_phyla$size_pct, na.rm = TRUE)

other_genome_pct <- max(0, 100 - top_genome_pct)
other_species_pct <- max(0, 100 - top_species_pct)
other_otu_pct <- max(0, 100 - top_otu_pct)
other_size_pct <- max(0, 100 - top_size_pct)

# Create long format data for alluvial plot
cat("Creating percentage-based alluvial data...\n")
long_data <- data.frame()

# Add top phyla data
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"),
    stratum = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"),
    percentage = c(
      top_phyla$genome_pct[i],    # Node 1: NCBI Genome %
      top_phyla$size_pct[i],      # Node 2: 18S Sequence %
      top_phyla$otu_pct[i],       # Node 3: 18S OTU %
      top_phyla$species_pct[i]    # Node 4: NCBI Species %
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"),
  stratum = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"),
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)

long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"))
long_data$stratum <- factor(long_data$stratum, levels = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Fix aesthetic issues
cat("Applying advanced alluvial preprocessing...\n")

# 1. Ensure each phylum has a row at every axis (fill missing with 0)
long_data_f <- long_data %>%
  complete(x, phylum, fill = list(percentage = 0)) %>%
  group_by(phylum) %>%
  mutate(alluvium = cur_group_id()) %>%  # stable id per phylum
  ungroup()

# 2. Order phyla by size at the first axis (prettier strata stacking)
first_axis <- sort(unique(long_data_f$x))[1]
sizes_first <- long_data_f %>%
  filter(x == first_axis) %>%
  arrange(desc(percentage)) %>%
  select(phylum) %>% pull()

long_data_f <- long_data_f %>%
  mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

cat("Advanced preprocessing complete - data optimized for clean alluvial flows\n")

# Create professional color palette using shared color system
phyla_names <- unique(long_data_f$phylum)
cat("Assigning colors for", length(phyla_names), "eukaryotic divisions...\n")

# Use shared color mapping system
colors <- get_domain_colors(phyla_names, "Eukaryota", color_config)

# Print color assignments for verification
print_color_summary(phyla_names, colors, "Eukaryota")

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
  scale_fill_manual(values = colors, name = "Division") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
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
    y = "Percentage (%)",
    x = NULL
  )

# Save advanced percentage alluvial plot with optimized dimensions
ggsave("alluvial_18s_pct_values_only.png", p_pct, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("alluvial_18s_pct_values_only.pdf", p_pct, width = 24, height = 10, dpi = 300, bg = "white")

# Create detailed flow annotations with percentages and absolute values
cat("Creating detailed flow annotations file...\n")

# Create detailed annotations for each taxon at each node (including "Other")
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

  # Node 3: 18S OTUs
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "18S_OTU_%",
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
other_genome_count <- sum(combined_data$ncbi_genome_count[!combined_data$phylum %in% top_phyla$phylum], na.rm = TRUE)
other_size_count <- sum(combined_data$census_size_count[!combined_data$phylum %in% top_phyla$phylum], na.rm = TRUE)
other_otu_count <- sum(combined_data$census_otu_count[!combined_data$phylum %in% top_phyla$phylum], na.rm = TRUE)
other_species_count <- sum(combined_data$ncbi_species_count[!combined_data$phylum %in% top_phyla$phylum], na.rm = TRUE)

# Other category flows
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
  Node = "18S_OTU_%",
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

write.table(flow_annotations, "alluvial_18s_pct_flow_annotations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    "Genbank Total Genomes (Eukaryota)",
    "IMG Genome Count (18S sequences)",
    "18S OTU Count",
    "Genbank Total Species (Eukaryota)"
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

write.table(node_descriptions, "alluvial_18s_pct_node_descriptions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary statistics file
cat("Creating summary statistics file...\n")
summary_stats <- data.frame(
  Metric = c(
    "Total_Phyla_Shown",
    "Top_Phyla_Count",
    "Other_Category_Included",
    "Total_Eukaryota_Genomes",
    "Total_18S_Sequences",
    "Total_18S_OTUs",
    "Total_Eukaryota_Species",
    "Filtering_Method",
    "Color_System"
  ),
  Value = c(
    nrow(top_phyla) + 1,  # +1 for Other category
    nrow(top_phyla),
    "Yes",
    scales::comma(total_genome_count),
    scales::comma(total_size_count),
    scales::comma(total_otu_count),
    scales::comma(total_species_count),
    "Top 8 by total representation",
    "Shared taxonomic color mapping"
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, "alluvial_18s_pct_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 18S Eukaryotic Percentage Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_18s_pct_values_only.png\n")
cat("  - alluvial_18s_pct_values_only.pdf\n")
cat("  - alluvial_18s_pct_flow_annotations.tsv\n")
cat("  - alluvial_18s_pct_node_descriptions.tsv\n")
cat("  - alluvial_18s_pct_summary.tsv\n")
cat("\n18S eukaryotic percentage alluvial plot generated with:\n")
cat("  - Clean merged data approach\n")
cat("  - Percentage normalization (0-100%)\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional eukaryotic color scheme\n")
cat("  - Thin nodes and elegant flows\n")
cat("  - Legend positioned far right\n")
cat("  - Node titles moved down\n\n")

cat("Percentage plot advantages:\n")
cat("  - Better comparison across different data types\n")
cat("  - Normalized scale (0-100%) for all nodes\n")
cat("  - Reveals relative proportions more clearly\n")
cat("  - Reduces dominance effects of large absolute numbers\n\n")
