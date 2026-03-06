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
library(yaml)

# Load shared color mapping functions and alluvial filtering functions
source("../../shared_config/color_mapping_functions.R")
source("../config/alluvial_filtering_functions.R")

cat("=== 18S rRNA Eukaryotic Alluvial Plot Generator (YAML Config-Based) ===\n\n")

# Process data using YAML configuration
cat("Processing 18S eukaryotic data using YAML configuration...\n")
result <- process_alluvial_data("eukaryota_18s")

processed_data <- result$data
config_info <- result$config
summary_info <- result$summary

# Load shared color configuration
cat("Loading shared color configuration...\n")
color_config <- load_taxonomic_colors("../../shared_config/taxonomic_color_mapping.yaml")

# Data is already processed by YAML configuration
# processed_data contains the filtered taxa based on the hybrid strategy
cat(paste("YAML filtering completed. Selected", nrow(processed_data), "taxa using",
          config_info$domain$filtering$strategy, "strategy\n"))

# Display summary of selected taxa
cat("\nSelected taxa summary:\n")
for (i in 1:min(10, nrow(processed_data))) {
  cat(paste("  ", i, ".", processed_data$phylum[i],
            "- OTUs:", processed_data$census_otu_count[i],
            "- Sequences:", processed_data$census_size_count[i], "\n"))
}

# Create long format data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Add processed data (filtered by YAML configuration)
for (i in 1:nrow(processed_data)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", processed_data$phylum[i]), 4),
    x = c("Genbank_Genome_Count", "IMG_Genome_Count", "18S_OTU_Count", "Genbank_Species_Count"),
    stratum = c("Genbank_Genome_Count", "IMG_Genome_Count", "18S_OTU_Count", "Genbank_Species_Count"),
    absolute_count = c(
      processed_data$ncbi_genome_count[i],        # Node 1: Genbank Genome Count
      processed_data$census_size_count[i],        # Node 2: IMG Genome Count (sequences)
      processed_data$census_otu_count[i],         # Node 3: 18S OTU Count
      processed_data$ncbi_species_count[i]        # Node 4: Genbank Species Count
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_Count", "IMG_Genome_Count", "18S_OTU_Count", "Genbank_Species_Count"))
long_data$stratum <- factor(long_data$stratum, levels = c("Genbank_Genome_Count", "IMG_Genome_Count", "18S_OTU_Count", "Genbank_Species_Count"))

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

# Create professional color palette using shared color system
phyla_names <- unique(long_data_f$phylum)
cat("Assigning colors for", length(phyla_names), "eukaryotic divisions...\n")

# Use shared color mapping system
colors <- get_domain_colors(phyla_names, "Eukaryota", color_config)

# Print color assignments for verification
print_color_summary(phyla_names, colors, "Eukaryota")

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

# Create detailed flow annotations with absolute values and percentages
cat("Creating detailed flow annotations file...\n")

# Calculate totals for percentage calculations
total_genomes <- sum(processed_data$ncbi_genome_count, na.rm = TRUE)
total_sequences <- sum(processed_data$census_size_count, na.rm = TRUE)
total_otus <- sum(processed_data$census_otu_count, na.rm = TRUE)
total_species <- sum(processed_data$ncbi_species_count, na.rm = TRUE)

# Create detailed annotations for each taxon at each node
flow_annotations <- data.frame()

for (i in 1:nrow(processed_data)) {
  taxon_name <- processed_data$phylum[i]

  # Node 1: Genbank Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Genome_Count",
    Node_Order = 1,
    Absolute_Count = processed_data$ncbi_genome_count[i],
    Percentage = round((processed_data$ncbi_genome_count[i] / total_genomes) * 100, 2),
    Flow_Width = processed_data$ncbi_genome_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 2: IMG Genomes
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "IMG_Genome_Count",
    Node_Order = 2,
    Absolute_Count = processed_data$census_size_count[i],
    Percentage = round((processed_data$census_size_count[i] / total_sequences) * 100, 2),
    Flow_Width = processed_data$census_size_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 3: 18S OTUs
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "18S_OTU_Count",
    Node_Order = 3,
    Absolute_Count = processed_data$census_otu_count[i],
    Percentage = round((processed_data$census_otu_count[i] / total_otus) * 100, 2),
    Flow_Width = processed_data$census_otu_count[i],
    stringsAsFactors = FALSE
  ))

  # Node 4: Genbank Species
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = taxon_name,
    Node = "Genbank_Species_Count",
    Node_Order = 4,
    Absolute_Count = processed_data$ncbi_species_count[i],
    Percentage = round((processed_data$ncbi_species_count[i] / total_species) * 100, 2),
    Flow_Width = processed_data$ncbi_species_count[i],
    stringsAsFactors = FALSE
  ))
}

# Sort by node order and then by flow width (descending)
flow_annotations <- flow_annotations %>%
  arrange(Node_Order, desc(Flow_Width))

write.table(flow_annotations, "alluvial_18s_abs_flow_annotations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_Count", "IMG_Genome_Count", "18S_OTU_Count", "Genbank_Species_Count"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    "Genbank Total Genomes (Eukaryota)",
    "IMG Genome Count (18S sequences)",
    "18S OTU Count",
    "Genbank Total Species (Eukaryota)"
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

write.table(node_descriptions, "alluvial_18s_abs_node_descriptions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary statistics file
cat("Creating summary statistics file...\n")
summary_stats <- data.frame(
  Metric = c(
    "Total_Taxa_Shown",
    "Filtering_Method",
    "Color_System",
    "YAML_Configuration"
  ),
  Value = c(
    nrow(processed_data),
    paste("YAML-based", config_info$domain$filtering$strategy, "strategy"),
    "Shared taxonomic color mapping",
    "Enabled"
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, "alluvial_18s_abs_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 18S Eukaryotic Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - alluvial_18s_abs_values_only.png\n")
cat("  - alluvial_18s_abs_values_only.pdf\n")
cat("  - alluvial_18s_abs_flow_annotations.tsv\n")
cat("  - alluvial_18s_abs_node_descriptions.tsv\n")
cat("  - alluvial_18s_abs_summary.tsv\n")
cat("\n18S eukaryotic alluvial plot generated with:\n")
cat("  - Clean merged data approach\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional eukaryotic color scheme\n")
cat("  - Thin nodes and elegant flows\n")
cat("  - Legend positioned far right\n")
cat("  - Node titles moved down\n\n")
