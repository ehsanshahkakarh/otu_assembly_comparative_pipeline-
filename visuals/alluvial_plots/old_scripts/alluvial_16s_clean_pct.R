#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - CLEAN PERCENTAGE VERSION
# ===================================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#
# For prokaryotic (16S) data only
# Shows top 10 shared phyla by total representation
# Uses PERCENTAGES for flow widths with "Other" calculated as 100 - sum of top 10 percentages
# Clean implementation with proper legend positioning
#
# Author: Clean Alluvial Team
# Date: 2025
#

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)

# Use default fonts to avoid font issues
default_font <- ""

cat("=== 16S Prokaryotic Clean Alluvial Plot Generator - PERCENTAGE VERSION ===\n\n")

# Load the 16S data
cat("Loading 16S data...\n")
data_16s <- read.csv("../../Eukcensus_merge/merged_output/16s_merged/results/16s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("16S data loaded:", nrow(data_16s), "rows\n")

# Filter for matched phyla only and remove N/A values
matched_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(coverage_percentage) & !is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum))

cat("Matched phyla found:", nrow(matched_data), "\n")

# Use the actual percentages from CSV for ranking (no recalculation needed)
# Calculate total representation using the CSV percentages and get top 10
matched_data$total_representation <- matched_data$genome_pct_db + matched_data$size_percentage +
                                    matched_data$otu_percentage + matched_data$species_pct

top_phyla <- matched_data %>%
  arrange(desc(total_representation)) %>%
  head(10)

cat("Top 10 phyla selected:", nrow(top_phyla), "\n")

# Print top 5 phyla for reference
cat("\nTop 5 prokaryotic phyla by total representation:\n")
for (i in 1:min(5, nrow(top_phyla))) {
  cat(sprintf("%d. %s: Genome=%.1f%%, Size=%.1f%%, OTU=%.1f%%, Species=%.1f%%\n",
              i, top_phyla$phylum[i], top_phyla$genome_pct_db[i],
              top_phyla$size_percentage[i], top_phyla$otu_percentage[i], top_phyla$species_pct[i]))
}

# Calculate "Other" percentages as 100 - sum of top 10 percentages for each node
other_genome_pct <- max(0, 100 - sum(top_phyla$genome_pct_db, na.rm = TRUE))
other_size_pct <- max(0, 100 - sum(top_phyla$size_percentage, na.rm = TRUE))
other_otu_pct <- max(0, 100 - sum(top_phyla$otu_percentage, na.rm = TRUE))
other_species_pct <- max(0, 100 - sum(top_phyla$species_pct, na.rm = TRUE))

cat("\nOther percentages (100 - sum of top 10):\n")
cat("  Other Genomes:", round(other_genome_pct, 1), "%\n")
cat("  Other Sequences:", round(other_size_pct, 1), "%\n")
cat("  Other OTUs:", round(other_otu_pct, 1), "%\n")
cat("  Other Species:", round(other_species_pct, 1), "%\n")

# Verification - Each node should sum to 100%
cat("\nVerification - Each node should sum to 100%:\n")
cat("  Genomes: Top10=", round(sum(top_phyla$genome_pct_db, na.rm = TRUE), 1), "% + Other=", round(other_genome_pct, 1), "% = ", round(sum(top_phyla$genome_pct_db, na.rm = TRUE) + other_genome_pct, 1), "%\n")
cat("  Sequences: Top10=", round(sum(top_phyla$size_percentage, na.rm = TRUE), 1), "% + Other=", round(other_size_pct, 1), "% = ", round(sum(top_phyla$size_percentage, na.rm = TRUE) + other_size_pct, 1), "%\n")
cat("  OTUs: Top10=", round(sum(top_phyla$otu_percentage, na.rm = TRUE), 1), "% + Other=", round(other_otu_pct, 1), "% = ", round(sum(top_phyla$otu_percentage, na.rm = TRUE) + other_otu_pct, 1), "%\n")
cat("  Species: Top10=", round(sum(top_phyla$species_pct, na.rm = TRUE), 1), "% + Other=", round(other_species_pct, 1), "% = ", round(sum(top_phyla$species_pct, na.rm = TRUE) + other_species_pct, 1), "%\n")

# Create long format data for alluvial plot using PERCENTAGES for flow widths
cat("\nPreparing data for 4-node alluvial plot with percentage-based flow widths...\n")
long_data <- data.frame()

# Add top 10 phyla data using percentages for flow widths
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    y = c(
      top_phyla$genome_pct_db[i],        # Node 1: NCBI Total Genomes (percentage)
      top_phyla$size_percentage[i],      # Node 2: 16S EukCensus Sequences (percentage)
      top_phyla$otu_percentage[i],       # Node 3: 16S EukCensus OTUs (percentage)
      top_phyla$species_pct[i]           # Node 4: NCBI Total Species (percentage)
    ),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],    # Keep absolute counts for labels
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    percentage = c(
      top_phyla$genome_pct_db[i],        # NCBI Total Genomes percentage
      top_phyla$size_percentage[i],      # 16S EukCensus Sequences percentage
      top_phyla$otu_percentage[i],       # 16S EukCensus OTUs percentage
      top_phyla$species_pct[i]           # NCBI Total Species percentage
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category using calculated percentages
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
  stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
  y = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  absolute_count = c(0, 0, 0, 0),  # Placeholder for absolute counts
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))

cat("Long data created with", nrow(long_data), "rows\n")

# Create node labels with COMPLETELY FIXED positioning
# ggalluvial stacks segments in the order they appear in the data
# We need to calculate positions exactly as ggalluvial does
node_labels <- long_data %>%
  group_by(x) %>%
  # CRITICAL: Sort by alluvium to match ggalluvial's internal ordering
  arrange(alluvium) %>%
  mutate(
    # Calculate the bottom position of each segment
    segment_bottom = lag(cumsum(y), default = 0),
    # Calculate the top position of each segment
    segment_top = cumsum(y),
    # Position label at the exact center of each segment
    label_y = segment_bottom + (y / 2)
  ) %>%
  ungroup() %>%
  filter(y >= 2) %>%  # Show labels for phyla >= 2%
  mutate(
    label_text = ifelse(phylum == "Other",
                       paste0("Other\n", round(y, 1), "%"),
                       paste0(substr(phylum, 4, 25), "\n", round(y, 1), "%"))
  )

# Debug: Print detailed label positions for first node to verify alignment
cat("DETAILED label positioning debug (NCBI_Total_Genomes node):\n")
first_node_debug <- node_labels %>%
  filter(x == "NCBI_Total_Genomes") %>%
  arrange(alluvium)

for(i in 1:nrow(first_node_debug)) {
  cat(sprintf("  %s: segment=%.1f%%-%.1f%%, center=%.1f%%, label_y=%.1f%% %s\n",
              first_node_debug$phylum[i],
              first_node_debug$segment_bottom[i],
              first_node_debug$segment_top[i],
              first_node_debug$segment_bottom[i] + first_node_debug$y[i]/2,
              first_node_debug$label_y[i],
              ifelse(abs((first_node_debug$segment_bottom[i] + first_node_debug$y[i]/2) - first_node_debug$label_y[i]) < 0.01, "✓", "✗")))
}

# Create professional color palette
n_colors <- length(unique(long_data$phylum))
colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#bebada",
  "#fb8072"
)[1:n_colors]
names(colors) <- unique(long_data$phylum)

cat("Creating alluvial plot...\n")
cat("Each node will show percentages with flow widths based on percentages\n")

# Create the main alluvial plot with legend on the right
p <- ggplot(long_data,
            aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with percentage widths
  geom_flow(alpha = 0.7, curve_type = "xspline", width = 0.05,
            color = "white", linewidth = 0.1) +
  # All nodes with consistent styling
  geom_stratum(alpha = 0.9, color = "white", linewidth = 0.8, width = 0.05) +
  # Add white outline by drawing text multiple times with slight offsets
  geom_text(data = node_labels,
            aes(x = x, y = label_y, label = label_text),
            size = 5, color = "white", fontface = "bold",
            nudge_x = 0.02, nudge_y = 0.5,
            inherit.aes = FALSE) +
  geom_text(data = node_labels,
            aes(x = x, y = label_y, label = label_text),
            size = 5, color = "white", fontface = "bold",
            nudge_x = -0.02, nudge_y = -0.5,
            inherit.aes = FALSE) +
  geom_text(data = node_labels,
            aes(x = x, y = label_y, label = label_text),
            size = 5, color = "white", fontface = "bold",
            nudge_x = 0.02, nudge_y = -0.5,
            inherit.aes = FALSE) +
  geom_text(data = node_labels,
            aes(x = x, y = label_y, label = label_text),
            size = 5, color = "white", fontface = "bold",
            nudge_x = -0.02, nudge_y = 0.5,
            inherit.aes = FALSE) +
  # Add black text on top (main text)
  geom_text(data = node_labels,
            aes(x = x, y = label_y, label = label_text),
            size = 5, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  # Professional color scheme
  scale_fill_manual(values = colors, name = "Phylum") +
  # Clean axis labels with expanded spacing
  scale_x_discrete(labels = c("NCBI\nTotal Genomes", "16S EukCensus\nSequences",
                             "16S EukCensus\nOTUs", "NCBI\nTotal Species"),
                   expand = expansion(mult = c(0.15, 0.15))) +
  # Format y-axis with percentage labels and ensure 0-100% scale
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  # Clean theme
  theme_minimal(base_family = default_font) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "cm"),
    legend.margin = margin(l = 20),
    legend.box.margin = margin(l = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 120, 20, 60)
  ) +
  labs(
    title = "16S Prokaryotic Data Flow: NCBI Genomes → EukCensus Sequences → EukCensus OTUs → NCBI Species",
    subtitle = "Top 10 Phyla by Total Representation (Percentage-based Flow Widths) • Other = 100% - Top 10%",
    y = "Percentage (%)"
  )

# Save the plot with extra width for better readability and legend
cat("Saving plot...\n")
ggsave("alluvial_16s_clean_pct.png", plot = p, width = 20, height = 12, dpi = 300, bg = "white")
ggsave("alluvial_16s_clean_pct.pdf", plot = p, width = 20, height = 12, bg = "white")

# Create summary statistics
summary_stats <- data.frame(
  Dataset = "16S_Prokaryotic_Percentage",
  Total_Phyla = nrow(matched_data),
  Top_10_Phyla = nrow(top_phyla),
  Top_10_Genome_Coverage = round(sum(top_phyla$genome_pct_db, na.rm = TRUE), 1),
  Top_10_Sequence_Coverage = round(sum(top_phyla$size_percentage, na.rm = TRUE), 1),
  Top_10_OTU_Coverage = round(sum(top_phyla$otu_percentage, na.rm = TRUE), 1),
  Top_10_Species_Coverage = round(sum(top_phyla$species_pct, na.rm = TRUE), 1),
  Other_Genome_Percentage = round(other_genome_pct, 1),
  Other_Sequence_Percentage = round(other_size_pct, 1),
  Other_OTU_Percentage = round(other_otu_pct, 1),
  Other_Species_Percentage = round(other_species_pct, 1)
)

write.csv(summary_stats, "alluvial_16s_clean_pct_summary.csv", row.names = FALSE)

cat("\n✅ 16S Clean Percentage Alluvial Plot completed successfully!\n")
cat("📁 Files saved:\n")
cat("   - alluvial_16s_clean_pct.png (high-res image)\n")
cat("   - alluvial_16s_clean_pct.pdf (vector format)\n")
cat("   - alluvial_16s_clean_pct_summary.csv (statistics)\n")
cat("\n📊 Summary:\n")
cat("   - Total phyla:", nrow(matched_data), "\n")
cat("   - Top 10 coverage: Genomes", summary_stats$Top_10_Genome_Coverage, "%, Sequences",
    summary_stats$Top_10_Sequence_Coverage, "%, OTUs", summary_stats$Top_10_OTU_Coverage,
    "%, Species", summary_stats$Top_10_Species_Coverage, "%\n")
cat("   - Other percentages: Genomes", summary_stats$Other_Genome_Percentage, "%, Sequences",
    summary_stats$Other_Sequence_Percentage, "%, OTUs", summary_stats$Other_OTU_Percentage,
    "%, Species", summary_stats$Other_Species_Percentage, "%\n")
