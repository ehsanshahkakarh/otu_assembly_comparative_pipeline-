#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - CLEAN ABSOLUTE VALUES VERSION
# =======================================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#
# For prokaryotic (16S) data only
# Shows top 10 shared phyla by total representation
# Uses ABSOLUTE VALUES for both node heights and flow widths
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

cat("=== 16S rRNA Prokaryotic Clean Alluvial Plot Generator ===\n\n")

# Load the 16S rRNA data
cat("Loading 16S rRNA data...\n")
data_16s <- read.csv("../../Eukcensus_merge/merged_output/16s_merged/results/16s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("16S rRNA data loaded:", nrow(data_16s), "rows\n")

# Filter for matched phyla only and remove N/A values
matched_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(coverage_percentage) & !is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum))

cat("Matched phyla found:", nrow(matched_data), "\n")

# ADDITION: Load original 16S rRNA census data to extract .U. entries
census_division_data <- read.csv("../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv", stringsAsFactors = FALSE)

# Extract .U. entries from original 16S rRNA census data for separate visualization
u_entries <- census_division_data %>%
  filter(grepl("\\.U\\.", Name_to_use)) %>%
  filter(otu_percentage >= 0.5) %>%  # Only include significant .U. entries (>= 0.5%)
  select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count,
         otu_percentage, size_percentage) %>%
  mutate(
    # Set NCBI values to 0 for .U. entries since they have no matches
    ncbi_genome_count = 0,
    ncbi_species_count = 0,
    isolate_count = 0,
    genome_pct_db = 0,
    species_pct = 0,
    isolate_percentage = 0,
    coverage_percentage = 0,
    match_status = 'unidentified',
    domain = ifelse(grepl("Bacteria", phylum), "Bacteria", "Archaea")
  )

cat("Significant .U. entries found:", nrow(u_entries), "\n")
if (nrow(u_entries) > 0) {
  cat("U entries:\n")
  print(u_entries[, c("phylum", "otu_percentage", "size_percentage")])
}

# Calculate total counts from ORIGINAL RAW ASSEMBLY files for proper "Other" calculation
# Load original NCBI phylum data for complete totals
ncbi_phylum_data <- read.csv("../../../ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv", stringsAsFactors = FALSE)

# Load original raw 16S assembly file to get true totals
cat("Loading original 16S assembly file for true totals...\n")
raw_16s_data <- read.csv("../../../16S_censusparse/metadata/eukcensus_16S.clusters.97.tsv",
                         sep = "\t", stringsAsFactors = FALSE)

# Use COMPLETE dataset totals from original assemblies (not parsed summaries)
total_genome_count <- sum(ncbi_phylum_data$phylum_genome_count, na.rm = TRUE)  # All NCBI genomes
total_species_count <- sum(ncbi_phylum_data$phylum_species_count, na.rm = TRUE)  # All NCBI species
total_size_count <- nrow(raw_16s_data)  # Total sequences from original 16S assembly
total_otu_count <- length(unique(raw_16s_data$OTU))  # Total unique OTUs from original 16S assembly

cat("Total counts from ORIGINAL RAW ASSEMBLY files:\n")
cat("  Genomes:", scales::comma(total_genome_count), "\n")
cat("  Sequences:", scales::comma(total_size_count), "\n")
cat("  OTUs:", scales::comma(total_otu_count), "\n")
cat("  Species:", scales::comma(total_species_count), "\n")

# Use the actual percentages from CSV for ranking (no recalculation needed)
# Calculate total representation using the CSV percentages and get top entries
matched_data$total_representation <- matched_data$genome_pct_db + matched_data$size_percentage +
                                    matched_data$otu_percentage + matched_data$species_pct

# UPDATED: Combine matched data with .U. entries for comprehensive visualization
combined_data <- bind_rows(matched_data, u_entries)

# Select top entries by size_percentage (since .U. entries have no genome data)
top_phyla <- combined_data %>%
  arrange(desc(size_percentage)) %>%
  head(12)  # Increased to 12 to accommodate .U. entries

cat("Top entries selected (matched + .U.):", nrow(top_phyla), "\n")

# Print top 5 phyla for reference
cat("\nTop 5 prokaryotic phyla by total representation:\n")
for (i in 1:min(5, nrow(top_phyla))) {
  cat(sprintf("%d. %s: Genome=%.1f%%, Size=%.1f%%, OTU=%.1f%%, Species=%.1f%%\n",
              i, top_phyla$phylum[i], top_phyla$genome_pct_db[i],
              top_phyla$size_percentage[i], top_phyla$otu_percentage[i], top_phyla$species_pct[i]))
}

# Calculate "Other" counts (everything not in top selected entries)
other_genome_count <- total_genome_count - sum(top_phyla$ncbi_genome_count, na.rm = TRUE)
other_size_count <- total_size_count - sum(top_phyla$census_size_count, na.rm = TRUE)
other_otu_count <- total_otu_count - sum(top_phyla$census_otu_count, na.rm = TRUE)
other_species_count <- total_species_count - sum(top_phyla$ncbi_species_count, na.rm = TRUE)

# Safety check for negative values
other_genome_count <- max(0, other_genome_count)
other_size_count <- max(0, other_size_count)
other_otu_count <- max(0, other_otu_count)
other_species_count <- max(0, other_species_count)

cat("\nOther counts (remaining taxa):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Verify that top 10 + Other = 100% for each node
top10_genome_sum <- sum(top_phyla$ncbi_genome_count)
top10_size_sum <- sum(top_phyla$census_size_count)
top10_otu_sum <- sum(top_phyla$census_otu_count)
top10_species_sum <- sum(top_phyla$ncbi_species_count)

cat("\nVerification - Each node should sum to 100% of its dataset:\n")
cat("  Genomes: Top10=", scales::comma(top10_genome_sum), " + Other=", scales::comma(other_genome_count), " = ", scales::comma(top10_genome_sum + other_genome_count), " (Total=", scales::comma(total_genome_count), ")\n")
cat("  Sequences: Top10=", scales::comma(top10_size_sum), " + Other=", scales::comma(other_size_count), " = ", scales::comma(top10_size_sum + other_size_count), " (Total=", scales::comma(total_size_count), ")\n")
cat("  OTUs: Top10=", scales::comma(top10_otu_sum), " + Other=", scales::comma(other_otu_count), " = ", scales::comma(top10_otu_sum + other_otu_count), " (Total=", scales::comma(total_otu_count), ")\n")
cat("  Species: Top10=", scales::comma(top10_species_sum), " + Other=", scales::comma(other_species_count), " = ", scales::comma(top10_species_sum + other_species_count), " (Total=", scales::comma(total_species_count), ")\n")

# Create long format data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Add top 10 phyla data
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    y = c(
      top_phyla$genome_pct_db[i],    # Node 1: NCBI Total Genomes percentage (for flow width)
      top_phyla$size_percentage[i],  # Node 2: 16S EukCensus Sequences percentage (for flow width)
      top_phyla$otu_percentage[i],   # Node 3: 16S EukCensus OTUs percentage (for flow width)
      top_phyla$species_pct[i]       # Node 4: NCBI Total Species percentage (for flow width)
    ),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    percentage = c(
      top_phyla$genome_pct_db[i],  # NCBI Total Genomes percentage
      top_phyla$size_percentage[i],  # 16S EukCensus Sequences percentage
      top_phyla$otu_percentage[i],   # 16S EukCensus OTUs percentage
      top_phyla$species_pct[i]       # NCBI Total Species percentage
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
  y = c(
    max(0, 100 - sum(top_phyla$genome_pct_db, na.rm = TRUE)),  # Remaining genome percentage (for flow width)
    max(0, 100 - sum(top_phyla$size_percentage, na.rm = TRUE)),  # Remaining sequences percentage (for flow width)
    max(0, 100 - sum(top_phyla$otu_percentage, na.rm = TRUE)),   # Remaining OTUs percentage (for flow width)
    max(0, 100 - sum(top_phyla$species_pct, na.rm = TRUE))       # Remaining species percentage (for flow width)
  ),
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
  percentage = c(
    max(0, 100 - sum(top_phyla$genome_pct_db, na.rm = TRUE)),  # Remaining genome percentage
    max(0, 100 - sum(top_phyla$size_percentage, na.rm = TRUE)),  # Remaining sequences percentage
    max(0, 100 - sum(top_phyla$otu_percentage, na.rm = TRUE)),   # Remaining OTUs percentage
    max(0, 100 - sum(top_phyla$species_pct, na.rm = TRUE))       # Remaining species percentage
  ),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))

cat("Long data created with", nrow(long_data), "rows\n")

# Create node labels for phyla >= 1.5% on each node
# UPDATED: Match annotation order to actual alluvium stacking order
node_labels <- long_data %>%
  group_by(x, phylum) %>%
  summarise(
    y = first(y),  # Use percentage values for positioning (matches plot y-axis)
    absolute_count = first(absolute_count),
    percentage = first(percentage),
    alluvium = first(alluvium),  # Get alluvium number for proper ordering
    .groups = "drop"
  ) %>%
  filter(percentage >= 1.5) %>%  # Lowered threshold to show .U. entries (Bacteria.U.phylum = 1.61%)
  group_by(x) %>%
  # Order by alluvium number to match actual flow stacking (alluvium 1 = bottom, higher = top)
  arrange(alluvium) %>%
  mutate(
    cumulative_y = cumsum(y),
    label_y = lag(cumulative_y, default = 0) + y/2,
    label_text = ifelse(phylum == "Other",
                       paste0("Other\n", scales::comma(absolute_count, accuracy = 1), "\n(", round(percentage, 1), "%)"),
                       paste0(substr(phylum, 4, 25), "\n", scales::comma(absolute_count, accuracy = 1), "\n(", round(percentage, 1), "%)"))
  ) %>%
  ungroup()

cat("Node labels created:", nrow(node_labels), "annotations\n")
if (nrow(node_labels) > 0) {
  cat("Sample annotations (ordered by alluvium):\n")
  print(head(node_labels[, c("x", "phylum", "alluvium", "percentage", "label_y")], 5))
}

# Create professional color palette with distinct colors for .U. entries and "Other"
# Ensure legend order matches flow order (largest to smallest, then "Other")
phyla_names <- unique(long_data$phylum)
n_colors <- length(phyla_names)

# Use distinct colors for top phyla, special colors for .U. entries, and gray for "Other"
base_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#8dd3c7", "#bebada", "#fb8072",
  "#ffd92f", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"
)

# Initialize colors
colors <- rep(NA, n_colors)
names(colors) <- phyla_names

# Assign special colors for .U. entries (distinctive patterns)
u_colors <- c("#ff6b6b", "#4ecdc4", "#45b7d1", "#96ceb4", "#feca57")  # Bright, distinctive colors
u_counter <- 1

for (i in 1:length(phyla_names)) {
  phylum <- phyla_names[i]
  if (grepl("\\.U\\.", phylum)) {
    # Assign distinctive colors to .U. entries
    colors[phylum] <- u_colors[min(u_counter, length(u_colors))]
    u_counter <- u_counter + 1
  } else if (phylum == "Other") {
    # Dark gray for "Other"
    colors[phylum] <- "#404040"
  } else {
    # Regular colors for matched phyla
    regular_counter <- sum(!grepl("\\.U\\.", phyla_names[1:i]) & phyla_names[1:i] != "Other")
    colors[phylum] <- base_colors[min(regular_counter, length(base_colors))]
  }
}

cat("Color assignments:\n")
for (i in 1:min(5, length(colors))) {
  cat(sprintf("  %s: %s\n", names(colors)[i], colors[i]))
}

cat("Creating alluvial plot...\n")
cat("Each node will show the complete dataset (Top phyla including .U. entries + Other = 100%)\n")

# Create the main alluvial plot with legend on the right
p <- ggplot(long_data,
            aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with percentage-based widths
  geom_flow(alpha = 0.7, curve_type = "xspline", width = 0.05,
            color = "white", size = 0.1) +
  # All nodes with consistent styling
  geom_stratum(alpha = 0.9, color = "white", size = 0.8, width = 0.05) +
  # Add node labels with ggrepel - larger for professional presentation
  geom_text_repel(data = node_labels,
                  aes(x = x, y = label_y, label = label_text),
                  size = 7, color = "black", fontface = "bold",  # Much larger node labels
                  bg.color = "white", bg.r = 0.2,  # Larger background
                  box.padding = 0.6, point.padding = 0.4,  # More padding
                  segment.color = "black", segment.size = 0.5,  # Thicker connector lines
                  max.overlaps = 20,  # Allow more overlapping labels
                  inherit.aes = FALSE) +
  # Professional color scheme
  scale_fill_manual(values = colors, name = "Phylum") +
  # Clean axis labels with expanded x-axis to reduce white space
  scale_x_discrete(labels = c("NCBI\nTotal Genomes", "16S rRNA EukCensus\nSequences",
                             "16S rRNA EukCensus\nOTUs", "NCBI\nTotal Species"),
                   expand = expansion(mult = c(0.05, 0.05))) +
  # Format y-axis for percentage display
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 20),
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  # Clean theme with professional sizing
  theme_minimal(base_family = default_font) +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),  # Even larger title
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "gray40"),  # Larger subtitle
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold"),  # Much larger y-axis title
    axis.text.x = element_text(size = 22, face = "bold"),  # Much larger node titles
    axis.text.y = element_text(size = 14),  # Larger y-axis text
    legend.position = "right",  # Position legend on the right
    legend.title = element_text(size = 20, face = "bold"),  # Moderately large legend title
    legend.text = element_text(size = 18),  # Larger phyla names for better readability
    legend.key.size = unit(2.2, "cm"),  # Moderately large legend symbols
    legend.key.height = unit(2.5, "cm"),  # Moderately large legend keys
    legend.spacing.y = unit(0.6, "cm"),  # Moderate spacing between legend items
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(30, 120, 30, 30)  # Adjusted margins for wider plot: top, right, bottom, left
  ) +
  labs(
    title = "16S rRNA Prokaryotic Data Flow: NCBI Genomes → EukCensus Sequences → EukCensus OTUs → NCBI Species",
    subtitle = paste0("Flow widths represent database percentages; Annotations show absolute counts • Total: ",
                     scales::comma(total_genome_count), " genomes, ",
                     scales::comma(total_size_count), " sequences, ",
                     scales::comma(total_otu_count), " OTUs, ",
                     scales::comma(total_species_count), " species"),
    y = "Database Percentage (%)"
  )

# Save the plot with larger dimensions for professional presentation
cat("Saving plot...\n")
ggsave("alluvial_16s_clean_pct.png", plot = p, width = 32, height = 16, dpi = 300, bg = "white")
ggsave("alluvial_16s_clean_pct.pdf", plot = p, width = 32, height = 16, dpi = 300, bg = "white")

# Create summary statistics
summary_stats <- data.frame(
  Dataset = "16S_Prokaryotic",
  Total_Phyla = nrow(matched_data),
  Top_10_Phyla = nrow(top_phyla),
  Total_Genomes = total_genome_count,
  Total_Sequences = total_size_count,
  Total_OTUs = total_otu_count,
  Total_Species = total_species_count,
  Top_10_Genome_Coverage = round((sum(top_phyla$ncbi_genome_count) / total_genome_count) * 100, 1),
  Top_10_Sequence_Coverage = round((sum(top_phyla$census_size_count) / total_size_count) * 100, 1),
  Top_10_OTU_Coverage = round((sum(top_phyla$census_otu_count) / total_otu_count) * 100, 1),
  Top_10_Species_Coverage = round((sum(top_phyla$ncbi_species_count) / total_species_count) * 100, 1)
)

write.csv(summary_stats, "alluvial_16s_clean_abs_summary.csv", row.names = FALSE)

cat("\n✅ 16S Clean Alluvial Plot completed successfully!\n")
cat("📁 Files saved:\n")
cat("   - alluvial_16s_clean_abs.png (high-res image)\n")
cat("   - alluvial_16s_clean_abs.pdf (vector format)\n")
cat("   - alluvial_16s_clean_abs_summary.csv (statistics)\n")
cat("\n📊 Summary:\n")
cat("   - Total phyla:", nrow(matched_data), "\n")
cat("   - Top 10 coverage: Genomes", summary_stats$Top_10_Genome_Coverage, "%, Sequences",
    summary_stats$Top_10_Sequence_Coverage, "%, OTUs", summary_stats$Top_10_OTU_Coverage,
    "%, Species", summary_stats$Top_10_Species_Coverage, "%\n")
