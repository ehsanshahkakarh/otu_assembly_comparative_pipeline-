#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - PERCENTAGE VERSION
# ============================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
# 
# For prokaryotic (16S) data only
# Shows top 15 shared phyla by total representation
# Uses PERCENTAGES for both node heights and flow widths (normalized view)
#
# Author: Enhanced EukCensus Analysis Team
# Date: 2024
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

cat("=== 16S Prokaryotic Alluvial Plot Generator - PERCENTAGE VERSION ===\n\n")

# Load the 16S data from new vectorized merger
cat("Loading 16S data...\n")
data_16s <- read.csv("../../Eukcensus_merge/merged_output/16s_merged/results/16s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("16S data loaded:", nrow(data_16s), "rows\n")

# Filter for matched phyla only and remove all N/A values
shared_data <- data_16s %>%
  filter(match_status == "matched") %>%
  filter(!is.na(coverage_percentage) & !is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum))

cat("Shared phyla found:", nrow(shared_data), "\n")

# Add total representation across all 4 nodes and get top 10
# Use correct column names from the CSV file
shared_data$total_representation <- shared_data$genome_pct_db +
                                   shared_data$size_percentage +
                                   shared_data$otu_percentage +
                                   shared_data$species_pct

top_phyla <- shared_data %>%
  arrange(desc(total_representation)) %>%
  head(10)

cat("Top phyla selected:", nrow(top_phyla), "\n")
cat("Top 10 total counts:\n")
cat("  Top 10 Genomes:", scales::comma(sum(top_phyla$NCBI_Genome_Count)), "\n")
cat("  Top 10 Sequences:", scales::comma(sum(top_phyla$Census_Size_Count)), "\n")
cat("  Top 10 OTUs:", scales::comma(sum(top_phyla$Census_OTU_Count)), "\n")
cat("  Top 10 Species:", scales::comma(sum(top_phyla$NCBI_Species_Count)), "\n")

# Print top 5 phyla for reference with all 4 nodes
cat("\nTop 5 prokaryotic phyla by total representation (4-node flow):\n")
for (i in 1:min(5, nrow(top_phyla))) {
  cat(sprintf("%d. %s:\n", i, top_phyla$Phylum[i]))
  cat(sprintf("   Genome=%.2f%%, Size=%.2f%%, OTU=%.2f%%, Species=%.2f%%\n",
              top_phyla$Genome_Percentage[i],
              top_phyla$Size_Percentage[i],
              top_phyla$OTU_Percentage[i],
              top_phyla$Species_Coverage[i]))
}

# Create long format data for 4-node alluvial plot using percentages
cat("\nPreparing data for 4-node alluvial plot with percentages...\n")
long_data <- data.frame()

# Calculate total counts for each node from ENTIRE DATASET to make nodes 100% normalized
# Use correct column names from the CSV file
total_genome_count <- sum(data_16s$ncbi_genome_count, na.rm = TRUE)
total_size_count <- sum(data_16s$census_size_count, na.rm = TRUE)
total_otu_count <- sum(data_16s$census_otu_count, na.rm = TRUE)
total_species_count <- sum(data_16s$ncbi_species_count, na.rm = TRUE)

cat("Total counts from ENTIRE dataset:\n")
cat("  Genomes:", scales::comma(total_genome_count), "\n")
cat("  Sequences:", scales::comma(total_size_count), "\n")
cat("  OTUs:", scales::comma(total_otu_count), "\n")
cat("  Species:", scales::comma(total_species_count), "\n")

# Calculate "Other" counts (everything not in top 10 shared phyla)
# Use correct column names from the CSV file
other_genome_count <- total_genome_count - sum(top_phyla$ncbi_genome_count)
other_size_count <- total_size_count - sum(top_phyla$census_size_count)
other_otu_count <- total_otu_count - sum(top_phyla$census_otu_count)
other_species_count <- total_species_count - sum(top_phyla$ncbi_species_count)

cat("Other counts (remaining phyla):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Calculate "Other" percentages
other_genome_pct <- (other_genome_count / total_genome_count) * 100
other_size_pct <- (other_size_count / total_size_count) * 100
other_otu_pct <- (other_otu_count / total_otu_count) * 100
other_species_pct <- (other_species_count / total_species_count) * 100

cat("Other percentages:\n")
cat("  Other Genomes:", round(other_genome_pct, 1), "%\n")
cat("  Other Sequences:", round(other_size_pct, 1), "%\n")
cat("  Other OTUs:", round(other_otu_pct, 1), "%\n")
cat("  Other Species:", round(other_species_pct, 1), "%\n")

for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),  # Use correct column name
    x = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
    y = c(
      (top_phyla$ncbi_genome_count[i] / total_genome_count) * 100,        # Use correct column names
      (top_phyla$census_size_count[i] / total_size_count) * 100,
      (top_phyla$census_otu_count[i] / total_otu_count) * 100,
      (top_phyla$ncbi_species_count[i] / total_species_count) * 100
    ),
    flow_width = c(
      (top_phyla$ncbi_genome_count[i] / total_genome_count) * 100,        # Use correct column names
      (top_phyla$census_size_count[i] / total_size_count) * 100,
      (top_phyla$census_otu_count[i] / total_otu_count) * 100,
      (top_phyla$ncbi_species_count[i] / total_species_count) * 100
    ),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],        # Use correct column names
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    percentage = c(
      (top_phyla$ncbi_genome_count[i] / total_genome_count) * 100,        # Use correct column names
      (top_phyla$census_size_count[i] / total_size_count) * 100,
      (top_phyla$census_otu_count[i] / total_otu_count) * 100,
      (top_phyla$ncbi_species_count[i] / total_species_count) * 100
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Add "Other" category to fill remaining node heights (normalized)
other_data <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
  stratum = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
  y = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),  # Use percentages for normalized heights
  flow_width = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),  # Use percentages for flow widths
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),  # Keep absolute counts for labels
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

cat("Long data created with", nrow(long_data), "rows\n")

# Create node labels for phyla >= 1% on each node (including "Other")
# Fix the positioning to match ggalluvial's stacking order
node_labels <- long_data %>%
  group_by(x, phylum) %>%
  summarise(
    y = first(percentage),  # Use percentage for positioning in normalized nodes
    absolute_count = first(absolute_count),  # This is absolute count (for labels)
    percentage = first(percentage),  # This is percentage (for filtering and labels)
    alluvium = first(alluvium),  # Keep alluvium for proper ordering
    .groups = "drop"
  ) %>%
  # Filter out annotations below 5% threshold
  filter(percentage >= 5) %>%
  # Calculate position for labels (center of each node segment) - FIXED positioning logic
  group_by(x) %>%
  # Sort by alluvium to match ggalluvial's internal stacking order
  arrange(alluvium) %>%
  mutate(
    # Calculate cumulative position from bottom up (ggalluvial stacks from bottom)
    cumulative_y = cumsum(y),
    # Position label at center of each segment
    label_y = lag(cumulative_y, default = 0) + y/2,
    # Create large, informative labels emphasizing percentages for normalized comparison
    label_text = ifelse(phylum == "Other",
                       paste0("Other\n", round(percentage, 1), "%\n(", scales::comma(absolute_count, accuracy = 1), ")"),
                       paste0(substr(phylum, 4, 25), "\n", round(percentage, 1), "%\n(", scales::comma(absolute_count, accuracy = 1), ")"))
  ) %>%
  ungroup()

# Add debugging output
cat("Node labels created:", nrow(node_labels), "labels (>= 5% threshold)\n")
cat("Average labels per node:", round(nrow(node_labels) / length(unique(node_labels$x)), 1), "\n")

# Debug label positioning for first node
cat("\nLabel positioning debug (first node):\n")
first_node_labels <- node_labels %>% filter(x == "NCBI_Total_Genomes") %>% arrange(alluvium)
for(i in 1:nrow(first_node_labels)) {
  cat(sprintf("  %s: y=%.1f%%, label_y=%.1f%%\n",
              first_node_labels$phylum[i],
              first_node_labels$y[i],
              first_node_labels$label_y[i]))
}

# Generate enhanced color palette
n_phyla <- length(unique(long_data$phylum))

cat("Number of unique phyla detected:", n_phyla, "\n")

# Check if we have valid data
if (n_phyla == 0) {
  stop("No phyla found in the data. Please check your input data.")
}

# Create a sophisticated color palette combining multiple schemes
if (n_phyla <= 3) {
  # Use minimum 3 colors from Set2 for very small numbers
  colors <- brewer.pal(3, "Set2")[1:n_phyla]
} else if (n_phyla <= 8) {
  # Use Set2 for small number of phyla (more distinct colors)
  colors <- brewer.pal(max(3, n_phyla), "Set2")[1:n_phyla]
} else if (n_phyla <= 12) {
  # Use Set3 for medium number of phyla
  colors <- brewer.pal(max(3, n_phyla), "Set3")[1:n_phyla]
} else {
  # For larger numbers, combine multiple palettes
  base_colors <- c(
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2")
  )

  # Add viridis colors if we need more
  if (n_phyla > 16) {
    additional_colors <- viridis(n_phyla - 16, option = "plasma")
    base_colors <- c(base_colors, additional_colors)
  }

  colors <- base_colors[1:n_phyla]
}

# Ensure colors are distinct and vibrant
colors <- adjustcolor(colors, alpha.f = 0.8)

# Set "Other" to a distinct color that doesn't conflict with phyla colors
other_index <- which(unique(long_data$phylum) == "Other")
if (length(other_index) > 0) {
  colors[other_index] <- "#8B4513"  # Saddle brown for "Other" - distinct from phyla colors
}

cat("Generated", length(colors), "enhanced colors for", n_phyla, "phyla (including Other)\n")

# Create enhanced alluvial plot
cat("Creating enhanced 16S alluvial plot with percentages...\n")

# Calculate total percentages for subtitle across all 4 nodes
# Use correct column names from the CSV file
total_genome_pct <- round(sum(top_phyla$genome_pct_db), 1)
total_size_pct <- round(sum(top_phyla$size_percentage), 1)
total_otu_pct <- round(sum(top_phyla$otu_percentage), 1)
total_species_pct <- round(sum(top_phyla$species_pct), 1)

p <- ggplot(long_data,
            aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with percentage widths (shows coverage within each dataset)
  geom_flow(alpha = 0.7, curve_type = "xspline", width = 0.05,
            color = "white", size = 0.1) +
  # All nodes (including "Other") with consistent styling - minimized width
  geom_stratum(alpha = 0.9, color = "white", size = 0.8, width = 0.05) +
  # Add IMPROVED node labels with ggrepel - phyla >= 5% on each node
  geom_text_repel(data = node_labels,
                  aes(x = x, y = label_y, label = label_text),
                  size = 8, color = "black", fontface = "bold",  # Reduced from 10 to 8 for better legibility
                  bg.color = "white", bg.r = 0.2,  # Increased background radius for better contrast
                  box.padding = 0.6, point.padding = 0.5,  # Increased padding to prevent overlap
                  segment.color = "black", segment.size = 0.5,  # Darker, thicker connector lines
                  max.overlaps = Inf, force = 3,  # Increased force for better separation
                  direction = "both",
                  min.segment.length = 0,  # Always show connector lines
                  seed = 42,  # Set seed for reproducible positioning
                  inherit.aes = FALSE) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(limits = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
                   labels = c("NCBI Total\nGenomes", "16S rRNA EukCensus\nSequences", "16S rRNA EukCensus\nOTUs", "NCBI Total\nSpecies")) +
  # Use percentage scale for normalized node heights
  scale_y_continuous(labels = percent_format(scale = 1, suffix = "%")) +
  labs(
    title = "Bacterial Representation Across NCBI Genomic and EukCensus Environmental Data",
    subtitle = paste0("Top 10 shared bacterial phyla + Other • Proportional Representation:\n",
                     "Total: ", scales::comma(total_genome_count), " Genomes → ",
                     scales::comma(total_size_count), " Sequences → ",
                     scales::comma(total_otu_count), " OTUs → ",
                     scales::comma(total_species_count), " Species"),
    fill = "Phylum",
    y = "Dataset Coverage (%)",
    caption = paste0("Each node: 100% of total dataset • Flow widths: percentage coverage • n = ", nrow(top_phyla),
                    " phyla + Other (", nrow(shared_data) - nrow(top_phyla), " remaining phyla) • Generated: ", Sys.Date())
  ) +
  theme_minimal() +
  theme(
    # Enhanced typography with optimal text sizes
    text = element_text(),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5,
                             margin = margin(b = 18), color = "#2c3e50"),
    plot.subtitle = element_text(size = 24, hjust = 0.5,
                                margin = margin(b = 25), color = "#34495e"),
    plot.caption = element_text(size = 18, hjust = 1,
                               margin = margin(t = 20), color = "#7f8c8d"),

    # Enhanced axis styling with optimal text sizes
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "#2c3e50"),
    axis.text.x = element_text(size = 20, face = "bold", color = "#2c3e50"),
    axis.text.y = element_text(size = 18, color = "#34495e"),

    # HUGE horizontal legend at bottom for presentation
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 32, face = "bold", color = "#2c3e50"),
    legend.text = element_text(size = 24, face = "bold", color = "#34495e"),
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.margin = margin(t = 10),  # Small margin from plot
    legend.box.margin = margin(t = 5),

    # Clean background
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#ecf0f1", size = 0.3),
    panel.background = element_rect(fill = "#fafafa", color = NA),
    plot.background = element_rect(fill = "white", color = NA),

    # Adjusted margins for bottom legend
    plot.margin = margin(20, 20, 10, 20)
  )

# Save enhanced plots in multiple formats
output_base <- "alluvial_16s_ncbi_pct"

# High-resolution PNG for publications (compact dimensions for clean look)
png_file <- paste0(output_base, ".png")
ggsave(png_file, plot = p, width = 24, height = 14, dpi = 300, bg = "white")

# PDF for vector graphics
pdf_file <- paste0(output_base, ".pdf")
ggsave(pdf_file, plot = p, width = 24, height = 14, bg = "white")

cat("Enhanced 16S alluvial plots saved:\n")
cat("  PNG:", png_file, "(", file.exists(png_file), ")\n")
cat("  PDF:", pdf_file, "(", file.exists(pdf_file), ")\n")

# Enhanced summary statistics with visual formatting
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("🧬 16S PROKARYOTIC ALLUVIAL PLOT ANALYSIS SUMMARY - PERCENTAGE VERSION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("📊 DATASET OVERVIEW:\n")
cat("   • Total prokaryotic phyla analyzed:", nrow(top_phyla), "\n")
cat("   • Coverage representation: Top 10 by total percentage\n")
cat("   • Data sources: NCBI Genomes → Sequences → OTUs → Species\n")
cat("   • Visualization: Normalized percentages for comparison\n\n")

cat("📈 COVERAGE STATISTICS (4-Node Flow):\n")
cat("   • Total genome coverage:", round(sum(top_phyla$genome_pct_db), 2), "%\n")
cat("   • Total sequence coverage:", round(sum(top_phyla$size_percentage), 2), "%\n")
cat("   • Total OTU coverage:", round(sum(top_phyla$otu_percentage), 2), "%\n")
cat("   • Total species coverage:", round(sum(top_phyla$species_pct), 2), "%\n\n")

cat("📋 AVERAGE PERCENTAGES (4-Node Flow):\n")
cat("   • Average genome percentage:", round(mean(top_phyla$genome_pct_db), 2), "%\n")
cat("   • Average sequence percentage:", round(mean(top_phyla$size_percentage), 2), "%\n")
cat("   • Average OTU percentage:", round(mean(top_phyla$otu_percentage), 2), "%\n")
cat("   • Average species percentage:", round(mean(top_phyla$species_pct), 2), "%\n\n")

cat("🎨 VISUALIZATION FEATURES:\n")
cat("   • Enhanced color palette with", length(colors), "distinct colors\n")
cat("   • Percentage labels on phyla >= 5% flows\n")
cat("   • Normalized nodes (100% height) for comparison\n")
cat("   • Multiple output formats (PNG, PDF)\n")
cat("   • High-resolution (300 DPI) for publication quality\n\n")

# Create and export comprehensive 4-node summary table
cat("📋 Creating comprehensive 4-node summary data table...\n")
summary_table <- top_phyla %>%
  select(phylum, genome_pct_db, size_percentage, otu_percentage, species_pct, total_representation) %>%
  arrange(desc(total_representation)) %>%
  mutate(
    Rank = row_number(),
    Genome_Percentage = round(genome_pct_db, 2),
    Size_Percentage = round(size_percentage, 2),
    OTU_Percentage = round(otu_percentage, 2),
    Species_Coverage = round(species_pct, 2),
    Total_Representation = round(total_representation, 2)
  ) %>%
  select(Rank, Phylum = phylum, Genome_Percentage, Size_Percentage, OTU_Percentage, Species_Coverage, Total_Representation)

# Export summary table
summary_file <- paste0(output_base, "_summary.csv")
write.csv(summary_table, summary_file, row.names = FALSE)
cat("Summary table exported:", summary_file, "\n")

# Display top 10 in console with all 4 nodes
cat("\n🏆 TOP 10 PROKARYOTIC PHYLA BY TOTAL REPRESENTATION (4-Node Flow):\n")
cat(paste(rep("-", 100), collapse = ""), "\n")
cat(sprintf("%-4s %-20s %8s %8s %8s %8s %8s\n", "Rank", "Phylum", "Genome%", "Size%", "OTU%", "Species%", "Total%"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in 1:min(10, nrow(summary_table))) {
  cat(sprintf("%-4d %-20s %7.2f%% %7.2f%% %7.2f%% %7.2f%% %7.2f%%\n",
              summary_table$Rank[i],
              substr(summary_table$Phylum[i], 1, 20),
              summary_table$Genome_Percentage[i],
              summary_table$Size_Percentage[i],
              summary_table$OTU_Percentage[i],
              summary_table$Species_Coverage[i],
              summary_table$Total_Representation[i]))
}
cat(paste(rep("-", 100), collapse = ""), "\n")

cat("✅ Enhanced 16S prokaryotic alluvial plot generation completed (PERCENTAGE VERSION)!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Display the plot (if running interactively)
if (interactive()) {
  print(p)
}
