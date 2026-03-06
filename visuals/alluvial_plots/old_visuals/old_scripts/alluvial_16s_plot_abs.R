#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - ABSOLUTE VALUES VERSION
# =================================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#
# For prokaryotic (16S) data only
# Shows top 15 shared phyla by total representation
# Uses ABSOLUTE VALUES for both node heights and flow widths
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

# Try to load extrafont for better typography (optional)
extrafont_available <- FALSE
tryCatch({
  library(extrafont)
  extrafont_available <- TRUE
  cat("✓ extrafont loaded - enhanced typography available\n")
}, error = function(e) {
  cat("⚠ extrafont not available - using default fonts\n")
})

# Set default font family
default_font <- if(extrafont_available) "Arial" else ""

cat("=== 16S Prokaryotic Alluvial Plot Generator ===\n\n")

# Load the 16S data
cat("Loading 16S data...\n")
data_16s <- read.csv("../../Eukcensus_merge/merged_output/16s_merged/results/16s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("16S data loaded:", nrow(data_16s), "rows\n")

# Filter for matched phyla only (match_status = 'matched') and remove all N/A values
shared_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(genome_pct_db) & !is.na(isolate_percentage) & !is.na(size_percentage) &
         !is.na(otu_percentage) & !is.na(coverage_percentage) & !is.na(species_pct) &
         !is.na(ncbi_genome_count) & !is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum))

cat("Shared phyla found:", nrow(shared_data), "\n")

# Calculate percentages for each node based on ENTIRE dataset (matched + unmatched)
total_genome_count <- sum(data_16s$ncbi_genome_count, na.rm = TRUE)
total_size_count <- sum(data_16s$census_size_count, na.rm = TRUE)
total_otu_count <- sum(data_16s$census_otu_count, na.rm = TRUE)
total_species_count <- sum(data_16s$ncbi_species_count, na.rm = TRUE)

shared_data$Genome_Percentage <- (shared_data$ncbi_genome_count / total_genome_count) * 100
shared_data$Size_Percentage <- (shared_data$census_size_count / total_size_count) * 100
shared_data$OTU_Percentage <- (shared_data$census_otu_count / total_otu_count) * 100
shared_data$Species_Coverage <- (shared_data$ncbi_species_count / total_species_count) * 100

# Add total representation across all 4 nodes and get top 10
shared_data$total_representation <- shared_data$Genome_Percentage +
                                   shared_data$Size_Percentage +
                                   shared_data$OTU_Percentage +
                                   shared_data$Species_Coverage

top_phyla <- shared_data %>%
  arrange(desc(total_representation)) %>%
  head(10)

cat("Top phyla selected:", nrow(top_phyla), "\n")
cat("Top 10 total counts:\n")
cat("  Top 10 Genomes:", scales::comma(sum(top_phyla$ncbi_genome_count)), "\n")
cat("  Top 10 Sequences:", scales::comma(sum(top_phyla$census_size_count)), "\n")
cat("  Top 10 OTUs:", scales::comma(sum(top_phyla$census_otu_count)), "\n")
cat("  Top 10 Species:", scales::comma(sum(top_phyla$ncbi_species_count)), "\n")

# Print top 5 phyla for reference with all 4 nodes
cat("\nTop 5 prokaryotic phyla by total representation (4-node flow):\n")
for (i in 1:min(5, nrow(top_phyla))) {
  cat(sprintf("%d. %s:\n", i, top_phyla$phylum[i]))
  cat(sprintf("   Genome=%.2f%%, Size=%.2f%%, OTU=%.2f%%, Species=%.2f%%\n",
              top_phyla$Genome_Percentage[i],
              top_phyla$Size_Percentage[i],
              top_phyla$OTU_Percentage[i],
              top_phyla$Species_Coverage[i]))
}

# Create long format data for 4-node alluvial plot using absolute counts
cat("\nPreparing data for 4-node alluvial plot with absolute counts...\n")
long_data <- data.frame()

# Calculate total counts for each node from ENTIRE DATASET to show full scope
total_genome_count <- sum(data_16s$NCBI_Genome_Count, na.rm = TRUE)
total_size_count <- sum(data_16s$Census_Size_Count, na.rm = TRUE)
total_otu_count <- sum(data_16s$Census_OTU_Count, na.rm = TRUE)
total_species_count <- sum(data_16s$NCBI_Species_Count, na.rm = TRUE)

cat("Total counts from ENTIRE dataset:\n")
cat("  Genomes:", scales::comma(total_genome_count), "\n")
cat("  Sequences:", scales::comma(total_size_count), "\n")
cat("  OTUs:", scales::comma(total_otu_count), "\n")
cat("  Species:", scales::comma(total_species_count), "\n")

# Calculate "Other" counts (everything not in top 10 shared phyla)
other_genome_count <- total_genome_count - sum(top_phyla$ncbi_genome_count)
other_size_count <- total_size_count - sum(top_phyla$census_size_count)
other_otu_count <- total_otu_count - sum(top_phyla$census_otu_count)
other_species_count <- total_species_count - sum(top_phyla$ncbi_species_count)

# Safety check for negative values
other_genome_count <- max(0, other_genome_count)
other_size_count <- max(0, other_size_count)
other_otu_count <- max(0, other_otu_count)
other_species_count <- max(0, other_species_count)

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
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),  # Add numbering to phylum names
    x = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
    y = c(
      top_phyla$ncbi_genome_count[i],        # Use absolute counts for node heights
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    flow_width = c(
      top_phyla$ncbi_genome_count[i],        # Use absolute counts for flow widths
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],        # Keep absolute counts for labels
      top_phyla$census_size_count[i],
      top_phyla$census_otu_count[i],
      top_phyla$ncbi_species_count[i]
    ),
    percentage = c(
      (top_phyla$ncbi_genome_count[i] / total_genome_count) * 100,        # Calculate percentage of TOTAL dataset for labels
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
  y = c(other_genome_count, other_size_count, other_otu_count, other_species_count),  # Use absolute counts for node heights
  flow_width = c(other_genome_count, other_size_count, other_otu_count, other_species_count),  # Use absolute counts for flow widths
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),  # Keep absolute counts for labels
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

cat("Long data created with", nrow(long_data), "rows\n")

# Create node labels for phyla >= 5% on each node (including "Other")
node_labels <- long_data %>%
  group_by(x, phylum) %>%
  summarise(
    y = first(absolute_count),  # Use absolute count for positioning in absolute nodes
    absolute_count = first(absolute_count),  # This is absolute count (for labels)
    percentage = first(percentage),  # This is percentage (for filtering and labels)
    .groups = "drop"
  ) %>%
  # Filter out annotations below 5% threshold
  filter(percentage >= 5) %>%
  # Calculate position for labels (center of each node segment) - fixed positioning logic
  group_by(x) %>%
  arrange(y) %>%  # Sort by absolute count ascending for proper stacking
  mutate(
    cumulative_y = cumsum(y),
    label_y = lag(cumulative_y, default = 0) + y/2,  # Fixed positioning calculation
    # Create large, informative labels emphasizing absolute counts for comparison
    label_text = ifelse(phylum == "Other",
                       paste0("Other\n", scales::comma(absolute_count, accuracy = 1), "\n(", round(percentage, 1), "%)"),
                       paste0(substr(phylum, 4, 25), "\n", scales::comma(absolute_count, accuracy = 1), "\n(", round(percentage, 1), "%)"))
  ) %>%
  ungroup()

# Add debugging output
cat("Node labels created:", nrow(node_labels), "labels (>= 5% threshold)\n")
cat("Average labels per node:", round(nrow(node_labels) / length(unique(node_labels$x)), 1), "\n")

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
cat("Creating enhanced 16S alluvial plot...\n")

# Calculate total percentages for subtitle across all 4 nodes
total_genome_pct <- round(sum(top_phyla$Genome_Percentage), 1)
total_size_pct <- round(sum(top_phyla$Size_Percentage), 1)
total_otu_pct <- round(sum(top_phyla$OTU_Percentage), 1)
total_species_pct <- round(sum(top_phyla$Species_Coverage), 1)

p <- ggplot(long_data,
            aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with absolute count widths (shows actual data volumes)
  geom_flow(alpha = 0.7, curve_type = "xspline", width = 0.05,
            color = "white", size = 0.1) +
  # All nodes (including "Other") with consistent styling - minimized width
  geom_stratum(alpha = 0.9, color = "white", size = 0.8, width = 0.05) +
  # Add large node labels with ggrepel - phyla >= 5% on each node
  geom_text_repel(data = node_labels,
                  aes(x = x, y = label_y, label = label_text),
                  size = 10, color = "black", fontface = "bold",  # HUGE size for presentation
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.4, point.padding = 0.3,  # Increased padding for larger text
                  segment.color = "grey40", segment.size = 0.4,
                  max.overlaps = Inf, force = 2,  # Increased force for larger labels
                  direction = "both",
                  inherit.aes = FALSE) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(limits = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
                   labels = c("NCBI Total\nGenomes", "16S rRNA EukCensus\nSequences", "16S rRNA EukCensus\nOTUs", "NCBI Total\nSpecies")) +
  # Use absolute count scale for node heights
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Bacterial Representation Across NCBI Genomic and EukCensus Environmental Data",
    subtitle = paste0("Top 10 shared bacterial phyla + Other • Absolute Count:\n",
                     "Total Genomes: ", scales::comma(sum(top_phyla$NCBI_Genome_Count)),
                     " → Sequences: ", scales::comma(sum(top_phyla$Census_Size_Count)),
                     " → OTUs: ", scales::comma(sum(top_phyla$Census_OTU_Count)),
                     " → Species: ", scales::comma(sum(top_phyla$NCBI_Species_Count))),
    fill = "Phylum",
    y = "Absolute Count",
    caption = paste0("Node heights and flow widths: absolute counts • n = ", nrow(top_phyla),
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
output_base <- "alluvial_16s_ncbi_abs"

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
cat("🧬 16S PROKARYOTIC ALLUVIAL PLOT ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("📊 DATASET OVERVIEW:\n")
cat("   • Total prokaryotic phyla analyzed:", nrow(top_phyla), "\n")
cat("   • Coverage representation: Top 10 by total percentage\n")
cat("   • Data sources: NCBI Genomes → Isolates → EukCensus OTUs\n\n")

cat("📈 COVERAGE STATISTICS (6-Node Flow):\n")
cat("   • Total genome coverage:", round(sum(top_phyla$Genome_Percentage), 2), "%\n")
cat("   • Total isolate coverage:", round(sum(top_phyla$Isolate_Percentage), 2), "%\n")
cat("   • Total strain sequence coverage:", round(sum(top_phyla$Size_Percentage), 2), "%\n")
cat("   • Total OTU coverage:", round(sum(top_phyla$OTU_Percentage), 2), "%\n")
cat("   • Total isolate species coverage:", round(sum(top_phyla$Isolate_Species_Coverage), 2), "%\n")
cat("   • Total species coverage:", round(sum(top_phyla$Species_Coverage), 2), "%\n\n")

cat("📋 AVERAGE PERCENTAGES (6-Node Flow):\n")
cat("   • Average genome percentage:", round(mean(top_phyla$Genome_Percentage), 2), "%\n")
cat("   • Average isolate percentage:", round(mean(top_phyla$Isolate_Percentage), 2), "%\n")
cat("   • Average strain sequence percentage:", round(mean(top_phyla$Size_Percentage), 2), "%\n")
cat("   • Average OTU percentage:", round(mean(top_phyla$OTU_Percentage), 2), "%\n")
cat("   • Average isolate species percentage:", round(mean(top_phyla$Isolate_Species_Coverage), 2), "%\n")
cat("   • Average species percentage:", round(mean(top_phyla$Species_Coverage), 2), "%\n\n")

cat("🎨 VISUALIZATION FEATURES:\n")
cat("   • Enhanced color palette with", length(colors), "distinct colors\n")
cat("   • Percentage labels on top 5 phyla flows\n")
cat("   • Multiple output formats (PNG, PDF, SVG)\n")
cat("   • High-resolution (300 DPI) for publication quality\n\n")

cat("✅ Enhanced 16S prokaryotic alluvial plot generation completed!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Create alternative dark theme visualization
cat("\n🌙 Creating alternative dark theme visualization...\n")

p_dark <- ggplot(long_data,
                aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with absolute count widths (shows actual data volumes)
  geom_flow(alpha = 0.8, curve_type = "xspline", width = 0.05,
            color = "#2c3e50", size = 0.2) +
  # All nodes (including "Other") with consistent styling - minimized width
  geom_stratum(alpha = 0.9, color = "#34495e", size = 0.8, width = 0.05) +
  # Add large node labels with ggrepel (dark theme) - phyla >= 5% on each node
  geom_text_repel(data = node_labels,
                  aes(x = x, y = label_y, label = label_text),
                  size = 6, color = "white", fontface = "bold",  # Increased from 4 to 6
                  bg.color = "#2c3e50", bg.r = 0.15,
                  box.padding = 0.4, point.padding = 0.3,  # Increased padding for larger text
                  segment.color = "grey70", segment.size = 0.4,
                  max.overlaps = Inf, force = 2,  # Increased force for larger labels
                  direction = "both",
                  inherit.aes = FALSE) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(limits = c("NCBI_Total_Genomes", "16S_rRNA_EukCensus_Sequences", "16S_rRNA_EukCensus_OTUs", "NCBI_Total_Species"),
                   labels = c("NCBI Total\nGenomes", "16S rRNA EukCensus\nSequences", "16S rRNA EukCensus\nOTUs", "NCBI Total\nSpecies")) +
  # Use absolute count scale for node heights
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Bacterial Representation Across NCBI Genomic and EukCensus Environmental Data",
    subtitle = paste0("Top 10 shared bacterial phyla + Other • Absolute Count:\n",
                     "Total: ", scales::comma(total_genome_count), " Genomes → ",
                     scales::comma(total_size_count), " Sequences → ",
                     scales::comma(total_otu_count), " OTUs → ",
                     scales::comma(total_species_count), " Species"),
    fill = "Phylum",
    y = "Absolute Count",
    caption = paste0("Node heights and flow widths: absolute counts • n = ", nrow(top_phyla),
                    " phyla + Other (", nrow(shared_data) - nrow(top_phyla), " remaining phyla) • Generated: ", Sys.Date())
  ) +
  theme_minimal() +
  theme(
    text = element_text(color = "white"),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5,
                             margin = margin(b = 18), color = "white"),
    plot.subtitle = element_text(size = 24, hjust = 0.5,
                                margin = margin(b = 25), color = "#bdc3c7"),
    plot.caption = element_text(size = 18, hjust = 1,
                               margin = margin(t = 20), color = "#95a5a6"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "white"),
    axis.text.x = element_text(size = 20, face = "bold", color = "white"),
    axis.text.y = element_text(size = 18, color = "#bdc3c7"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 32, face = "bold", color = "white"),
    legend.text = element_text(size = 24, face = "bold", color = "#bdc3c7"),
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.margin = margin(t = 10),  # Small margin from plot
    legend.box.margin = margin(t = 5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#34495e", size = 0.3),
    panel.background = element_rect(fill = "#2c3e50", color = NA),
    plot.background = element_rect(fill = "#2c3e50", color = NA),
    plot.margin = margin(20, 20, 10, 20)  # Adjusted margins for bottom legend
  )

# Save dark theme version (compact dimensions for clean look)
dark_file <- paste0(output_base, "_dark.png")
ggsave(dark_file, plot = p_dark, width = 24, height = 14, dpi = 300, bg = "#2c3e50")
cat("Dark theme version saved:", dark_file, "\n")

# Create and export comprehensive 4-node summary table
cat("\n📋 Creating comprehensive 4-node summary data table...\n")
summary_table <- top_phyla %>%
  select(Phylum, Genome_Percentage, Size_Percentage, OTU_Percentage, Species_Coverage, total_representation) %>%
  arrange(desc(total_representation)) %>%
  mutate(
    Rank = row_number(),
    Genome_Percentage = round(Genome_Percentage, 2),
    Size_Percentage = round(Size_Percentage, 2),
    OTU_Percentage = round(OTU_Percentage, 2),
    Species_Coverage = round(Species_Coverage, 2),
    Total_Representation = round(total_representation, 2)
  ) %>%
  select(Rank, Phylum, Genome_Percentage, Size_Percentage, OTU_Percentage, Species_Coverage, Total_Representation)

# Export summary table
summary_file <- paste0(output_base, "_summary.csv")
write.csv(summary_table, summary_file, row.names = FALSE)
cat("Summary table exported:", summary_file, "\n")

# Display top 15 in console with all 4 nodes
cat("\n🏆 TOP 15 PROKARYOTIC PHYLA BY TOTAL REPRESENTATION (4-Node Flow):\n")
cat(paste(rep("-", 100), collapse = ""), "\n")
cat(sprintf("%-4s %-20s %8s %8s %8s %8s %8s\n", "Rank", "Phylum", "Genome%", "Size%", "OTU%", "Species%", "Total%"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in 1:min(15, nrow(summary_table))) {
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

# Display the plots (if running interactively)
if (interactive()) {
  print(p)
  cat("\nPress Enter to view dark theme version...")
  readline()
  print(p_dark)
}
