#!/usr/bin/env Rscript
#
# 18S Eukaryotic Alluvial Plot Generator - CLEAN ABSOLUTE VALUES VERSION
# ======================================================================
#
# Creates an alluvial plot showing the flow from:
# NCBI Total Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Total Species
#
# For eukaryotic (18S) data only
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

cat("=== 18S rRNA Eukaryotic Clean Alluvial Plot Generator ===\n\n")

# Load the 18S rRNA data
cat("Loading 18S rRNA data...\n")
data_18s <- read.csv("../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv", stringsAsFactors = FALSE)
cat("18S rRNA data loaded:", nrow(data_18s), "rows\n")

# Filter for matched phyla only
matched_data <- data_18s %>%
  filter(match_status == 'matched')

cat("Matched phyla found:", nrow(matched_data), "\n")

# Debug: Check the data structure
cat("Sample matched data:\n")
if (nrow(matched_data) > 0) {
  print(head(matched_data[, c("phylum", "census_otu_count", "census_size_count", "ncbi_genome_count", "ncbi_species_count", "match_status")]))
} else {
  cat("No matched data found! Checking all data...\n")
  print(table(data_18s$match_status))
}

# Calculate total counts from ORIGINAL RAW ASSEMBLY files for proper "Other" calculation
# Load original NCBI phylum data for complete totals
ncbi_phylum_data <- read.csv("../../../ncbi_parse/csv_ncbi/ncbi_phylum_counts.csv", stringsAsFactors = FALSE)

# FILTER NCBI data by Eukaryota domain only for 18S analysis
ncbi_eukaryota_data <- ncbi_phylum_data %>%
  filter(domain == "Eukaryota")

cat("Total NCBI phyla:", nrow(ncbi_phylum_data), "\n")
cat("Eukaryota NCBI phyla:", nrow(ncbi_eukaryota_data), "\n")

# Load original raw 18S assembly file to get true totals
cat("Loading original 18S assembly file for true totals...\n")
raw_18s_data <- read.csv("../../../18S_censusparse/metadata/eukcensus_18S.clusters.97.tsv",
                         sep = "\t", stringsAsFactors = FALSE)

# Keep parsed division data for .U. entry extraction
census_division_data <- read.csv("../../../18S_censusparse/csv_outputs/eukcensus_18S_by_division.csv", stringsAsFactors = FALSE)

# ADDITION: Extract .U. entries from original 18S rRNA EukCensus data for separate visualization
u_entries <- census_division_data %>%
  filter(grepl("\\.U\\.", Name_to_use)) %>%
  filter(otu_count >= 50) %>%  # Only include significant .U. entries (>= 50 OTUs)
  select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
  mutate(
    # Set NCBI values to 0 for .U. entries since they have no matches
    ncbi_genome_count = 0,
    ncbi_species_count = 0,
    isolate_count = 0,
    genome_pct_db = 0,
    species_pct = 0,
    isolate_percentage = 0,
    otu_percentage = 0,
    size_percentage = 0,
    coverage_percentage = 0,
    match_status = 'unidentified',
    domain = 'Eukaryota'
  )

cat("Significant .U. entries found:", nrow(u_entries), "\n")
if (nrow(u_entries) > 0) {
  cat("U entries:\n")
  print(u_entries[, c("phylum", "otu_percentage", "size_percentage")])
}

# Use COMPLETE dataset totals from original assemblies (not parsed summaries)
# UPDATED: Use only Eukaryota NCBI data for totals
total_genome_count <- sum(ncbi_eukaryota_data$phylum_genome_count, na.rm = TRUE)  # Eukaryota NCBI genomes only
total_species_count <- sum(ncbi_eukaryota_data$phylum_species_count, na.rm = TRUE)  # Eukaryota NCBI species only
total_size_count <- nrow(raw_18s_data)  # Total sequences from original 18S assembly
total_otu_count <- length(unique(raw_18s_data$centroid))  # Total unique OTUs from original 18S assembly

cat("Total counts from ORIGINAL RAW ASSEMBLY files (Eukaryota only):\n")
cat("  Eukaryota Genomes:", scales::comma(total_genome_count), "\n")
cat("  18S Sequences:", scales::comma(total_size_count), "\n")
cat("  18S OTUs:", scales::comma(total_otu_count), "\n")
cat("  Eukaryota Species:", scales::comma(total_species_count), "\n")

# Calculate percentages from scratch using raw counts
cat("=== PERCENTAGE CALCULATION DETAILS ===\n")

# Calculate totals from census division data for proper 18S percentages
total_census_otu_count <- sum(census_division_data$otu_count, na.rm = TRUE)
total_census_size_count <- sum(census_division_data$size_count, na.rm = TRUE)

cat("📊 CENSUS TOTALS (from census division data):\n")
cat("   - Total OTUs:", scales::comma(total_census_otu_count), "\n")
cat("   - Total Sequences:", scales::comma(total_census_size_count), "\n")

cat("\n📊 NCBI TOTALS (Eukaryota only):\n")
cat("   - Total Genomes:", scales::comma(total_genome_count), "\n")
cat("   - Total Species:", scales::comma(total_species_count), "\n")

cat("\n🧮 CALCULATING CENSUS PERCENTAGES:\n")
cat("   Formula: (phylum_count / total_count) * 100\n")

# Calculate fresh percentages for all data
matched_data <- matched_data %>%
  mutate(
    # Calculate 18S percentages using census division totals
    otu_percentage = (census_otu_count / total_census_otu_count) * 100,
    size_percentage = (census_size_count / total_census_size_count) * 100,
    # RECALCULATE NCBI percentages using EUKARYOTA-ONLY totals
    genome_pct_db = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    total_representation = genome_pct_db + size_percentage + otu_percentage + species_pct
  )

cat("\n🧮 NCBI PERCENTAGES (RECALCULATED for Eukaryota only):\n")
cat("   - genome_pct_db: (phylum_genomes / eukaryota_total_genomes) * 100\n")
cat("   - species_pct: (phylum_species / eukaryota_total_species) * 100\n")

cat("\n📋 EXAMPLE CALCULATIONS FOR TOP 3 PHYLA:\n")
for (i in 1:min(3, nrow(matched_data))) {
  phylum <- matched_data$phylum[i]
  otu_count <- matched_data$census_otu_count[i]
  size_count <- matched_data$census_size_count[i]
  genome_count <- matched_data$ncbi_genome_count[i]
  species_count <- matched_data$ncbi_species_count[i]

  otu_pct <- matched_data$otu_percentage[i]
  size_pct <- matched_data$size_percentage[i]
  genome_pct <- matched_data$genome_pct_db[i]
  species_pct <- matched_data$species_pct[i]

  cat(sprintf("\n%d. %s:\n", i, phylum))
  cat(sprintf("   🧬 Census OTUs: %s / %s = %.2f%%\n",
              scales::comma(otu_count), scales::comma(total_census_otu_count), otu_pct))
  cat(sprintf("   🧬 Census Sequences: %s / %s = %.2f%%\n",
              scales::comma(size_count), scales::comma(total_census_size_count), size_pct))
  cat(sprintf("   🏛️ NCBI Genomes: %s / %s = %.2f%%\n",
              scales::comma(genome_count), scales::comma(total_genome_count), genome_pct))
  cat(sprintf("   🏛️ NCBI Species: %s / %s = %.2f%%\n",
              scales::comma(species_count), scales::comma(total_species_count), species_pct))
}

# Also calculate percentages for .U. entries using the same totals
if (nrow(u_entries) > 0) {
  u_entries <- u_entries %>%
    mutate(
      otu_percentage = (census_otu_count / total_census_otu_count) * 100,
      size_percentage = (census_size_count / total_census_size_count) * 100,
      # .U. entries have 0 NCBI counts, so percentages are 0
      genome_pct_db = (ncbi_genome_count / total_genome_count) * 100,  # = 0
      species_pct = (ncbi_species_count / total_species_count) * 100,  # = 0
      total_representation = genome_pct_db + size_percentage + otu_percentage + species_pct
    )

  cat("\n📋 .U. ENTRIES CALCULATIONS:\n")
  for (i in 1:nrow(u_entries)) {
    phylum <- u_entries$phylum[i]
    otu_count <- u_entries$census_otu_count[i]
    size_count <- u_entries$census_size_count[i]
    otu_pct <- u_entries$otu_percentage[i]
    size_pct <- u_entries$size_percentage[i]

    cat(sprintf("\n%d. %s:\n", i, phylum))
    cat(sprintf("   🧬 Census OTUs: %s / %s = %.2f%%\n",
                scales::comma(otu_count), scales::comma(total_census_otu_count), otu_pct))
    cat(sprintf("   🧬 Census Sequences: %s / %s = %.2f%%\n",
                scales::comma(size_count), scales::comma(total_census_size_count), size_pct))
    cat(sprintf("   🏛️ NCBI Genomes: 0 / %s = 0.00%%\n", scales::comma(total_genome_count)))
    cat(sprintf("   🏛️ NCBI Species: 0 / %s = 0.00%%\n", scales::comma(total_species_count)))
  }
}

# UPDATED: Combine matched data with .U. entries for comprehensive visualization
combined_data <- bind_rows(matched_data, u_entries)

# Select top entries by size_percentage (since .U. entries have no genome data)
top_phyla <- combined_data %>%
  arrange(desc(size_percentage)) %>%
  head(12)  # Increased to 12 to accommodate .U. entries

cat("Top entries selected (matched + .U.):", nrow(top_phyla), "\n")

# Print top 5 phyla for reference
cat("\nTop 5 eukaryotic phyla by total representation:\n")
for (i in 1:min(5, nrow(top_phyla))) {
  cat(sprintf("%d. %s: Genome=%.1f%%, Size=%.1f%%, OTU=%.1f%%, Species=%.1f%%\n",
              i, top_phyla$phylum[i], top_phyla$genome_pct_db[i],
              top_phyla$size_percentage[i], top_phyla$otu_percentage[i], top_phyla$species_pct[i]))
}

# Calculate "Other" counts (everything not in top selected entries from COMPLETE datasets)
# UPDATED: Account for both matched phyla and .U. entries in top selection
top_genome_count <- sum(top_phyla$ncbi_genome_count, na.rm = TRUE)
top_size_count <- sum(top_phyla$census_size_count, na.rm = TRUE)
top_otu_count <- sum(top_phyla$census_otu_count, na.rm = TRUE)
top_species_count <- sum(top_phyla$ncbi_species_count, na.rm = TRUE)

other_genome_count <- total_genome_count - top_genome_count
other_size_count <- total_size_count - top_size_count
other_otu_count <- total_otu_count - top_otu_count
other_species_count <- total_species_count - top_species_count

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

# Verify that top entries + Other = 100% for each node
top_genome_sum <- sum(top_phyla$ncbi_genome_count, na.rm = TRUE)
top_size_sum <- sum(top_phyla$census_size_count, na.rm = TRUE)
top_otu_sum <- sum(top_phyla$census_otu_count, na.rm = TRUE)
top_species_sum <- sum(top_phyla$ncbi_species_count, na.rm = TRUE)

cat("\nVerification - Each node should sum to 100% of its dataset:\n")
cat("  Genomes: Top=", scales::comma(top_genome_sum), " + Other=", scales::comma(other_genome_count), " = ", scales::comma(top_genome_sum + other_genome_count), " (Total=", scales::comma(total_genome_count), ")\n")
cat("  Sequences: Top=", scales::comma(top_size_sum), " + Other=", scales::comma(other_size_count), " = ", scales::comma(top_size_sum + other_size_count), " (Total=", scales::comma(total_size_count), ")\n")
cat("  OTUs: Top=", scales::comma(top_otu_sum), " + Other=", scales::comma(other_otu_count), " = ", scales::comma(top_otu_sum + other_otu_count), " (Total=", scales::comma(total_otu_count), ")\n")
cat("  Species: Top=", scales::comma(top_species_sum), " + Other=", scales::comma(other_species_count), " = ", scales::comma(top_species_sum + other_species_count), " (Total=", scales::comma(total_species_count), ")\n")

# Create long format data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Add top phyla data (including .U. entries)
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
    y = c(
      top_phyla$genome_pct_db[i],    # Node 1: NCBI Total Genomes percentage (for flow width)
      top_phyla$size_percentage[i],  # Node 2: 18S EukCensus Sequences percentage (for flow width)
      top_phyla$otu_percentage[i],   # Node 3: 18S EukCensus OTUs percentage (for flow width)
      top_phyla$species_pct[i]       # Node 4: NCBI Total Species percentage (for flow width)
    ),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],        # Node 1: NCBI Total Genomes
      top_phyla$census_size_count[i],        # Node 2: 18S EukCensus Sequences
      top_phyla$census_otu_count[i],         # Node 3: 18S EukCensus OTUs
      top_phyla$ncbi_species_count[i]        # Node 4: NCBI Total Species
    ),
    percentage = c(
      top_phyla$genome_pct_db[i],  # NCBI Total Genomes percentage
      top_phyla$size_percentage[i],  # 18S EukCensus Sequences percentage
      top_phyla$otu_percentage[i],   # 18S EukCensus OTUs percentage
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
  x = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
  stratum = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"),
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
long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "18S_EukCensus_Sequences", "18S_EukCensus_OTUs", "NCBI_Total_Species"))

cat("Long data created with", nrow(long_data), "rows\n")

# Debug: Check long_data structure
cat("Sample long_data:\n")
if (nrow(long_data) > 0) {
  print(head(long_data[, c("phylum", "x", "y", "absolute_count", "percentage")]))
  cat("Y value range:", range(long_data$y, na.rm = TRUE), "\n")
  cat("Any NA values in y:", any(is.na(long_data$y)), "\n")
  cat("Any infinite values in y:", any(is.infinite(long_data$y)), "\n")
} else {
  cat("ERROR: No long_data created!\n")
}

# Create node labels for phyla >= 2% on each node
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
  filter(percentage >= 2) %>%  # Lowered threshold for more annotations
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

# Debug: Print node labels info
cat("Node labels created:", nrow(node_labels), "labels\n")
if (nrow(node_labels) > 0) {
  cat("Sample node labels:\n")
  print(head(node_labels))
} else {
  cat("WARNING: No node labels created! Check percentage threshold.\n")
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

# Create node annotations CSV file
cat("Creating node annotations CSV file...\n")
node_annotations_csv <- node_labels %>%
  select(node = x, phylum, absolute_count, percentage, alluvium) %>%
  arrange(node, desc(percentage)) %>%
  mutate(
    node = case_when(
      node == "NCBI_Total_Genomes" ~ "NCBI_Eukaryota_Genomes",
      node == "18S_EukCensus_Sequences" ~ "18S_EukCensus_Sequences",
      node == "18S_EukCensus_OTUs" ~ "18S_EukCensus_OTUs",
      node == "NCBI_Total_Species" ~ "NCBI_Eukaryota_Species",
      TRUE ~ as.character(node)
    )
  )

write.csv(node_annotations_csv, "alluvial_18s_node_annotations.csv", row.names = FALSE)
cat("Node annotations saved to: alluvial_18s_node_annotations.csv\n")

# Create the main alluvial plot WITHOUT annotations and legend
p_base <- ggplot(long_data,
            aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = phylum)) +
  # Flows with percentage-based widths
  geom_flow(alpha = 0.7, curve_type = "xspline", width = 0.05,
            color = "white", size = 0.1) +
  # All nodes with consistent styling
  geom_stratum(alpha = 0.9, color = "white", size = 0.8, width = 0.05) +
  # Professional color scheme
  scale_fill_manual(values = colors, name = "Phylum") +
  # Clean axis labels with reduced white space
  scale_x_discrete(labels = c("NCBI\nEukaryota Genomes", "18S rRNA EukCensus\nSequences",
                             "18S rRNA EukCensus\nOTUs", "NCBI\nEukaryota Species"),
                   expand = expansion(mult = c(0.05, 0.05))) +  # Reduce x-axis white space
  # Format y-axis for percentage display with 10% intervals
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 10),  # Every 10%
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  # Clean theme with professional sizing - NO LEGEND
  theme_minimal(base_family = default_font) +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),  # Even larger title
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "gray40"),  # Larger subtitle
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold"),  # Much larger y-axis title
    axis.text.x = element_text(size = 22, face = "bold"),  # Much larger node titles
    axis.text.y = element_text(size = 24, face = "bold"),  # Much larger y-axis percentage text
    legend.position = "none",  # REMOVE LEGEND
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(30, 40, 30, 40)  # Adjusted margins without legend space
  ) +
  labs(
    title = "18S rRNA Eukaryotic Data Flow: NCBI Eukaryota Genomes → EukCensus Sequences → EukCensus OTUs → NCBI Eukaryota Species",
    subtitle = paste0("Flow widths represent database percentages • Eukaryota Total: ",
                     scales::comma(total_genome_count), " genomes, ",
                     scales::comma(total_size_count), " sequences, ",
                     scales::comma(total_otu_count), " OTUs, ",
                     scales::comma(total_species_count), " species"),
    y = "Database Percentage (%)"
  )

# Save the plot
cat("Saving plot...\n")
ggsave("alluvial_18s_clean_eukaryota.png", plot = p_base, width = 32, height = 16, dpi = 300, bg = "white")
ggsave("alluvial_18s_clean_eukaryota.pdf", plot = p_base, width = 32, height = 16, bg = "white")

# Create separate legend plot with same colors
cat("Creating separate phyla legend...\n")

# Create legend data frame with the same colors used in the plot (excluding "Other")
legend_data <- data.frame(
  phylum = names(colors)[names(colors) != "Other"],
  color = colors[names(colors) != "Other"],
  stringsAsFactors = FALSE
)

# Create a vertical legend with larger phylum names
legend_plot <- ggplot(legend_data, aes(x = 1, y = seq_along(phylum) * 3, fill = phylum)) +
  geom_tile(width = 2.5, height = 2.5, color = "white", linewidth = 3) +
  geom_text(aes(label = phylum), x = 3.5, hjust = 0, vjust = 0.5, size = 12, fontface = "bold") +
  scale_fill_manual(values = legend_data$color) +
  scale_x_continuous(limits = c(-0.5, 8)) +
  scale_y_continuous(limits = c(0.5, length(legend_data$phylum) * 3 + 0.5)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 60, 20, 20),  # Right margin for text labels
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save the vertical legend with larger tiles as a separate PNG
ggsave("alluvial_18s_phyla_legend.png", plot = legend_plot, width = 8, height = 20, dpi = 300, bg = "white")

cat("Creating absolute count version of alluvial plot...\n")

# Create absolute count version of the data
long_data_abs <- long_data %>%
  mutate(y_abs = absolute_count)

# Create absolute count alluvial plot
p_abs <- ggplot(long_data_abs, aes(x = x, stratum = phylum, alluvium = alluvium, y = y_abs, fill = phylum)) +
  geom_alluvium(alpha = 0.8, decreasing = FALSE) +
  geom_stratum(alpha = 0.9, decreasing = FALSE, color = "white", linewidth = 0.8) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 24, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5, color = "gray40"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "18S rRNA Eukaryotic Data Flow (Absolute Counts)",
    subtitle = paste("NCBI Eukaryota Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Eukaryota Species"),
    y = "Absolute Count"
  )

# Save absolute count plot
ggsave("alluvial_18s_clean_eukaryota_abs_counts.png", plot = p_abs, width = 32, height = 16, dpi = 300, bg = "white")
ggsave("alluvial_18s_clean_eukaryota_abs_counts.pdf", plot = p_abs, width = 32, height = 16, dpi = 300, bg = "white")

cat("Absolute count alluvial plot saved as: alluvial_18s_clean_eukaryota_abs_counts.png\n")
cat("Legend saved as: alluvial_18s_phyla_legend.png\n")

# Create summary statistics
summary_stats <- data.frame(
  Dataset = "18S_Eukaryotic",
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

write.csv(summary_stats, "alluvial_18s_clean_abs_summary.csv", row.names = FALSE)

cat("\n✅ 18S Clean Alluvial Plot completed successfully!\n")
cat("📁 Files saved:\n")
cat("   - alluvial_18s_clean_eukaryota.png (high-res image)\n")
cat("   - alluvial_18s_clean_eukaryota.pdf (vector format)\n")
cat("   - alluvial_18s_phyla_legend.png (separate legend)\n")
cat("   - alluvial_18s_clean_abs_summary.csv (statistics)\n")
cat("   - alluvial_18s_node_annotations.csv (node data with percentages)\n")
cat("\n📊 Summary (Eukaryota-only NCBI data):\n")
cat("   - Total eukaryotic phyla:", nrow(matched_data), "\n")
cat("   - Node annotations:", nrow(node_annotations_csv), "entries\n")
cat("   - Top 10 coverage: Eukaryota Genomes", summary_stats$Top_10_Genome_Coverage, "%, 18S Sequences",
    summary_stats$Top_10_Sequence_Coverage, "%, 18S OTUs", summary_stats$Top_10_OTU_Coverage,
    "%, Eukaryota Species", summary_stats$Top_10_Species_Coverage, "%\n")
