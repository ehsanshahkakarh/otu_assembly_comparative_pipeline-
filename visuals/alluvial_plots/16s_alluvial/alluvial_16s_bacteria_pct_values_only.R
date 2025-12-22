#!/usr/bin/env Rscript

# ============================================================================
# 16S Bacterial Alluvial Plot - PERCENTAGE VALUES (Bacteria-Specific Approach)
# ============================================================================
# This script creates percentage-based alluvial plots for 16S bacterial data
# using ONLY bacteria-specific data filtering and percentage calculations.
# 
# Data Flow: Genbank Genomes → IMG Sequences → 16S OTUs → Genbank Species
# 
# Key Features:
# - Uses only bacteria domain data for percentage calculations
# - Advanced alluvial preprocessing for clean flows
# - Professional bacteria color schemes from shared config
# - Optimized aesthetics (thin nodes, elegant flows)
# - Percentage normalization based on bacteria-only totals
# - Bacteria-specific phyla organization and filtering
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)

cat("=== 16S Bacteria Alluvial Plot (Percentage Values) ===\n")
cat("Using bacteria-specific data approach for reliable visualization\n\n")

# Create output directories if they don't exist
if (!dir.exists("figures")) {
  dir.create("figures")
  cat("Created 'figures/' directory for plot outputs\n")
}

if (!dir.exists("annotations")) {
  dir.create("annotations")
  cat("Created 'annotations/' directory for TSV outputs\n")
}

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

# Load shared color configuration
source("../../shared_config/color_mapping_functions.R")

# Load colors from shared configuration
cat("Loading colors from shared configuration...\n")
color_config <- load_taxonomic_colors("../../shared_config/taxonomic_color_mapping.yaml")

# Get bacterial colors from shared configuration
get_bacterial_colors <- function() {
  bacteria_list <- color_config$bacteria_colors
  bacteria_colors <- character(length(bacteria_list))
  names(bacteria_colors) <- names(bacteria_list)
  
  for (i in 1:length(bacteria_list)) {
    bacteria_colors[names(bacteria_list)[i]] <- bacteria_list[[i]]
  }
  
  return(bacteria_colors)
}

bacterial_colors <- get_bacterial_colors()
cat("Loaded", length(bacterial_colors), "bacterial colors from shared configuration\n")

# Extended color palette for unmapped phyla
extended_colors <- c(
  "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
  "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
  "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9"
)

# Process bacteria-specific data
process_bacteria_data <- function() {
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Bacteria")  # BACTERIA ONLY
  
  return(matched_data)
}

# Get bacteria .U. entries
get_bacteria_u_entries <- function(census_data) {
  # Filter for bacterial .U. entries with specific patterns
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 50) %>%
    filter(grepl("Bacteria", Name_to_use) |
           grepl("Proteobacteria", Name_to_use) |
           grepl("Firmicutes", Name_to_use) |
           grepl("Bacteroidetes", Name_to_use) |
           grepl("Actinobacteria", Name_to_use) |
           (!grepl("Eukaryota", Name_to_use) & !grepl("Archaea", Name_to_use))) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = "Bacteria",
      match_status = "census_only"
    )
  
  return(u_entries)
}

# Process the data
matched_data <- process_bacteria_data()
u_entries <- get_bacteria_u_entries(census_division_data)

cat("Bacteria data processing results:\n")
cat("  - Matched bacteria phyla:", nrow(matched_data), "\n")
cat("  - Bacteria .U. entries:", nrow(u_entries), "\n")

# Display bacteria phyla found
if (nrow(matched_data) > 0) {
  cat("\nBacteria phyla in dataset:\n")
  bacteria_phyla <- unique(matched_data$phylum)
  for (i in 1:length(bacteria_phyla)) {
    phylum_data <- matched_data[matched_data$phylum == bacteria_phyla[i], ]
    total_genomes <- sum(phylum_data$ncbi_genome_count, na.rm = TRUE)
    total_species <- sum(phylum_data$ncbi_species_count, na.rm = TRUE)
    total_otus <- sum(phylum_data$census_otu_count, na.rm = TRUE)
    total_sequences <- sum(phylum_data$census_size_count, na.rm = TRUE)

    cat(sprintf("  %d. %-25s | Genomes: %8s | Species: %8s | OTUs: %8s | Sequences: %8s\n",
                i, bacteria_phyla[i],
                format(total_genomes, big.mark = ","),
                format(total_species, big.mark = ","),
                format(total_otus, big.mark = ","),
                format(total_sequences, big.mark = ",")))
  }
}

# Combine matched bacteria data with bacteria .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals from BACTERIA-ONLY data for proper percentage calculations
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)

# For sequences and OTUs, use COMBINED bacteria data (matched + .U. entries)
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)

cat(paste("\nTotal counts (BACTERIA-ONLY calculations):\n"))
cat(paste("  Bacteria Genomes:", scales::comma(total_genome_count), "\n"))
cat("  Bacteria 16S Sequences:", scales::comma(total_size_count), "\n")
cat("  Bacteria 16S OTUs:", scales::comma(total_otu_count), "\n")
cat(paste("  Bacteria Species:", scales::comma(total_species_count), "\n\n"))

# Calculate percentages for matched data
matched_data <- matched_data %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    size_pct = (census_size_count / total_size_count) * 100,
    otu_pct = (census_otu_count / total_otu_count) * 100
  )

# Calculate percentages for .U. entries
if (nrow(u_entries) > 0) {
  u_entries <- u_entries %>%
    mutate(
      genome_pct = 0,  # No genomes for .U. entries
      species_pct = 0,  # No species for .U. entries
      size_pct = (census_size_count / total_size_count) * 100,
      otu_pct = (census_otu_count / total_otu_count) * 100
    )
}

# Combine matched and .U. data
combined_data <- bind_rows(matched_data, u_entries)

# Calculate total representation for each phylum (sum of all 4 percentages)
combined_data <- combined_data %>%
  mutate(total_representation = genome_pct + species_pct + size_pct + otu_pct)

# Select top 12 bacteria phyla by total representation
n_top <- 12
top_phyla <- combined_data %>%
  arrange(desc(total_representation)) %>%
  head(n_top)

cat(paste("Top", n_top, "bacteria phyla selected by total representation\n\n"))

# Display selected phyla
cat("=== Top", n_top, "Bacteria Phyla ===\n")
for (i in 1:nrow(top_phyla)) {
  entry_type <- if(grepl("\\.U\\.", top_phyla$phylum[i])) "(.U.)" else "(matched)"
  cat(sprintf("%2d. %-25s %10s | Total Rep: %6.2f%% | Genomes: %6.2f%% | Species: %6.2f%% | OTUs: %6.2f%% | Seqs: %6.2f%%\n",
              i,
              top_phyla$phylum[i],
              entry_type,
              top_phyla$total_representation[i],
              top_phyla$genome_pct[i],
              top_phyla$species_pct[i],
              top_phyla$otu_pct[i],
              top_phyla$size_pct[i]))
}

# Calculate "Other" category
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

cat(sprintf("\nOther category: %.2f%% genomes, %.2f%% species, %.2f%% OTUs, %.2f%% sequences\n\n",
            other_genome_pct, other_species_pct, other_otu_pct, other_size_pct))

# Create long format data for alluvial plot
long_data <- data.frame()

for (i in 1:nrow(top_phyla)) {
  phylum_name <- top_phyla$phylum[i]
  is_u_entry <- grepl("\\.U\\.", phylum_name)

  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(phylum_name, 4),
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
other_data_long <- data.frame(
  alluvium = rep(nrow(top_phyla) + 1, 4),
  phylum = rep("Other", 4),
  x = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  stratum = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)

long_data <- rbind(long_data, other_data_long)

# Fix x-axis ordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))
long_data$stratum <- factor(long_data$stratum, levels = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Keep original alluvium IDs intact
cat("Applying advanced alluvial preprocessing...\n")

# Keep original data structure - DO NOT reassign alluvium IDs!
long_data_f <- long_data

# Replace zero values with minimal visible values (0.1%) for .U. entries only
long_data_f <- long_data_f %>%
  mutate(percentage = ifelse(percentage == 0 & grepl("\\.U\\.", phylum), 0.1, percentage))

# Order phyla by size at the first axis (prettier strata stacking)
first_axis <- sort(unique(long_data_f$x))[1]
sizes_first <- long_data_f %>%
  filter(x == first_axis) %>%
  arrange(desc(percentage)) %>%
  select(phylum) %>% pull()

long_data_f <- long_data_f %>%
  mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

cat("Advanced preprocessing complete - preserved original alluvium IDs and data structure\n")

# Assign colors to phyla using shared configuration
cat("\nAssigning colors to phyla from shared configuration...\n")

phyla_in_plot <- unique(long_data_f$phylum)
colors <- character(length(phyla_in_plot))
names(colors) <- phyla_in_plot

for (i in 1:length(phyla_in_plot)) {
  phylum_name <- as.character(phyla_in_plot[i])

  # Handle "Other" category
  if (phylum_name == "Other") {
    colors[i] <- "#CCCCCC"  # Light gray for "Other"
    next
  }

  # Clean phylum name for matching (remove .U. suffix if present)
  clean_phylum <- gsub("\\.U\\..*$", "", phylum_name)

  # Try to match with bacterial color palette from shared config
  if (clean_phylum %in% names(bacterial_colors)) {
    colors[i] <- bacterial_colors[clean_phylum]
  } else {
    # Use extended fallback colors for unmapped phyla
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
  }
}

cat("Colors assigned to", length(colors), "phyla\n")

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

# Create ADVANCED alluvial plot with optimized aesthetics
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
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 110),  # Increased to accommodate node labels
    expand = expansion(mult = c(0, 0))
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
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Percentage (%)")

# Save plot to figures/ directory
cat("\nSaving plot...\n")
ggsave("figures/alluvial_16s_bacteria_pct_values_only.png", p_pct, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("figures/alluvial_16s_bacteria_pct_values_only.pdf", p_pct, width = 24, height = 10, dpi = 300, bg = "white")

cat("Plot saved successfully to figures/ directory\n")

# Create detailed flow annotations file
cat("\nCreating flow annotations file...\n")

flow_annotations <- data.frame(
  Taxon = character(),
  Node = character(),
  Node_Order = integer(),
  Absolute_Count = numeric(),
  Percentage = numeric(),
  Flow_Width = numeric(),
  stringsAsFactors = FALSE
)

# Add annotations for each phylum at each node
for (i in 1:nrow(top_phyla)) {
  phylum_name <- top_phyla$phylum[i]

  # Node 1: Genbank Genome %
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = phylum_name,
    Node = "Genbank_Genome_%",
    Node_Order = 1,
    Absolute_Count = top_phyla$ncbi_genome_count[i],
    Percentage = round(top_phyla$genome_pct[i], 2),
    Flow_Width = top_phyla$genome_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 2: IMG Genome %
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = phylum_name,
    Node = "IMG_Genome_%",
    Node_Order = 2,
    Absolute_Count = top_phyla$census_size_count[i],
    Percentage = round(top_phyla$size_pct[i], 2),
    Flow_Width = top_phyla$size_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 3: 16S OTU %
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = phylum_name,
    Node = "16S_OTU_%",
    Node_Order = 3,
    Absolute_Count = top_phyla$census_otu_count[i],
    Percentage = round(top_phyla$otu_pct[i], 2),
    Flow_Width = top_phyla$otu_pct[i],
    stringsAsFactors = FALSE
  ))

  # Node 4: Genbank Species %
  flow_annotations <- rbind(flow_annotations, data.frame(
    Taxon = phylum_name,
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

write.table(flow_annotations, "annotations/alluvial_16s_bacteria_pct_flow_annotations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    "Genbank Total Genomes (Bacteria)",
    "IMG Genome Count (16S sequences)",
    "16S OTU Count",
    "Genbank Total Species (Bacteria)"
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

write.table(node_descriptions, "annotations/alluvial_16s_bacteria_pct_node_descriptions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary statistics file
summary_stats <- data.frame(
  Metric = c(
    "Total_Taxa_Shown",
    "Top_Taxa_Count",
    "Other_Category_Included",
    "Total_Bacteria_Genomes",
    "Total_16S_Sequences",
    "Total_16S_OTUs",
    "Total_Bacteria_Species",
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
    "Bacteria color mapping from shared config",
    "Bacteria-only totals"
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, "annotations/alluvial_16s_bacteria_pct_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 16S Bacteria Percentage Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  📊 Figures:\n")
cat("     - figures/alluvial_16s_bacteria_pct_values_only.png\n")
cat("     - figures/alluvial_16s_bacteria_pct_values_only.pdf\n")
cat("  📋 Annotations:\n")
cat("     - annotations/alluvial_16s_bacteria_pct_flow_annotations.tsv\n")
cat("     - annotations/alluvial_16s_bacteria_pct_node_descriptions.tsv\n")
cat("     - annotations/alluvial_16s_bacteria_pct_summary.tsv\n")

cat("\n16S Bacteria percentage alluvial plot generated with:\n")
cat("  - Bacteria-specific data filtering and percentage calculations\n")
cat("  - Advanced alluvial aesthetics with optimized flow guidance\n")
cat("  - Professional bacteria color scheme from shared config\n")
cat("  - Detailed flow annotations (TSV format)\n")
cat("  - Node descriptions and summary statistics\n")
cat("  - Percentage normalization based on bacteria-only totals\n")
cat("  - Clean minimal appearance with thin nodes\n")
cat("  - Preserved alluvium IDs for proper flow continuity\n")

