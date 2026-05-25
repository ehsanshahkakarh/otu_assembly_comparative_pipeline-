#!/usr/bin/env Rscript

# ============================================================================
# 16S Bacterial/Archaea Alluvial Plot - PERCENTAGE VALUES (Clean Merged Data Approach)
# ============================================================================
# This script creates percentage-based alluvial plots for 16S bacterial/archaea data
# using ONLY the clean merged data files for consistency and reliability.
# 
# Data Flow: NCBI Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Species
# 
# Key Features:
# - Uses only merged CSV files (no raw assembly dependencies)
# - Advanced alluvial preprocessing for clean flows
# - Professional bacterial/archaea color schemes
# - Optimized aesthetics (thin nodes, elegant flows)
# - Percentage normalization for better comparison
# - Configurable for both Bacteria and Archaea domains
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)

# CONFIGURATION: Change this to switch between domains
PROCESS_DOMAIN <- "Bacteria"  # Change to "Archaea" for archaea plot

cat("=== 16S", PROCESS_DOMAIN, "Alluvial Plot (Percentage Values) ===\n")
cat("Using clean merged data approach for reliable visualization\n\n")

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../../misc/Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../../misc/16S/csv_16S/eukcensus16S_by_division.csv"

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
cat("  - 16S merged entries:", nrow(data_16s), "\n")
cat("  - Census division entries:", nrow(census_division_data), "\n\n")

# Process domain-specific data
process_domain_data <- function(domain_name) {
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == domain_name)
  
  return(matched_data)
}

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
  for (name in names(bacteria_list)) {
    bacteria_colors[name] <- as.character(bacteria_list[[name]])
  }
  return(bacteria_colors)
}

# Extended fallback colors for bacteria phyla not in main palette
get_extended_bacterial_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

# Excluded phyla that use fallback colors (same as mega 16S script)
get_excluded_phyla <- function() {
  c("Campylobacterota", "Mycoplasmatota", "Thermotogota")
}

# Master archaea color palette from mega 16S script
get_archaea_colors <- function() {
  archaea_color_map <- c(
    "Euryarchaeota" = "#55ee79",        # Most common archaea phylum - bright green
    "Nitrososphaerota" = "#ef9e17",     # Important ammonia-oxidizing archaea
    "Thermoproteota" = "#f61e5d",       # Hyperthermophilic archaea
    "Nanoarchaeota" = "#4ad9d5"         # Highest novelty factor (117.7×)
  )
  return(archaea_color_map)
}

# Extended fallback colors for archaea phyla not in main palette
get_extended_archaea_colors <- function() {
  c("#8B4513", "#2F4F4F", "#800080", "#008B8B", "#B22222", "#FF4500", "#32CD32", "#9932CC")
}

# Extract domain-specific .U. entries from census data
get_domain_u_entries <- function(census_data, domain_name) {
  if (domain_name == "Bacteria") {
    pattern <- "Bacteria\\.U\\."
  } else if (domain_name == "Archaea") {
    pattern <- "Archaea\\.U\\."
  } else {
    return(data.frame())  # Return empty for other domains
  }
  
  u_entries <- census_data %>%
    filter(grepl(pattern, Name_to_use)) %>%
    filter(otu_count >= 10) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = domain_name,
      match_status = "census_only"
    )
  
  return(u_entries)
}

# Process data for the specified domain
matched_data <- process_domain_data(PROCESS_DOMAIN)
u_entries <- get_domain_u_entries(census_division_data, PROCESS_DOMAIN)

cat("Processed", PROCESS_DOMAIN, "data:\n")
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
cat("  - Total", PROCESS_DOMAIN, "Genomes:", scales::comma(total_genome_count), "\n")
cat("  - Total", PROCESS_DOMAIN, "Species:", scales::comma(total_species_count), "\n")
cat("  - Total 16S OTUs:", scales::comma(total_otu_count), "\n")
cat("  - Total 16S Sequences:", scales::comma(total_size_count), "\n\n")

# Select top entries by total representation (domain-appropriate number)
top_n <- if (PROCESS_DOMAIN == "Bacteria") 12 else 8

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

cat("Selected top", top_n, PROCESS_DOMAIN, "phyla by total representation\n")

cat("Top", top_n, PROCESS_DOMAIN, "phyla selected for visualization\n\n")

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

# Add top phyla data (handling .U. entries correctly)
for (i in 1:nrow(top_phyla)) {
  # Check if this is a .U. entry (census_only)
  is_u_entry <- top_phyla$match_status[i] == "census_only"

  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),  # Restore numbered prefixes for color matching
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

# 1. Use original data structure (don't create artificial zero entries)
# The complete() function was creating too many artificial entries that interfered with alluvium structure
# IMPORTANT: Keep the original alluvium IDs - DO NOT reassign them!
long_data_f <- long_data

# 2. Only replace zero values for .U. entries (not all zeros)
# This prevents ggalluvial issues while preserving real data structure
long_data_f <- long_data_f %>%
  mutate(percentage = ifelse(percentage == 0 & grepl("\\.U\\.", phylum), 0.1, percentage))

# 3. Order phyla by size at the first axis (prettier strata stacking)
first_axis <- sort(unique(long_data_f$x))[1]
sizes_first <- long_data_f %>%
  filter(x == first_axis) %>%
  arrange(desc(percentage)) %>%
  select(phylum) %>% pull()

long_data_f <- long_data_f %>%
  mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

cat("Advanced preprocessing complete - preserved original alluvium IDs and data structure\n")

# Create professional color palette using bacterial colors from mega 16S script
phyla_names <- unique(long_data_f$phylum)
n_colors <- length(phyla_names)

# Get bacterial color mappings
bacterial_colors <- get_bacterial_colors()
extended_colors <- get_extended_bacterial_colors()
excluded_phyla <- get_excluded_phyla()

# Assign colors to phyla (matching mega 16S script logic)
colors <- character(n_colors)
names(colors) <- phyla_names

for (i in seq_along(phyla_names)) {
  phylum_name <- phyla_names[i]

  # Remove numbering prefix for color matching (e.g., "1. Pseudomonadota" -> "Pseudomonadota")
  clean_phylum <- gsub("^[0-9]+\\. ", "", phylum_name)

  if (clean_phylum == "Other") {
    colors[i] <- "#808080"  # Gray for "Other"
  } else if (grepl("\\.U\\.", clean_phylum)) {
    # Distinctive yellow-green color for all .U. entries - reduces purple and distinguishes from other greens
    colors[i] <- "#9ACD32"  # Yellow-green - distinctive for unclassified entries
  } else if (clean_phylum %in% excluded_phyla) {
    # Excluded phyla use fallback colors (same as mega 16S script)
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
  } else if (clean_phylum %in% names(bacterial_colors)) {
    colors[i] <- bacterial_colors[clean_phylum]
  } else {
    # Use extended fallback colors for unmapped phyla
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
  }
}

cat("Color mapping complete for", length(colors), PROCESS_DOMAIN, "phyla\n")

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

# DEBUG: Print color assignments
cat("\n=== DEBUG: COLOR ASSIGNMENTS ===\n")
for (i in seq_along(phyla_names)) {
  cat(sprintf("%2d. %-30s -> %s\n", i, phyla_names[i], colors[i]))
}

# DEBUG: Check for Pseudomonadota
pseudo_check <- grep("Pseudomonadota", phyla_names, value = TRUE)
if (length(pseudo_check) > 0) {
  cat("\n✅ Pseudomonadota found in phyla_names:", pseudo_check, "\n")
  cat("   Color:", colors[pseudo_check], "\n")
} else {
  cat("\n❌ Pseudomonadota NOT found in phyla_names!\n")
}

# DEBUG: Check data structure
cat("\n=== DEBUG: DATA STRUCTURE ===\n")
cat("Total rows in long_data_f:", nrow(long_data_f), "\n")
cat("Unique phyla:", length(unique(long_data_f$phylum)), "\n")
cat("Unique alluvium IDs:", length(unique(long_data_f$alluvium)), "\n")

pseudo_rows <- long_data_f[grepl("Pseudomonadota", long_data_f$phylum), ]
if (nrow(pseudo_rows) > 0) {
  cat("\n✅ Pseudomonadota data rows:", nrow(pseudo_rows), "\n")
  cat("Pseudomonadota percentages:\n")
  print(pseudo_rows[, c("phylum", "x", "percentage", "alluvium")])
} else {
  cat("\n❌ NO Pseudomonadota rows in long_data_f!\n")
}
cat("=== END DEBUG ===\n\n")

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

# Create domain-specific output filename
domain_suffix <- tolower(PROCESS_DOMAIN)
output_png <- paste0("alluvial_16s_", domain_suffix, "_pct_values_only.png")
output_pdf <- paste0("alluvial_16s_", domain_suffix, "_pct_values_only.pdf")

# Save advanced percentage alluvial plot with optimized dimensions
ggsave(output_png, p_pct, width = 24, height = 10, dpi = 300, bg = "white")
ggsave(output_pdf, p_pct, width = 24, height = 10, dpi = 300, bg = "white")

# Create detailed flow annotations with percentages and absolute values
cat("Creating detailed flow annotations file...\n")

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

write.table(flow_annotations, paste0("alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_flow_annotations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Create simple node descriptions file
node_descriptions <- data.frame(
  Node = c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"),
  Node_Order = c(1, 2, 3, 4),
  Description = c(
    paste("Genbank Total Genomes (", PROCESS_DOMAIN, ")", sep = ""),
    "IMG Genome Count (16S sequences)",
    "16S OTU Count",
    paste("Genbank Total Species (", PROCESS_DOMAIN, ")", sep = "")
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

write.table(node_descriptions, paste0("alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_node_descriptions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Create summary statistics file
summary_stats <- data.frame(
  Metric = c(
    "Total_Taxa_Shown",
    "Top_Taxa_Count",
    "Other_Category_Included",
    paste("Total_", PROCESS_DOMAIN, "_Genomes", sep = ""),
    "Total_16S_Sequences",
    "Total_16S_OTUs",
    paste("Total_", PROCESS_DOMAIN, "_Species", sep = ""),
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
    paste("Top", nrow(top_phyla), "by total representation"),
    paste(PROCESS_DOMAIN, "color mapping")
  ),
  stringsAsFactors = FALSE
)

write.table(summary_stats, paste0("alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== 16S", PROCESS_DOMAIN, "Percentage Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  -", output_png, "\n")
cat("  -", output_pdf, "\n")
cat(paste("  - alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_flow_annotations.tsv\n", sep = ""))
cat(paste("  - alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_node_descriptions.tsv\n", sep = ""))
cat(paste("  - alluvial_16s_", tolower(PROCESS_DOMAIN), "_pct_summary.tsv\n", sep = ""))

cat("\n16S", PROCESS_DOMAIN, "percentage alluvial plot generated with:\n")
cat("  - Clean merged data approach\n")
cat("  - Percentage normalization (0-100%)\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional", PROCESS_DOMAIN, "color scheme\n")
cat("  - Thin nodes and elegant flows\n")
cat("  - Legend positioned far right\n")
cat("  - Node titles moved down\n")
cat("  - Detailed flow annotations (TSV format)\n")
cat("  - Node descriptions and summary statistics\n\n")

cat("Percentage plot advantages:\n")
cat("  - Better comparison across different data types\n")
cat("  - Normalized scale (0-100%) for all nodes\n")
cat("  - Reveals relative proportions more clearly\n")
cat("  - Reduces dominance effects of large absolute numbers\n\n")

cat("To generate archaea plot, change PROCESS_DOMAIN to 'Archaea' at line 21\n")
