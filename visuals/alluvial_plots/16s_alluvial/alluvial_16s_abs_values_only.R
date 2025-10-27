#!/usr/bin/env Rscript
#
# 16S Prokaryotic Alluvial Plot Generator - BACTERIA & ARCHAEA
# ============================================================
#
# Creates alluvial plots showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#
# Generates separate plots for both Bacteria and Archaea domains
# Shows top phyla by total representation for each domain
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

cat("=== 16S rRNA Prokaryotic Alluvial Plot Generator (Bacteria & Archaea) ===\n\n")

# Configuration: Set which domain to process
PROCESS_DOMAIN <- "Bacteria"  # Change to "Archaea" for archaea plot
cat(paste("Processing domain:", PROCESS_DOMAIN, "\n\n"))

# Robust path handling - works from script location
# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
census_data_path <- "../../../16S_censusparse/csv_16S/eukcensus16S_by_division.csv"

# Load the 16S rRNA merged data
cat("Loading 16S rRNA merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 16S merged data file at: ", merged_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
cat("16S rRNA data loaded:", nrow(data_16s), "rows\n")

# Load 16S division data for .U. entries (bacterial domain)
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 16S census data file at: ", census_data_path,
       "\nPlease ensure you're running this script from the correct directory.")
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

# Function to process domain-specific data
process_domain_data <- function(domain_name) {
  cat(paste("\n=== Processing", domain_name, "Domain ===\n"))

  # Filter for matched phyla for specific domain
  matched_data <- data_16s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == domain_name)

  cat("Matched", tolower(domain_name), "phyla found:", nrow(matched_data), "\n")

  return(matched_data)
}

# Process the configured domain
matched_data <- process_domain_data(PROCESS_DOMAIN)

# Extract domain-specific .U. entries from 16S census data
get_domain_u_entries <- function(census_data, domain_name) {
  if (domain_name == "Bacteria") {
    # Filter for bacterial .U. entries
    u_entries <- census_data %>%
      filter(grepl("\\.U\\.", Name_to_use)) %>%
      filter(otu_count >= 50) %>%  # Only significant .U. entries
      # Filter for bacterial .U. entries (exclude eukaryotic/archaeal .U.)
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
  } else if (domain_name == "Archaea") {
    # Filter for archaea .U. entries
    u_entries <- census_data %>%
      filter(grepl("\\.U\\.", Name_to_use)) %>%
      filter(otu_count >= 10) %>%  # Lower threshold for archaea (less abundant)
      filter(grepl("Archaea", Name_to_use)) %>%
      select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
      mutate(
        ncbi_genome_count = 0,
        ncbi_species_count = 0,
        domain = "Archaea",
        match_status = "census_only"
      )
  } else {
    u_entries <- data.frame()
  }

  return(u_entries)
}

u_entries <- get_domain_u_entries(census_division_data, PROCESS_DOMAIN)
cat(paste(PROCESS_DOMAIN, ".U. entries found:", nrow(u_entries), "\n"))

# Calculate totals from merged file + census totals for proper percentages
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)
total_size_count <- sum(census_division_data$size_count, na.rm = TRUE)  # Use full census for sequences
total_otu_count <- sum(census_division_data$otu_count, na.rm = TRUE)    # Use full census for OTUs

cat(paste("Total counts (", PROCESS_DOMAIN, "matched + full census for sequences/OTUs):\n"))
cat(paste("  ", PROCESS_DOMAIN, "Genomes:", scales::comma(total_genome_count), "\n"))
cat("  16S Sequences:", scales::comma(total_size_count), "\n")
cat("  16S OTUs:", scales::comma(total_otu_count), "\n")
cat(paste("  ", PROCESS_DOMAIN, "Species:", scales::comma(total_species_count), "\n"))

# Calculate percentages for matched data
matched_data <- matched_data %>%
  mutate(
    otu_percentage = (census_otu_count / total_otu_count) * 100,
    size_percentage = (census_size_count / total_size_count) * 100,
    genome_pct_db = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100
  )

# Calculate percentages for .U. entries
if (nrow(u_entries) > 0) {
  u_entries <- u_entries %>%
    mutate(
      otu_percentage = (census_otu_count / total_otu_count) * 100,
      size_percentage = (census_size_count / total_size_count) * 100,
      genome_pct_db = 0,  # No genomes for .U. entries
      species_pct = 0     # No species for .U. entries
    )
}

# Combine matched domain data with domain .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Select top N entries by size_percentage (domain-appropriate)
n_top <- if (PROCESS_DOMAIN == "Bacteria") 12 else min(8, nrow(combined_data))
top_phyla <- combined_data %>%
  arrange(desc(size_percentage)) %>%
  head(n_top)

cat(paste("Top", n_top, "entries selected (", tolower(PROCESS_DOMAIN), "phyla +", tolower(PROCESS_DOMAIN), ".U.):", nrow(top_phyla), "\n"))

# Debug: Print all selected entries with their counts
cat("\n=== DEBUG: Top 12 entries selected (bacterial + .U.) ===\n")
for (i in 1:nrow(top_phyla)) {
  entry_type <- if(grepl("\\.U\\.", top_phyla$phylum[i])) "(.U.)" else "(matched)"
  cat(sprintf("%2d. %-25s %6s | Genomes: %8s | Species: %8s | OTUs: %8s | Sequences: %8s\n",
              i,
              top_phyla$phylum[i],
              entry_type,
              format(top_phyla$ncbi_genome_count[i], big.mark = ","),
              format(top_phyla$ncbi_species_count[i], big.mark = ","),
              format(top_phyla$census_otu_count[i], big.mark = ","),
              format(top_phyla$census_size_count[i], big.mark = ",")))
}
cat("=== END DEBUG ===\n\n")

# Calculate "Other" category from remaining entries not in top 12
other_data <- combined_data %>%
  filter(!phylum %in% top_phyla$phylum)

other_genome_count <- sum(other_data$ncbi_genome_count, na.rm = TRUE)
other_size_count <- sum(other_data$census_size_count, na.rm = TRUE)
other_otu_count <- sum(other_data$census_otu_count, na.rm = TRUE)
other_species_count <- sum(other_data$ncbi_species_count, na.rm = TRUE)

cat("\nOther category (remaining", nrow(other_data), "bacterial entries):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Create long format data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Add top phyla data (including .U. entries)
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    absolute_count = c(
      top_phyla$ncbi_genome_count[i],        # Node 1: NCBI Total Genomes
      top_phyla$census_size_count[i],        # Node 2: 16S EukCensus Sequences
      top_phyla$census_otu_count[i],         # Node 3: 16S EukCensus OTUs
      top_phyla$ncbi_species_count[i]        # Node 4: NCBI Total Species
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
  absolute_count = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
  stringsAsFactors = FALSE
)
long_data <- rbind(long_data, other_data)

# Fix x-axis ordering to prevent alphabetical reordering
long_data$x <- factor(long_data$x, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"))

cat("Long data created with", nrow(long_data), "rows\n")

# ADVANCED ALLUVIAL PREPROCESSING - Fix aesthetic issues
cat("Applying advanced alluvial preprocessing...\n")

# 1. Ensure each phylum has a row at every axis (fill missing with 0)
library(tidyr)
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

# Debug: Check the data structure and counts for NCBI nodes
cat("\n=== DEBUG: NCBI Node Data Check ===\n")
ncbi_genome_data <- long_data[long_data$x == "NCBI_Total_Genomes", ]
ncbi_species_data <- long_data[long_data$x == "NCBI_Total_Species", ]

cat("NCBI Genome Node flows:\n")
for (i in 1:nrow(ncbi_genome_data)) {
  if (ncbi_genome_data$absolute_count[i] > 0) {
    cat(sprintf("  %-25s: %8s genomes\n",
                ncbi_genome_data$phylum[i],
                format(ncbi_genome_data$absolute_count[i], big.mark = ",")))
  }
}

cat("\nNCBI Species Node flows:\n")
for (i in 1:nrow(ncbi_species_data)) {
  if (ncbi_species_data$absolute_count[i] > 0) {
    cat(sprintf("  %-25s: %8s species\n",
                ncbi_species_data$phylum[i],
                format(ncbi_species_data$absolute_count[i], big.mark = ",")))
  }
}
cat("=== END DEBUG ===\n\n")

# Master bacterial color palette from mega 16S script
get_bacterial_colors <- function() {
  bacteria_color_map <- c(
    # Top novel phyla (highest factors)
    "Abditibacteriota" = "#87bd4b",        # 36.71× novel - highest factor
    "Calditrichota" = "#e4a98f",           # 35.00× novel
    "Planctomycetota" = "#923052",         # 26.34× novel - unique cell wall-less
    "Thermomicrobiota" = "#6cd490",        # 25.55× novel
    "Acidobacteriota" = "#985739",         # 24.66× novel - important soil bacteria
    "Gemmatimonadota" = "#9d9719",         # 19.97× novel
    "Chloroflexota" = "#b2f3fd",           # 18.45× novel - green non-sulfur bacteria
    "Verrucomicrobiota" = "#8d0571",       # 16.85× novel - environmental bacteria
    "Armatimonadota" = "#e99953",          # 15.34× novel
    "Myxococcota" = "#fae692",             # 14.31× novel
    "Bdellovibrionota" = "#164c3f",        # 10.71× novel
    "Thermodesulfobacteriota" = "#01414d", # 6.88× novel
    "Spirochaetota" = "#dcb94b",           # 6.82× novel - spiral-shaped bacteria
    "Bacteroidota" = "#de724d",            # 4.25× novel - important gut bacteria
    # Top overrepresented phyla
    "Chlorobiota" = "#458f9c",             # 18.67× overrepresented
    "Chlamydiota" = "#5b1f06",             # 9.50× overrepresented - dark brown for obligate parasites
    "Ignavibacteriota" = "#f16127"         # 1063.64× overrepresented - blood orange for distinction
  )
  return(bacteria_color_map)
}

# Extended fallback colors for bacteria phyla not in main palette
get_extended_bacterial_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

# Excluded phyla that use fallback colors (same as mega 16S script)
get_excluded_phyla <- function() {
  c("Campylobacterota", "Pseudomonadota", "Mycoplasmatota", "Thermotogota", "Bacillota")
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

# Create professional color palette using bacterial colors from mega 16S script
phyla_names <- unique(long_data$phylum)
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
    # Special colors for .U. entries
    u_colors <- c("#ff9999", "#ffcc99", "#99ccff", "#ccffcc")
    u_index <- sum(grepl("\\.U\\.", phyla_names[1:i]))
    colors[i] <- u_colors[((u_index - 1) %% length(u_colors)) + 1]
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

cat("Colors assigned to", length(colors), "phyla\n")

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
  scale_fill_manual(values = colors, name = "Phylum") +
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
ggsave("alluvial_16s_abs_values_only.png", p_abs, width = 24, height = 10, dpi = 300, bg = "white")
ggsave("alluvial_16s_abs_values_only.pdf", p_abs, width = 24, height = 10, dpi = 300, bg = "white")

cat("📁 Files saved:\n")
cat("   - alluvial_16s_abs_values_only.png (high-res image)\n")
cat("   - alluvial_16s_abs_values_only.pdf (vector format)\n")

cat("\n🔍 FINAL CHECK: Phyla that should be visible in NCBI nodes:\n")
ncbi_visible <- top_phyla[top_phyla$ncbi_genome_count > 0, ]
cat(sprintf("   Expected %d phyla with genomes in NCBI nodes:\n", nrow(ncbi_visible)))
for (i in 1:nrow(ncbi_visible)) {
  cat(sprintf("   %2d. %-20s: %8s genomes\n",
              i, ncbi_visible$phylum[i],
              format(ncbi_visible$ncbi_genome_count[i], big.mark = ",")))
}

cat("\n📊 Summary (Bacteria NCBI data):\n")
cat("   - Total bacterial phyla:", nrow(matched_data), "\n")
cat("   - Top entries displayed:", nrow(top_phyla), "\n")
cat("   - Phyla with NCBI genomes:", nrow(ncbi_visible), "\n")
cat("   - Total Bacterial Genomes:", scales::comma(total_genome_count), "\n")
cat("   - Total 16S Sequences:", scales::comma(total_size_count), "\n")
cat("   - Total 16S OTUs:", scales::comma(total_otu_count), "\n")
cat("   - Total Bacterial Species:", scales::comma(total_species_count), "\n")

cat("\n⚠️  If you only see 2 phyla in the plot, check:\n")
cat("   1. File timestamps - make sure you're viewing the latest files\n")
cat("   2. Plot zoom level - smaller flows might be hard to see\n")
cat("   3. Color contrast - similar colors might blend together\n")

# Create CSV table with figure annotations in same order as alluvial plot
cat("\nCreating CSV table with figure data...\n")

# Prepare data for CSV - same order as appears in the alluvial plot
csv_data <- data.frame(
  Phylum = character(),
  NCBI_Total_Genomes = integer(),
  EukCensus_16S_Sequences = integer(),
  EukCensus_16S_OTUs = integer(),
  NCBI_Total_Species = integer(),
  stringsAsFactors = FALSE
)

# Add top phyla data in order
for (i in 1:nrow(top_phyla)) {
  csv_data <- rbind(csv_data, data.frame(
    Phylum = paste0(i, ". ", top_phyla$phylum[i]),
    NCBI_Total_Genomes = top_phyla$ncbi_genome_count[i],
    EukCensus_16S_Sequences = top_phyla$census_size_count[i],
    EukCensus_16S_OTUs = top_phyla$census_otu_count[i],
    NCBI_Total_Species = top_phyla$ncbi_species_count[i],
    stringsAsFactors = FALSE
  ))
}

# Add "Other" category
csv_data <- rbind(csv_data, data.frame(
  Phylum = "Other",
  NCBI_Total_Genomes = other_genome_count,
  EukCensus_16S_Sequences = other_size_count,
  EukCensus_16S_OTUs = other_otu_count,
  NCBI_Total_Species = other_species_count,
  stringsAsFactors = FALSE
))

# Add totals row
csv_data <- rbind(csv_data, data.frame(
  Phylum = "TOTAL",
  NCBI_Total_Genomes = total_genome_count,
  EukCensus_16S_Sequences = total_size_count,
  EukCensus_16S_OTUs = total_otu_count,
  NCBI_Total_Species = total_species_count,
  stringsAsFactors = FALSE
))

# Save CSV table
csv_filename <- "16s_bacterial_alluvial_data_table.csv"
write.csv(csv_data, csv_filename, row.names = FALSE)

cat("CSV table saved:", csv_filename, "\n")
cat("Table contains", nrow(csv_data), "rows (top 12 + Other + Total)\n")

# Print preview of CSV table
cat("\nCSV Table Preview:\n")
cat("Phylum                        | NCBI Genomes | 16S Sequences | 16S OTUs | NCBI Species\n")
cat("------------------------------|--------------|---------------|----------|-------------\n")
for (i in 1:min(5, nrow(csv_data))) {
  cat(sprintf("%-29s | %12s | %13s | %8s | %11s\n",
              csv_data$Phylum[i],
              format(csv_data$NCBI_Total_Genomes[i], big.mark = ","),
              format(csv_data$EukCensus_16S_Sequences[i], big.mark = ","),
              format(csv_data$EukCensus_16S_OTUs[i], big.mark = ","),
              format(csv_data$NCBI_Total_Species[i], big.mark = ",")))
}
if (nrow(csv_data) > 5) {
  cat("... (", nrow(csv_data) - 5, " more rows)\n")
}

cat("\n✅ 16S Bacterial alluvial plot generation complete!\n")
cat("✅ CSV data table created with figure annotations!\n")
