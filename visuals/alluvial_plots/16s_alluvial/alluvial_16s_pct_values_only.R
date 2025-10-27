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

# CONFIGURATION: Change this to switch between domains
PROCESS_DOMAIN <- "Bacteria"  # Change to "Archaea" for archaea plot

cat("=== 16S", PROCESS_DOMAIN, "Alluvial Plot (Percentage Values) ===\n")
cat("Using clean merged data approach for reliable visualization\n\n")

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

# Get bacterial color mappings (from mega 16S script)
get_bacteria_colors <- function() {
  bacteria_color_map <- c(
    "Pseudomonadota" = "#1f77b4",    # Blue
    "Bacillota" = "#ff7f0e",         # Orange  
    "Actinomycetota" = "#2ca02c",    # Green
    "Bacteroidota" = "#d62728",      # Red
    "Cyanobacteriota" = "#9467bd",   # Purple
    "Chloroflexota" = "#b2f3fd",     # Light blue
    "Verrucomicrobiota" = "#8d0571", # Dark purple
    "Planctomycetota" = "#923052"    # Dark red
  )
  return(bacteria_color_map)
}

# Get archaea color mappings (from mega 16S script)
get_archaea_colors <- function() {
  archaea_color_map <- c(
    "Euryarchaeota" = "#55ee79",     # Light green
    "Nitrososphaerota" = "#ef9e17",  # Orange
    "Thermoproteota" = "#f61e5d",    # Pink/red
    "Nanoarchaeota" = "#4ad9d5"      # Cyan
  )
  return(archaea_color_map)
}

# Extended color palette for unmapped phyla
get_extended_colors <- function() {
  extended_colors <- c(
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5"
  )
  return(extended_colors)
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

# Add top phyla data
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", top_phyla$phylum[i]), 4),
    x = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"),
    stratum = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"),
    percentage = c(
      top_phyla$genome_pct[i],    # Node 1: NCBI Genome %
      top_phyla$size_pct[i],      # Node 2: 16S Sequence %
      top_phyla$otu_pct[i],       # Node 3: 16S OTU %
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
  x = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"),
  stratum = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"),
  percentage = c(other_genome_pct, other_size_pct, other_otu_pct, other_species_pct),
  stringsAsFactors = FALSE
)

long_data <- rbind(long_data, other_data)

# Fix x-axis ordering
long_data$x <- factor(long_data$x, levels = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"))
long_data$stratum <- factor(long_data$stratum, levels = c("NCBI_Genomes_%", "16S_Sequences_%", "16S_OTUs_%", "NCBI_Species_%"))

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

# Create professional color palette using domain-specific colors
phyla_names <- unique(long_data_f$phylum)
n_colors <- length(phyla_names)

# Get domain-specific color mappings
if (PROCESS_DOMAIN == "Bacteria") {
  domain_colors <- get_bacteria_colors()
  u_color <- "#ff9999"  # Light red for bacterial .U. entries
} else {
  domain_colors <- get_archaea_colors()
  u_color <- "#ffcc99"  # Light orange for archaea .U. entries
}

extended_colors <- get_extended_colors()

# Assign colors to phyla
colors <- character(n_colors)
names(colors) <- phyla_names

for (i in seq_along(phyla_names)) {
  phylum_name <- phyla_names[i]
  clean_phylum <- gsub("^\\d+\\. ", "", phylum_name)
  
  if (phylum_name == "Other") {
    colors[i] <- "#808080"  # Gray for Other
  } else if (grepl("\\.U\\.", clean_phylum)) {
    colors[i] <- u_color  # Domain-specific color for .U. entries
  } else if (clean_phylum %in% names(domain_colors)) {
    colors[i] <- domain_colors[clean_phylum]
  } else {
    # Use extended fallback colors for unmapped phyla
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
  }
}

cat("Color mapping complete for", length(colors), "phyla\n")

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
  scale_fill_manual(values = colors, name = "Phylum") +
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

# Create domain-specific output filename
domain_suffix <- tolower(PROCESS_DOMAIN)
output_png <- paste0("alluvial_16s_", domain_suffix, "_pct_values_only.png")
output_pdf <- paste0("alluvial_16s_", domain_suffix, "_pct_values_only.pdf")

# Save advanced percentage alluvial plot with optimized dimensions
ggsave(output_png, p_pct, width = 24, height = 10, dpi = 300, bg = "white")
ggsave(output_pdf, p_pct, width = 24, height = 10, dpi = 300, bg = "white")

cat("\n=== 16S", PROCESS_DOMAIN, "Percentage Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  -", output_png, "\n")
cat("  -", output_pdf, "\n")
cat("\n16S", PROCESS_DOMAIN, "percentage alluvial plot generated with:\n")
cat("  - Clean merged data approach\n")
cat("  - Percentage normalization (0-100%)\n")
cat("  - Advanced alluvial aesthetics\n")
cat("  - Optimized flow guidance\n")
cat("  - Professional", PROCESS_DOMAIN, "color scheme\n")
cat("  - Thin nodes and elegant flows\n")
cat("  - Legend positioned far right\n")
cat("  - Node titles moved down\n\n")

cat("Percentage plot advantages:\n")
cat("  - Better comparison across different data types\n")
cat("  - Normalized scale (0-100%) for all nodes\n")
cat("  - Reveals relative proportions more clearly\n")
cat("  - Reduces dominance effects of large absolute numbers\n\n")

cat("To generate archaea plot, change PROCESS_DOMAIN to 'Archaea' at line 21\n")
