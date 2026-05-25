#!/usr/bin/env Rscript

# ============================================================================
# 18S Eukaryotic STACKED Alluvial Plot - ABSOLUTE + PERCENTAGE
# ============================================================================
# This script creates a stacked alluvial plot for 18S eukaryotic data
# by combining ABSOLUTE value plot (top) and PERCENTAGE plot (bottom)
# vertically with a single shared legend.
#
# Data Flow: NCBI Eukaryota Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Eukaryota Species
#
# Key Features:
# - Absolute values plot on TOP
# - Percentage values plot on BOTTOM
# - Single shared legend for both plots
# - Uses patchwork for plot composition
# - Professional eukaryotic color scheme
# ============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)
library(patchwork)  # For combining plots
library(ggrepel)    # For non-overlapping stratum labels

# Hide labels for strata below this percentage of their column total.
# Set to 0 to keep all labels. Tune up to reduce label overlap further.
LABEL_MIN_PCT <- 3.0

# Anchor stratum labels at their true y-centers so the visible top-to-bottom
# order of labels strictly matches the per-axis descending magnitude of the
# strata. ggrepel's y-direction repulsion otherwise swaps adjacent labels.
# Mirrors decreasing=FALSE so the stat re-computes the same per-axis ordering
# that geom_stratum uses (otherwise the embedded StatStratum defaults to
# factor-level order and produces a different y-center for each label).
geom_text_repel <- function(..., direction = NULL, max.overlaps = NULL,
                            force = NULL, force_pull = NULL,
                            box.padding = NULL, point.padding = NULL,
                            min.segment.length = NULL, segment.color = NULL,
                            segment.size = NULL, segment.alpha = NULL,
                            seed = NULL, max.time = NULL, max.iter = NULL,
                            xlim = NULL, ylim = NULL, arrow = NULL,
                            verbose = NULL) {
  args <- list(...)
  if (identical(args$stat, "stratum") && is.null(args$decreasing)) {
    args$decreasing <- FALSE
  }
  do.call(ggplot2::geom_text, args)
}

# Load shared color mapping functions
source("../../../shared_config/color_mapping_functions.R")

cat("=== 18S Eukaryotic STACKED Alluvial Plot (Absolute + Percentage) ===\n")
cat("Creating stacked visualization: Absolute (top) + Percentage (bottom)\n")
cat("Single shared legend for both plots\n\n")

# Define relative paths from alluvial script location to data directories
merged_data_path <- "../../../../final_merger/outputs/18s_ncbi_merged_division.csv"
census_data_path <- "../../../../eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_division.csv"

# Load merged 18S data
cat("Loading 18S merged data...\n")
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 18S merged data file at: ", merged_data_path)
}
data_18s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
if ("division" %in% colnames(data_18s) && !"phylum" %in% colnames(data_18s)) {
  data_18s <- data_18s %>% rename(phylum = division)
}

# Load census division data for .U. entries
if (!file.exists(census_data_path)) {
  stop("ERROR: Cannot find 18S census data file at: ", census_data_path)
}
census_division_data <- read.csv(census_data_path, stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("  - 18S merged entries:", nrow(data_18s), "\n")
cat("  - Census division entries:", nrow(census_division_data), "\n\n")

# Process eukaryotic data (filter for Eukaryota domain only)
process_eukaryotic_data <- function() {
  matched_data <- data_18s %>%
    filter(match_status == 'matched') %>%
    filter(!is.na(census_size_count) & !is.na(census_otu_count) &
           !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
    filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
    filter(domain == "Eukaryota")
  
  return(matched_data)
}

# Load shared color configuration
cat("Loading shared color configuration...\n")
color_config <- load_taxonomic_colors("../../../shared_config/taxonomic_color_mapping.yaml")

# Extract eukaryotic .U. entries from census data
get_eukaryotic_u_entries <- function(census_data) {
  u_entries <- census_data %>%
    filter(grepl("\\.U\\.", Name_to_use)) %>%
    filter(otu_count >= 10) %>%
    filter(grepl("Eukaryota", Name_to_use) |
           (!grepl("Bacteria", Name_to_use) & !grepl("Archaea", Name_to_use))) %>%
    select(phylum = Name_to_use, census_otu_count = otu_count, census_size_count = size_count) %>%
    mutate(
      ncbi_genome_count = 0,
      ncbi_species_count = 0,
      domain = "Eukaryota",
      match_status = "census_only"
    )
  
  return(u_entries)
}

# Process data
matched_data <- process_eukaryotic_data()
u_entries <- get_eukaryotic_u_entries(census_division_data)

cat("Processed eukaryotic data:\n")
cat("  - Matched entries:", nrow(matched_data), "\n")
cat("  - .U. entries:", nrow(u_entries), "\n\n")

# Combine matched data with .U. entries
combined_data <- bind_rows(matched_data, u_entries)

# Calculate totals for percentage calculations
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)
total_otu_count <- sum(combined_data$census_otu_count, na.rm = TRUE)
total_size_count <- sum(combined_data$census_size_count, na.rm = TRUE)

cat("Total counts (Eukaryota-only):\n")
cat("  - Genomes:", scales::comma(total_genome_count), "\n")
cat("  - Species:", scales::comma(total_species_count), "\n")
cat("  - OTUs:", scales::comma(total_otu_count), "\n")
cat("  - Sequences:", scales::comma(total_size_count), "\n\n")

# Calculate percentages
combined_data <- combined_data %>%
  mutate(
    genome_pct = (ncbi_genome_count / total_genome_count) * 100,
    species_pct = (ncbi_species_count / total_species_count) * 100,
    size_pct = (census_size_count / total_size_count) * 100,
    otu_pct = (census_otu_count / total_otu_count) * 100,
    total_representation = genome_pct + species_pct + size_pct + otu_pct
  )

# Select top 8 eukaryotic divisions by total representation
n_top <- 8
top_phyla <- combined_data %>%
  arrange(desc(total_representation)) %>%
  head(n_top)

cat("Top", n_top, "eukaryotic divisions selected for visualization\n\n")

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
  long_data <- rbind(long_data, data.frame(
    phylum = top_phyla$phylum[i],
    x = "Genbank_Genome_%",
    percentage = top_phyla$genome_pct[i],
    alluvium = i,
    stringsAsFactors = FALSE
  ))
  long_data <- rbind(long_data, data.frame(
    phylum = top_phyla$phylum[i],
    x = "IMG_Genome_%",
    percentage = top_phyla$size_pct[i],
    alluvium = i,
    stringsAsFactors = FALSE
  ))
  long_data <- rbind(long_data, data.frame(
    phylum = top_phyla$phylum[i],
    x = "18S_OTU_%",
    percentage = top_phyla$otu_pct[i],
    alluvium = i,
    stringsAsFactors = FALSE
  ))
  long_data <- rbind(long_data, data.frame(
    phylum = top_phyla$phylum[i],
    x = "Genbank_Species_%",
    percentage = top_phyla$species_pct[i],
    alluvium = i,
    stringsAsFactors = FALSE
  ))
}

# Add "Other" category
other_alluvium <- nrow(top_phyla) + 1
long_data <- rbind(long_data, data.frame(
  phylum = "Other",
  x = "Genbank_Genome_%",
  percentage = other_genome_pct,
  alluvium = other_alluvium,
  stringsAsFactors = FALSE
))
long_data <- rbind(long_data, data.frame(
  phylum = "Other",
  x = "IMG_Genome_%",
  percentage = other_size_pct,
  alluvium = other_alluvium,
  stringsAsFactors = FALSE
))
long_data <- rbind(long_data, data.frame(
  phylum = "Other",
  x = "18S_OTU_%",
  percentage = other_otu_pct,
  alluvium = other_alluvium,
  stringsAsFactors = FALSE
))
long_data <- rbind(long_data, data.frame(
  phylum = "Other",
  x = "Genbank_Species_%",
  percentage = other_species_pct,
  alluvium = other_alluvium,
  stringsAsFactors = FALSE
))

# Set factor levels for proper ordering
long_data$x <- factor(long_data$x, levels = c("Genbank_Genome_%", "IMG_Genome_%", "18S_OTU_%", "Genbank_Species_%"))

# Order phylum factor by total percentage (descending), keep Other last
phylum_totals <- long_data %>%
  group_by(phylum) %>%
  summarise(total = sum(percentage, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total))
ordered_phyla <- c(setdiff(phylum_totals$phylum, "Other"), "Other")
long_data$phylum <- factor(long_data$phylum, levels = ordered_phyla)

# Advanced preprocessing to preserve alluvium IDs
long_data_f <- long_data %>%
  arrange(x, desc(percentage)) %>%
  mutate(label_pct = ifelse(percentage >= LABEL_MIN_PCT,
                            paste0(phylum, "\n", sprintf("%.1f%%", percentage)),
                            ""))

# Create professional color palette using shared color system
phyla_names <- unique(long_data_f$phylum)
cat("Assigning colors for", length(phyla_names), "eukaryotic divisions...\n")

# Use shared color mapping system
colors <- get_domain_colors(phyla_names, "Eukaryota", color_config)

# Print color assignments for verification
print_color_summary(phyla_names, colors, "Eukaryota")

# Prepare absolute values data
cat("\nPreparing absolute values data...\n")
long_data_abs <- long_data_f %>%
  mutate(
    absolute_value = case_when(
      x == "Genbank_Genome_%" ~ percentage / 100 * total_genome_count,
      x == "IMG_Genome_%" ~ percentage / 100 * total_size_count,
      x == "18S_OTU_%" ~ percentage / 100 * total_otu_count,
      x == "Genbank_Species_%" ~ percentage / 100 * total_species_count,
      TRUE ~ 0
    ),
    label_abs = ifelse(percentage >= LABEL_MIN_PCT,
                       paste0(phylum, "\n", scales::comma(round(absolute_value))),
                       "")
  )

# Create STACKED alluvial plots with shared legend
cat("\nCreating stacked alluvial plots (Absolute + Percentage)...\n")

# Plot 1: ABSOLUTE VALUES (TOP)
p_abs <- ggplot(
  long_data_abs,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_value, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "Genbank_Genome_%"),
    stat = "stratum",
    aes(label = label_abs),
    hjust = 0, nudge_x = 0.16,
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "Genbank_Species_%"),
    stat = "stratum",
    aes(label = label_abs),
    hjust = 1, nudge_x = -0.16,
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "IMG_Genome_%"),
    stat = "stratum",
    aes(label = label_abs),
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "18S_OTU_%"),
    stat = "stratum",
    aes(label = label_abs),
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  scale_fill_manual(values = colors, name = "Division") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "18S_OTU_%" = "18S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(1 - 0.09, 4 + 0.09), clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 22, face = "bold", angle = 0, hjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 22, face = "bold", colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold", colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold", colour = "black"),
    legend.text = element_text(size = 16, face = "bold", colour = "black"),
    legend.key.size = unit(1.8, "cm"),
    legend.margin = margin(l = 10),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, colour = "black"),
    plot.margin = margin(10, 20, 5, 20)
  ) +
  labs(
    title = "Absolute Values",
    y = "Count",
    x = NULL
  )

# Plot 2: PERCENTAGE VALUES (BOTTOM)
p_pct <- ggplot(
  long_data_f,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, decreasing = FALSE, width = 0.18,
                lode.guidance = "forward", knot.pos = 0.35) +
  geom_stratum(alpha = 0.95, decreasing = FALSE, color = "white",
               linewidth = 0.2, width = 0.02) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "Genbank_Genome_%"),
    stat = "stratum",
    aes(label = label_pct),
    hjust = 0, nudge_x = 0.16,
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "Genbank_Species_%"),
    stat = "stratum",
    aes(label = label_pct),
    hjust = 1, nudge_x = -0.16,
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "IMG_Genome_%"),
    stat = "stratum",
    aes(label = label_pct),
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  geom_text_repel(
    data = function(d) d %>% filter(x == "18S_OTU_%"),
    stat = "stratum",
    aes(label = label_pct),
    size = 6.0, fontface = "bold", lineheight = 0.95,
    direction = "y", max.overlaps = Inf, force = 0.5,
    box.padding = 0.25, min.segment.length = 0,
    segment.color = NA, seed = 42
  ) +
  scale_fill_manual(values = colors, name = "Division") +
  scale_x_discrete(
    expand = expansion(mult = 0, add = 0),
    labels = c("Genbank_Genome_%" = "Genbank\nGenomes",
               "IMG_Genome_%" = "IMG/M\nSequences",
               "18S_OTU_%" = "18S\nOTUs",
               "Genbank_Species_%" = "Genbank\nSpecies")
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(1 - 0.09, 4 + 0.09), clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 22, face = "bold", angle = 0, hjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 22, face = "bold", colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold", colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold", colour = "black"),
    legend.text = element_text(size = 16, face = "bold", colour = "black"),
    legend.key.size = unit(1.8, "cm"),
    legend.margin = margin(l = 10),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, colour = "black"),
    plot.margin = margin(5, 20, 10, 20)
  ) +
  labs(
    title = "Percentage Values",
    y = "Percentage (%)",
    x = NULL
  )

# Stack plots vertically with shared legend using patchwork
cat("Combining plots with patchwork...\n")
p_stacked <- p_abs / p_pct +
  plot_layout(guides = "collect") &
  theme(legend.box.spacing = unit(0.3, "cm"))

# Save the stacked plot
cat("\nSaving stacked alluvial plot...\n")
ggsave("../figures/alluvial_18s_stacked_abs_pct.png", p_stacked, width = 24, height = 16, dpi = 300, bg = "white")
ggsave("../figures/alluvial_18s_stacked_abs_pct.pdf", p_stacked, width = 24, height = 16, dpi = 300, bg = "white")

cat("\n=== 18S Eukaryotic STACKED Alluvial Plot Created Successfully ===\n")
cat("Files saved:\n")
cat("  - ../figures/alluvial_18s_stacked_abs_pct.png\n")
cat("  - ../figures/alluvial_18s_stacked_abs_pct.pdf\n")
cat("\n18S eukaryotic stacked alluvial plot generated with:\n")
cat("  - Absolute values plot (TOP)\n")
cat("  - Percentage values plot (BOTTOM)\n")
cat("  - Single shared legend\n")
cat("  - Professional eukaryotic color scheme\n")
cat("  - Advanced alluvial aesthetics\n\n")

