#!/usr/bin/env Rscript
# ============================================================================
# 16S Bacterial Alluvial Plot - Configuration-Based Version
# ============================================================================
# Creates alluvial plots for 16S bacterial data using centralized configuration
# 
# Features:
# - Uses alluvial_filtering_config.yaml for all filtering logic
# - Supports multiple filtering strategies (top_abundance, threshold_based, hybrid)
# - Integrates with shared color mapping system
# - Clean, maintainable code with minimal hardcoded values
# - Generates both absolute and percentage-based plots
#
# Author: Config-Based Alluvial Team
# Date: 2025-01-17
# ============================================================================

# Load required libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(scales)
library(tidyr)
library(yaml)

# Load configuration and helper functions
source("../config/alluvial_filtering_functions.R")
source("../../shared_config/color_mapping_functions.R")

cat("=== 16S Bacterial Alluvial Plot (Configuration-Based) ===\n")
cat("Using centralized configuration for clean, maintainable code\n\n")

# Process data using configuration
cat("Processing 16S bacterial data...\n")
result <- process_alluvial_data("bacteria_16s")

processed_data <- result$data
config_info <- result$config
summary_info <- result$summary

# Load shared color configuration
cat("Loading shared color configuration...\n")
color_config <- load_taxonomic_colors()

# Prepare data for alluvial plot
cat("\nPreparing data for 4-node alluvial plot...\n")
long_data <- data.frame()

# Create long format data
for (i in 1:nrow(processed_data)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 4),
    phylum = rep(paste0(i, ". ", processed_data$phylum[i]), 4),
    x = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    stratum = c("NCBI_Total_Genomes", "16S_EukCensus_Sequences", "16S_EukCensus_OTUs", "NCBI_Total_Species"),
    absolute_count = c(
      processed_data$ncbi_genome_count[i],        # Node 1: NCBI Total Genomes
      processed_data$census_size_count[i],        # Node 2: 16S EukCensus Sequences
      processed_data$census_otu_count[i],         # Node 3: 16S EukCensus OTUs
      processed_data$ncbi_species_count[i]        # Node 4: NCBI Total Species
    ),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

# Clean and prepare data for plotting
long_data_f <- long_data %>%
  filter(absolute_count > 0) %>%
  group_by(phylum, x) %>%
  summarise(absolute_count = sum(absolute_count, na.rm = TRUE), .groups = 'drop') %>%
  ungroup()

# Order phyla by size at first axis for better stacking
first_axis <- sort(unique(long_data_f$x))[1]
sizes_first <- long_data_f %>%
  filter(x == first_axis) %>%
  arrange(desc(absolute_count)) %>%
  select(phylum) %>% pull()

long_data_f <- long_data_f %>%
  mutate(phylum = factor(phylum, levels = unique(c(sizes_first, setdiff(phylum, sizes_first)))))

# Get colors using shared color system
phyla_names <- unique(long_data_f$phylum)
cat("Assigning colors for", length(phyla_names), "bacterial phyla...\n")
colors <- get_domain_colors(phyla_names, "Bacteria", color_config)

# Print color assignments
print_color_summary(phyla_names, colors, "Bacteria")

# Create absolute values alluvial plot
cat("\nCreating absolute values alluvial plot...\n")
plot_abs <- ggplot(long_data_f, aes(x = x, stratum = phylum, alluvium = phylum, y = absolute_count, fill = phylum)) +
  geom_flow(alpha = config_info$global$plot_aesthetics$flow_alpha %||% 0.7, 
            color = "white", size = 0.2, curve_type = "xspline") +
  geom_stratum(alpha = 0.9, color = "white", size = 0.3, 
               width = config_info$global$plot_aesthetics$node_width %||% 0.1) +
  scale_fill_manual(values = colors, name = "Bacterial Phylum") +
  scale_y_continuous(labels = comma_format(), expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.1, 0)) +
  labs(
    title = "16S Bacterial Data Flow: NCBI Genomes → EukCensus Sequences → EukCensus OTUs → NCBI Species",
    subtitle = paste("Showing", length(phyla_names), "phyla using", 
                     config_info$domain$filtering$strategy, "filtering strategy"),
    x = NULL,
    y = "Absolute Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 20)),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save absolute values plot
output_abs_png <- "alluvial_16s_bacteria_config_based_abs.png"
output_abs_pdf <- "alluvial_16s_bacteria_config_based_abs.pdf"

plot_config <- config_info$global$plot_aesthetics %||% list(plot_width = 16, plot_height = 10, dpi = 300)

ggsave(output_abs_png, plot_abs, 
       width = plot_config$plot_width %||% 16, 
       height = plot_config$plot_height %||% 10,
       dpi = plot_config$dpi %||% 300, 
       bg = "white")

ggsave(output_abs_pdf, plot_abs,
       width = plot_config$plot_width %||% 16, 
       height = plot_config$plot_height %||% 10,
       bg = "white")

cat(paste("Absolute values plot saved as:", output_abs_png, "and", output_abs_pdf, "\n"))

# Create percentage-based plot
cat("Creating percentage-based alluvial plot...\n")

# Calculate percentages for each axis
long_data_pct <- long_data_f %>%
  group_by(x) %>%
  mutate(percentage = (absolute_count / sum(absolute_count)) * 100) %>%
  ungroup()

plot_pct <- ggplot(long_data_pct, aes(x = x, stratum = phylum, alluvium = phylum, y = percentage, fill = phylum)) +
  geom_flow(alpha = 0.7, color = "white", size = 0.2, curve_type = "xspline") +
  geom_stratum(alpha = 0.9, color = "white", size = 0.3, width = 0.1) +
  scale_fill_manual(values = colors, name = "Bacterial Phylum") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.1, 0)) +
  labs(
    title = "16S Bacterial Data Flow: Percentage Distribution",
    subtitle = paste("Showing", length(phyla_names), "phyla using", 
                     config_info$domain$filtering$strategy, "filtering strategy"),
    x = NULL,
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 20)),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save percentage plot
output_pct_png <- "alluvial_16s_bacteria_config_based_pct.png"
output_pct_pdf <- "alluvial_16s_bacteria_config_based_pct.pdf"

ggsave(output_pct_png, plot_pct,
       width = plot_config$plot_width %||% 16,
       height = plot_config$plot_height %||% 10,
       dpi = plot_config$dpi %||% 300,
       bg = "white")

ggsave(output_pct_pdf, plot_pct,
       width = plot_config$plot_width %||% 16,
       height = plot_config$plot_height %||% 10,
       bg = "white")

cat(paste("Percentage plot saved as:", output_pct_png, "and", output_pct_pdf, "\n"))

# Create summary CSV
cat("Creating summary CSV...\n")
summary_data <- data.frame(
  Dataset = "16S_Bacterial_Config",
  Filtering_Strategy = config_info$domain$filtering$strategy,
  Total_Phyla = summary_info$total_taxa,
  Selected_Phyla = summary_info$total_taxa - ifelse(summary_info$has_other, 1, 0),
  Has_Other_Category = summary_info$has_other,
  Total_Genomes = summary_info$total_genomes,
  Total_Sequences = summary_info$total_sequences,
  Total_OTUs = summary_info$total_otus,
  Total_Species = summary_info$total_species,
  stringsAsFactors = FALSE
)

summary_file <- "alluvial_16s_bacteria_config_based_summary.csv"
write.csv(summary_data, summary_file, row.names = FALSE)
cat(paste("Summary saved as:", summary_file, "\n"))

# Create node annotations CSV
cat("Creating node annotations CSV...\n")
node_annotations <- long_data_f %>%
  mutate(percentage = round((absolute_count / sum(absolute_count)) * 100, 2)) %>%
  select(node = x, phylum, absolute_count, percentage) %>%
  arrange(node, desc(absolute_count))

annotations_file <- "alluvial_16s_bacteria_config_based_node_annotations.csv"
write.csv(node_annotations, annotations_file, row.names = FALSE)
cat(paste("Node annotations saved as:", annotations_file, "\n"))

cat("\n✅ 16S Bacterial alluvial plots completed successfully!\n")
cat("📊 Configuration-based approach used for clean, maintainable code\n")
cat("🎨 Colors synchronized with comprehensive scatter plot system\n\n")
