#!/usr/bin/env Rscript
#
# Static Sankey-Style Plot Generator
# ==================================
#
# Creates a static PNG version of the Sankey diagram using ggplot2
# Shows the flow from: NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
#

# Load required libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)

cat("=== Static Sankey-Style Plot Generator ===\n\n")

# Load the 16S data
cat("Loading 16S data...\n")
data_16s <- read.csv("merged_output/16s_ncbi_merged_phylum.csv", stringsAsFactors = FALSE)

# Filter for shared phyla only (in_both = 1)
shared_data <- data_16s %>%
  filter(in_both == 1) %>%
  filter(!is.na(Genome_Percentage) & !is.na(Size_Percentage) &
         !is.na(OTU_Percentage) & !is.na(Species_Coverage))

# Add total representation and get top 15
shared_data$total_representation <- shared_data$Genome_Percentage +
                                   shared_data$Size_Percentage +
                                   shared_data$OTU_Percentage +
                                   shared_data$Species_Coverage

top_phyla <- shared_data %>%
  arrange(desc(total_representation)) %>%
  head(15)

# Calculate totals and "Other"
total_genome_count <- sum(shared_data$NCBI_Genome_Count, na.rm = TRUE)
total_size_count <- sum(shared_data$Census_Size_Count, na.rm = TRUE)
total_otu_count <- sum(shared_data$Census_OTU_Count, na.rm = TRUE)
total_species_count <- sum(shared_data$NCBI_Species_Count, na.rm = TRUE)

other_genome_count <- total_genome_count - sum(top_phyla$NCBI_Genome_Count)
other_size_count <- total_size_count - sum(top_phyla$Census_Size_Count)
other_otu_count <- total_otu_count - sum(top_phyla$Census_OTU_Count)
other_species_count <- total_species_count - sum(top_phyla$NCBI_Species_Count)

# Create data for plotting
plot_data <- data.frame()

# Add top 15 phyla
for (i in 1:nrow(top_phyla)) {
  phylum_data <- data.frame(
    phylum = rep(paste0(i, ". ", substr(top_phyla$Phylum[i], 1, 15)), 4),
    stage = c("Genomes", "Sequences", "OTUs", "Species"),
    x = c(1, 2, 3, 4),
    value = c(
      top_phyla$NCBI_Genome_Count[i],
      top_phyla$Census_Size_Count[i],
      top_phyla$Census_OTU_Count[i],
      top_phyla$NCBI_Species_Count[i]
    ),
    percentage = c(
      top_phyla$Genome_Percentage[i],
      top_phyla$Size_Percentage[i],
      top_phyla$OTU_Percentage[i],
      top_phyla$Species_Coverage[i]
    ),
    stringsAsFactors = FALSE
  )
  plot_data <- rbind(plot_data, phylum_data)
}

# Add "Other" category
other_data <- data.frame(
  phylum = rep("Other", 4),
  stage = c("Genomes", "Sequences", "OTUs", "Species"),
  x = c(1, 2, 3, 4),
  value = c(other_genome_count, other_size_count, other_otu_count, other_species_count),
  percentage = c(
    (other_genome_count / total_genome_count) * 100,
    (other_size_count / total_size_count) * 100,
    (other_otu_count / total_otu_count) * 100,
    (other_species_count / total_species_count) * 100
  ),
  stringsAsFactors = FALSE
)
plot_data <- rbind(plot_data, other_data)

# Generate colors
phyla_names <- c(paste0(1:15, ". ", substr(top_phyla$Phylum, 1, 15)), "Other")
n_phyla <- length(phyla_names)

if (n_phyla <= 8) {
  colors <- brewer.pal(max(3, n_phyla), "Set2")
} else if (n_phyla <= 12) {
  colors <- brewer.pal(max(3, n_phyla), "Set3")
} else {
  base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
  colors <- base_colors[1:n_phyla]
}

# Set "Other" to gray
colors[length(colors)] <- "#CCCCCC"

# Create the plot
cat("Creating static Sankey-style plot...\n")

p <- ggplot(plot_data, aes(x = x, y = value, fill = phylum)) +
  geom_col(position = "stack", width = 0.6, color = "white", size = 0.5) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c("NCBI Total\nGenomes", "16S EukCensus\nSequences", 
                               "16S EukCensus\nOTUs", "NCBI Total\nSpecies")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Prokaryotic Taxonomic Flow Analysis: Static Sankey-Style View",
    subtitle = paste0("Top 15 shared prokaryotic phyla + Other • 4-Stage Data Flow:\n",
                     "Genomes → Sequences → OTUs → Species"),
    x = "Data Stage",
    y = "Absolute Count",
    caption = "Stacked bars show absolute counts • Colors represent phyla • Generated for static viewing"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1))

# Save the plot
output_file <- "static_sankey_16s_ncbi_abs.png"
ggsave(output_file, plot = p, width = 16, height = 10, dpi = 300, bg = "white")

cat("Static Sankey-style plot saved:", output_file, "\n")

# Create a horizontal version showing flows
cat("Creating horizontal flow version...\n")

p_horizontal <- ggplot(plot_data, aes(x = stage, y = value, fill = phylum)) +
  geom_col(position = "stack", width = 0.7, color = "white", size = 0.5) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Prokaryotic Data Flow: Horizontal Stacked View",
    subtitle = "Top 15 shared prokaryotic phyla + Other",
    x = "Data Type",
    y = "Absolute Count",
    caption = "Each bar represents 100% of that data type • Flow shows data reduction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1))

# Save horizontal version
horizontal_file <- "horizontal_flow_16s_ncbi_abs.png"
ggsave(horizontal_file, plot = p_horizontal, width = 14, height = 10, dpi = 300, bg = "white")

cat("Horizontal flow plot saved:", horizontal_file, "\n")

# Print summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("📊 STATIC SANKEY-STYLE PLOTS GENERATED\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("📁 OUTPUT FILES:\n")
cat("   •", output_file, "- Vertical stacked bars\n")
cat("   •", horizontal_file, "- Horizontal flow view\n\n")

cat("📈 DATA SUMMARY:\n")
cat("   • Total genomes:", scales::comma(total_genome_count), "\n")
cat("   • Total sequences:", scales::comma(total_size_count), "\n")
cat("   • Total OTUs:", scales::comma(total_otu_count), "\n")
cat("   • Total species:", scales::comma(total_species_count), "\n\n")

cat("✅ Static plots generated successfully!\n")
cat("   These PNG files can be viewed directly on the system.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
