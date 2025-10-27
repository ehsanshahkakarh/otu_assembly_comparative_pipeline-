#!/usr/bin/env Rscript
# Debug script to test alluvial plot with simplified data
# This will help identify if the issue is with data processing or plot rendering

library(ggplot2)
library(ggalluvial)
library(dplyr)

cat("=== DEBUGGING ALLUVIAL PLOT ISSUE ===\n")

# Create simplified test data with the 10 phyla that should show in NCBI nodes
test_data <- data.frame(
  phylum = c("Opisthokonta", "Streptophyta", "Stramenopiles", "Alveolata", 
             "Chlorophyta", "Discoba", "Rhizaria", "Metamonada", "Evosea", "Tubulinea"),
  genomes = c(41803, 5885, 684, 570, 440, 328, 68, 66, 46, 2),
  species = c(17372, 2349, 330, 175, 223, 121, 15, 22, 31, 1),
  stringsAsFactors = FALSE
)

cat("Test data created with", nrow(test_data), "phyla:\n")
for (i in 1:nrow(test_data)) {
  cat(sprintf("  %s: %s genomes\n", test_data$phylum[i], format(test_data$genomes[i], big.mark = ",")))
}

# Create long format for alluvial plot (simplified 2-node version)
long_data <- data.frame()
for (i in 1:nrow(test_data)) {
  phylum_data <- data.frame(
    alluvium = rep(i, 2),
    phylum = rep(paste0(i, ". ", test_data$phylum[i]), 2),
    x = c("NCBI_Genomes", "NCBI_Species"),
    stratum = c("NCBI_Genomes", "NCBI_Species"),
    absolute_count = c(test_data$genomes[i], test_data$species[i]),
    stringsAsFactors = FALSE
  )
  long_data <- rbind(long_data, phylum_data)
}

cat("\nLong data created with", nrow(long_data), "rows\n")

# Check what's in the long data
cat("\nGenome node data:\n")
genome_data <- long_data[long_data$x == "NCBI_Genomes", ]
for (i in 1:nrow(genome_data)) {
  cat(sprintf("  %s: %s\n", genome_data$phylum[i], format(genome_data$absolute_count[i], big.mark = ",")))
}

# Create colors
colors <- rainbow(nrow(test_data))
names(colors) <- unique(long_data$phylum)

cat("\nColors assigned:", length(colors), "\n")

# Create simple alluvial plot
cat("\nCreating test alluvial plot...\n")
p_test <- ggplot(long_data, aes(x = x, stratum = phylum, alluvium = alluvium, y = absolute_count, fill = phylum)) +
  geom_alluvium(alpha = 0.8, decreasing = FALSE, width = 0.8) +
  geom_stratum(alpha = 0.9, decreasing = FALSE, color = "white", linewidth = 0.8, width = 0.8) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "DEBUG: Simple 2-Node Alluvial Test",
    y = "Count"
  )

# Save test plot
ggsave("debug_alluvial_test.png", plot = p_test, width = 16, height = 8, dpi = 300, bg = "white")
ggsave("debug_alluvial_test.pdf", plot = p_test, width = 16, height = 8, dpi = 300, bg = "white")

cat("Test plot saved as debug_alluvial_test.png and debug_alluvial_test.pdf\n")
cat("Check this plot to see if all 10 phyla are visible!\n")
