#!/usr/bin/env Rscript

# Minimal test script to isolate Pseudomonadota visibility issue
# This script creates a simple alluvial plot with just a few phyla including Pseudomonadota

library(ggplot2)
library(ggalluvial)

cat("=== Testing Pseudomonadota Visibility ===\n")

# Create minimal test data with known values
test_data <- data.frame(
  phylum = rep(c("Pseudomonadota", "Bacillota", "Other"), each = 4),
  x = rep(c("Genbank_Genome_%", "IMG_Genome_%", "16S_OTU_%", "Genbank_Species_%"), 3),
  percentage = c(
    # Pseudomonadota (should be highly visible)
    57.83, 28.68, 22.93, 39.53,
    # Bacillota (for comparison)
    23.65, 16.11, 21.7, 17.0,
    # Other (small values)
    1.66, 8.27, 5.0, 3.0
  ),
  alluvium = rep(1:3, each = 4),
  stringsAsFactors = FALSE
)

# Simple color mapping
colors <- c(
  "Pseudomonadota" = "#1f77b4",  # Blue - should be very visible
  "Bacillota" = "#ff7f0e",       # Orange
  "Other" = "#808080"            # Gray
)

cat("Test data created:\n")
print(test_data)
cat("\nColors assigned:\n")
print(colors)

# Create simple alluvial plot
p <- ggplot(test_data, aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)) +
  geom_alluvium(alpha = 0.8, width = 0.2) +
  geom_stratum(alpha = 0.9, color = "white", linewidth = 0.3, width = 0.1) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0.1)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Test: Pseudomonadota Visibility",
    subtitle = "Blue flow should be dominant at nodes 1 and 4",
    x = "Data Source",
    y = "Percentage"
  )

# Save the test plot
output_file <- "test_pseudomonadota_visibility.png"
ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300, bg = "white")

cat("Test plot saved as:", output_file, "\n")
cat("Check if Pseudomonadota (blue) is visible at nodes 1 and 4\n")
cat("=== Test Complete ===\n")
