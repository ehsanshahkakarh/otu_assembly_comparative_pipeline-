#!/usr/bin/env Rscript

# Minimal test to verify Pseudomonadota appears in a simple alluvial plot

library(ggplot2)
library(ggalluvial)

# Create minimal test data with just 3 phyla
test_data <- data.frame(
  phylum = rep(c("1. Pseudomonadota", "2. Bacillota", "3. Other"), each = 4),
  x = rep(c("Node1", "Node2", "Node3", "Node4"), 3),
  percentage = c(
    # Pseudomonadota - should be DOMINANT and BLUE
    57.83, 28.68, 22.93, 39.53,
    # Bacillota - should be green
    23.65, 16.11, 21.31, 17.06,
    # Other - should be gray
    18.52, 55.21, 55.76, 43.41
  ),
  alluvium = rep(1:3, each = 4),
  stringsAsFactors = FALSE
)

# Set factor levels
test_data$x <- factor(test_data$x, levels = c("Node1", "Node2", "Node3", "Node4"))
test_data$phylum <- factor(test_data$phylum, levels = c("1. Pseudomonadota", "2. Bacillota", "3. Other"))

# Define colors
colors <- c(
  "1. Pseudomonadota" = "#1f77b4",  # BLUE
  "2. Bacillota" = "#228B22",        # GREEN
  "3. Other" = "#808080"             # GRAY
)

cat("=== MINIMAL TEST PLOT ===\n")
cat("Creating test plot with 3 phyla\n")
cat("Pseudomonadota should be BLUE and DOMINANT\n\n")

# Create plot
p <- ggplot(
  test_data,
  aes(x = x, stratum = phylum, alluvium = alluvium, y = percentage, fill = phylum)
) +
  geom_alluvium(alpha = 0.85, width = 0.18) +
  geom_stratum(alpha = 0.95, color = "white", linewidth = 0.2, width = 0.02) +
  scale_fill_manual(values = colors, name = "Phylum") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 16),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "MINIMAL TEST: Pseudomonadota Visibility Check",
    subtitle = "Blue flow = Pseudomonadota (should be dominant)",
    y = "Percentage (%)",
    x = NULL
  )

# Save
ggsave("TEST_pseudomonadota_minimal.png", p, width = 16, height = 8, dpi = 300, bg = "white")
ggsave("TEST_pseudomonadota_minimal.pdf", p, width = 16, height = 8, bg = "white")

cat("✅ Test plot saved as TEST_pseudomonadota_minimal.png\n")
cat("   If Pseudomonadota appears in this test plot but not in the main plot,\n")
cat("   then there's an issue with the main script's data processing.\n")

