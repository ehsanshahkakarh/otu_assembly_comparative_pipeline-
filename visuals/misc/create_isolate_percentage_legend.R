#!/usr/bin/env Rscript
# Create Isolate Percentage Circle Size Legend
# Generated: 2025-01-20
# Purpose: Create standalone PNG legend showing inverted isolate percentage circle sizing

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Create legend data (vertical layout)
legend_data <- data.frame(
  isolate_percentage = c(0, 50, 100),
  circle_size = c(30, 9, 4),  # Updated larger sizes
  label = c("Completely\nmetagenomic",
            "Moderately\ncultured",
            "Fully\ncultured"),
  x = rep(1, 3),     # All circles in same column
  y = c(3, 2, 1),    # Vertical arrangement: 0% at top, 100% at bottom
  stringsAsFactors = FALSE
)

# Create the legend plot
legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
  # Add circles with different sizes
  geom_point(aes(size = circle_size),
             color = "black",
             fill = "black",
             shape = 21,  # Circles with slightly thicker outline
             stroke = 0.25,
             alpha = 0.8) +

  # Add percentage labels to the left of circles
  geom_text(aes(label = paste0(isolate_percentage, "%")),
            hjust = 1.5,
            size = 6,
            fontface = "bold",
            color = "black") +

  # Add description labels to the right of circles
  geom_text(aes(label = label),
            hjust = -0.2,
            size = 4.5,
            lineheight = 0.9,
            color = "black") +

  # Set manual size scale to match actual circle sizes
  scale_size_identity() +

  # Clean theme
  theme_void() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 20), color = "gray30"),
    plot.margin = margin(25, 25, 25, 25),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +

  # Set plot limits for vertical layout
  xlim(-0.5, 3) +
  ylim(0.5, 3.5) +

  # Add title and subtitle
  labs(
    title = "Circle Size Legend: Isolate Representation",
    subtitle = "Larger circles highlight poorly cultured taxa (cultivation gaps)"
  )

# Save the legend as PNG
output_file <- "isolate_percentage_circle_legend.png"

# Create high-quality PNG (adjusted for vertical layout)
ggsave(
  filename = output_file,
  plot = legend_plot,
  width = 8,
  height = 10,  # Taller for vertical layout
  dpi = 300,
  bg = "white"
)

cat("Legend saved as:", output_file, "\n")
cat("Legend shows inverted isolate percentage circle sizing (vertical layout):\n")
cat("  - 0% isolates = Size 30 (largest circle)\n")
cat("  - 50% isolates = Size 9 (medium circle)\n")
cat("  - 100% isolates = Size 4 (smallest circle)\n")
