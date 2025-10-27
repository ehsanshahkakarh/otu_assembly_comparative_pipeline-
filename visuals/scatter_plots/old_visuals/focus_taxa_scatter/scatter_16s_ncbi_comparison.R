#!/usr/bin/env Rscript

# 16S EukCensus vs NCBI Species Count Scatter Plot Analysis - COMPREHENSIVE
# Optimized for all taxonomic levels with bacteria and archaea analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(),
                                          file.path(getwd(), "scatter_16s_ncbi_comparison.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Comprehensive configuration for 16S analysis with robust paths
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "merged_output", "16s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  ncbi_data_dir = file.path("..", "..", "ncbi_parse", "csv_ncbi"),
  plot_width = 40, plot_height = 18, dpi = 300,  # Extra wide dimensions for better layout
  top_n = 10, text_size = 8,  # Optimal text size for annotations
  size_range = c(8, 35)  # Large circle size range for better visibility
)

# Verify critical paths exist
if (!dir.exists(config$data_dir)) {
  stop(paste("Data directory not found:", config$data_dir))
}
if (!dir.exists(config$ncbi_data_dir)) {
  stop(paste("NCBI data directory not found:", config$ncbi_data_dir))
}
cat(paste("Data directory verified:", config$data_dir, "\n"))
cat(paste("NCBI data directory verified:", config$ncbi_data_dir, "\n"))

# Fixed phylum mapping with robust path handling
get_phylum_mapping <- function(level) {
  if (level == "phylum") return(NULL)

  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) {
    cat(paste("Warning: NCBI", level, "file not found:", file_path, "\n"))
    return(NULL)
  }

  data <- read.csv(file_path, stringsAsFactors = FALSE)
  mapping <- setNames(
    sapply(1:nrow(data), function(i) {
      lineage_parts <- strsplit(data$lineage[i], ";")[[1]]
      rank_parts <- strsplit(data$lineage_ranks[i], ";")[[1]]
      phylum_idx <- which(rank_parts == "phylum")[1]
      if (!is.na(phylum_idx) && phylum_idx <= length(lineage_parts)) {
        trimws(lineage_parts[phylum_idx])
      } else {
        "Other"
      }
    }),
    data[[level]]
  )
  return(mapping)
}

# Function to identify archaea vs bacteria
is_archaea <- function(taxon_name) {
  archaea_patterns <- c(
    "archaeota", "archaeia", "archaea", "methanobacteriota", "thermoplasmatota",
    "euryarchaeota", "crenarchaeota", "thaumarchaeota", "korarchaeota",
    "nanoarchaeota", "promethearchaeota", "odinarchaeia", "heimdallarchaeia",
    "halobacteriota", "asgardarchaeota", "lokiarchaeia", "thorarchaeia"
  )
  return(any(sapply(archaea_patterns, function(pattern) grepl(pattern, tolower(taxon_name)))))
}

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Streamlined data processing with robust path handling
load_data <- function(level, domain = "bacteria") {
  filepath <- file.path(config$data_dir, paste0("16s_ncbi_merged_clean_", level, ".csv"))
  if (!file.exists(filepath)) {
    cat(paste("File not found:", filepath, "\n"))
    cat("Available files in data directory:\n")
    print(list.files(config$data_dir, pattern = "*.csv"))
    stop(paste("16S NCBI merged file not found:", filepath))
  }

  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop("No data found")

  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  names(data)[names(data) == "census_otu_count"] <- "Census_OTU_Count"
  names(data)[names(data) == "ncbi_species_count"] <- "NCBI_Species_Count"
  names(data)[names(data) == "ncbi_genome_count"] <- "NCBI_Genome_Count"
  names(data)[names(data) == "isolate_count"] <- "Isolate_Count"

  # Filter by domain first, then apply appropriate data filters
  if (domain == "archaea") {
    archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
    data <- data %>% filter(Taxon %in% archaea_taxa)
    cat(paste("Found", nrow(data), "total archaea taxa\n"))

    # For archaea, be more lenient - allow taxa with either Census OR NCBI data
    data <- data %>% filter(Census_OTU_Count > 0 | NCBI_Species_Count > 0)
    cat(paste("After lenient filtering:", nrow(data), "archaea taxa\n"))
  } else {
    # For bacteria, use the standard strict filter (both Census AND NCBI data)
    data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)
    cat(paste("Found", nrow(data), "total taxa with both Census and NCBI data\n"))

    # Then filter out archaea to keep only bacteria
    archaea_taxa <- data$Taxon[sapply(data$Taxon, is_archaea)]
    data <- data %>% filter(!Taxon %in% archaea_taxa)
    cat(paste("After removing archaea:", nrow(data), "bacteria taxa\n"))
  }

  if (nrow(data) == 0) stop(paste("No", domain, "data found after filtering"))

  # Add phylum mapping for family/genus
  if (level != "phylum") {
    phylum_map <- get_phylum_mapping(level)
    data$Phylum <- ifelse(data$Taxon %in% names(phylum_map), phylum_map[data$Taxon], "Other")
  }

  # Calculate key metrics with safe division
  data$Novelty_Ratio <- ifelse(data$NCBI_Species_Count > 0,
                               data$Census_OTU_Count / data$NCBI_Species_Count,
                               Inf)
  data$Coverage_Factor <- ifelse(data$Census_OTU_Count > 0,
                                data$NCBI_Species_Count / data$Census_OTU_Count,
                                Inf)

  # Frederik Schulz suggestion: Circle size based on genome-to-isolate ratio
  # Larger circles = fewer isolates relative to genomes (higher ratio)
  # Handle division by zero: if no isolates, use maximum ratio
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))

  # Scale the ratio for circle sizing (log transform to handle wide range)
  data$Circle_Size <- log10(data$Genome_Isolate_Ratio + 1)  # +1 to avoid log(0)

  # Identify top taxa
  data$Is_Top_Novelty <- rank(-data$Novelty_Ratio) <= config$top_n
  data$Is_Top_Coverage <- rank(-data$Coverage_Factor) <= config$top_n

  return(data)
}

# Generate summary tables
generate_summary <- function(data, type = "coverage") {
  col <- if (type == "coverage") "Coverage_Factor" else "Novelty_Ratio"
  data %>% arrange(desc(.data[[col]])) %>% head(config$top_n) %>%
    select(Taxon, Coverage_Factor, Novelty_Ratio, Census_OTU_Count, NCBI_Species_Count) %>%
    mutate(across(c(Coverage_Factor, Novelty_Ratio), ~ round(.x, 2)))
}

# Ultra-streamlined plotting with proper circle size legend
create_scatter_plot <- function(data, level, domain = "bacteria") {
  # Base plot with expanded x-axis to spread out clustered points horizontally
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  expand = expansion(mult = c(0.25, 0.25))) +  # 25% expansion on x-axis for horizontal spreading
    scale_y_log10(labels = comma_format(),
                  expand = expansion(mult = c(0.05, 0.05)))    # Minimal y-axis expansion

  # Background points (grey) - sized by isolate percentage
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "grey70", alpha = 0.6)
  }

  # Highlighted points with colors by phyla
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  if (nrow(top_data) > 0) {
    # Generate enough colors for all categories
    get_colors <- function(n) {
      if (n <= 3) return(c("#1f77b4", "#ff7f0e", "#2ca02c")[1:n])
      if (n <= 8) return(RColorBrewer::brewer.pal(n, "Set2"))
      if (n <= 12) return(RColorBrewer::brewer.pal(n, "Set3"))
      # For more than 12, use rainbow colors
      return(rainbow(n, alpha = 0.8))
    }

    if (level == "phylum") {
      n_taxa <- length(unique(top_data$Taxon))
      colors <- get_colors(n_taxa)

      # Add black outline layer first
      p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                         color = "black", alpha = 0.9, stroke = 0.5)

      # Add colored fill layer on top
      p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Taxon)),
                         alpha = 0.8) +
        scale_color_manual(values = colors, name = "Top Taxa",
                          guide = guide_legend(override.aes = list(size = 8, alpha = 1, shape = 15),
                                             keywidth = unit(1.5, "cm"), keyheight = unit(1.2, "cm")))
    } else {
      # Show only phyla of the annotated (top) taxa in the legend
      top_phyla <- unique(top_data$Phylum[top_data$Phylum != "Other"])
      top_phyla <- sort(top_phyla)  # Sort alphabetically for consistent ordering
      n_phyla <- length(top_phyla)
      colors <- get_colors(n_phyla)

      # Add black outline layer first
      p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                         color = "black", alpha = 0.9, stroke = 0.5)

      # Add colored fill layer on top
      p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Phylum, levels = top_phyla)),
                         alpha = 0.8) +
        scale_color_manual(values = setNames(colors, top_phyla), name = "Phyla",
                          guide = guide_legend(override.aes = list(size = 8, alpha = 1, shape = 15),
                                             keywidth = unit(1.5, "cm"), keyheight = unit(1.2, "cm")))
    }
  }

  # Hide size legend and add proportionality line
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "darkblue", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "right",
          # Dramatically increased title and axis sizes
          plot.title = element_text(size = 48, hjust = 0.5, face = "bold", color = "#2c3e50"),
          axis.title = element_text(size = 42, face = "bold", color = "#34495e"),
          axis.text = element_text(size = 28, color = "#34495e"),
          # Legend styling - dramatically increased to fill right side
          legend.text = element_text(size = 36),      # Massive legend text to fill side
          legend.title = element_text(size = 42, face = "bold"),  # Massive legend title
          legend.key.size = unit(3.0, "cm"),         # Much larger legend keys (squares)
          legend.key.width = unit(3.5, "cm"),        # Wider legend keys
          legend.key.height = unit(2.5, "cm"),       # Taller legend keys
          legend.spacing.y = unit(2.0, "cm"),        # Much more spacing between items
          legend.spacing.x = unit(1.5, "cm"),        # More horizontal spacing
          legend.margin = margin(l = 50, r = 50, t = 50, b = 50)) +  # Much larger margins
    labs(title = paste("16S EukCensus vs NCBI:", tools::toTitleCase(domain), tools::toTitleCase(level)),
         x = "16S EukCensus OTU Count", y = "NCBI Species Count",
         caption = paste(nrow(data), domain, "taxa"))

  # Add superimposed circle size legend in top-left quadrant
  x_range <- range(data$Census_OTU_Count, na.rm = TRUE)
  y_range <- range(data$NCBI_Species_Count, na.rm = TRUE)

  # Position legend in top-left area (log scale coordinates)
  legend_x <- 10^(log10(x_range[1]) + 0.05 * diff(log10(x_range)))
  legend_y_top <- 10^(log10(y_range[2]) - 0.05 * diff(log10(y_range)))

  # Add legend title - much larger
  p <- p + annotate("text", x = legend_x, y = legend_y_top,
                   label = "Circle Size", hjust = 0, vjust = 1,
                   size = 7, fontface = "bold", color = "black")

  # Add subtitle - Frederik Schulz suggestion - much larger
  p <- p + annotate("text", x = legend_x, y = 10^(log10(legend_y_top) - 0.05 * diff(log10(y_range))),
                   label = "Genomes/Isolates Ratio", hjust = 0, vjust = 1,
                   size = 5.5, fontface = "italic", color = "black")

  # Add example circles with labels based on genome-to-isolate ratio ranges
  # Calculate representative ratios from the actual data
  ratio_range <- range(data$Genome_Isolate_Ratio, na.rm = TRUE)
  circle_ratios <- c(
    quantile(data$Genome_Isolate_Ratio, 0.25, na.rm = TRUE),  # Low ratio
    quantile(data$Genome_Isolate_Ratio, 0.5, na.rm = TRUE),   # Medium ratio
    quantile(data$Genome_Isolate_Ratio, 0.75, na.rm = TRUE)   # High ratio
  )
  circle_labels <- c(
    paste0("Few isolates (", round(circle_ratios[1], 1), "x)"),
    paste0("Moderate (", round(circle_ratios[2], 1), "x)"),
    paste0("Many isolates (", round(circle_ratios[3], 1), "x)")
  )

  for (i in seq_along(circle_ratios)) {
    y_pos <- 10^(log10(legend_y_top) - ((i + 1) * 0.08 * diff(log10(y_range))))

    # Add circle - scale the size to match the actual plot scaling using log-transformed ratios
    log_ratio <- log10(circle_ratios[i] + 1)  # Same transformation as in data
    circle_size_range <- range(data$Circle_Size, na.rm = TRUE)

    # Proper scaling to match ggplot2's scale_size_continuous behavior
    if (diff(circle_size_range) > 0) {
      normalized_size <- (log_ratio - circle_size_range[1]) / diff(circle_size_range)
    } else {
      normalized_size <- 0.5  # Default to middle if no variation
    }
    scaled_size <- normalized_size * (config$size_range[2] - config$size_range[1]) + config$size_range[1]

    # Ensure minimum size for visibility
    scaled_size <- max(scaled_size, config$size_range[1])

    p <- p + annotate("point", x = legend_x, y = y_pos,
                     size = scaled_size, color = "darkgrey", alpha = 0.8)

    # Add label with much larger text
    label_x <- 10^(log10(legend_x) + 0.05 * diff(log10(x_range)))
    p <- p + annotate("text", x = label_x, y = y_pos,
                     label = circle_labels[i], hjust = 0, vjust = 0.5,
                     size = 6, color = "black", fontface = "bold")  # Much larger legend text
  }

  # Add all text labels in black
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_data) > 0) {
    novelty_data <- top_data[top_data$Is_Top_Novelty, ]
    coverage_data <- top_data[top_data$Is_Top_Coverage & !top_data$Is_Top_Novelty, ]

    if (nrow(novelty_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = novelty_data,
        aes(label = paste0(Taxon, " (", round(Novelty_Ratio, 1), "x)")),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        # Enhanced ggrepel settings for better spacing
        box.padding = 1.5,          # Much more padding around text boxes
        point.padding = 1.0,        # More padding around points
        force = 5,                  # Stronger repulsion force
        force_pull = 2,             # Pull labels toward their points
        max.overlaps = Inf,         # Allow all labels to show
        min.segment.length = 0,     # Always show connector lines
        segment.color = "grey50",   # Visible connector lines
        segment.size = 0.5,         # Thicker connector lines
        segment.alpha = 0.8,        # Semi-transparent connectors
        nudge_x = 0.1,              # Slight nudge to spread labels
        nudge_y = 0.1,
        seed = 42                   # Reproducible label placement
      )
    }
    if (nrow(coverage_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = coverage_data,
        aes(label = paste0(Taxon, " (", round(Coverage_Factor, 1), "x)")),
        color = "black",
        size = config$text_size,
        fontface = "bold",
        # Enhanced ggrepel settings for better spacing
        box.padding = 1.5,          # Much more padding around text boxes
        point.padding = 1.0,        # More padding around points
        force = 5,                  # Stronger repulsion force
        force_pull = 2,             # Pull labels toward their points
        max.overlaps = Inf,         # Allow all labels to show
        min.segment.length = 0,     # Always show connector lines
        segment.color = "grey50",   # Visible connector lines
        segment.size = 0.5,         # Thicker connector lines
        segment.alpha = 0.8,        # Semi-transparent connectors
        nudge_x = -0.1,             # Slight nudge in opposite direction to novelty
        nudge_y = -0.1,
        seed = 123                  # Different seed for different placement
      )
    }
  }

  return(p)
}



# Comprehensive main function for bacteria and archaea
main <- function() {
  cat("16S EukCensus vs NCBI scatter plot analysis - COMPREHENSIVE\n")

  for (domain in c("bacteria", "archaea")) {
    cat(paste("\n=== Processing", toupper(domain), "===\n"))

    for (level in c("phylum", "family", "genus")) {
      cat(paste("Processing", domain, level, "...\n"))

      tryCatch({
        # Load, process, and plot
        data <- load_data(level, domain)
        plot <- create_scatter_plot(data, level, domain)

        # Save plot and summaries with domain-specific naming
        ggsave(file.path(config$output_dir, paste0("16s_ncbi_scatter_", domain, "_", level, ".png")), plot,
               width = config$plot_width, height = config$plot_height, dpi = config$dpi, bg = "white")

        write.csv(generate_summary(data, "coverage"),
                  file.path(config$output_dir, paste0("16s_ncbi_top_coverage_", domain, "_", level, ".csv")), row.names = FALSE)
        write.csv(generate_summary(data, "novelty"),
                  file.path(config$output_dir, paste0("16s_ncbi_top_novelty_", domain, "_", level, ".csv")), row.names = FALSE)

        cat(paste("  Saved:", nrow(data), domain, "taxa,", sum(data$Is_Top_Novelty), "novel,", sum(data$Is_Top_Coverage), "coverage\n"))
      }, error = function(e) {
        cat(paste("  Error processing", domain, level, ":", e$message, "\n"))
      })
    }
  }

  cat("\n16S comprehensive analysis complete!\n")
  cat("Generated files for both bacteria and archaea:\n")
  for (domain in c("bacteria", "archaea")) {
    for (level in c("phylum", "family", "genus")) {
      cat(paste("  - 16s_ncbi_scatter_", domain, "_", level, ".png\n", sep = ""))
      cat(paste("  - 16s_ncbi_top_coverage_", domain, "_", level, ".csv\n", sep = ""))
      cat(paste("  - 16s_ncbi_top_novelty_", domain, "_", level, ".csv\n", sep = ""))
    }
  }
}

# Run analysis
if (!interactive()) main()
