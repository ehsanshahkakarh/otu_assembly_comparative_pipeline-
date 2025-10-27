#!/usr/bin/env Rscript

# 18S EukCensus vs NCBI Species Count Scatter Plot Analysis - COMPREHENSIVE
# Optimized for all taxonomic levels (phylum, family, genus) with enhanced visualization

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
})

# Robust path handling - detect script location and set paths relative to it
script_dir <- dirname(normalizePath(ifelse(interactive(),
                                          file.path(getwd(), "scatter_18s_ncbi_comparison.R"),
                                          commandArgs(trailingOnly = FALSE)[4])))
cat(paste("Script directory:", script_dir, "\n"))

# Set working directory to script location for consistent path resolution
setwd(script_dir)
cat(paste("Working directory set to:", getwd(), "\n"))

# Comprehensive configuration for 18S analysis with robust paths
config <- list(
  data_dir = file.path("..", "Eukcensus_merge", "merged_output", "18s_merged", "results"),
  output_dir = file.path("final_visualizations"),
  census_data_dir = file.path("..", "..", "18S_censusparse", "csv_outputs"),
  plot_width = 40, plot_height = 18, dpi = 300,  # Extra wide dimensions for better layout
  top_n = 10, text_size = 8,  # Optimal text size for annotations
  size_range = c(8, 40)  # Large circle size range to fill white space
)

# Verify critical paths exist
if (!dir.exists(config$data_dir)) {
  stop(paste("Data directory not found:", config$data_dir))
}
if (!dir.exists(config$census_data_dir)) {
  stop(paste("Census data directory not found:", config$census_data_dir))
}
cat(paste("Data directory verified:", config$data_dir, "\n"))
cat(paste("Census data directory verified:", config$census_data_dir, "\n"))

# Enhanced phylum mapping for 18S data using lineage parsing from EukCensus data
get_18s_phylum_from_eukcensus <- function(level) {
  if (level == "phylum") return(NULL)

  # Use the original EukCensus parsed files which have lineage information
  file_path <- file.path(config$census_data_dir, paste0("eukcensus_18S_by_", level, ".csv"))
  if (!file.exists(file_path)) {
    cat(paste("Warning: EukCensus", level, "file not found:", file_path, "\n"))
    return(NULL)
  }

  cat(paste("Loading EukCensus lineage data from:", file_path, "\n"))
  data <- read.csv(file_path, stringsAsFactors = FALSE)

  # Check if we have the required columns
  if (!("lineage" %in% colnames(data)) || !("lineage_ranks" %in% colnames(data))) {
    cat("Warning: lineage or lineage_ranks columns not found in EukCensus data\n")
    return(NULL)
  }

  mapping <- setNames(
    sapply(1:nrow(data), function(i) {
      lineage_parts <- strsplit(data$lineage[i], ";")[[1]]
      rank_parts <- strsplit(data$lineage_ranks[i], ";")[[1]]
      phylum_idx <- which(rank_parts == "phylum")[1]
      if (!is.na(phylum_idx) && phylum_idx <= length(lineage_parts)) {
        phylum_name <- trimws(lineage_parts[phylum_idx])
        # Clean up any trailing characters or suffixes
        phylum_name <- gsub("_.*$", "", phylum_name)  # Remove _XX suffixes
        return(phylum_name)
      } else {
        return("Other")
      }
    }),
    data[[1]]  # First column should be the taxon name
  )

  cat(paste("Created phylum mapping for", length(mapping), level, "entries\n"))
  return(mapping)
}

# Comprehensive 18S phylum pattern matching function
get_18s_phylum_by_pattern <- function(taxon_name) {
  # Clean the taxon name by removing suffixes like _XXXX, _XX, _X
  taxon_clean <- gsub("_X+$", "", taxon_name)
  taxon_lower <- tolower(taxon_clean)

  # Streptophyta patterns (land plants) - check first as it's very specific
  if (grepl("^streptophyta", taxon_lower)) return("Streptophyta")
  if (grepl("plant|embryo|charophy|viridiplant", taxon_lower)) return("Streptophyta")

  # Specific fungal phyla (more specific patterns first)
  if (grepl("chytridiomyc|chytrid|neocallimastig|rhizophyd|polychytri", taxon_lower)) return("Chytridiomycota")
  if (grepl("mucoromyc|zoopagomyc", taxon_lower)) return("Mucoromycota")
  if (grepl("microsporidia", taxon_lower)) return("Microsporidia")
  if (grepl("cryptomyc", taxon_lower)) return("Cryptomycota")
  if (grepl("sanchytriomyc", taxon_lower)) return("Sanchytriomycota")
  if (grepl("olpidiomyc", taxon_lower)) return("Olpidiomycota")
  if (grepl("dacrymyc", taxon_lower)) return("Basidiomycota")
  if (grepl("schizosaccharomyc", taxon_lower)) return("Ascomycota")
  if (grepl("lipomyc", taxon_lower)) return("Ascomycota")

  # Specific animal phyla (vertebrates)
  if (grepl("chordata|vertebrat|teleostei|lepidosauria|amphibia", taxon_lower)) return("Chordata")
  if (grepl("arthropod|panarthropod|insect", taxon_lower)) return("Panarthropoda")
  if (grepl("mollusc", taxon_lower)) return("Mollusca")
  if (grepl("nematod", taxon_lower)) return("Nematoda")
  if (grepl("cnidari", taxon_lower)) return("Cnidaria")
  if (grepl("platyhelminth", taxon_lower)) return("Platyhelminthes")
  if (grepl("porifera|sponge", taxon_lower)) return("Porifera")
  if (grepl("nemertea", taxon_lower)) return("Nemertea")
  if (grepl("rotifera", taxon_lower)) return("Rotifera")
  if (grepl("tardigrada", taxon_lower)) return("Tardigrada")
  if (grepl("placozoa", taxon_lower)) return("Placozoa")
  if (grepl("bryozoa", taxon_lower)) return("Bryozoa")
  if (grepl("brachiopod", taxon_lower)) return("Brachiopoda")
  if (grepl("chaetognath", taxon_lower)) return("Chaetognatha")
  if (grepl("priapulid", taxon_lower)) return("Priapulida")
  if (grepl("phoronid", taxon_lower)) return("Phoronida")
  if (grepl("orthonectid|mesozoa|dicyemid", taxon_lower)) return("Mesozoa")

  # Chlorophyta patterns (green algae)
  if (grepl("chloro|trebouxio|volvoc|chlamydom|scenedes|selenastra", taxon_lower)) return("Chlorophyta")

  # Stramenopiles patterns (diatoms, brown algae, oomycetes)
  if (grepl("skeletonem|diatom|phaeophyc|oomyc|chrysophyc|bacillario|stram|bicoecea|limnofilidae", taxon_lower)) return("Stramenopiles")
  if (grepl("ochrophyt", taxon_lower)) return("Stramenopiles")

  # Alveolata patterns (ciliates, dinoflagellates, apicomplexans)
  if (grepl("parameci|oxytrich|dysterii|cerati|dinophyc|theileri|plasmod|stylocephal", taxon_lower)) return("Alveolata")
  if (grepl("ciliophora", taxon_lower)) return("Alveolata")
  if (grepl("perkinsozoa", taxon_lower)) return("Alveolata")

  # Rhizaria patterns (foraminifera, radiolaria, cercozoans)
  if (grepl("foram|radiol|cercoz|glaucoc|rhizar|endomyxa|ascetosporea", taxon_lower)) return("Rhizaria")
  if (grepl("retaria", taxon_lower)) return("Rhizaria")

  # Evosea patterns (amoebozoans)
  if (grepl("amoeb|dictyo|physar|evose", taxon_lower)) return("Evosea")

  # Discoba patterns (euglenids, kinetoplastids)
  if (grepl("euglena|trypano|leishman|kineto|discob", taxon_lower)) return("Discoba")

  # Discosea patterns
  if (grepl("discose", taxon_lower)) return("Discosea")

  # Tubulinea patterns
  if (grepl("tubulin", taxon_lower)) return("Tubulinea")

  # Metamonada patterns
  if (grepl("metamona|giardia|trichomonas", taxon_lower)) return("Metamonada")
  if (grepl("parabasalia", taxon_lower)) return("Parabasalia")
  if (grepl("preaxostyla", taxon_lower)) return("Preaxostyla")

  # Centroplasthelida patterns
  if (grepl("centroplast", taxon_lower)) return("Centroplasthelida")

  # Rhodophyta patterns (red algae)
  if (grepl("rhodophyt|red.algae|cyanidiales", taxon_lower)) return("Rhodophyta")

  # Hemimastigophora patterns
  if (grepl("hemimastig", taxon_lower)) return("Hemimastigophora")

  # Other specific phyla
  if (grepl("picozoa", taxon_lower)) return("Picozoa")
  if (grepl("nibbleridia", taxon_lower)) return("Nibbleridia")
  if (grepl("nebulidia", taxon_lower)) return("Nebulidia")

  # Ichthyophonida - a specific group
  if (grepl("ichthyophon", taxon_lower)) return("Ichthyophonida")

  # General fungal patterns (broader catch-all for fungi)
  if (grepl("saccharomyc|mycet|fungi", taxon_lower)) return("Ascomycota")

  # General animal patterns (broader catch-all for animals)
  if (grepl("animal|metazoa", taxon_lower)) return("Chordata")

  # Choanoflagellates
  if (grepl("choano", taxon_lower)) return("Choanoflagellata")

  # Default to Other for truly unmatched taxa
  return("Other")
}

# Create output directory
if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)

# Comprehensive data loading function for all 18S taxonomic levels
load_18s_data <- function(level) {
  filepath <- file.path(config$data_dir, paste0("18s_ncbi_merged_clean_", level, ".csv"))
  if (!file.exists(filepath)) {
    cat(paste("File not found:", filepath, "\n"))
    cat("Available files in data directory:\n")
    print(list.files(config$data_dir, pattern = "*.csv"))
    stop(paste("18S NCBI merged file not found:", filepath))
  }

  data <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(data) == 0) stop(paste("No 18S", level, "data found"))

  # Rename columns to match expected format
  colnames(data)[1] <- "Taxon"
  names(data)[names(data) == "census_otu_count"] <- "Census_OTU_Count"
  names(data)[names(data) == "ncbi_species_count"] <- "NCBI_Species_Count"
  names(data)[names(data) == "ncbi_genome_count"] <- "NCBI_Genome_Count"
  names(data)[names(data) == "isolate_count"] <- "Isolate_Count"

  data <- data %>% filter(Census_OTU_Count > 0, NCBI_Species_Count > 0)

  # Add phylum mapping for family/genus levels
  if (level != "phylum") {
    # First try to get phylum mapping from EukCensus lineage data
    phylum_map <- get_18s_phylum_from_eukcensus(level)
    if (!is.null(phylum_map)) {
      cat(paste("Using EukCensus lineage data for", level, "phylum assignment\n"))
      data$Phylum <- ifelse(data$Taxon %in% names(phylum_map), phylum_map[data$Taxon],
                           sapply(data$Taxon, get_18s_phylum_by_pattern))

      # Count how many were mapped vs pattern matched
      mapped_count <- sum(data$Taxon %in% names(phylum_map))
      pattern_count <- nrow(data) - mapped_count
      cat(paste("Mapped from lineage:", mapped_count, "| Pattern matched:", pattern_count, "\n"))
    } else {
      # Fallback to pattern matching
      cat(paste("Using pattern matching for", level, "phylum assignment\n"))
      data$Phylum <- sapply(data$Taxon, get_18s_phylum_by_pattern)
    }

    # Report phylum distribution
    phylum_counts <- table(data$Phylum)
    cat(paste("Phylum distribution in 18S", level, ":\n"))
    print(phylum_counts)
    cat(paste("Total unique phyla found:", length(unique(data$Phylum)), "\n"))
  }

  # Calculate key metrics
  data$Novelty_Ratio <- data$Census_OTU_Count / data$NCBI_Species_Count
  data$Coverage_Factor <- data$NCBI_Species_Count / data$Census_OTU_Count

  # Frederik Schulz suggestion: Circle size based on genome-to-isolate ratio
  # Larger circles = fewer isolates relative to genomes (higher ratio)
  # Handle division by zero: if no isolates, use maximum ratio
  data$Genome_Isolate_Ratio <- ifelse(data$Isolate_Count > 0,
                                      data$NCBI_Genome_Count / data$Isolate_Count,
                                      max(data$NCBI_Genome_Count / pmax(data$Isolate_Count, 1), na.rm = TRUE))

  # Large circle sizing to fill white space effectively
  # Use square root transformation for balanced size differences
  data$Circle_Size_Raw <- sqrt(data$Genome_Isolate_Ratio)

  # Normalize to a larger range (3-20) to fill more white space
  min_size <- min(data$Circle_Size_Raw, na.rm = TRUE)
  max_size <- max(data$Circle_Size_Raw, na.rm = TRUE)
  data$Circle_Size <- 3 + 17 * (data$Circle_Size_Raw - min_size) / (max_size - min_size)

  # Identify top taxa
  data$Is_Top_Novelty <- rank(-data$Novelty_Ratio) <= config$top_n
  data$Is_Top_Coverage <- rank(-data$Coverage_Factor) <= config$top_n

  return(data)
}

# Generate summary tables
generate_18s_summary <- function(data, type = "coverage") {
  col <- if (type == "coverage") "Coverage_Factor" else "Novelty_Ratio"
  data %>% arrange(desc(.data[[col]])) %>% head(config$top_n) %>%
    select(Taxon, Coverage_Factor, Novelty_Ratio, Census_OTU_Count, NCBI_Species_Count) %>%
    mutate(across(c(Coverage_Factor, Novelty_Ratio), ~ round(.x, 2)))
}

# Comprehensive 18S scatter plot creation function
create_18s_scatter_plot <- function(data, level) {
  # Base plot with reduced expansion for tighter layout
  p <- ggplot(data, aes(x = Census_OTU_Count, y = NCBI_Species_Count)) +
    scale_x_log10(labels = comma_format(),
                  expand = expansion(mult = c(0.1, 0.1))) +  # Reduced x-axis expansion for tighter layout
    scale_y_log10(labels = comma_format(),
                  expand = expansion(mult = c(0.1, 0.1)))    # Reduced y-axis expansion for tighter layout

  # Background points (grey) - large to fill white space
  bg_data <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                       color = "grey70", alpha = 0.6)
  }

  # Highlighted points with colors by phyla
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  if (nrow(top_data) > 0) {
    # Generate pastel colors for all categories
    get_pastel_colors <- function(n) {
      if (n <= 8) {
        # Use Set2 which has nice pastel colors
        return(RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n])
      } else if (n <= 12) {
        # Use Pastel1 and Pastel2 for more pastel options
        pastel1 <- RColorBrewer::brewer.pal(9, "Pastel1")
        pastel2 <- RColorBrewer::brewer.pal(8, "Pastel2")
        return(c(pastel1, pastel2)[1:n])
      } else {
        # For many categories, create custom pastel palette
        pastel_colors <- c(
          "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF",  # Light pastels
          "#E6B3FF", "#FFB3E6", "#B3FFB3", "#B3E6FF", "#FFE6B3",  # More pastels
          "#D4B3FF", "#B3FFD4", "#FFD4B3", "#B3D4FF", "#D4FFB3",  # Even more
          "#FFCCCB", "#E0BBE4", "#C7CEEA", "#B5EAD7", "#FFDAC1"   # Additional soft colors
        )
        return(pastel_colors[1:min(n, length(pastel_colors))])
      }
    }

    if (level == "phylum") {
      n_taxa <- length(unique(top_data$Taxon))
      colors <- get_pastel_colors(n_taxa)

      # Add black outline layer first
      p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                         color = "black", alpha = 0.9, stroke = 0.5)

      # Add colored fill layer on top
      p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Taxon)),
                         alpha = 0.9) +
        scale_color_manual(values = colors, name = "Top Taxa",
                          guide = guide_legend(override.aes = list(size = 12, alpha = 1, shape = 15),
                                             keywidth = unit(3.0, "cm"), keyheight = unit(2.5, "cm")))
    } else {
      # Show only phyla of the annotated (top) taxa in the legend
      top_phyla <- unique(top_data$Phylum[top_data$Phylum != "Other"])
      top_phyla <- sort(top_phyla)  # Sort alphabetically for consistent ordering
      n_phyla <- length(top_phyla)
      colors <- get_pastel_colors(n_phyla)

      # Add black outline layer first
      p <- p + geom_point(data = top_data, aes(size = Circle_Size),
                         color = "black", alpha = 0.9, stroke = 0.5)

      # Add colored fill layer on top
      p <- p + geom_point(data = top_data, aes(size = Circle_Size * 0.8, color = factor(Phylum, levels = top_phyla)),
                         alpha = 0.9) +
        scale_color_manual(values = setNames(colors, top_phyla), name = "Phyla",
                          guide = guide_legend(override.aes = list(size = 12, alpha = 1, shape = 15),
                                             keywidth = unit(3.0, "cm"), keyheight = unit(2.5, "cm")))
    }
  }

  # Hide size legend and add proportionality line
  p <- p +
    scale_size_continuous(range = config$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "darkblue", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "right",
          # Professional styling with larger, darker axes
          plot.title = element_text(size = 48, hjust = 0.5, face = "bold", color = "black"),
          axis.title = element_text(size = 46, face = "bold", color = "black"),
          axis.text = element_text(size = 32, face = "bold", color = "black"),
          # Legend styling - extra large to dominate right side
          legend.text = element_text(size = 48),      # Even more massive legend text
          legend.title = element_text(size = 56, face = "bold"),  # Huge legend title
          legend.key.size = unit(4.0, "cm"),         # Enormous legend keys (squares)
          legend.key.width = unit(4.5, "cm"),        # Extra wide legend keys
          legend.key.height = unit(3.5, "cm"),       # Extra tall legend keys
          legend.spacing.y = unit(3.0, "cm"),        # Massive spacing between items
          legend.spacing.x = unit(2.0, "cm"),        # More horizontal spacing
          legend.margin = margin(l = 80, r = 80, t = 80, b = 80),  # Huge margins
          # Reduced plot margins for tighter layout
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
    labs(title = paste("18S EukCensus vs NCBI:", tools::toTitleCase(level)),
         x = "18S EukCensus OTU Count", y = "NCBI Species Count",
         caption = paste(nrow(data), "taxa"))

  # Add superimposed circle size legend in bottom-left quadrant
  x_range <- range(data$Census_OTU_Count, na.rm = TRUE)
  y_range <- range(data$NCBI_Species_Count, na.rm = TRUE)

  # Position legend well within the plot area in bottom-left
  legend_x <- 10^(log10(x_range[1]) + 0.05 * diff(log10(x_range)))
  legend_y_bottom <- 10^(log10(y_range[1]) + 0.25 * diff(log10(y_range)))

  # Add legend title - much larger
  p <- p + annotate("text", x = legend_x, y = legend_y_bottom,
                   label = "18S rRNA Circle Size", hjust = 0, vjust = 0,
                   size = 7, fontface = "bold", color = "black")

  # Add subtitle - Frederik Schulz suggestion - much larger
  p <- p + annotate("text", x = legend_x, y = 10^(log10(legend_y_bottom) + 0.08 * diff(log10(y_range))),
                   label = "Genomes/Isolates Ratio", hjust = 0, vjust = 0,
                   size = 5.5, fontface = "italic", color = "black")

  # Add example circles with labels based on ACTUAL genome-to-isolate ratio ranges from the data
  # Use actual data quantiles for realistic scaling
  actual_ratios <- data$Genome_Isolate_Ratio[is.finite(data$Genome_Isolate_Ratio)]
  if (length(actual_ratios) > 0) {
    circle_ratios <- c(
      quantile(actual_ratios, 0.25, na.rm = TRUE),  # 25th percentile
      quantile(actual_ratios, 0.5, na.rm = TRUE),   # Median
      quantile(actual_ratios, 0.75, na.rm = TRUE)   # 75th percentile
    )
    circle_labels <- c(
      paste0("Many isolates (", round(circle_ratios[1], 1), "x)"),
      paste0("Moderate (", round(circle_ratios[2], 1), "x)"),
      paste0("Few isolates (", round(circle_ratios[3], 1), "x)")
    )
  } else {
    # Fallback if no valid ratios
    circle_ratios <- c(1, 10, 100)
    circle_labels <- c("Many isolates (1x)", "Moderate (10x)", "Few isolates (100x)")
  }

  for (i in seq_along(circle_ratios)) {
    # Position circles going DOWNWARD from the starting position
    y_pos <- 10^(log10(legend_y_bottom) - (i * 0.12 * diff(log10(y_range))))

    # Large scaling to match the main plot circles
    if (i == 1) {
      scaled_size <- 8   # Small circle
    } else if (i == 2) {
      scaled_size <- 20  # Medium circle
    } else {
      scaled_size <- 35  # Large circle
    }

    p <- p + annotate("point", x = legend_x, y = y_pos,
                     size = scaled_size, color = "darkgrey", alpha = 0.8, stroke = 1)

    # Add label with much larger text, positioned to the right
    label_x <- 10^(log10(legend_x) + 0.12 * diff(log10(x_range)))
    p <- p + annotate("text", x = label_x, y = y_pos,
                     label = circle_labels[i], hjust = 0, vjust = 0.5,
                     size = 6, color = "black", fontface = "bold")
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
        # Tighter ggrepel settings for closer annotations
        box.padding = 0.5,          # Reduced padding around text boxes
        point.padding = 0.3,        # Reduced padding around points
        force = 2,                  # Reduced repulsion force for closer labels
        force_pull = 3,             # Stronger pull toward points
        max.overlaps = Inf,         # Allow all labels to show
        min.segment.length = 0,     # Always show connector lines
        segment.color = "grey50",   # Visible connector lines
        segment.size = 0.3,         # Thinner connector lines
        segment.alpha = 0.7,        # Semi-transparent connectors
        nudge_x = 0.05,             # Smaller nudge for tighter layout
        nudge_y = 0.05,
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
        # Tighter ggrepel settings for closer annotations
        box.padding = 0.5,          # Reduced padding around text boxes
        point.padding = 0.3,        # Reduced padding around points
        force = 2,                  # Reduced repulsion force for closer labels
        force_pull = 3,             # Stronger pull toward points
        max.overlaps = Inf,         # Allow all labels to show
        min.segment.length = 0,     # Always show connector lines
        segment.color = "grey50",   # Visible connector lines
        segment.size = 0.3,         # Thinner connector lines
        segment.alpha = 0.7,        # Semi-transparent connectors
        nudge_x = -0.05,            # Smaller nudge in opposite direction
        nudge_y = -0.05,
        seed = 123                  # Different seed for different placement
      )
    }
  }

  return(p)
}

# Comprehensive main function for all 18S taxonomic levels
main <- function() {
  cat("18S EukCensus vs NCBI scatter plot analysis - COMPREHENSIVE\n")

  for (level in c("phylum", "family", "genus")) {
    cat(paste("Processing 18S", level, "...\n"))

    # Load, process, and plot
    data <- load_18s_data(level)
    plot <- create_18s_scatter_plot(data, level)

    # Save plot and summaries with 18S-specific naming
    ggsave(file.path(config$output_dir, paste0("18s_ncbi_scatter_", level, ".png")), plot,
           width = config$plot_width, height = config$plot_height, dpi = config$dpi, bg = "white")

    write.csv(generate_18s_summary(data, "coverage"),
              file.path(config$output_dir, paste0("18s_ncbi_top_coverage_", level, ".csv")), row.names = FALSE)
    write.csv(generate_18s_summary(data, "novelty"),
              file.path(config$output_dir, paste0("18s_ncbi_top_novelty_", level, ".csv")), row.names = FALSE)

    cat(paste("  Saved:", nrow(data), "taxa,", sum(data$Is_Top_Novelty), "novel,", sum(data$Is_Top_Coverage), "coverage\n"))
  }

  cat("18S comprehensive analysis complete!\n")
  cat("Generated files:\n")
  for (level in c("phylum", "family", "genus")) {
    cat(paste("  - 18s_ncbi_scatter_", level, ".png\n", sep = ""))
    cat(paste("  - 18s_ncbi_top_coverage_", level, ".csv\n", sep = ""))
    cat(paste("  - 18s_ncbi_top_novelty_", level, ".csv\n", sep = ""))
  }
}

# Run comprehensive 18S analysis
if (!interactive()) main()
