#!/usr/bin/env Rscript
# ============================================================================
# Standalone Phyla Legend Generator
# ============================================================================
# This script generates a combined horizontal phyla legend showing all taxa
# that appear in the comprehensive scatter plots (top 10 novelty + top 10 overrepresentation)
#
# Output: PNG and PDF files of the combined phyla legend
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

cat("🎨 STANDALONE PHYLA LEGEND GENERATOR\n")
cat("====================================\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Simple path configuration
config <- list(
  data_dir_16s = "../../../Eukcensus_merge/16s_merged/csv_results",
  data_dir_18s = "../../../Eukcensus_merge/18s_merged/csv_results",
  output_dir = ".",  # Output to current directory (mega_scrip)
  ncbi_data_dir = "../../../ncbi_parse/csv_ncbi"
)

# Legend output configuration
LEGEND_CONFIG <- list(
  width = 20,        # Width in inches
  height = 4,        # Height in inches
  dpi = 300,         # Resolution
  output_png = "combined_phyla_legend.png",
  output_pdf = "combined_phyla_legend.pdf"
)

# ============================================================================
# LOAD SHARED COLOR CONFIGURATION
# ============================================================================

load_shared_color_config <- function() {
  config_path <- file.path("..", "..", "shared_config", "taxonomic_color_mapping.yaml")
  if (!file.exists(config_path)) {
    stop("Shared color config not found at: ", config_path)
  }
  yaml::read_yaml(config_path)
}

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

# Load 16S data function
load_16s_data <- function(level, domain) {
  filename <- paste0("16s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_16s, filename)

  if (!file.exists(filepath)) stop(paste("File not found:", filepath))

  data <- read.csv(filepath, stringsAsFactors = FALSE) %>%
    filter(domain == !!domain, census_otu_count > 0, ncbi_species_count > 0)

  if (nrow(data) == 0) {
    return(data.frame(Taxon = character(), novelty_factor = numeric(),
                     overrepresentation_factor = numeric(), Phylum = character()))
  }

  data <- data %>%
    mutate(
      novelty_factor = census_otu_count / ncbi_species_count,
      overrepresentation_factor = ncbi_species_count / census_otu_count,
      Taxon = if (level == "phylum") phylum else if (level == "family") family else genus
    )

  return(data)
}

# Load 18S data function
load_18s_data <- function(level) {
  filename <- paste0("18s_ncbi_merged_clean_", level, ".csv")
  filepath <- file.path(config$data_dir_18s, filename)

  if (!file.exists(filepath)) stop(paste("File not found:", filepath))

  data <- read.csv(filepath, stringsAsFactors = FALSE) %>%
    filter(census_otu_count > 0, ncbi_species_count > 0)

  if (nrow(data) == 0) {
    return(data.frame(Taxon = character(), novelty_factor = numeric(),
                     overrepresentation_factor = numeric(), Division = character()))
  }

  # Calculate factors and extract taxon name
  data <- data %>%
    mutate(
      novelty_factor = census_otu_count / ncbi_species_count,
      overrepresentation_factor = ncbi_species_count / census_otu_count
    )

  # Add Taxon column based on level
  if (level == "phylum") {
    data$Taxon <- data$phylum
    data$Division <- data$phylum
  } else if (level == "family") {
    data$Taxon <- data$family
    # For family level, we'll need to get the phylum/division later
  } else {
    data$Taxon <- data$genus
  }

  return(data)
}

# Get phylum mapping for family level
get_phylum_for_taxa <- function(taxa, level) {
  file_path <- file.path(config$ncbi_data_dir, paste0("ncbi_", level, "_counts.csv"))
  if (!file.exists(file_path)) return(rep("Unknown", length(taxa)))

  ncbi_data <- read.csv(file_path, stringsAsFactors = FALSE)
  taxon_col <- intersect(c("family_name", "genus_name", "family", "genus", "taxon", "Taxon"),
                        colnames(ncbi_data))[1]

  if (is.na(taxon_col)) return(rep("Unknown", length(taxa)))

  phylum_col <- intersect(c("phylum_name", "phylum", "Phylum"), colnames(ncbi_data))[1]
  if (is.na(phylum_col)) return(rep("Unknown", length(taxa)))

  phylum_map <- setNames(ncbi_data[[phylum_col]], ncbi_data[[taxon_col]])
  result <- phylum_map[taxa]
  result[is.na(result)] <- "Unknown"

  return(as.character(result))
}

# ============================================================================
# LEGEND EXTRACTION FUNCTIONS
# ============================================================================

# Helper function to extract phyla/divisions that appear in domain data (top taxa only)
extract_domain_phyla <- function(all_data, domain) {
  domain_phyla <- c()
  threshold <- 1.0  # Same threshold as used in plots
  top_n <- 10      # Same top_n as used in plots

  if (domain == "Eukaryota") {
    # For Eukaryota: only divisions that are in top 10 novelty or top 10 overrepresentation
    data_key <- paste("18S", domain, "phylum", sep = "_")
    if (data_key %in% names(all_data) && nrow(all_data[[data_key]]) > 0) {
      data <- all_data[[data_key]]

      # Get top novelty taxa (same logic as in plots)
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      # Get top overrepresentation taxa
      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      # Combine top taxa from both categories
      domain_phyla <- unique(c(top_novelty_taxa, top_coverage_taxa))
    }
  } else if (domain == "Bacteria") {
    # For Bacteria: only phyla that are in top 10 novelty or top 10 overrepresentation
    data_key <- paste("16S", domain, "phylum", sep = "_")
    if (data_key %in% names(all_data) && nrow(all_data[[data_key]]) > 0) {
      data <- all_data[[data_key]]

      # Get top novelty taxa
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      # Get top overrepresentation taxa
      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      # Combine top taxa from both categories
      domain_phyla <- unique(c(top_novelty_taxa, top_coverage_taxa))
    }
  } else if (domain == "Archaea") {
    # For Archaea: top phyla from phylum row + additional top phyla from family row

    # First get top phyla from phylum row
    phylum_data_key <- paste("16S", domain, "phylum", sep = "_")
    if (phylum_data_key %in% names(all_data) && nrow(all_data[[phylum_data_key]]) > 0) {
      data <- all_data[[phylum_data_key]]

      # Get top novelty and overrepresentation taxa from phylum level
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      top_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        top_novelty_taxa <- head(novelty_candidates$Taxon, top_n)
      }

      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      top_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        top_coverage_taxa <- head(coverage_candidates$Taxon, top_n)
      }

      domain_phyla <- c(domain_phyla, unique(c(top_novelty_taxa, top_coverage_taxa)))
    }

    # Then add any additional top phyla from family row that aren't already included
    family_data_key <- paste("16S", domain, "family", sep = "_")
    if (family_data_key %in% names(all_data) && nrow(all_data[[family_data_key]]) > 0) {
      data <- all_data[[family_data_key]]

      # Get top taxa from family level and extract their phyla
      novelty_candidates <- data[data$novelty_factor > threshold, ]
      family_novelty_taxa <- c()
      if (nrow(novelty_candidates) > 0) {
        novelty_candidates <- novelty_candidates[order(-novelty_candidates$novelty_factor), ]
        family_novelty_taxa <- head(novelty_candidates$Phylum, top_n)
      }

      coverage_candidates <- data[data$overrepresentation_factor > threshold, ]
      family_coverage_taxa <- c()
      if (nrow(coverage_candidates) > 0) {
        coverage_candidates <- coverage_candidates[order(-coverage_candidates$overrepresentation_factor), ]
        family_coverage_taxa <- head(coverage_candidates$Phylum, top_n)
      }

      family_phyla <- unique(c(family_novelty_taxa, family_coverage_taxa))

      # Only add phyla that aren't already in the list and are actual phyla (not family names)
      color_config <- load_shared_color_config()
      additional_phyla <- family_phyla[!family_phyla %in% domain_phyla &
                                     family_phyla %in% names(color_config$archaea_colors)]
      domain_phyla <- c(domain_phyla, additional_phyla)
    }
  }

  # Remove unknowns and empty values
  domain_phyla <- unique(domain_phyla[!domain_phyla %in% c("Unknown", "", "Other", NA)])

  cat(paste("🔍", domain, "legend will show", length(domain_phyla), "top phyla:", paste(domain_phyla, collapse = ", "), "\n"))

  return(sort(domain_phyla))
}


# ============================================================================
# LEGEND CREATION FUNCTIONS
# ============================================================================

# Function to create one long combined phyla legend for all domains (only phyla that appear in data)
create_combined_phyla_legend <- function(all_data) {
  color_config <- load_shared_color_config()

  # Extract phyla that actually appear in the data for each domain
  bacteria_phyla <- extract_domain_phyla(all_data, "Bacteria")
  archaea_phyla <- extract_domain_phyla(all_data, "Archaea")
  eukaryota_divisions <- extract_domain_phyla(all_data, "Eukaryota")

  # Combine all taxa with domain labels
  all_taxa <- c()
  all_colors <- c()
  all_domains <- c()

  # Add bacteria (only those that appear in data)
  if (length(bacteria_phyla) > 0) {
    bacteria_colors <- sapply(bacteria_phyla, function(phylum) {
      if (phylum %in% names(color_config$bacteria_colors)) {
        return(color_config$bacteria_colors[[phylum]])
      } else {
        # CROSS-DOMAIN RECYCLING: Bacteria uses Eukaryota fallback colors
        shared_pool <- unique(unlist(color_config$fallback_colors$eukaryota))
        pool_index <- ((match(phylum, bacteria_phyla) - 1) %% length(shared_pool)) + 1
        return(shared_pool[pool_index])
      }
    })
    all_taxa <- c(all_taxa, bacteria_phyla)
    all_colors <- c(all_colors, unname(bacteria_colors))
    all_domains <- c(all_domains, rep("Bacteria", length(bacteria_phyla)))
  }

  # Add archaea (only those that appear in data)
  if (length(archaea_phyla) > 0) {
    archaea_colors <- sapply(archaea_phyla, function(phylum) {
      if (phylum %in% names(color_config$archaea_colors)) {
        return(color_config$archaea_colors[[phylum]])
      } else {
        # CROSS-DOMAIN RECYCLING: Archaea uses combined Bacteria + Eukaryota fallback colors
        shared_pool <- unique(c(
          unlist(color_config$fallback_colors$bacteria),
          unlist(color_config$fallback_colors$eukaryota)
        ))
        pool_index <- ((match(phylum, archaea_phyla) - 1) %% length(shared_pool)) + 1
        return(shared_pool[pool_index])
      }
    })
    all_taxa <- c(all_taxa, archaea_phyla)
    all_colors <- c(all_colors, unname(archaea_colors))
    all_domains <- c(all_domains, rep("Archaea", length(archaea_phyla)))
  }

  # Add eukaryota (only those that appear in data)
  if (length(eukaryota_divisions) > 0) {
    euk_colors <- sapply(eukaryota_divisions, function(div) {
      if (div %in% names(color_config$eukaryota_colors)) {
        return(color_config$eukaryota_colors[[div]])
      } else {
        # CROSS-DOMAIN RECYCLING: Eukaryota uses Bacteria fallback colors
        shared_pool <- unique(unlist(color_config$fallback_colors$bacteria))
        pool_index <- ((match(div, eukaryota_divisions) - 1) %% length(shared_pool)) + 1
        return(shared_pool[pool_index])
      }
    })
    all_taxa <- c(all_taxa, eukaryota_divisions)
    all_colors <- c(all_colors, unname(euk_colors))
    all_domains <- c(all_domains, rep("Eukaryota", length(eukaryota_divisions)))
  }

  # Create the combined legend
  combined_legend <- create_long_horizontal_legend(all_taxa, all_colors, all_domains)

  return(combined_legend)
}

# Function to create one long horizontal legend with domain separators
create_long_horizontal_legend <- function(all_taxa, all_colors, all_domains) {
  if (length(all_taxa) == 0) {
    return(ggplot() + theme_void())
  }

  # Create legend data
  legend_data <- data.frame(
    taxon = all_taxa,
    color = all_colors,
    domain = all_domains,
    x = seq_along(all_taxa),
    y = rep(1, length(all_taxa)),
    stringsAsFactors = FALSE
  )

  # Create domain boundary positions for separators
  domain_changes <- which(c(TRUE, all_domains[-1] != all_domains[-length(all_domains)]))
  domain_labels <- data.frame(
    domain = unique(all_domains),
    x_pos = domain_changes + (c(domain_changes[-1], length(all_taxa) + 1) - domain_changes) / 2,
    y_pos = rep(1.8, length(unique(all_domains))),
    stringsAsFactors = FALSE
  )

  cat(paste("🎨 Creating combined legend with", nrow(legend_data), "total taxa\n"))
  cat(paste("   Bacteria:", sum(all_domains == "Bacteria"), "taxa\n"))
  cat(paste("   Archaea:", sum(all_domains == "Archaea"), "taxa\n"))
  cat(paste("   Eukaryota:", sum(all_domains == "Eukaryota"), "taxa\n"))

  # Create the plot with clean format matching the eukaryota legend
  p <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_tile(aes(fill = color), color = "black", width = 0.8, height = 0.6, size = 0.5) +
    geom_text(aes(label = taxon), y = 0.3, angle = 45, hjust = 1, vjust = 1, size = 3) +
    scale_fill_identity() +
    theme_void() +
    ggtitle("Combined Phyla Legend: Bacteria | Archaea | Eukaryota") +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    xlim(0.5, length(all_taxa) + 0.5) + ylim(0, 1.5)

  # Add subtle domain separators (minimal visual impact)
  if (length(domain_changes) > 1) {
    separator_x <- domain_changes[-1] - 0.5
    for (sep_x in separator_x) {
      p <- p + geom_vline(xintercept = sep_x, color = "gray80", size = 0.5, linetype = "dotted")
    }
  }

  return(p)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

main <- function() {
  cat("📂 Loading data...\n")

  # Define levels and domains
  levels <- c("phylum", "family")
  domains_16s <- c("Bacteria", "Archaea")
  domain_18s <- "Eukaryota"

  # Collect all data
  all_data <- list()

  cat("Loading 16S data...\n")
  for (level in levels) {
    for (domain in domains_16s) {
      cat(paste("  Loading", domain, level, "data...\n"))
      data <- load_16s_data(level, domain)
      all_data[[paste("16S", domain, level, sep = "_")]] <- data
    }
  }

  cat("Loading 18S data...\n")
  for (level in levels) {
    cat(paste("  Loading eukaryote", level, "data...\n"))
    data <- load_18s_data(level)
    all_data[[paste("18S", domain_18s, level, sep = "_")]] <- data
  }

  cat("\n🎨 Generating combined phyla legend...\n")
  combined_phyla_legend <- create_combined_phyla_legend(all_data)

  # Save PNG
  output_png <- file.path(config$output_dir, LEGEND_CONFIG$output_png)
  cat(paste("\n💾 Saving PNG:", output_png, "\n"))
  ggsave(output_png, combined_phyla_legend,
         width = LEGEND_CONFIG$width,
         height = LEGEND_CONFIG$height,
         dpi = LEGEND_CONFIG$dpi,
         bg = "white")

  # Save PDF
  output_pdf <- file.path(config$output_dir, LEGEND_CONFIG$output_pdf)
  cat(paste("💾 Saving PDF:", output_pdf, "\n"))
  ggsave(output_pdf, combined_phyla_legend,
         width = LEGEND_CONFIG$width,
         height = LEGEND_CONFIG$height,
         bg = "white")

  cat("\n✅ LEGEND GENERATION COMPLETE!\n")
  cat(paste("   PNG output:", output_png, "\n"))
  cat(paste("   PDF output:", output_pdf, "\n"))
  cat(paste("   Dimensions:", LEGEND_CONFIG$width, "x", LEGEND_CONFIG$height, "inches\n"))
  cat(paste("   Resolution:", LEGEND_CONFIG$dpi, "DPI\n"))
  cat(paste("   Format: Bacteria | Archaea | Eukaryota in one long horizontal strip\n"))
}

# Run main function
main()

