#!/usr/bin/env Rscript
# Cultivation Success Spectrum Visualization
# Created: 2025-11-13
# Purpose: Create an innovative "spectrum" visualization showing the cultivation success
#          across all domains of life, from completely uncultured to well-represented taxa

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(cowplot)
  library(RColorBrewer)
})

# Set working directory to script location
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- getwd()
}
setwd(script_dir)

cat("🎨 Creating Cultivation Success Spectrum Visualization...\n")
cat("📁 Working directory:", getwd(), "\n")

# Configuration
CONFIG <- list(
  # Data paths (relative to visuals directory)
  prokaryotic_16s_gtdb = "../Eukcensus_merge/16s_merged/csv_results/16s_gtdb_merged_clean_phylum.csv",
  prokaryotic_16s_ncbi = "../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv", 
  eukaryotic_18s_eukprot = "../Eukcensus_merge/18s_merged/csv_results/18s_eukprot_merged_division.csv",
  eukaryotic_18s_ncbi = "../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv",
  
  # Visual parameters
  plot_width = 16,
  plot_height = 12,
  dpi = 300,
  text_size = 12,
  
  # Spectrum categories
  spectrum_breaks = c(0, 0.1, 1, 10, 50, 100),
  spectrum_labels = c("Uncultured\n(0%)", "Critically\nUndercultured\n(<0.1%)", 
                     "Severely\nUndercultured\n(0.1-1%)", "Moderately\nUndercultured\n(1-10%)",
                     "Well\nCultured\n(10-50%)", "Extensively\nCultured\n(>50%)")
)

# Function to load and process data
load_and_process_data <- function() {
  cat("📊 Loading merger data files...\n")
  
  # Load prokaryotic data (16S)
  if (file.exists(CONFIG$prokaryotic_16s_gtdb)) {
    prok_gtdb <- read.csv(CONFIG$prokaryotic_16s_gtdb, stringsAsFactors = FALSE)
    prok_gtdb$database <- "GTDB"
    prok_gtdb$domain_type <- "Prokaryotic"
    prok_gtdb$coverage_pct <- prok_gtdb$gtdb_species_count / prok_gtdb$census_otu_count * 100
    prok_gtdb$coverage_pct[is.infinite(prok_gtdb$coverage_pct)] <- 0
    prok_gtdb$coverage_pct[is.na(prok_gtdb$coverage_pct)] <- 0
    cat("✅ Loaded GTDB prokaryotic data:", nrow(prok_gtdb), "phyla\n")
  } else {
    cat("❌ GTDB prokaryotic file not found\n")
    prok_gtdb <- data.frame()
  }
  
  # Load eukaryotic data (18S)
  if (file.exists(CONFIG$eukaryotic_18s_eukprot)) {
    euk_eukprot <- read.csv(CONFIG$eukaryotic_18s_eukprot, stringsAsFactors = FALSE)
    euk_eukprot$database <- "EukProt"
    euk_eukprot$domain_type <- "Eukaryotic"
    euk_eukprot$phylum <- euk_eukprot$division  # Standardize column name
    euk_eukprot$coverage_pct <- euk_eukprot$coverage_percentage
    cat("✅ Loaded EukProt eukaryotic data:", nrow(euk_eukprot), "divisions\n")
  } else {
    cat("❌ EukProt eukaryotic file not found\n")
    euk_eukprot <- data.frame()
  }
  
  # Combine datasets
  if (nrow(prok_gtdb) > 0 && nrow(euk_eukprot) > 0) {
    # Select common columns
    prok_subset <- prok_gtdb %>%
      select(phylum, domain_type, database, coverage_pct, census_otu_count, novelty_factor) %>%
      filter(coverage_pct >= 0)
    
    euk_subset <- euk_eukprot %>%
      select(phylum, domain_type, database, coverage_pct, census_otu_count, novelty_factor) %>%
      filter(coverage_pct >= 0)
    
    combined_data <- rbind(prok_subset, euk_subset)
    
    # Add spectrum categories
    combined_data$spectrum_category <- cut(combined_data$coverage_pct, 
                                         breaks = CONFIG$spectrum_breaks,
                                         labels = CONFIG$spectrum_labels,
                                         include.lowest = TRUE, right = FALSE)
    
    cat("🔬 Combined dataset:", nrow(combined_data), "taxa across domains\n")
    return(combined_data)
  } else {
    cat("❌ Insufficient data for visualization\n")
    return(data.frame())
  }
}

# Function to create spectrum colors
create_spectrum_colors <- function() {
  # Create a cultivation success spectrum from red (uncultured) to green (well-cultured)
  spectrum_colors <- c(
    "#8B0000",  # Dark red - Uncultured
    "#CD5C5C",  # Indian red - Critically undercultured  
    "#FF6347",  # Tomato - Severely undercultured
    "#FFD700",  # Gold - Moderately undercultured
    "#9ACD32",  # Yellow green - Well cultured
    "#228B22"   # Forest green - Extensively cultured
  )
  names(spectrum_colors) <- CONFIG$spectrum_labels
  return(spectrum_colors)
}

# Function to create the main spectrum plot
create_spectrum_plot <- function(data) {
  cat("🎨 Creating cultivation success spectrum plot...\n")
  
  spectrum_colors <- create_spectrum_colors()
  
  # Create the main plot
  p <- ggplot(data, aes(x = log10(census_otu_count + 1), y = log10(novelty_factor + 0.1))) +
    geom_point(aes(fill = spectrum_category, size = coverage_pct), 
               shape = 21, alpha = 0.8, stroke = 0.5) +
    scale_fill_manual(values = spectrum_colors, name = "Cultivation\nSuccess") +
    scale_size_continuous(range = c(3, 12), name = "Coverage %",
                         breaks = c(0, 1, 10, 50, 100),
                         labels = c("0%", "1%", "10%", "50%", "100%")) +
    facet_wrap(~ domain_type, scales = "free", ncol = 2) +
    labs(
      title = "Cultivation Success Spectrum Across Domains of Life",
      subtitle = "Environmental diversity (OTU count) vs Cultivation success (novelty factor)\nColor represents cultivation success from uncultured (red) to well-cultured (green)",
      x = "Environmental Diversity (log10 OTU count)",
      y = "Cultivation Gap (log10 novelty factor)",
      caption = "Data: EukCensus environmental surveys vs GTDB/EukProt genomic databases"
    ) +
    theme_minimal(base_size = CONFIG$text_size) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

# Main execution
main <- function() {
  cat("🚀 Starting Cultivation Success Spectrum Analysis...\n")
  
  # Load and process data
  data <- load_and_process_data()
  
  if (nrow(data) == 0) {
    cat("❌ No data available for visualization\n")
    return()
  }
  
  # Create visualization
  spectrum_plot <- create_spectrum_plot(data)
  
  # Save the plot
  output_file <- "cultivation_success_spectrum.png"
  cat("💾 Saving plot to:", output_file, "\n")
  
  ggsave(output_file, spectrum_plot, 
         width = CONFIG$plot_width, height = CONFIG$plot_height, 
         dpi = CONFIG$dpi, bg = "white")
  
  cat("✅ Cultivation Success Spectrum visualization complete!\n")
  cat("📊 Summary statistics:\n")
  
  # Print summary statistics
  summary_stats <- data %>%
    group_by(domain_type, spectrum_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(domain_type, spectrum_category)
  
  print(summary_stats)
}

# Execute main function
main()
