#!/usr/bin/env Rscript

# Extract taxonomic combinations for top novelty genera
# This script identifies all taxonomic combinations (OTU clusters) that contribute 
# to genera marked as Is_Top_Novelty = TRUE

library(dplyr)

# Function to load top novelty genera from source data
load_top_novelty_genera <- function() {
  cat("Loading top novelty genera from source data files...\n")
  
  # Load all three domain files
  bacteria_file <- "source_data/comprehensive_genus_bacteria_top_taxa.tsv"
  archaea_file <- "source_data/comprehensive_genus_archaea_top_taxa.tsv"
  eukaryota_file <- "source_data/comprehensive_genus_eukaryota_top_taxa.tsv"
  
  all_genera <- data.frame()
  
  # Define common columns to keep
  common_cols <- c("genus", "family", "census_otu_count", "ncbi_genome_count", "ncbi_species_count",
                   "isolate_count", "isolate_percentage", "novelty_factor", "overrepresentation_factor",
                   "domain", "Is_Top_Novelty", "Is_Top_Coverage", "Domain")

  # Load bacteria
  if (file.exists(bacteria_file)) {
    bacteria <- read.table(bacteria_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    bacteria$data_source <- "16S_Bacteria"
    bacteria$Division <- NA  # Add Division column for consistency
    bacteria <- bacteria[, c(common_cols, "data_source", "Division")]
    all_genera <- rbind(all_genera, bacteria)
  }

  # Load archaea
  if (file.exists(archaea_file)) {
    archaea <- read.table(archaea_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    archaea$data_source <- "16S_Archaea"
    archaea$Division <- NA  # Add Division column for consistency
    archaea <- archaea[, c(common_cols, "data_source", "Division")]
    all_genera <- rbind(all_genera, archaea)
  }

  # Load eukaryota
  if (file.exists(eukaryota_file)) {
    eukaryota <- read.table(eukaryota_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    eukaryota$data_source <- "18S_Eukaryota"
    # Keep existing Division column
    eukaryota <- eukaryota[, c(common_cols, "data_source", "Division")]
    all_genera <- rbind(all_genera, eukaryota)
  }
  
  # Filter for top novelty genera only
  top_novelty <- all_genera[all_genera$Is_Top_Novelty == TRUE, ]
  
  cat("Found", nrow(top_novelty), "top novelty genera across all domains\n")
  cat("Breakdown by domain:\n")
  print(table(top_novelty$Domain))
  
  return(top_novelty)
}

# Function to load taxonomic combinations
load_taxonomic_combinations <- function() {
  cat("\nLoading taxonomic combinations from census data...\n")
  
  combinations_file <- "../18S_censusparse/metadata/sanity_check/taxonomic_combinations_detailed.csv"
  
  if (!file.exists(combinations_file)) {
    stop("Taxonomic combinations file not found: ", combinations_file)
  }
  
  # Read the file, skipping comment lines
  combinations <- read.csv(combinations_file, comment.char="#", stringsAsFactors=FALSE)
  
  cat("Loaded", nrow(combinations), "taxonomic combinations\n")
  cat("Columns:", paste(colnames(combinations), collapse=", "), "\n")
  
  return(combinations)
}

# Function to extract combinations for top novelty genera
extract_novelty_combinations <- function(top_novelty_genera, taxonomic_combinations) {
  cat("\nExtracting taxonomic combinations for top novelty genera...\n")
  
  # Get unique top novelty genera (focusing on 18S data since that's where we have detailed combinations)
  euk_novelty <- top_novelty_genera[top_novelty_genera$data_source == "18S_Eukaryota", ]
  
  if (nrow(euk_novelty) == 0) {
    cat("No eukaryotic top novelty genera found\n")
    return(data.frame())
  }
  
  cat("Processing", nrow(euk_novelty), "eukaryotic top novelty genera:\n")
  print(euk_novelty$genus)
  
  # Match genera to taxonomic combinations
  novelty_combinations <- data.frame()
  
  for (i in 1:nrow(euk_novelty)) {
    genus_name <- euk_novelty$genus[i]
    novelty_factor <- euk_novelty$novelty_factor[i]
    census_count <- euk_novelty$census_otu_count[i]
    
    cat("\nProcessing genus:", genus_name, "(novelty factor:", novelty_factor, ", census OTUs:", census_count, ")\n")
    
    # Find all combinations for this genus
    genus_combinations <- taxonomic_combinations[taxonomic_combinations$genus == genus_name, ]
    
    if (nrow(genus_combinations) > 0) {
      # Add metadata
      genus_combinations$target_genus <- genus_name
      genus_combinations$novelty_factor <- novelty_factor
      genus_combinations$census_otu_count <- census_count
      
      cat("  Found", nrow(genus_combinations), "taxonomic combinations\n")
      cat("  Total OTU clusters:", sum(genus_combinations$row_count), "\n")
      cat("  Sample combinations:\n")
      for (j in 1:min(3, nrow(genus_combinations))) {
        cat("    ", genus_combinations$taxonomic_combination[j], "\n")
      }
      
      novelty_combinations <- rbind(novelty_combinations, genus_combinations)
    } else {
      cat("  No combinations found for", genus_name, "\n")
    }
  }
  
  return(novelty_combinations)
}

# Main function
main <- function() {
  cat("=== Extracting Taxonomic Combinations for Top Novelty Genera ===\n\n")
  
  # Load top novelty genera
  top_novelty <- load_top_novelty_genera()
  
  # Load taxonomic combinations
  combinations <- load_taxonomic_combinations()
  
  # Extract combinations for novelty genera
  novelty_combinations <- extract_novelty_combinations(top_novelty, combinations)
  
  if (nrow(novelty_combinations) > 0) {
    # Export results
    output_file <- "source_data/top_novelty_taxonomic_combinations.tsv"
    write.table(novelty_combinations, output_file, sep="\t", row.names=FALSE, quote=FALSE)
    cat("\n=== RESULTS EXPORTED ===\n")
    cat("Taxonomic combinations for top novelty genera exported to:", output_file, "\n")
    cat("Total combinations:", nrow(novelty_combinations), "\n")
    cat("Total OTU clusters represented:", sum(novelty_combinations$row_count), "\n")
    
    # Summary by genus
    cat("\nSummary by genus:\n")
    summary_by_genus <- novelty_combinations %>%
      group_by(target_genus, novelty_factor, census_otu_count) %>%
      summarise(
        num_combinations = n(),
        total_clusters = sum(row_count),
        total_size = sum(total_size),
        avg_cluster_size = round(mean(avg_size_per_cluster), 2),
        .groups = 'drop'
      ) %>%
      arrange(desc(novelty_factor))
    
    print(summary_by_genus)
    
    # Export summary
    summary_file <- "source_data/top_novelty_genus_summary.tsv"
    write.table(summary_by_genus, summary_file, sep="\t", row.names=FALSE, quote=FALSE)
    cat("\nSummary exported to:", summary_file, "\n")
    
  } else {
    cat("\nNo taxonomic combinations found for top novelty genera\n")
  }
  
  cat("\n=== EXTRACTION COMPLETE ===\n")
}

# Run if called directly
if (!interactive()) {
  main()
}
