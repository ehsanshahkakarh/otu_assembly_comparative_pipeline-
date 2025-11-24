#!/usr/bin/env Rscript

# Debug Pseudomonadota Color Assignment
# =====================================

# Load required libraries
library(yaml)

# Load the shared config
config_file <- "../../shared_config/taxonomic_color_mapping.yaml"
if (file.exists(config_file)) {
  color_config <- yaml.load_file(config_file)
  cat("✅ Loaded shared config from:", config_file, "\n")
} else {
  stop("❌ Cannot find shared config file")
}

# Test the hardcoded functions from the script
get_bacterial_colors <- function() {
  bacteria_list <- color_config$bacteria_colors
  bacteria_colors <- character(length(bacteria_list))
  names(bacteria_colors) <- names(bacteria_list)
  for (name in names(bacteria_list)) {
    bacteria_colors[name] <- as.character(bacteria_list[[name]])
  }
  return(bacteria_colors)
}

get_excluded_phyla <- function() {
  c("Campylobacterota", "Mycoplasmatota", "Thermotogota")
}

get_extended_bacterial_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

# Test color assignment
bacterial_colors <- get_bacterial_colors()
extended_colors <- get_extended_bacterial_colors()
excluded_phyla <- get_excluded_phyla()

cat("\n=== DEBUGGING PSEUDOMONADOTA COLOR ASSIGNMENT ===\n")

# Test phylum names
test_phyla <- c("1. Pseudomonadota", "2. Bacillota", "3. Bacteroidota", "4. Bacteria.U.phylum")

for (i in seq_along(test_phyla)) {
  phylum_name <- test_phyla[i]
  clean_phylum <- gsub("^[0-9]+\\. ", "", phylum_name)
  
  cat(sprintf("\n--- Testing: %s ---\n", phylum_name))
  cat(sprintf("Clean name: %s\n", clean_phylum))
  cat(sprintf("In excluded_phyla: %s\n", clean_phylum %in% excluded_phyla))
  cat(sprintf("In bacterial_colors: %s\n", clean_phylum %in% names(bacterial_colors)))
  
  if (clean_phylum %in% names(bacterial_colors)) {
    cat(sprintf("Shared config color: %s\n", bacterial_colors[clean_phylum]))
  }
  
  # Apply the same logic as the script
  assigned_color <- NULL
  if (clean_phylum == "Other") {
    assigned_color <- "#808080"
    cat("Assigned: Other color\n")
  } else if (grepl("\\.U\\.", clean_phylum)) {
    assigned_color <- "#9ACD32"
    cat("Assigned: .U. entry color\n")
  } else if (clean_phylum %in% excluded_phyla) {
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    assigned_color <- extended_colors[fallback_index]
    cat(sprintf("Assigned: Excluded phyla fallback color [%d]: %s\n", fallback_index, assigned_color))
  } else if (clean_phylum %in% names(bacterial_colors)) {
    assigned_color <- bacterial_colors[clean_phylum]
    cat(sprintf("Assigned: Shared config color: %s\n", assigned_color))
  } else {
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    assigned_color <- extended_colors[fallback_index]
    cat(sprintf("Assigned: General fallback color [%d]: %s\n", fallback_index, assigned_color))
  }
  
  cat(sprintf("FINAL COLOR: %s\n", assigned_color))
}

cat("\n=== SHARED CONFIG CONTENTS ===\n")
cat("Pseudomonadota in bacteria_colors:\n")
if ("Pseudomonadota" %in% names(bacterial_colors)) {
  cat(sprintf("  Pseudomonadota: %s\n", bacterial_colors["Pseudomonadota"]))
} else {
  cat("  ❌ Pseudomonadota NOT FOUND in bacteria_colors\n")
  cat("Available bacteria colors:\n")
  for (name in names(bacterial_colors)) {
    cat(sprintf("  %s: %s\n", name, bacterial_colors[name]))
  }
}
