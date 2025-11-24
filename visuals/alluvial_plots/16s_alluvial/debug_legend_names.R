#!/usr/bin/env Rscript

# Debug what phylum names are actually in the legend
library(yaml)

# Load the shared config
config_file <- "../../shared_config/taxonomic_color_mapping.yaml"
color_config <- yaml.load_file(config_file)

# Simulate the phylum names as they appear in the plot
# Based on line 202: paste0(i, ". ", top_phyla$phylum[i])
simulated_phyla <- c(
  "1. Pseudomonadota",
  "2. Bacillota", 
  "3. Campylobacterota",
  "4. Bacteroidota",
  "5. Actinomycetota",
  "6. Chloroflexota",
  "7. Acidobacteriota",
  "8. Verrucomicrobiota",
  "9. Planctomycetota",
  "10. Thermodesulfobacteriota",
  "11. Cyanobacteriota",
  "12. Bacteria.U.phylum",
  "Other"
)

cat("\n=== SIMULATED LEGEND NAMES ===\n")
for (i in seq_along(simulated_phyla)) {
  cat(sprintf("%2d. %s\n", i, simulated_phyla[i]))
}

# Now test color assignment
get_bacterial_colors <- function() {
  bacteria_list <- color_config$bacteria_colors
  bacteria_colors <- character(length(bacteria_list))
  names(bacteria_colors) <- names(bacteria_list)
  for (name in names(bacteria_list)) {
    bacteria_colors[name] <- as.character(bacteria_list[[name]])
  }
  return(bacteria_colors)
}

get_extended_bacterial_colors <- function() {
  c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
    "#F1948A", "#C39BD3", "#D7BDE2", "#A9DFBF", "#F9E79F")
}

get_excluded_phyla <- function() {
  c("Campylobacterota", "Mycoplasmatota", "Thermotogota")
}

bacterial_colors <- get_bacterial_colors()
extended_colors <- get_extended_bacterial_colors()
excluded_phyla <- get_excluded_phyla()

cat("\n=== COLOR ASSIGNMENTS ===\n")
colors <- character(length(simulated_phyla))
names(colors) <- simulated_phyla

for (i in seq_along(simulated_phyla)) {
  phylum_name <- simulated_phyla[i]
  clean_phylum <- gsub("^[0-9]+\\. ", "", phylum_name)
  
  if (clean_phylum == "Other") {
    colors[i] <- "#808080"
    assignment <- "Other (gray)"
  } else if (grepl("\\.U\\.", clean_phylum)) {
    colors[i] <- "#9ACD32"
    assignment <- ".U. entry (yellow-green)"
  } else if (clean_phylum %in% excluded_phyla) {
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
    assignment <- sprintf("Excluded phyla fallback [%d]", fallback_index)
  } else if (clean_phylum %in% names(bacterial_colors)) {
    colors[i] <- bacterial_colors[clean_phylum]
    assignment <- "Shared config"
  } else {
    fallback_index <- ((i - 1) %% length(extended_colors)) + 1
    colors[i] <- extended_colors[fallback_index]
    assignment <- sprintf("General fallback [%d]", fallback_index)
  }
  
  cat(sprintf("%2d. %-30s -> %-10s (%s)\n", i, phylum_name, colors[i], assignment))
}

cat("\n=== CHECKING FOR BLUE COLOR (#1f77b4) ===\n")
blue_phyla <- names(colors)[colors == "#1f77b4"]
if (length(blue_phyla) > 0) {
  cat("✅ Blue color assigned to:\n")
  for (p in blue_phyla) {
    cat(sprintf("   - %s\n", p))
  }
} else {
  cat("❌ Blue color (#1f77b4) NOT assigned to any phylum!\n")
}

