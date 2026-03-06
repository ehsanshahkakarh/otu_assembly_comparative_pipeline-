#!/usr/bin/env Rscript
#
# 18S Eukaryotic Sankey Diagram with COLORED FLOWS
# =================================================
#
# Creates a Sankey diagram showing the flow from:
# NCBI Total Genomes → 18S EukCensus Sequences → 18S EukCensus OTUs → NCBI Total Species
# 
# For eukaryotic (18S) data
# Uses shared taxonomic color mapping for FLOWS (not just nodes)
# Shows top divisions by total representation + Other
#
# Author: Enhanced Sankey Team
# Date: 2026-01-28
#

# Load required libraries
library(dplyr)
library(networkD3)
library(htmlwidgets)
library(scales)
library(yaml)

cat("=== 18S Eukaryotic Sankey Diagram with Colored Flows ===\n\n")

# Source shared color mapping functions
color_functions_path <- "../../shared_config/color_mapping_functions.R"
if (!file.exists(color_functions_path)) {
  stop("ERROR: Cannot find color mapping functions at: ", color_functions_path)
}
source(color_functions_path)

# Load color configuration
cat("Loading taxonomic color configuration...\n")
color_config <- load_taxonomic_colors()

# Load alluvial filtering configuration
config_path <- "../config/alluvial_filtering_config.yaml"
if (!file.exists(config_path)) {
  stop("ERROR: Cannot find alluvial config at: ", config_path)
}
alluvial_config <- yaml.load_file(config_path)

# Get eukaryota-specific configuration
domain_config <- alluvial_config$eukaryota_18s
top_n <- domain_config$filtering$hybrid$max_total_taxa

# Load the 18S data
cat("Loading 18S data...\n")
merged_data_path <- "../../../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_phylum.csv"
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 18S merged data file at: ", merged_data_path)
}
data_18s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
cat("18S data loaded:", nrow(data_18s), "rows\n")

# Filter for matched divisions for Eukaryota
matched_data <- data_18s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
  filter(domain == "Eukaryota")

cat("Matched eukaryotic divisions found:", nrow(matched_data), "\n")

# Add total representation across all 4 nodes and get top N
matched_data$total_representation <- matched_data$ncbi_genome_count +
                                     matched_data$census_size_count +
                                     matched_data$census_otu_count +
                                     matched_data$ncbi_species_count

top_divisions <- matched_data %>%
  arrange(desc(total_representation)) %>%
  head(top_n)

cat("Top", top_n, "divisions selected:", nrow(top_divisions), "\n")

# Calculate total counts for each node from MATCHED DATA ONLY
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_size_count <- sum(matched_data$census_size_count, na.rm = TRUE)
total_otu_count <- sum(matched_data$census_otu_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)

# Calculate "Other" counts (remaining divisions not in top N)
other_genome_count <- total_genome_count - sum(top_divisions$ncbi_genome_count)
other_size_count <- total_size_count - sum(top_divisions$census_size_count)
other_otu_count <- total_otu_count - sum(top_divisions$census_otu_count)
other_species_count <- total_species_count - sum(top_divisions$ncbi_species_count)

cat("Other counts (remaining divisions):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Get division names and assign colors using shared color mapping
division_names <- c(top_divisions$phylum, "Other")
cat("\nAssigning colors from shared taxonomic color mapping...\n")
division_colors <- get_domain_colors(division_names, "Eukaryota", color_config)

# Print color mapping
print_color_summary(division_names, division_colors, "Eukaryota")

# Create nodes dataframe
# Node structure: Each division appears in each stage
stages <- c("Genomes", "Sequences", "OTUs", "Species")

# Create all combinations of divisions and stages
nodes <- data.frame()
node_id <- 0

for (stage in stages) {
  for (i in seq_along(division_names)) {
    division <- division_names[i]
    nodes <- rbind(nodes, data.frame(
      id = node_id,
      name = paste0(division, "_", stage),
      division = division,
      stage = stage,
      color = division_colors[i],
      stringsAsFactors = FALSE
    ))
    node_id <- node_id + 1
  }
}

cat("Created", nrow(nodes), "nodes\n")

# Create links dataframe with color information
links <- data.frame()

# Helper function to get node ID
get_node_id <- function(division, stage) {
  return(nodes$id[nodes$division == division & nodes$stage == stage])
}

# Create links between consecutive stages for top divisions
for (i in 1:nrow(top_divisions)) {
  division <- top_divisions$phylum[i]
  link_color <- division_colors[i]
  
  # Genomes -> Sequences
  links <- rbind(links, data.frame(
    source = get_node_id(division, "Genomes"),
    target = get_node_id(division, "Sequences"),
    value = top_divisions$census_size_count[i],
    division = division,
    color = link_color,
    stringsAsFactors = FALSE
  ))

  # Sequences -> OTUs
  links <- rbind(links, data.frame(
    source = get_node_id(division, "Sequences"),
    target = get_node_id(division, "OTUs"),
    value = top_divisions$census_otu_count[i],
    division = division,
    color = link_color,
    stringsAsFactors = FALSE
  ))

  # OTUs -> Species
  links <- rbind(links, data.frame(
    source = get_node_id(division, "OTUs"),
    target = get_node_id(division, "Species"),
    value = top_divisions$ncbi_species_count[i],
    division = division,
    color = link_color,
    stringsAsFactors = FALSE
  ))
}

# Add "Other" links with gray color
other_color <- color_config$special_colors$Other
links <- rbind(links, data.frame(
  source = get_node_id("Other", "Genomes"),
  target = get_node_id("Other", "Sequences"),
  value = other_size_count,
  division = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "Sequences"),
  target = get_node_id("Other", "OTUs"),
  value = other_otu_count,
  division = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "OTUs"),
  target = get_node_id("Other", "Species"),
  value = other_species_count,
  division = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

cat("Created", nrow(links), "links with individual colors\n")

# Create the Sankey diagram with COLORED FLOWS
cat("Creating Sankey diagram with colored flows...\n")

# Prepare color scale for nodes (by division)
node_color_scale <- paste0('d3.scaleOrdinal().domain([',
                           paste0('"', division_names, '"', collapse = ','),
                           ']).range([',
                           paste0('"', division_colors, '"', collapse = ','),
                           '])')

# Create custom JavaScript to color links
links$link_group <- links$division

link_color_scale <- paste0('d3.scaleOrdinal().domain([',
                           paste0('"', division_names, '"', collapse = ','),
                           ']).range([',
                           paste0('"', division_colors, '"', collapse = ','),
                           '])')

sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "division",
  LinkGroup = "link_group",  # This colors the flows!
  colourScale = node_color_scale,
  linkColour = link_color_scale,  # This applies colors to links
  fontSize = 12,
  fontFamily = "Arial",
  nodeWidth = 20,
  nodePadding = 10,
  margin = list(top = 50, right = 50, bottom = 50, left = 50),
  height = 800,
  width = 1200,
  iterations = 100,
  sinksRight = FALSE
)

# Save the Sankey diagram
output_file <- "sankey_18s_eukaryota_colored_flows.html"
saveWidget(sankey_plot, file = output_file, selfcontained = TRUE)

cat("Sankey diagram saved:", output_file, "\n")

# Try to save a static PNG version using webshot (if available)
tryCatch({
  library(webshot)
  png_file <- "sankey_18s_eukaryota_colored_flows.png"
  webshot(output_file, png_file, vwidth = 1200, vheight = 800)
  cat("Static PNG version saved:", png_file, "\n")
}, error = function(e) {
  cat("Note: webshot package not available for PNG export\n")
  cat("HTML file can be transferred to local machine for viewing\n")
})

# Display summary statistics
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("🧬 EUKARYOTA SANKEY DIAGRAM WITH COLORED FLOWS - SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("📊 DATASET OVERVIEW:\n")
cat("   • Domain: Eukaryota\n")
cat("   • Total divisions analyzed:", nrow(top_divisions), "+ Other\n")
cat("   • Data flow: Genomes → Sequences → OTUs → Species\n")
cat("   • Interactive HTML output with colored flows\n")
cat("   • Colors from shared taxonomic color mapping\n\n")

cat("📈 TOTAL COUNTS:\n")
cat("   • Total genomes:", scales::comma(total_genome_count), "\n")
cat("   • Total sequences:", scales::comma(total_size_count), "\n")
cat("   • Total OTUs:", scales::comma(total_otu_count), "\n")
cat("   • Total species:", scales::comma(total_species_count), "\n\n")

cat("🎨 VISUALIZATION FEATURES:\n")
cat("   • Interactive Sankey diagram with hover details\n")
cat("   • COLORED FLOWS matching division colors\n")
cat("   • Color-coded divisions with", length(division_colors), "distinct colors\n")
cat("   • Flow widths represent absolute counts\n")
cat("   • HTML output for web viewing and sharing\n\n")

# Create a text-based flow summary
cat("📋 FLOW SUMMARY:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
            "Division", "Genomes", "Sequences", "OTUs", "Species", "Color"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:nrow(top_divisions)) {
  cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
              substr(top_divisions$phylum[i], 1, 25),
              scales::comma(top_divisions$ncbi_genome_count[i]),
              scales::comma(top_divisions$census_size_count[i]),
              scales::comma(top_divisions$census_otu_count[i]),
              scales::comma(top_divisions$ncbi_species_count[i]),
              division_colors[i]))
}

cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
            "Other",
            scales::comma(other_genome_count),
            scales::comma(other_size_count),
            scales::comma(other_otu_count),
            scales::comma(other_species_count),
            other_color))

cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("✅ EUKARYOTA Sankey diagram with colored flows completed!\n")
cat("   📁 HTML file:", output_file, "\n")
cat("   💡 Transfer to local machine to view interactive diagram\n")
cat("   🎨 Flows are colored to match their division colors\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

