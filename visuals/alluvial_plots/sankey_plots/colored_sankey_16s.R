#!/usr/bin/env Rscript
#
# 16S Prokaryotic Sankey Diagram with COLORED FLOWS
# ==================================================
#
# Creates a Sankey diagram showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
# 
# For prokaryotic (16S) data - Bacteria and Archaea
# Uses shared taxonomic color mapping for FLOWS (not just nodes)
# Shows top phyla by total representation + Other
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

cat("=== 16S Prokaryotic Sankey Diagram with Colored Flows ===\n\n")

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

# Configuration: Set which domain to process
PROCESS_DOMAIN <- "Bacteria"  # Change to "Archaea" for archaea plot
cat(paste("Processing domain:", PROCESS_DOMAIN, "\n\n"))

# Get domain-specific configuration
if (PROCESS_DOMAIN == "Bacteria") {
  domain_config <- alluvial_config$bacteria_16s
  top_n <- domain_config$filtering$hybrid$max_total_taxa
} else if (PROCESS_DOMAIN == "Archaea") {
  domain_config <- alluvial_config$archaea_16s
  top_n <- domain_config$filtering$top_n
} else {
  stop("Invalid domain. Must be 'Bacteria' or 'Archaea'")
}

# Load the 16S data
cat("Loading 16S data...\n")
merged_data_path <- "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
if (!file.exists(merged_data_path)) {
  stop("ERROR: Cannot find 16S merged data file at: ", merged_data_path)
}
data_16s <- read.csv(merged_data_path, stringsAsFactors = FALSE)
cat("16S data loaded:", nrow(data_16s), "rows\n")

# Filter for matched phyla for specific domain
matched_data <- data_16s %>%
  filter(match_status == 'matched') %>%
  filter(!is.na(census_size_count) & !is.na(census_otu_count) &
         !is.na(ncbi_genome_count) & !is.na(ncbi_species_count) & !is.na(phylum)) %>%
  filter(phylum != "N/A" & phylum != "" & !is.null(phylum)) %>%
  filter(domain == PROCESS_DOMAIN)

cat("Matched", tolower(PROCESS_DOMAIN), "phyla found:", nrow(matched_data), "\n")

# Add total representation across all 4 nodes and get top N
matched_data$total_representation <- matched_data$ncbi_genome_count +
                                     matched_data$census_size_count +
                                     matched_data$census_otu_count +
                                     matched_data$ncbi_species_count

top_phyla <- matched_data %>%
  arrange(desc(total_representation)) %>%
  head(top_n)

cat("Top", top_n, "phyla selected:", nrow(top_phyla), "\n")

# Calculate total counts for each node from MATCHED DATA ONLY
total_genome_count <- sum(matched_data$ncbi_genome_count, na.rm = TRUE)
total_size_count <- sum(matched_data$census_size_count, na.rm = TRUE)
total_otu_count <- sum(matched_data$census_otu_count, na.rm = TRUE)
total_species_count <- sum(matched_data$ncbi_species_count, na.rm = TRUE)

# Calculate "Other" counts (remaining phyla not in top N)
other_genome_count <- total_genome_count - sum(top_phyla$ncbi_genome_count)
other_size_count <- total_size_count - sum(top_phyla$census_size_count)
other_otu_count <- total_otu_count - sum(top_phyla$census_otu_count)
other_species_count <- total_species_count - sum(top_phyla$ncbi_species_count)

cat("Other counts (remaining phyla):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Get phyla names and assign colors using shared color mapping
phyla_names <- c(top_phyla$phylum, "Other")
cat("\nAssigning colors from shared taxonomic color mapping...\n")
phyla_colors <- get_domain_colors(phyla_names, PROCESS_DOMAIN, color_config)

# Print color mapping
print_color_summary(phyla_names, phyla_colors, PROCESS_DOMAIN)

# Create nodes dataframe
# Node structure: Each phylum appears in each stage
stages <- c("Genomes", "Sequences", "OTUs", "Species")

# Create all combinations of phyla and stages
nodes <- data.frame()
node_id <- 0

for (stage in stages) {
  for (i in seq_along(phyla_names)) {
    phylum <- phyla_names[i]
    nodes <- rbind(nodes, data.frame(
      id = node_id,
      name = paste0(phylum, "_", stage),
      phylum = phylum,
      stage = stage,
      color = phyla_colors[i],
      stringsAsFactors = FALSE
    ))
    node_id <- node_id + 1
  }
}

cat("Created", nrow(nodes), "nodes\n")

# Create links dataframe with color information
links <- data.frame()

# Helper function to get node ID
get_node_id <- function(phylum, stage) {
  return(nodes$id[nodes$phylum == phylum & nodes$stage == stage])
}

# Create links between consecutive stages for top phyla
for (i in 1:nrow(top_phyla)) {
  phylum <- top_phyla$phylum[i]
  link_color <- phyla_colors[i]

  # Genomes -> Sequences
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "Genomes"),
    target = get_node_id(phylum, "Sequences"),
    value = top_phyla$census_size_count[i],
    phylum = phylum,
    color = link_color,
    stringsAsFactors = FALSE
  ))

  # Sequences -> OTUs
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "Sequences"),
    target = get_node_id(phylum, "OTUs"),
    value = top_phyla$census_otu_count[i],
    phylum = phylum,
    color = link_color,
    stringsAsFactors = FALSE
  ))

  # OTUs -> Species
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "OTUs"),
    target = get_node_id(phylum, "Species"),
    value = top_phyla$ncbi_species_count[i],
    phylum = phylum,
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
  phylum = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "Sequences"),
  target = get_node_id("Other", "OTUs"),
  value = other_otu_count,
  phylum = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "OTUs"),
  target = get_node_id("Other", "Species"),
  value = other_species_count,
  phylum = "Other",
  color = other_color,
  stringsAsFactors = FALSE
))

cat("Created", nrow(links), "links with individual colors\n")

# Create the Sankey diagram with COLORED FLOWS
cat("Creating Sankey diagram with colored flows...\n")

# Prepare color scale for nodes (by phylum)
node_color_scale <- paste0('d3.scaleOrdinal().domain([',
                           paste0('"', phyla_names, '"', collapse = ','),
                           ']).range([',
                           paste0('"', phyla_colors, '"', collapse = ','),
                           '])')

# Create custom JavaScript to color links
# We'll use the LinkGroup parameter and create a color scale for it
links$link_group <- links$phylum

link_color_scale <- paste0('d3.scaleOrdinal().domain([',
                           paste0('"', phyla_names, '"', collapse = ','),
                           ']).range([',
                           paste0('"', phyla_colors, '"', collapse = ','),
                           '])')

sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "phylum",
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
domain_lower <- tolower(PROCESS_DOMAIN)
output_file <- paste0("sankey_16s_", domain_lower, "_colored_flows.html")
saveWidget(sankey_plot, file = output_file, selfcontained = TRUE)

cat("Sankey diagram saved:", output_file, "\n")

# Try to save a static PNG version using webshot (if available)
tryCatch({
  library(webshot)
  png_file <- paste0("sankey_16s_", domain_lower, "_colored_flows.png")
  webshot(output_file, png_file, vwidth = 1200, vheight = 800)
  cat("Static PNG version saved:", png_file, "\n")
}, error = function(e) {
  cat("Note: webshot package not available for PNG export\n")
  cat("HTML file can be transferred to local machine for viewing\n")
})

# Display summary statistics
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("🧬", toupper(PROCESS_DOMAIN), "SANKEY DIAGRAM WITH COLORED FLOWS - SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("📊 DATASET OVERVIEW:\n")
cat("   • Domain:", PROCESS_DOMAIN, "\n")
cat("   • Total phyla analyzed:", nrow(top_phyla), "+ Other\n")
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
cat("   • COLORED FLOWS matching phylum colors\n")
cat("   • Color-coded phyla with", length(phyla_colors), "distinct colors\n")
cat("   • Flow widths represent absolute counts\n")
cat("   • HTML output for web viewing and sharing\n\n")

# Create a text-based flow summary
cat("📋 FLOW SUMMARY:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
            "Phylum", "Genomes", "Sequences", "OTUs", "Species", "Color"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:nrow(top_phyla)) {
  cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
              substr(top_phyla$phylum[i], 1, 25),
              scales::comma(top_phyla$ncbi_genome_count[i]),
              scales::comma(top_phyla$census_size_count[i]),
              scales::comma(top_phyla$census_otu_count[i]),
              scales::comma(top_phyla$ncbi_species_count[i]),
              phyla_colors[i]))
}

cat(sprintf("%-25s %12s %12s %12s %12s %10s\n",
            "Other",
            scales::comma(other_genome_count),
            scales::comma(other_size_count),
            scales::comma(other_otu_count),
            scales::comma(other_species_count),
            other_color))

cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("✅", toupper(PROCESS_DOMAIN), "Sankey diagram with colored flows completed!\n")
cat("   📁 HTML file:", output_file, "\n")
cat("   💡 Transfer to local machine to view interactive diagram\n")
cat("   🎨 Flows are colored to match their phylum colors\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

