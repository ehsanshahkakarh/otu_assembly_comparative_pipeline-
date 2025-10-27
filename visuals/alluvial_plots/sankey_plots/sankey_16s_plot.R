#!/usr/bin/env Rscript
#
# 16S Prokaryotic Sankey Diagram Generator
# =========================================
#
# Creates a Sankey diagram showing the flow from:
# NCBI Total Genomes → 16S EukCensus Sequences → 16S EukCensus OTUs → NCBI Total Species
# 
# For prokaryotic (16S) data only
# Shows top 15 shared phyla by total representation + Other
#
# Author: Enhanced EukCensus Analysis Team
# Date: 2024
#

# Load required libraries
library(dplyr)
library(networkD3)
library(RColorBrewer)
library(viridis)
library(scales)
library(htmlwidgets)

cat("=== 16S Prokaryotic Sankey Diagram Generator ===\n\n")

# Load the 16S data
cat("Loading 16S data...\n")
data_16s <- read.csv("merged_output/16s_ncbi_merged_phylum.csv", stringsAsFactors = FALSE)
cat("16S data loaded:", nrow(data_16s), "rows\n")

# Filter for shared phyla only (in_both = 1) and check all required columns
shared_data <- data_16s %>%
  filter(in_both == 1) %>%
  filter(!is.na(Genome_Percentage) & !is.na(Size_Percentage) &
         !is.na(OTU_Percentage) & !is.na(Species_Coverage))

cat("Shared phyla found:", nrow(shared_data), "\n")

# Add total representation across all 4 nodes and get top 15
shared_data$total_representation <- shared_data$Genome_Percentage +
                                   shared_data$Size_Percentage +
                                   shared_data$OTU_Percentage +
                                   shared_data$Species_Coverage

top_phyla <- shared_data %>%
  arrange(desc(total_representation)) %>%
  head(15)

cat("Top 15 phyla selected:", nrow(top_phyla), "\n")

# Calculate total counts for each node from SHARED DATA ONLY
total_genome_count <- sum(shared_data$NCBI_Genome_Count, na.rm = TRUE)
total_size_count <- sum(shared_data$Census_Size_Count, na.rm = TRUE)
total_otu_count <- sum(shared_data$Census_OTU_Count, na.rm = TRUE)
total_species_count <- sum(shared_data$NCBI_Species_Count, na.rm = TRUE)

# Calculate "Other" counts (remaining phyla not in top 15)
other_genome_count <- total_genome_count - sum(top_phyla$NCBI_Genome_Count)
other_size_count <- total_size_count - sum(top_phyla$Census_Size_Count)
other_otu_count <- total_otu_count - sum(top_phyla$Census_OTU_Count)
other_species_count <- total_species_count - sum(top_phyla$NCBI_Species_Count)

cat("Other counts (remaining phyla):\n")
cat("  Other Genomes:", scales::comma(other_genome_count), "\n")
cat("  Other Sequences:", scales::comma(other_size_count), "\n")
cat("  Other OTUs:", scales::comma(other_otu_count), "\n")
cat("  Other Species:", scales::comma(other_species_count), "\n")

# Create nodes dataframe
# Node structure: Each phylum appears in each stage
phyla_names <- c(top_phyla$Phylum, "Other")
stages <- c("Genomes", "Sequences", "OTUs", "Species")

# Create all combinations of phyla and stages
nodes <- data.frame()
node_id <- 0

for (stage in stages) {
  for (phylum in phyla_names) {
    nodes <- rbind(nodes, data.frame(
      id = node_id,
      name = paste0(phylum, "_", stage),
      phylum = phylum,
      stage = stage,
      stringsAsFactors = FALSE
    ))
    node_id <- node_id + 1
  }
}

cat("Created", nrow(nodes), "nodes\n")

# Create links dataframe
links <- data.frame()

# Helper function to get node ID
get_node_id <- function(phylum, stage) {
  return(nodes$id[nodes$phylum == phylum & nodes$stage == stage])
}

# Create links between consecutive stages
for (i in 1:nrow(top_phyla)) {
  phylum <- top_phyla$Phylum[i]
  
  # Genomes -> Sequences
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "Genomes"),
    target = get_node_id(phylum, "Sequences"),
    value = top_phyla$Census_Size_Count[i],
    phylum = phylum,
    stringsAsFactors = FALSE
  ))
  
  # Sequences -> OTUs
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "Sequences"),
    target = get_node_id(phylum, "OTUs"),
    value = top_phyla$Census_OTU_Count[i],
    phylum = phylum,
    stringsAsFactors = FALSE
  ))
  
  # OTUs -> Species
  links <- rbind(links, data.frame(
    source = get_node_id(phylum, "OTUs"),
    target = get_node_id(phylum, "Species"),
    value = top_phyla$NCBI_Species_Count[i],
    phylum = phylum,
    stringsAsFactors = FALSE
  ))
}

# Add "Other" links
links <- rbind(links, data.frame(
  source = get_node_id("Other", "Genomes"),
  target = get_node_id("Other", "Sequences"),
  value = other_size_count,
  phylum = "Other",
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "Sequences"),
  target = get_node_id("Other", "OTUs"),
  value = other_otu_count,
  phylum = "Other",
  stringsAsFactors = FALSE
))

links <- rbind(links, data.frame(
  source = get_node_id("Other", "OTUs"),
  target = get_node_id("Other", "Species"),
  value = other_species_count,
  phylum = "Other",
  stringsAsFactors = FALSE
))

cat("Created", nrow(links), "links\n")

# Generate color palette
n_phyla <- length(phyla_names)
if (n_phyla <= 8) {
  colors <- brewer.pal(max(3, n_phyla), "Set2")
} else if (n_phyla <= 12) {
  colors <- brewer.pal(max(3, n_phyla), "Set3")
} else {
  base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
  if (n_phyla > 16) {
    additional_colors <- viridis(n_phyla - 16, option = "plasma")
    base_colors <- c(base_colors, additional_colors)
  }
  colors <- base_colors[1:n_phyla]
}

# Set "Other" to gray
other_index <- which(phyla_names == "Other")
if (length(other_index) > 0) {
  colors[other_index] <- "#CCCCCC"
}

# Assign colors to nodes
nodes$color <- colors[match(nodes$phylum, phyla_names)]

cat("Generated", length(colors), "colors for", n_phyla, "phyla\n")

# Create the Sankey diagram
cat("Creating Sankey diagram...\n")

sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "phylum",
  colourScale = paste0('d3.scaleOrdinal().domain([', 
                      paste0('"', phyla_names, '"', collapse = ','), 
                      ']).range([', 
                      paste0('"', colors, '"', collapse = ','), 
                      '])'),
  fontSize = 12,
  fontFamily = "Arial",
  nodeWidth = 20,
  nodePadding = 10,
  margin = list(top = 50, right = 50, bottom = 50, left = 50),
  height = 800,
  width = 1200,
  iterations = 100
)

# Save the Sankey diagram
output_file <- "sankey_16s_ncbi_abs.html"
saveWidget(sankey_plot, file = output_file, selfcontained = TRUE)

cat("Sankey diagram saved:", output_file, "\n")

# Try to save a static PNG version using webshot (if available)
tryCatch({
  library(webshot)
  png_file <- "sankey_16s_ncbi_abs.png"
  webshot(output_file, png_file, vwidth = 1200, vheight = 800)
  cat("Static PNG version saved:", png_file, "\n")
}, error = function(e) {
  cat("Note: webshot package not available for PNG export\n")
  cat("HTML file can be transferred to local machine for viewing\n")
})

# Display summary statistics
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("🧬 16S PROKARYOTIC SANKEY DIAGRAM ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("📊 DATASET OVERVIEW:\n")
cat("   • Total prokaryotic phyla analyzed:", nrow(top_phyla), "+ Other\n")
cat("   • Data flow: Genomes → Sequences → OTUs → Species\n")
cat("   • Interactive HTML output with hover details\n\n")

cat("📈 TOTAL COUNTS:\n")
cat("   • Total genomes:", scales::comma(total_genome_count), "\n")
cat("   • Total sequences:", scales::comma(total_size_count), "\n")
cat("   • Total OTUs:", scales::comma(total_otu_count), "\n")
cat("   • Total species:", scales::comma(total_species_count), "\n\n")

cat("🎨 VISUALIZATION FEATURES:\n")
cat("   • Interactive Sankey diagram with hover details\n")
cat("   • Color-coded phyla with", length(colors), "distinct colors\n")
cat("   • Flow widths represent absolute counts\n")
cat("   • HTML output for web viewing and sharing\n\n")

# Create a text-based flow summary
cat("📋 FLOW SUMMARY (Top 10 Phyla):\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-20s %12s %12s %12s %12s\n", "Phylum", "Genomes", "Sequences", "OTUs", "Species"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:min(10, nrow(top_phyla))) {
  cat(sprintf("%-20s %12s %12s %12s %12s\n",
              substr(top_phyla$Phylum[i], 1, 20),
              scales::comma(top_phyla$NCBI_Genome_Count[i]),
              scales::comma(top_phyla$Census_Size_Count[i]),
              scales::comma(top_phyla$Census_OTU_Count[i]),
              scales::comma(top_phyla$NCBI_Species_Count[i])))
}

cat(sprintf("%-20s %12s %12s %12s %12s\n",
            "Other",
            scales::comma(other_genome_count),
            scales::comma(other_size_count),
            scales::comma(other_otu_count),
            scales::comma(other_species_count)))

cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("✅ 16S prokaryotic Sankey diagram generation completed!\n")
cat("   📁 HTML file:", output_file, "\n")
cat("   💡 Transfer to local machine to view interactive diagram\n")
cat("   📊 Text summary above shows the data flows\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
