# Test column filtering
cat("Testing column filtering...\n")

# Read the actual data
data <- read.csv("../Eukcensus_merge/18s_merged/csv_results/18s_ncbi_merged_clean_genus.csv")
cat("Original data columns:", paste(colnames(data), collapse=", "), "\n")

# Define essential columns
essential_cols <- c("genus", "census_otu_count", "ncbi_genome_count", "ncbi_species_count", 
                   "isolate_count", "isolate_percentage", "novelty_factor", 
                   "overrepresentation_factor", "domain")

# Test filtering
available_cols <- intersect(essential_cols, colnames(data))
cat("Available essential columns:", paste(available_cols, collapse=", "), "\n")

# Filter data
filtered_data <- data[1:3, available_cols]
cat("Filtered data columns:", paste(colnames(filtered_data), collapse=", "), "\n")

# Write test file
write.table(filtered_data, "test_filtered.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("Test file written: test_filtered.tsv\n")
