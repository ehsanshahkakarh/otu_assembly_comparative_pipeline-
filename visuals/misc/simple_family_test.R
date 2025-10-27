# Simple family test
cat("Simple family test...\n")

# Create test data that mimics eukaryota structure
test_data <- data.frame(
  genus = c("Saccharomyces", "Aspergillus", "Penicillium"),
  census_otu_count = c(13, 48, 100),
  ncbi_genome_count = c(1936, 1433, 559),
  ncbi_species_count = c(21, 263, 174),
  isolate_count = c(1929, 1431, 558),
  isolate_percentage = c(99.64, 99.86, 99.82),
  novelty_factor = c(0.619, 0.183, 0.575),
  overrepresentation_factor = c(1.615, 5.479, 1.74),
  domain = c("Eukaryota", "Eukaryota", "Eukaryota"),
  Division = c("Opisthokonta", "Opisthokonta", "Opisthokonta"),
  Is_Top_Novelty = c(FALSE, FALSE, FALSE),
  Is_Top_Coverage = c(TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Add family column
test_data$family <- c("Saccharomycetales", "Aspergillaceae", "Trichocomaceae")
test_data$Domain <- "Eukaryota"

cat("Test data columns:", paste(colnames(test_data), collapse=", "), "\n")

# Define essential columns (same as in script)
essential_cols <- c("genus", "family", "census_otu_count", "ncbi_genome_count", "ncbi_species_count", 
                   "isolate_count", "isolate_percentage", "novelty_factor", 
                   "overrepresentation_factor", "domain", "Division", "Is_Top_Novelty", 
                   "Is_Top_Coverage", "Domain")

# Filter columns
available_cols <- intersect(essential_cols, colnames(test_data))
cat("Available columns:", paste(available_cols, collapse=", "), "\n")

test_export <- test_data[, available_cols]
cat("Export columns:", paste(colnames(test_export), collapse=", "), "\n")

# Write test file
write.table(test_export, "simple_family_test.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("Test file written: simple_family_test.tsv\n")
