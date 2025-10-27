# Test family column issue
cat("Testing family column issue...\n")

# Create test data
test_data <- data.frame(
  genus = c("Saccharomyces", "Aspergillus"),
  count = c(10, 20),
  stringsAsFactors = FALSE
)

# Add family column
test_data$family <- c("Saccharomycetales", "Aspergillaceae")

cat("Test data columns:", paste(colnames(test_data), collapse=", "), "\n")

# Define essential columns
essential_cols <- c("genus", "family", "count")

# Test intersect
available_cols <- intersect(essential_cols, colnames(test_data))
cat("Available essential columns:", paste(available_cols, collapse=", "), "\n")

# Filter data
filtered_data <- test_data[, available_cols]
cat("Filtered data columns:", paste(colnames(filtered_data), collapse=", "), "\n")

# Write test file
write.table(filtered_data, "test_family_output.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("Test file written with family column\n")
