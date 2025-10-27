# Test eukaryota export with family
cat("Testing eukaryota export with family...\n")

source('mega_genus_comprehensive_visual.R')

# Load data
eukaryota_data <- load_18s_data('genus')
eukaryota_top <- eukaryota_data[eukaryota_data$Is_Top_Novelty | eukaryota_data$Is_Top_Coverage, ]

# Load family mapping and add family info
family_mapping <- load_family_mapping()
eukaryota_top <- add_family_info(eukaryota_top, family_mapping, '18s')
eukaryota_top$Domain <- 'Eukaryota'

cat("Columns after adding family:", paste(colnames(eukaryota_top), collapse=", "), "\n")
cat("Sample family values:", paste(head(eukaryota_top$family, 5), collapse=", "), "\n")

# Define essential columns
essential_cols <- c("genus", "family", "census_otu_count", "ncbi_genome_count", "ncbi_species_count", 
                   "isolate_count", "isolate_percentage", "novelty_factor", 
                   "overrepresentation_factor", "domain", "Division", "Is_Top_Novelty", 
                   "Is_Top_Coverage", "Domain")

# Filter columns
available_cols <- intersect(essential_cols, colnames(eukaryota_top))
cat("Available columns:", paste(available_cols, collapse=", "), "\n")

eukaryota_export <- eukaryota_top[, available_cols]
cat("Export columns:", paste(colnames(eukaryota_export), collapse=", "), "\n")

# Write test file
write.table(eukaryota_export, "test_euk_with_family.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("Test file written: test_euk_with_family.tsv\n")
