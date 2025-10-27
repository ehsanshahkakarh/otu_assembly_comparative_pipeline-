#!/usr/bin/env Rscript

# Install required R packages for comprehensive scatter visualizer

cat("Installing required R packages...\n")

# List of required packages
required_packages <- c(
  "ggplot2", 
  "dplyr", 
  "readr", 
  "ggrepel", 
  "RColorBrewer", 
  "viridis", 
  "scales", 
  "stringr", 
  "patchwork"
)

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org/", dependencies = TRUE)
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("✓", pkg, "installed successfully\n")
    } else {
      cat("✗ Failed to install", pkg, "\n")
    }
  } else {
    cat("✓", pkg, "already installed\n")
  }
}

# Install all required packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("\nPackage installation complete!\n")

# Test loading all packages
cat("\nTesting package loading...\n")
all_loaded <- TRUE
for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat("✓", pkg, "loaded successfully\n")
  }, error = function(e) {
    cat("✗ Failed to load", pkg, ":", e$message, "\n")
    all_loaded <<- FALSE
  })
}

if (all_loaded) {
  cat("\n🎉 All packages installed and loaded successfully!\n")
  cat("You can now run the comprehensive_scatter_visualizer.R script\n")
} else {
  cat("\n❌ Some packages failed to load. Please check the error messages above.\n")
}
