#!/usr/bin/env Rscript

# Install required packages for scatter plot analysis

# List of required packages
required_packages <- c(
  "ggplot2",
  "dplyr",
  "scales",
  "viridis",
  "gridExtra",
  "RColorBrewer",
  "ggrepel",
  "ggtext"
)

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("Installing package:", pkg, "\n"))
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
    } else {
      cat(paste("Package already installed:", pkg, "\n"))
    }
  }
}

# Install packages
cat("Checking and installing required R packages...\n")
install_if_missing(required_packages)

# Test loading all packages
cat("\nTesting package loading...\n")
success <- TRUE
for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(paste("✓", pkg, "loaded successfully\n"))
  }, error = function(e) {
    cat(paste("✗", pkg, "failed to load:", e$message, "\n"))
    success <<- FALSE
  })
}

if (success) {
  cat("\n✓ All packages installed and loaded successfully!\n")
  cat("You can now run the scatter plot analysis.\n")
} else {
  cat("\n✗ Some packages failed to install or load.\n")
  cat("Please check the error messages above.\n")
}
