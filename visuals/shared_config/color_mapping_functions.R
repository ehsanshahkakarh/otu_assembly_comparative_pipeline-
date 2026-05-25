# Shared Color Mapping Functions
# ===============================
# Functions to load and apply taxonomic color mappings consistently
# across all visualization scripts (alluvial, scatter, etc.)

# Load required libraries
if (!require(yaml, quietly = TRUE)) {
  install.packages("yaml")
  library(yaml)
}

# Load colors from comprehensive visual config or standalone color mapping
load_taxonomic_colors <- function(config_path = NULL) {
  # Try multiple config file options
  config_options <- c(
    "scatter_plots/config/comprehensive_visual_config.yaml",
    "../scatter_plots/config/comprehensive_visual_config.yaml",
    "../../scatter_plots/config/comprehensive_visual_config.yaml",
    "shared_config/taxonomic_color_mapping.yaml",
    "../shared_config/taxonomic_color_mapping.yaml",
    "../../shared_config/taxonomic_color_mapping.yaml"
  )

  if (!is.null(config_path)) {
    config_options <- c(config_path, config_options)
  }

  config_file <- NULL
  for (path in config_options) {
    if (file.exists(path)) {
      config_file <- path
      break
    }
  }

  if (is.null(config_file)) {
    stop(paste("Could not find color configuration file. Tried paths:",
               paste(config_options, collapse = ", ")))
  }

  cat(paste("Loading colors from:", config_file, "\n"))
  config <- yaml.load_file(config_file)

  # Handle comprehensive visual config format
  if ("color_config" %in% names(config) && config$color_config$use_shared_colors) {
    # Load the referenced shared color file
    shared_file <- config$color_config$shared_color_file

    # Try multiple possible paths for the shared file
    shared_paths <- c(
      file.path(dirname(config_file), shared_file),
      shared_file,
      file.path("visuals", shared_file),
      file.path("..", shared_file),
      file.path("../..", shared_file)
    )

    for (shared_path in shared_paths) {
      if (file.exists(shared_path)) {
        cat(paste("Loading shared colors from:", shared_path, "\n"))
        return(yaml.load_file(shared_path))
      }
    }

    cat("Warning: Could not find shared color file, using fallback colors\n")
  }

  return(config)
}

# Deterministic fallback index based on taxon name.
# This keeps unmapped taxa color-consistent across different plots/scripts.
get_fallback_index <- function(taxon_name, pool_size) {
  if (is.na(taxon_name) || !nzchar(taxon_name) || pool_size <= 0) {
    return(1)
  }

  key <- tolower(trimws(taxon_name))
  chars <- utf8ToInt(key)

  if (length(chars) == 0) {
    return(1)
  }

  # Weighted checksum to reduce collisions for similar names.
  checksum <- sum(chars * seq_along(chars))
  ((checksum - 1) %% pool_size) + 1
}

# Get colors for a specific domain
get_domain_colors <- function(taxa_names, domain, color_config = NULL) {
  if (is.null(color_config)) {
    color_config <- load_taxonomic_colors()
  }
  
  # Get domain-specific color mappings
  domain_key <- paste0(tolower(domain), "_colors")
  domain_colors <- color_config[[domain_key]]
  fallback_colors <- color_config$fallback_colors[[tolower(domain)]]
  special_colors <- color_config$special_colors
  
  if (is.null(domain_colors)) {
    stop(paste("No color mapping found for domain:", domain))
  }
  
  # Initialize result vector
  result_colors <- character(length(taxa_names))
  names(result_colors) <- taxa_names
  
  for (i in seq_along(taxa_names)) {
    taxon <- taxa_names[i]
    
    # Clean taxon name (remove numbering, etc.)
    clean_taxon <- gsub("^[0-9]+\\. ", "", taxon)  # Remove "1. " prefixes
    clean_taxon <- trimws(clean_taxon)
    
    # Apply color assignment rules
    assigned_color <- NULL
    
    # 1. Check special cases first
    if (taxon == "Other" || clean_taxon == "Other") {
      assigned_color <- special_colors$Other
    } else if (taxon == "Unknown" || clean_taxon == "Unknown") {
      assigned_color <- special_colors$Unknown
    } else if (grepl("\\.U\\.", clean_taxon)) {
      # Unclassified entries: try full taxon key first (e.g. "Bacteria.U.phylum",
      # "Amoebozoa.U.division") so each .U. bin can have its own muted tone,
      # then fall back to the per-domain key.
      if (clean_taxon %in% names(special_colors)) {
        assigned_color <- special_colors[[clean_taxon]]
      } else {
        special_key <- paste0("unclassified_", tolower(domain))
        assigned_color <- special_colors[[special_key]]
      }
    } else {
      # 2. Try exact match
      if (clean_taxon %in% names(domain_colors)) {
        assigned_color <- domain_colors[[clean_taxon]]
      } else {
        # 3. Try case-insensitive match
        match_idx <- which(tolower(names(domain_colors)) == tolower(clean_taxon))
        if (length(match_idx) > 0) {
          assigned_color <- domain_colors[[match_idx[1]]]
        }
      }
    }
    
    # 4. Use fallback colors if no match found
    if (is.null(assigned_color) && !is.null(fallback_colors)) {
      fallback_idx <- get_fallback_index(clean_taxon, length(fallback_colors))
      assigned_color <- fallback_colors[fallback_idx]
    }
    
    # 5. Final fallback to gray
    if (is.null(assigned_color)) {
      assigned_color <- "#808080"
    }
    
    result_colors[i] <- assigned_color
  }
  
  return(result_colors)
}

# Convenience function to get all colors for mixed domain data
get_mixed_domain_colors <- function(taxa_data, taxon_col = "taxon", domain_col = "domain", 
                                   color_config = NULL) {
  if (is.null(color_config)) {
    color_config <- load_taxonomic_colors()
  }
  
  result_colors <- character(nrow(taxa_data))
  
  for (domain in unique(taxa_data[[domain_col]])) {
    domain_rows <- taxa_data[[domain_col]] == domain
    domain_taxa <- taxa_data[[taxon_col]][domain_rows]
    domain_colors <- get_domain_colors(domain_taxa, domain, color_config)
    result_colors[domain_rows] <- domain_colors
  }
  
  return(result_colors)
}

# Print color mapping summary
print_color_summary <- function(taxa_names, colors, domain = NULL) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  if (!is.null(domain)) {
    cat(paste("COLOR MAPPING SUMMARY -", toupper(domain), "\n"))
  } else {
    cat("COLOR MAPPING SUMMARY\n")
  }
  cat(paste(rep("=", 50), collapse = ""), "\n")

  for (i in seq_along(taxa_names)) {
    cat(sprintf("%2d. %-25s %s\n", i, taxa_names[i], colors[i]))
  }
  cat("\n")
}
