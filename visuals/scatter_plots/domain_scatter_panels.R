#!/usr/bin/env Rscript
# Domain-Specific Scatter Panels
# ===============================
# Generates three identically formatted scatter plots (one per domain)
# designed for stitching together in Keynote/PowerPoint.
#
# Each output: Phylum row (top) + Family row (bottom) + color legend
# Format matches mega_comprehensive_stacked_visual_v2.R exactly:
#   - Same axis limits (1–10000 log-log), dashed diagonal
#   - Same circle sizing (25/20/15/10 by isolate %), same shape/stroke/alpha
#   - Same label format (phylum: factor only; family: name + factor)
#   - Same text repel parameters, same black panel borders
#   - Same top-10 filtration (threshold > 1.0)
#
# Run with: conda activate r_env && Rscript domain_scatter_panels.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggrepel)
  library(yaml)
  library(grid)
  library(cowplot)
})

# ── Configuration (matches v2 exactly) ────────────────────────────────────────
SOURCE_DIR   <- file.path("mega_scrip", "source_data")
COLOR_CONFIG <- file.path("..", "shared_config", "taxonomic_color_mapping.yaml")
OUTPUT_DIR   <- "domain_panels"

PANEL_W <- 18   # inches width per domain figure
PANEL_H <- 20   # inches height per domain figure
DPI     <- 300

TOP_N     <- 10
THRESHOLD <- 1.0

# Match v2 plot config exactly
PLOT_CFG <- list(
  text_size    = 11.5,
  size_range   = c(10, 22),
  circle_shape = 21,
  circle_stroke = 0.6,
  circle_alpha = 0.9,
  bg_alpha     = 0.3
)

# ── Color assignment (matches v2 global registry) ─────────────────────────────
color_config <- yaml::read_yaml(COLOR_CONFIG)
TAXON_COLORS <- list()

assign_taxon_color <- function(taxon, domain) {
  if (taxon %in% names(TAXON_COLORS)) return(TAXON_COLORS[[taxon]])

  color <- NULL
  if (domain == "Bacteria"  && taxon %in% names(color_config$bacteria_colors))
    color <- color_config$bacteria_colors[[taxon]]
  if (domain == "Archaea"   && taxon %in% names(color_config$archaea_colors))
    color <- color_config$archaea_colors[[taxon]]
  if (domain == "Eukaryota" && taxon %in% names(color_config$eukaryota_colors))
    color <- color_config$eukaryota_colors[[taxon]]

  if (is.null(color)) {
    # Cross-domain fallback pool (same logic as v2)
    if (domain == "Bacteria") {
      pool <- c(unlist(color_config$archaea_colors),
                unlist(color_config$eukaryota_colors),
                unlist(color_config$fallback_colors$archaea),
                unlist(color_config$fallback_colors$eukaryota))
    } else if (domain == "Eukaryota") {
      pool <- c(unlist(color_config$bacteria_colors),
                unlist(color_config$archaea_colors),
                unlist(color_config$fallback_colors$bacteria),
                unlist(color_config$fallback_colors$archaea))
    } else {
      pool <- c(unlist(color_config$bacteria_colors),
                unlist(color_config$eukaryota_colors),
                unlist(color_config$fallback_colors$bacteria),
                unlist(color_config$fallback_colors$eukaryota))
    }
    pool <- unique(pool)
    used <- unlist(TAXON_COLORS)
    avail <- setdiff(pool, used)
    if (length(avail) > 0) {
      color <- avail[1]
    } else {
      color <- pool[(length(TAXON_COLORS) %% length(pool)) + 1]
    }
  }

  TAXON_COLORS[[taxon]] <<- color
  return(color)
}

init_color_registry <- function() { TAXON_COLORS <<- list() }

# ── Read source data ──────────────────────────────────────────────────────────
read_source <- function(domain, level) {
  fname <- paste0(domain, "_", level, "_source_data.csv")
  fpath <- file.path(SOURCE_DIR, fname)
  if (!file.exists(fpath)) stop("Missing: ", fpath)
  df <- read.csv(fpath, stringsAsFactors = FALSE)
  # Ensure Circle_Size exists (same formula as v2)
  if (!"Circle_Size" %in% names(df)) {
    iso_pct <- ifelse(df$ncbi_genome_count > 0,
                      (df$isolate_count / df$ncbi_genome_count) * 100, 0)
    df$Circle_Size <- case_when(
      iso_pct == 0  ~ 25,
      iso_pct < 10  ~ 20,
      iso_pct < 50  ~ 15,
      TRUE          ~ 10
    )
  }
  df
}

# ── Scatter plot (matches v2 create_individual_scatter exactly) ───────────────
create_individual_scatter <- function(data, level, domain) {
  top_data <- data[data$Is_Top_Novelty | data$Is_Top_Coverage, ]
  bg_data  <- data[!data$Is_Top_Novelty & !data$Is_Top_Coverage, ]

  # Base plot — same axis limits as v2
  p <- ggplot(data, aes(x = census_otu_count, y = ncbi_species_count)) +
    scale_x_log10(labels = comma_format(), limits = c(1, 10000)) +
    scale_y_log10(labels = comma_format(), limits = c(1, 10000))

  # Background points
  if (nrow(bg_data) > 0) {
    p <- p + geom_point(data = bg_data, aes(size = Circle_Size),
                        color = "lightgray", fill = "lightgray",
                        shape = PLOT_CFG$circle_shape, alpha = PLOT_CFG$bg_alpha,
                        stroke = PLOT_CFG$circle_stroke)
  }

  # Highlighted points with color by Phylum/Division
  if (nrow(top_data) > 0) {
    if (domain == "Eukaryota") {
      color_col <- "Division"
      plot_groups <- sort(unique(top_data$Division[!top_data$Division %in% c("Unknown", "", "Other", NA)]))
    } else {
      color_col <- "Phylum"
      plot_groups <- sort(unique(top_data$Phylum[!top_data$Phylum %in% c("Unknown", "", "Other", NA)]))
    }

    group_colors <- sapply(plot_groups, function(taxon) assign_taxon_color(taxon, domain))
    names(group_colors) <- plot_groups

    p <- p + geom_point(data = top_data,
                        aes_string(size = "Circle_Size",
                                   fill = paste0("factor(", color_col, ", levels = plot_groups)")),
                        color = "black",
                        shape = PLOT_CFG$circle_shape, alpha = PLOT_CFG$circle_alpha,
                        stroke = PLOT_CFG$circle_stroke) +
      scale_fill_manual(values = group_colors, guide = "none")
  }

  # Styling — matches v2 exactly
  p <- p +
    scale_size_continuous(range = PLOT_CFG$size_range, guide = "none") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, color = "grey50"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )

  # Labels — same repel params as v2
  if (nrow(top_data) > 0) {
    novelty_data <- top_data[top_data$Is_Top_Novelty == TRUE, ]
    overrep_data <- top_data[top_data$Is_Top_Coverage == TRUE, ]

    if (nrow(novelty_data) > 0) {
      novelty_data <- novelty_data[order(-novelty_data$novelty_factor), ]
      if (level == "phylum") {
        novelty_data$label <- paste0(sprintf("%.1f", novelty_data$novelty_factor), "\u00d7")
      } else {
        novelty_data$label <- paste0(novelty_data$Taxon, " (", sprintf("%.1f", novelty_data$novelty_factor), "\u00d7)")
      }
      p <- p + ggrepel::geom_text_repel(
        data = novelty_data, aes(label = label),
        size = PLOT_CFG$text_size * 0.9, fontface = "bold", color = "black",
        max.overlaps = Inf, force = 4, force_pull = 0.1,
        box.padding = 1.0, point.padding = 0.5,
        segment.color = "gray40", segment.size = 0.3, segment.alpha = 0.8,
        min.segment.length = 0, direction = "both",
        nudge_y = -0.5, nudge_x = 0,
        xlim = c(NA, NA), ylim = c(NA, NA), seed = 42
      )
    }

    if (nrow(overrep_data) > 0) {
      overrep_data <- overrep_data[order(-overrep_data$overrepresentation_factor), ]
      if (level == "phylum") {
        overrep_data$label <- paste0(sprintf("%.1f", overrep_data$overrepresentation_factor), "\u00d7")
      } else {
        overrep_data$label <- paste0(overrep_data$Taxon, " (", sprintf("%.1f", overrep_data$overrepresentation_factor), "\u00d7)")
      }
      p <- p + ggrepel::geom_text_repel(
        data = overrep_data, aes(label = label),
        size = PLOT_CFG$text_size * 0.9, fontface = "bold", color = "black",
        max.overlaps = Inf, force = 4, force_pull = 0.1,
        box.padding = 1.0, point.padding = 0.5,
        segment.color = "gray40", segment.size = 0.3, segment.alpha = 0.8,
        min.segment.length = 0, direction = "both",
        nudge_y = 0.5, nudge_x = 0,
        xlim = c(NA, NA), ylim = c(NA, NA), seed = 123
      )
    }
  }

  return(p)
}

# ── Domain-specific legend (grid of colored squares) ──────────────────────────
create_domain_legend <- function(all_data, domain) {
  # Only collect Phylum/Division values from the phylum-level data (not family)
  phylum_df <- all_data[["phylum"]]
  groups <- c()
  if (!is.null(phylum_df)) {
    top <- phylum_df[phylum_df$Is_Top_Novelty | phylum_df$Is_Top_Coverage, ]
    if (domain == "Eukaryota") {
      groups <- c(groups, top$Division)
    } else {
      groups <- c(groups, top$Phylum)
    }
  }
  groups <- sort(unique(groups[!groups %in% c("Unknown", "", "Other", NA)]))
  if (length(groups) == 0) return(ggplot() + theme_void())

  # Get colors (already assigned during scatter creation)
  group_colors <- sapply(groups, function(g) {
    if (g %in% names(TAXON_COLORS)) TAXON_COLORS[[g]]
    else assign_taxon_color(g, domain)
  })

  n <- length(groups)
  ncols <- min(n, 4)
  nrows <- ceiling(n / ncols)

  legend_df <- data.frame(
    taxon = groups,
    color = group_colors,
    col   = rep(1:ncols, length.out = n),
    row   = rep(1:nrows, each = ncols, length.out = n),
    stringsAsFactors = FALSE
  )

  ggplot(legend_df, aes(x = col, y = -row)) +
    geom_point(aes(fill = color), shape = 22, size = 14, color = "black", stroke = 0.8) +
    geom_text(aes(label = taxon), hjust = 0, nudge_x = 0.35, size = 7,
              color = "grey20", fontface = "bold") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.5, ncols + 5), ylim = c(-nrows - 0.8, 0.8)) +
    theme_void() +
    theme(plot.margin = margin(15, 20, 15, 20))
}



# ── Main: generate one figure per domain ──────────────────────────────────────
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

domains <- c("Bacteria", "Archaea", "Eukaryota")
levels  <- c("phylum", "family")

for (domain in domains) {
  cat(paste0("\n═══ ", domain, " ═══\n"))
  init_color_registry()

  # Load data for both levels
  domain_data <- list()
  for (level in levels) {
    domain_data[[level]] <- read_source(domain, level)
    cat(paste("  Loaded", level, ":", nrow(domain_data[[level]]), "taxa\n"))
  }

  # Create scatter plots (this also populates TAXON_COLORS)
  plots <- list()
  for (level in levels) {
    plots[[level]] <- create_individual_scatter(domain_data[[level]], level, domain)
  }

  # Frame each plot with black border (same as v2)
  framed_plots <- lapply(plots, function(p) {
    ggdraw(p) + theme(plot.background = element_rect(color = "black", fill = NA, size = 3))
  })

  # Row labels
  level_labels <- list(
    phylum = ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Phyla",
               size = 12, fontface = "bold", color = "grey30", angle = 90) +
      theme_void(),
    family = ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Family",
               size = 12, fontface = "bold", color = "grey30", angle = 90) +
      theme_void()
  )

  # Domain title header
  title_row <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = domain,
             size = 14, fontface = "bold", color = "grey50") +
    theme_void()

  # Build data rows
  phylum_row <- plot_grid(level_labels$phylum, framed_plots$phylum,
                          ncol = 2, rel_widths = c(0.08, 1))
  family_row <- plot_grid(level_labels$family, framed_plots$family,
                          ncol = 2, rel_widths = c(0.08, 1))

  # Legend
  legend_panel <- create_domain_legend(domain_data, domain)

  # Combine vertically: title / phylum / gap / family / legend
  complete <- plot_grid(
    title_row,
    phylum_row,
    ggplot() + theme_void(),  # spacer
    family_row,
    ggplot() + theme_void(),  # spacer
    legend_panel,
    ncol = 1,
    rel_heights = c(0.06, 1, 0.03, 1, 0.03, 0.50)
  ) + theme(plot.margin = margin(10, 40, 10, 10))  # extra right-side padding

  # Save PNG + PDF
  out_png <- file.path(OUTPUT_DIR, paste0(tolower(domain), "_scatter_panel.png"))
  out_pdf <- file.path(OUTPUT_DIR, paste0(tolower(domain), "_scatter_panel.pdf"))

  ggsave(out_png, complete, width = PANEL_W, height = PANEL_H,
         dpi = DPI, bg = "white", limitsize = FALSE)
  ggsave(out_pdf, complete, width = PANEL_W, height = PANEL_H,
         bg = "white", limitsize = FALSE)

  cat(paste("  ✓ Saved:", out_png, "\n"))
  cat(paste("  ✓ Saved:", out_pdf, "\n"))
}

cat("\n✅ All three domain panels generated in:", OUTPUT_DIR, "\n")
