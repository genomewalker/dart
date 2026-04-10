#!/usr/bin/env Rscript
# AGP damage estimates vs metaDMG ground truth
# Reproduces damage_validation_scatter.png from the DART wiki.
#
# Input:
#   data/damage_comparison/*.json        — DART JSON output per sample
#   data/damage_comparison/comparison_table.tsv  — pre-computed merge table
#
# Usage:
#   Rscript scripts/plot_damage_comparison.R
#
# Output:
#   figures/damage_validation_scatter.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(jsonlite)
})

script_dir  <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _), mustWork = FALSE))
if (length(script_dir) == 0 || is.na(script_dir)) script_dir <- "."

data_dir <- file.path(script_dir, "..", "data", "damage_comparison")
out_dir  <- file.path(script_dir, "..", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load pre-computed comparison table (generated from all 31 samples, including those
# whose JSON files are not redistributed due to size).
# To recompute from raw JSON: re-run with all JSON files present.
tsv_path <- file.path(data_dir, "agp_vs_metadmg_results.tsv")
dt_raw <- fread(tsv_path)
setnames(dt_raw,
    old = c("metaDMG_pct", "AGP_d_max", "AGP_validated"),
    new = c("metadmg_dmax", "d_max", "validated"))
dt       <- dt_raw
dt_valid <- dt[validated == TRUE]

r_val <- cor(dt_valid$d_max, dt_valid$metadmg_dmax)
cat("Samples:", nrow(dt), "| Validated:", nrow(dt_valid), "\n")
cat("Pearson r (d_max vs metaDMG):", round(r_val, 3), "\n")
cat("Mean bias (AGP - metaDMG):", round(mean(dt_valid$d_max - dt_valid$metadmg_dmax), 2), "%\n")

# Scatter: d_max vs metaDMG
p <- ggplot(dt, aes(x = metadmg_dmax, y = d_max, shape = validated)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(size = 3, aes(color = validated), alpha = 0.8) +
    geom_smooth(data = dt_valid, method = "lm", se = TRUE, color = "#E69F00", linewidth = 0.9) +
    annotate("text", x = 5, y = 55,
             label = sprintf("r = %.2f (n = %d validated)", r_val, nrow(dt_valid)),
             hjust = 0, size = 4.5) +
    scale_color_manual(values = c("FALSE" = "#999999", "TRUE" = "#0072B2"),
                       labels = c("FALSE" = "Not validated (Channel B)", "TRUE" = "Validated")) +
    scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16),
                       labels = c("FALSE" = "Not validated (Channel B)", "TRUE" = "Validated")) +
    labs(
        x      = "metaDMG d_max (%)",
        y      = "DART d_max (%)",
        color  = NULL,
        shape  = NULL,
        title  = "DART damage vs metaDMG (31 ancient metagenomes)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom") +
    coord_fixed(xlim = c(0, 60), ylim = c(0, 60))

out_pdf <- file.path(out_dir, "damage_validation_scatter.pdf")
out_png <- file.path(out_dir, "damage_validation_scatter.png")
ggsave(out_pdf, p, width = 6, height = 6)
ggsave(out_png, p, width = 6, height = 6, dpi = 300)
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
