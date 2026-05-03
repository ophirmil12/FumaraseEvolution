#!/usr/bin/env Rscript

# Usage:
# Rscript tree_plot_circular.R <metadata.txt> <tree.nwk> <output.png> [dpi]
# Example:
# Rscript tree_plot_circular.R data/processed/metadata.txt results/trees/tree.nwk results/figures/phylogeny_tree.png 600

suppressPackageStartupMessages({
  library(ggnewscale)
  library(phangorn)
  library(ggtree)
  library(ggtreeExtra)
  library(tidyverse)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(3, 4))) {
  stop("Usage: Rscript tree_plot_circular.R <metadata.txt> <tree.nwk> <output.png> [dpi]")
}

metadata_file <- args[1]
tree_file     <- args[2]
output_file   <- args[3]
# Default DPI set to 600
dpi <- ifelse(length(args) == 4, as.numeric(args[4]), 600)

if (is.na(dpi) || dpi < 72) {
  stop("DPI must be a valid number (e.g., >= 72)")
}

# --- Tree ---
tree <- read.tree(tree_file)
tree$tip.label <- sapply(strsplit(tree$tip.label, "\\|"), `[`, 1)       # Remove any suffix after '|' in tip labels
tree <- midpoint(tree)                                                  # Midpoint root the tree
tree_ultra <- compute.brlen(tree, method = "Grafen")                # Compute branch lengths using Grafen's method for better visualization

# --- Metadata ---
# Read the TSV file
raw_metadata <- read.delim(metadata_file, sep = "\t", stringsAsFactors = FALSE)

# Initialize the metadata dataframe that ggtree expects
metadata <- data.frame(
  label = raw_metadata$proteome_id,
  Group = NA_character_,
  stringsAsFactors = FALSE
)

# Extract the group from the taxonomic_lineage column using keyword matching
lineage <- raw_metadata$taxonomic_lineage
metadata$Group[grepl("Eukaryota", lineage, ignore.case = TRUE)] <- "Eukaryote"
metadata$Group[grepl("Archaea", lineage, ignore.case = TRUE)] <- "Archaea"
metadata$Group[grepl("Bacteria", lineage, ignore.case = TRUE)] <- "Bacteria"
# Apply Alphaproteobacteria AFTER Bacteria so it overwrites the general Bacteria label
metadata$Group[grepl("Alphaproteobacteria", lineage, ignore.case = TRUE)] <- "Alphaproteobacteria"

# Palette
group_colors <- c(
  "Eukaryote"           = "#1f77b4",
  "Alphaproteobacteria" = "#ff7f0e",
  "Bacteria"            = "#A6CEE3",
  "Archaea"             = "#d62728"
)

metadata$Group <- factor(metadata$Group, levels = names(group_colors))

# Diagnostic Check: See how many tips are actually missing from the TSV
missing_tips <- sum(!(tree$tip.label %in% metadata$label))
if(missing_tips > 0) {
  cat("    WARNING:", missing_tips, "tips on the tree were not found in the TSV file.\n")
}

# --- Plotting ---
p <- ggtree(tree_ultra, layout = "circular") %<+% metadata +
  geom_tree(aes(color = Group), linewidth = 0.6) +
  scale_color_manual(values = group_colors, guide = "none", na.translate = FALSE) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Group),
    width = 0.5,
    offset = -0.25,
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors, na.translate = FALSE) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Group),
    width = 0.05,
    offset = 0.3
  ) +
  scale_fill_manual(values = group_colors, name = "Taxonomic Group", na.translate = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.margin = margin(10, 10, 10, 10)
  )

# Ensure the output directory exists
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

ggsave(output_file, plot = p, width = 12, height = 12, dpi = dpi)
cat("Done! Tree for: ", tree_file, " saved to:", output_file, "\n")