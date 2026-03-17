#!/usr/bin/env Rscript

# Usage:
# Rscript tree_plot_unrooted.R <metadata.tsv> <tree.nwk> <output.png> [dpi]
# Example:
# Rscript tree_plot_unrooted.R proteome_metadata.tsv phylogenetic_tree.nwk output_tree_unrooted.png 600

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
  stop("Usage: Rscript tree_plot_unrooted.R <metadata.tsv> <tree.nwk> <output.png> [dpi]")
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
# Remove any suffix after '|' in tip labels to match proteome_ids
tree$tip.label <- sapply(strsplit(tree$tip.label, "\\|"), `[`, 1)       

# NOTE: For unrooted trees, we DO NOT use midpoint() or compute.brlen().
# We want to keep the raw, unaltered evolutionary branch lengths.

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

# Fill in the "holes" with an "Other" category for missing classifications
metadata$Group[is.na(metadata$Group)] <- "Other / Unknown"

# Palette
group_colors <- c(
  "Eukaryote"           = "#1f77b4",
  "Alphaproteobacteria" = "#ff7f0e",
  "Bacteria"            = "#A6CEE3",
  "Archaea"             = "#d62728",
  "Other / Unknown"     = "#e0e0e0"
)

metadata$Group <- factor(metadata$Group, levels = names(group_colors))

# Diagnostic Check: See how many tips are actually missing from the TSV
missing_tips <- sum(!(tree$tip.label %in% metadata$label))
if(missing_tips > 0) {
  cat("    WARNING:", missing_tips, "tips on the tree were not found in the TSV file.\n")
}

# --- Plotting ---
# Change layout to "daylight" and set a static neutral color/thickness for all branches
p <- ggtree(tree, layout = "daylight", color = "black", size = 0.3) %<+% metadata +
  
  # Draw colored dots at the ends of the branches
  geom_tippoint(aes(color = Group), size = 2) + 
  
  # Apply the colors ONLY to the tips now
  scale_color_manual(values = group_colors, name = "Taxonomic Group", na.translate = FALSE) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(output_file, plot = p, width = 12, height = 12, dpi = dpi)
cat("Done! Unrooted tree for: ", tree_file, " saved to:", output_file, "\n")