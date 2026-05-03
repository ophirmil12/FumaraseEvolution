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
  library(plotly)
  library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(3, 4))) {
  stop("Usage: Rscript tree_plot_unrooted.R <metadata.tsv> <tree.nwk> <output.png> [dpi]")
}

metadata_file <- args[1]
tree_file     <- args[2]
output_file   <- args[3]
dpi <- ifelse(length(args) == 4, as.numeric(args[4]), 600)

if (is.na(dpi) || dpi < 72) {
  stop("DPI must be a valid number (e.g., >= 72)")
}

# --- Tree ---
tree <- read.tree(tree_file)
tree$tip.label <- sapply(strsplit(tree$tip.label, "\\|"), `[`, 1)

# --- Metadata ---
raw_metadata <- read.delim(metadata_file, sep = "\t", stringsAsFactors = FALSE)

metadata <- data.frame(
  label = raw_metadata$proteome_id,
  organism = raw_metadata$organism,
  Group = NA_character_,
  stringsAsFactors = FALSE
)

lineage <- raw_metadata$taxonomic_lineage
metadata$Group[grepl("Eukaryota",          lineage, ignore.case = TRUE)] <- "Eukaryote"
metadata$Group[grepl("Archaea",            lineage, ignore.case = TRUE)] <- "Archaea"
metadata$Group[grepl("Bacteria",           lineage, ignore.case = TRUE)] <- "Bacteria"
metadata$Group[grepl("Alphaproteobacteria",lineage, ignore.case = TRUE)] <- "Alphaproteobacteria"
metadata$Group[is.na(metadata$Group)] <- "Other / Unknown"

group_colors <- c(
  "Eukaryote"           = "#1f77b4",
  "Alphaproteobacteria" = "#ff7f0e",
  "Bacteria"            = "#A6CEE3",
  "Archaea"             = "#d62728",
  "Other / Unknown"     = "#e0e0e0"
)

metadata$Group <- factor(metadata$Group, levels = names(group_colors))

missing_tips <- sum(!(tree$tip.label %in% metadata$label))
if (missing_tips > 0) {
  cat("    WARNING:", missing_tips, "tips on the tree were not found in the TSV file.\n")
}

# --- Plotting ---
p <- ggtree(tree, layout = "daylight", color = "black", linewidth = 0.3) %<+% metadata +
  geom_tippoint(aes(color = Group), size = 2) +
  scale_color_manual(values = group_colors, name = "Taxonomic Group", na.translate = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    plot.margin  = margin(10, 10, 10, 10)
  )

# --- Static PNG ---
ggsave(output_file, plot = p, width = 12, height = 12, dpi = dpi)
cat("Done! Unrooted tree saved to:", output_file, "\n")

# --- Interactive HTML: manual coordinate extraction ---
# ggplotly() fails on ggtree "daylight" layout because branch segments are
# built internally and come through with mode="" (invisible). Instead, we
# pull the coordinates directly from the ggtree object and build Plotly traces.

html_file <- sub("\\.[^.]+$", "_interactive.html", output_file)

tryCatch({
  library(plotly)
  library(htmlwidgets)

  # Extract the layout data from the ggtree object
  td <- p$data  # ggtree stores node/tip coordinates here after %<+% join

  # --- Branch segments ---
  # ggtree stores parent x/y in xend/yend; build a line per edge
  edges <- td[!is.na(td$parent) & td$node != td$parent, ]

  # Build one trace per edge (or a single trace with NAs to separate lines)
  edge_x <- c()
  edge_y <- c()
  for (i in seq_len(nrow(edges))) {
    parent_row <- td[td$node == edges$parent[i], ]
    if (nrow(parent_row) == 0) next
    edge_x <- c(edge_x, edges$x[i], parent_row$x[1], NA)
    edge_y <- c(edge_y, edges$y[i], parent_row$y[1], NA)
  }

  branch_trace <- list(
    x = edge_x, y = edge_y,
    type = "scatter", mode = "lines",
    line = list(color = "black", width = 0.8),
    hoverinfo = "none",
    showlegend = FALSE,
    name = "branches"
  )

  # --- Tip points (colored by Group) ---
  tips <- td[td$isTip == TRUE & !is.na(td$Group), ]

  tip_traces <- lapply(levels(tips$Group), function(grp) {
    sub_tips <- tips[tips$Group == grp, ]
    list(
      x = sub_tips$x,
      y = sub_tips$y,
      type = "scatter", mode = "markers",
      marker = list(color = group_colors[grp], size = 6),
      text = paste0("ID: ", sub_tips$label, "<br>Organism: ", sub_tips$organism, "<br>Group: ", grp),
      hoverinfo = "text",
      name = grp
    )
  })

  all_traces <- c(list(branch_trace), tip_traces)

  fig <- plot_ly() 
  for (tr in all_traces) {
    fig <- add_trace(fig,
      x = tr$x, y = tr$y,
      type = tr$type, mode = tr$mode,
      marker = tr$marker %||% NULL,
      line = tr$line %||% NULL,
      text = tr$text %||% NULL,
      hoverinfo = tr$hoverinfo,
      name = tr$name,
      showlegend = tr$showlegend %||% TRUE
    )
  }

  fig <- fig %>% layout(
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE, scaleanchor = "x"),
    plot_bgcolor  = "white",
    paper_bgcolor = "white",
    legend = list(title = list(text = "<b>Taxonomic Group</b>"))
  )

  my_local_tmp <- "R_tmp_libs"
  
  saveWidget(fig, file = basename(html_file), selfcontained = TRUE, libdir = my_local_tmp)
  file.rename(basename(html_file), html_file)
  
  unlink(my_local_tmp, recursive = TRUE)
  
  cat("Interactive HTML saved to:", html_file, "\n")

}, error = function(e) {
  cat("WARNING: Could not save interactive HTML.\n  Reason:", conditionMessage(e), "\n")
})




