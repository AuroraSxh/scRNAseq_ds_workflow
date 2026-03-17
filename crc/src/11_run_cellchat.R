#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
lib_dir <- file.path(root_dir, ".r_libs")
if (dir.exists(lib_dir)) {
  .libPaths(c(lib_dir, .libPaths()))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(future)
  library(patchwork)
  library(ggalluvial)
})

suppressPackageStartupMessages(library(CellChat))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4 * 1024^3)

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
worker_count <- min(8L, as.integer(detected_cores))
future::plan("multisession", workers = worker_count)
on.exit(future::plan("sequential"), add = TRUE)

input_dir <- file.path(root_dir, "results", "communication", "cellchat", "input")
output_dir <- file.path(root_dir, "results", "communication", "cellchat")
checkpoint_path <- file.path(output_dir, "cellchat_checkpoint.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading CellChat input matrix")
expr <- readMM(file.path(input_dir, "expression.mtx"))
expr <- as(expr, "dgCMatrix")
genes <- read.delim(file.path(input_dir, "genes.tsv"), header = FALSE, stringsAsFactors = FALSE)[, 1]
meta <- read.delim(file.path(input_dir, "cells.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(length(genes) == nrow(expr))
stopifnot(nrow(meta) == ncol(expr))
rownames(expr) <- genes
colnames(expr) <- meta$cell_id
rownames(meta) <- meta$cell_id

group_levels <- unique(meta$cellchat_group)
meta$cellchat_group <- factor(meta$cellchat_group, levels = group_levels)

sanitize_name <- function(x) {
  cleaned <- gsub("[^A-Za-z0-9]+", "_", x)
  cleaned <- gsub("^_+|_+$", "", cleaned)
  tolower(cleaned)
}

pretty_group_name <- function(x) {
  mapping <- c(
    "C1QC_TAM" = "C1QC TAM",
    "FOLR2_TAM" = "FOLR2 TAM",
    "SPP1_TAM" = "SPP1 TAM",
    "Inflammatory_Macro" = "Inflam Macro",
    "Tpex_like" = "Tpex-like",
    "Early_Tex" = "Early Tex",
    "Effector_memory" = "Effector Memory",
    "Terminal_Tex" = "Terminal Tex",
    "Stress_program" = "Stress Program"
  )
  unname(ifelse(x %in% names(mapping), mapping[x], gsub("_", " ", x)))
}

wrap_labels <- function(values, width = 10L) {
  vapply(
    values,
    function(x) paste(strwrap(as.character(x), width = width), collapse = "\n"),
    character(1)
  )
}

identify_pattern_with_fallback <- function(cellchat, pattern_name, k_candidates) {
  for (pattern_k in k_candidates) {
    attempt <- tryCatch(
      {
        tmp <- identifyCommunicationPatterns(cellchat, pattern = pattern_name, k = pattern_k)
        list(ok = TRUE, value = tmp, k = pattern_k)
      },
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )
    if (isTRUE(attempt$ok)) {
      return(list(cellchat = attempt$value, k = attempt$k))
    }
  }
  stop(sprintf("CellChat pattern identification failed for pattern '%s' across tested k values.", pattern_name))
}

save_pattern_overview <- function(cellchat, pattern_name, k_candidates, output_path, title_text) {
  pattern_result <- identify_pattern_with_fallback(cellchat, pattern_name, k_candidates)
  cellchat <- pattern_result$cellchat
  chosen_k <- pattern_result$k
  message(sprintf("Rendering %s communication pattern plot with k=%d", pattern_name, chosen_k))

  river_plot <- netAnalysis_river(cellchat, pattern = pattern_name)
  dot_plot <- netAnalysis_dot(cellchat, pattern = pattern_name)

  if (length(dot_plot$scales$scales) > 0) {
    dot_plot$scales$scales <- Filter(
      function(scale_obj) !("x" %in% scale_obj$aesthetics),
      dot_plot$scales$scales
    )
  }

  river_plot <- river_plot +
    ggplot2::theme(
      plot.margin = ggplot2::margin(12, 14, 8, 12),
      axis.text = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_blank(),
      legend.position = "none"
    )

  dot_plot <- dot_plot +
    ggplot2::scale_x_discrete(labels = function(x) wrap_labels(x, width = 9L)) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(8, 18, 22, 12),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  combined_plot <- river_plot / dot_plot +
    plot_layout(heights = c(1.0, 1.15)) +
    plot_annotation(
      title = title_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0,
          size = 13,
          face = "bold",
          margin = ggplot2::margin(0, 0, 8, 0)
        )
      )
    )
  ggplot2::ggsave(
    filename = output_path,
    plot = combined_plot,
    width = 12.6,
    height = 12.2,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )

  list(cellchat = cellchat, k = chosen_k)
}

extract_pathway_prob_sum <- function(cellchat, signaling, thresh = 0.05) {
  pairLR <- searchPair(
    signaling = signaling,
    pairLR.use = cellchat@LR$LRsig,
    key = "pathway_name",
    matching.exact = TRUE,
    pair.only = TRUE
  )
  net <- cellchat@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0

  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name, drop = FALSE], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name, drop = FALSE]) != 0]
  }
  if (!length(pairLR.name.use)) {
    stop(sprintf("There is no significant communication for pathway '%s'.", signaling))
  }

  prob <- prob[, , pairLR.name.use, drop = FALSE]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
  }
  apply(prob, c(1, 2), sum)
}

netVisual_hierarchy2_fixed <- function(
  net,
  vertex.receiver,
  color.use = NULL,
  title.name = NULL,
  sources.use = NULL,
  targets.use = NULL,
  remove.isolate = FALSE,
  top = 1,
  weight.scale = FALSE,
  vertex.weight = 20,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  edge.weight.max = NULL,
  edge.width.max = 8,
  alpha.edge = 0.6,
  label.dist = 2.8,
  space.v = 1.5,
  space.h = 1.6,
  shape = NULL,
  label.edge = FALSE,
  edge.curved = 0,
  margin = 0.2,
  vertex.label.cex = 0.6,
  vertex.label.color = "black",
  arrow.width = 1,
  arrow.size = 0.2,
  edge.label.color = "black",
  edge.label.cex = 0.5,
  vertex.size = NULL
) {
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    vertex.size.max <- if (length(unique(vertex.weight)) == 1) 5 else 15
  }
  options(warn = -1)
  thresh_value <- stats::quantile(net, probs = 1 - top)
  net[net < thresh_value] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    cells.level <- rownames(net)
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }

  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, , drop = FALSE]
    net <- net[, -idx, drop = FALSE]
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight / vertex.weight.max * vertex.size.max + 6

  m <- length(vertex.receiver)
  m0 <- nrow(net) - m
  net2 <- net
  reorder.row <- c(setdiff(seq_len(nrow(net)), vertex.receiver), vertex.receiver)
  net2 <- net2[reorder.row, vertex.receiver, drop = FALSE]
  m1 <- nrow(net2)
  n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[setdiff(seq_len(m1), vertex.receiver)], row.names(net)[vertex.receiver], rep("", m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(seq_len(m1), vertex.receiver)], color.use[vertex.receiver], rep("#FFFFFF", m))
  color.use3.frame <- c(color.use[setdiff(seq_len(m1), vertex.receiver)], color.use[vertex.receiver], color.use[vertex.receiver])

  if (length(vertex.weight) != 1) {
    vertex.weight <- c(vertex.weight[setdiff(seq_len(m1), vertex.receiver)], vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle", nrow(net3))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = TRUE)
  edge.start <- ends(g, es = igraph::E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0, 1] <- 0
  coords[(m0 + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h / 2
  coords[1:m0, 2] <- seq(space.v, 0, by = -space.v / (m0 - 1))
  coords[(m0 + 1):m1, 2] <- seq(space.v, 0, by = -space.v / (m1 - m0 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v / (n1 - 1))

  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- round(igraph::E(g)$weight, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight / edge.weight.max * edge.width.max
  } else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[, 1]], alpha.edge)

  label.distances <- c(rep(space.h * label.dist, m), rep(space.h * label.dist, m1 - m), rep(0, nrow(net3) - m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1 - m0), rep(-pi, nrow(net3) - m1))
  text.pos <- cbind(c(-space.h / 1.5, space.h / 22, space.h / 1.5), space.v - space.v / 7)
  cellchat_mycircle <- get("mycircle", envir = asNamespace("CellChat"))
  igraph::add.vertex.shape(
    "fcircle",
    clip = igraph::igraph.shape.noclip,
    plot = cellchat_mycircle,
    parameters = list(vertex.frame.color = 1, vertex.frame.width = 1)
  )
  plot(
    g,
    edge.curved = edge.curved,
    layout = coords,
    margin = margin,
    rescale = TRUE,
    vertex.shape = "fcircle",
    vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
    vertex.label.degree = label.locs,
    vertex.label.dist = label.distances,
    vertex.label.family = "Helvetica"
  )
  text(text.pos, c("Source", "Target", "Source"), cex = 0.8, col = c("#c51b7d", "#2f6661", "#2f6661"))
  shape::Arrows(-space.h / 1.5, space.v - space.v / 4, space.h / 1e+05, space.v - space.v / 4, col = "#c51b7d", arr.lwd = 1e-04, arr.length = 0.2, lwd = 0.8, arr.type = "triangle")
  shape::Arrows(space.h / 1.5, space.v - space.v / 4, space.h / 20, space.v - space.v / 4, col = "#2f6661", arr.lwd = 1e-04, arr.length = 0.2, lwd = 0.8, arr.type = "triangle")
  if (!is.null(title.name)) {
    text(space.h / 8, space.v, paste0(title.name, " signaling network"), cex = 1)
  }
  recordPlot()
}

rank_directional_pathways <- function(interactions, source_groups, target_groups, top_n = 3L) {
  selected <- interactions[
    interactions$source %in% source_groups &
      interactions$target %in% target_groups &
      !is.na(interactions$pathway_name) &
      nzchar(interactions$pathway_name),
    ,
    drop = FALSE
  ]
  if (!nrow(selected)) {
    return(data.frame())
  }

  ranked <- aggregate(prob ~ pathway_name, data = selected, FUN = sum)
  ranked <- ranked[order(-ranked$prob, ranked$pathway_name), , drop = FALSE]
  rownames(ranked) <- NULL
  head(ranked, top_n)
}

save_pathway_network_plot <- function(
  cellchat,
  signaling,
  output_path,
  title_text,
  layout_name,
  source_groups,
  target_groups
) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = output_path, width = 3600, height = 1800, res = 300, bg = "white")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mar = c(0.5, 0.5, 3.5, 0.5), xpd = TRUE)
  if (identical(layout_name, "hierarchy")) {
    prob.sum <- extract_pathway_prob_sum(cellchat, signaling)
    original_groups <- rownames(prob.sum)
    pretty_names <- pretty_group_name(original_groups)
    rownames(prob.sum) <- pretty_names
    colnames(prob.sum) <- pretty_names
    source_groups_pretty <- pretty_group_name(source_groups)
    target_groups_pretty <- pretty_group_name(target_groups)
    receiver_idx <- which(pretty_names %in% target_groups_pretty)
    graphics::par(mfrow = c(1, 2), oma = c(0, 0, 2.5, 0), ps = 18)
    netVisual_hierarchy1(
      net = prob.sum,
      vertex.receiver = receiver_idx,
      sources.use = source_groups_pretty,
      targets.use = target_groups_pretty,
      remove.isolate = FALSE,
      edge.width.max = 7,
      vertex.size.max = 18,
      vertex.label.cex = 0.75,
      title.name = NULL
    )
    netVisual_hierarchy2_fixed(
      net = prob.sum,
      vertex.receiver = setdiff(seq_len(nrow(prob.sum)), receiver_idx),
      sources.use = source_groups_pretty,
      targets.use = target_groups_pretty,
      remove.isolate = FALSE,
      edge.width.max = 7,
      vertex.size.max = 18,
      vertex.label.cex = 0.75,
      title.name = NULL
    )
    graphics::mtext(title_text, side = 3, outer = TRUE, cex = 1.05, line = 0.2, font = 2)
  } else {
    netVisual_aggregate(
      object = cellchat,
      signaling = signaling,
      layout = "circle",
      sources.use = source_groups,
      targets.use = target_groups,
      remove.isolate = TRUE,
      edge.width.max = 10,
      vertex.size.max = 18,
      vertex.label.cex = 1.15,
      pt.title = 18,
      title.space = 2
    )
  }
}

message(sprintf("Running CellChat on %d cells, %d genes, %d workers", ncol(expr), nrow(expr), worker_count))
if (file.exists(checkpoint_path)) {
  message("Loading existing CellChat checkpoint")
  cellchat <- readRDS(checkpoint_path)
} else {
  cellchat <- createCellChat(object = expr, meta = meta, group.by = "cellchat_group")
  data(CellChatDB.human)
  cellchat@DB <- CellChatDB.human

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 20)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  saveRDS(cellchat, checkpoint_path)
}

pattern_k_candidates <- seq.int(min(5L, length(group_levels) - 1L), 2L, by = -1L)
incoming_result <- save_pattern_overview(
  cellchat = cellchat,
  pattern_name = "incoming",
  k_candidates = pattern_k_candidates,
  output_path = file.path(output_dir, "cellchat.png"),
  title_text = "Incoming communication overview"
)
cellchat <- incoming_result$cellchat

outgoing_result <- save_pattern_overview(
  cellchat = cellchat,
  pattern_name = "outgoing",
  k_candidates = pattern_k_candidates,
  output_path = file.path(output_dir, "cellchat_outgoing.png"),
  title_text = "Outgoing communication overview"
)
cellchat <- outgoing_result$cellchat

interactions <- subsetCommunication(cellchat)
group_lineage <- unique(meta[, c("cellchat_group", "lineage")])
lineage_map <- setNames(group_lineage$lineage, group_lineage$cellchat_group)
interactions$source_lineage <- unname(lineage_map[interactions$source])
interactions$target_lineage <- unname(lineage_map[interactions$target])
cross_lineage <- interactions[interactions$source_lineage != interactions$target_lineage, , drop = FALSE]
write.csv(interactions, file.path(output_dir, "cellchat_all_interactions.csv"), row.names = FALSE)
write.csv(cross_lineage, file.path(output_dir, "cellchat_myeloid_cd8_interactions.csv"), row.names = FALSE)

interaction_summary <- aggregate(
  prob ~ source + target + source_lineage + target_lineage,
  data = cross_lineage,
  FUN = sum
)
colnames(interaction_summary)[colnames(interaction_summary) == "prob"] <- "summed_prob"
interaction_summary <- interaction_summary[order(-interaction_summary$summed_prob), , drop = FALSE]
write.csv(interaction_summary, file.path(output_dir, "cellchat_myeloid_cd8_summary.csv"), row.names = FALSE)

macrophage_groups <- group_lineage$cellchat_group[group_lineage$lineage == "Myeloid"]
cd8_groups <- group_lineage$cellchat_group[group_lineage$lineage == "CD8"]
mac_t_dir <- file.path(root_dir, "results", "communication", "mac_t")
mac_t_network_dir <- file.path(mac_t_dir, "pathway_networks")
dir.create(mac_t_network_dir, recursive = TRUE, showWarnings = FALSE)

mac_to_cd8_pathways <- rank_directional_pathways(interactions, macrophage_groups, cd8_groups, top_n = 3L)
cd8_to_mac_pathways <- rank_directional_pathways(interactions, cd8_groups, macrophage_groups, top_n = 3L)
write.csv(mac_to_cd8_pathways, file.path(mac_t_network_dir, "macrophage_to_cd8_top_network_pathways.csv"), row.names = FALSE)
write.csv(cd8_to_mac_pathways, file.path(mac_t_network_dir, "cd8_to_macrophage_top_network_pathways.csv"), row.names = FALSE)

plot_specs <- list(
  list(
    prefix = "macrophage_to_cd8",
    title_prefix = "Macrophage -> CD8",
    source_groups = macrophage_groups,
    target_groups = cd8_groups,
    ranked = mac_to_cd8_pathways
  ),
  list(
    prefix = "cd8_to_macrophage",
    title_prefix = "CD8 -> Macrophage",
    source_groups = cd8_groups,
    target_groups = macrophage_groups,
    ranked = cd8_to_mac_pathways
  )
)

for (spec in plot_specs) {
  if (!nrow(spec$ranked)) {
    next
  }
  for (i in seq_len(nrow(spec$ranked))) {
    pathway_name <- spec$ranked$pathway_name[[i]]
    safe_name <- sanitize_name(pathway_name)
    hierarchy_path <- file.path(mac_t_network_dir, sprintf("%s_%s_hierarchy.png", spec$prefix, safe_name))
    circle_path <- file.path(mac_t_network_dir, sprintf("%s_%s_circle.png", spec$prefix, safe_name))
    save_pathway_network_plot(
      cellchat = cellchat,
      signaling = pathway_name,
      output_path = hierarchy_path,
      title_text = sprintf("%s %s hierarchy network", spec$title_prefix, pathway_name),
      layout_name = "hierarchy",
      source_groups = spec$source_groups,
      target_groups = spec$target_groups
    )
    save_pathway_network_plot(
      cellchat = cellchat,
      signaling = pathway_name,
      output_path = circle_path,
      title_text = sprintf("%s %s circle network", spec$title_prefix, pathway_name),
      layout_name = "circle",
      source_groups = spec$source_groups,
      target_groups = spec$target_groups
    )
  }
}

file.copy(file.path(output_dir, "cellchat.png"), file.path(mac_t_dir, "cellchat_overview.png"), overwrite = TRUE)
file.copy(file.path(output_dir, "cellchat_outgoing.png"), file.path(mac_t_dir, "cellchat_outgoing_overview.png"), overwrite = TRUE)
saveRDS(cellchat, file.path(output_dir, "cellchat.rds"))

message("CellChat run finished")
