#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)     # SaveSeuratRds / LoadSeuratRds
  library(sctransform)
  library(magrittr)
  library(ggplot2)
  library(MASS)            # kde2d
  library(harmony)         # RunHarmony
  library(ComplexHeatmap)
  library(circlize)
  library(tidyr)
})

# ============================================================
# Utilities
# ============================================================

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n")
`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

# IMPORTANT: Seurat subset() does not like get("colname") inside subset=
# This helper subsets by meta.data column using explicit cell names.
subset_by_meta <- function(obj, col, values, invert = FALSE) {
  md <- obj@meta.data
  if (!col %in% colnames(md)) {
    stop("Metadata column not found: ", col, "\nAvailable columns: ",
         paste(colnames(md), collapse = ", "))
  }
  keep <- md[[col]] %in% values
  if (invert) keep <- !keep
  cells_keep <- rownames(md)[which(keep)]
  subset(obj, cells = cells_keep)
}

filter_min_cells <- function(obj, label_col, min_cells = 10) {
  md <- obj@meta.data
  if (!label_col %in% colnames(md)) {
    stop("Label column not found in meta.data: ", label_col, "\nAvailable columns: ",
         paste(colnames(md), collapse = ", "))
  }
  ct <- table(md[[label_col]])
  keep_types <- names(ct[ct > min_cells])
  subset_by_meta(obj, label_col, keep_types, invert = FALSE)
}

get_palette <- function(n) {
  if (requireNamespace("scCustomize", quietly = TRUE)) {
    cols <- scCustomize::DiscretePalette_scCustomize(num_colors = max(36, n), palette = "polychrome")
  } else {
    cols <- grDevices::hcl.colors(max(36, n), "Dynamic")
  }
  cols[seq_len(n)]
}

read_seurat_rds <- function(path, assay = "RNA") {
  if (!file.exists(path)) stop("File not found: ", path)
  obj <- LoadSeuratRds(path)
  if (assay %in% Assays(obj)) DefaultAssay(obj) <- assay
  obj
}

# ============================================================
# Core STEVE: integration + KDE Bayes model + prediction
# ============================================================

steve_integrate <- function(ref_obj,
                            query_obj,
                            integration = "harmony",
                            dims = 1:20,
                            cluster_resolution = 0.8,
                            seed = 123,
                            reduction_name = "integrated.result",
                            umap_name = "umap.result",
                            harmony_group_col = "orig.ident") {
  
  set.seed(seed)
  
  # define the ONLY batch variable we want Harmony to correct: reference vs query
  ref_obj$orig.ident <- "reference"
  query_obj$orig.ident <- "query"
  
  combined <- merge(query_obj, ref_obj)
  assay_use <- DefaultAssay(combined)
  
  # If this is a v5 assay with multiple counts.* layers, join them to ONE counts matrix.
  # This prevents accidental "layer = one subset wins" behavior.
  if (assay_use %in% Assays(combined)) {
    lay <- tryCatch(Layers(combined[[assay_use]]), error = function(e) character(0))
    if (length(lay) > 0) {
      # join any counts.* layers into a single "counts"
      combined <- tryCatch(
        JoinLayers(combined, assay = assay_use, layers = "counts"),
        error = function(e) combined
      )
    }
  }
  
  # Ensure the harmony grouping column exists and is clean
  if (!harmony_group_col %in% colnames(combined@meta.data)) {
    stop("Harmony group column not found: ", harmony_group_col)
  }
  combined@meta.data[[harmony_group_col]] <- factor(combined@meta.data[[harmony_group_col]])
  
  # Standard preprocessing
  combined <- NormalizeData(combined, verbose = FALSE)
  combined <- FindVariableFeatures(combined, verbose = FALSE)
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, verbose = FALSE)
  
  integration_lower <- tolower(integration)
  
  # ---- Integration choice ----
  # We DO NOT call IntegrateLayers() here, because it requires multiple normalized layers (data.*)
  # and easily breaks depending on SeuratObject/Seurat versions.
  # RunHarmony is stable and uses the grouping column explicitly.
  if (integration_lower %in% c("harmony", "harmonyintegration", "runharmony")) {
    combined <- RunHarmony(
      object = combined,
      group.by.vars = harmony_group_col,
      reduction = "pca",
      dims.use = dims,
      assay.use = assay_use,
      reduction.save = reduction_name,
      verbose = FALSE
    )
    reduction_for_umap <- reduction_name
    
  } else if (integration_lower %in% c("none", "no", "pca")) {
    reduction_for_umap <- "pca"
    
  } else {
    # If user passes something else, we just fall back to no integration (PCA)
    msg("Unknown integration='%s' -> using PCA (no integration).", integration)
    reduction_for_umap <- "pca"
  }
  
  combined <- RunUMAP(
    combined,
    reduction = reduction_for_umap,
    dims = dims,
    reduction.name = umap_name,
    verbose = FALSE
  )
  combined <- FindNeighbors(combined, reduction = reduction_for_umap, dims = dims, verbose = FALSE)
  combined <- FindClusters(combined, resolution = cluster_resolution, verbose = FALSE)
  
  combined
}

steve_build_kde_model <- function(combined,
                                  ref_ident = "reference",
                                  ref_label_col,
                                  umap_name = "umap.result",
                                  grid_n = 200,
                                  laplace = 1) {
  
  if (!ref_label_col %in% colnames(combined@meta.data)) {
    stop("Missing ref_label_col in meta.data: ", ref_label_col, "\nAvailable columns: ",
         paste(colnames(combined@meta.data), collapse = ", "))
  }
  if (!umap_name %in% names(combined@reductions)) stop("Missing reduction: ", umap_name)
  
  ref <- subset(combined, subset = orig.ident == ref_ident)
  umap_all <- combined@reductions[[umap_name]]@cell.embeddings
  
  celltypes <- levels(droplevels(as.factor(ref@meta.data[[ref_label_col]])))
  if (length(celltypes) < 2) stop("Need at least 2 reference cell types after filtering.")
  
  lims <- c(range(umap_all[, 1]), range(umap_all[, 2]))
  
  kde_list <- vector("list", length(celltypes))
  names(kde_list) <- celltypes
  
  for (ct in celltypes) {
    ref_ct <- subset_by_meta(ref, ref_label_col, ct, invert = FALSE)
    coords <- ref_ct@reductions[[umap_name]]@cell.embeddings
    
    if (nrow(coords) == 1) {
      coords <- rbind(coords, coords)
      coords[2, ] <- coords[2, ] + 0.1
    }
    
    kd <- kde2d(coords[, 1], coords[, 2], n = grid_n, lims = lims)
    lik <- kd$z / sum(kd$z)
    kde_list[[ct]] <- list(x = kd$x, y = kd$y, likelihood = lik)
  }
  
  ref_labels <- ref@meta.data[[ref_label_col]]
  counts <- table(ref_labels)
  priors <- (as.numeric(counts[celltypes]) + laplace) / (length(ref_labels) + laplace * length(celltypes))
  names(priors) <- celltypes
  
  list(
    celltypes = celltypes,
    priors = priors,
    kde = kde_list,
    umap_name = umap_name,
    grid_n = grid_n,
    lims = lims
  )
}

steve_predict <- function(combined,
                          model,
                          odds_cutoff = 2,
                          unknown_label = "unknown") {
  
  umap <- combined@reductions[[model$umap_name]]@cell.embeddings
  cell_ids <- rownames(umap)
  celltypes <- model$celltypes
  priors <- model$priors
  kde_list <- model$kde
  
  nearest_idx <- function(grid, value) which.min(abs(grid - value))
  
  pred_label <- character(nrow(umap))
  odds <- numeric(nrow(umap))
  top_p <- numeric(nrow(umap))
  second_p <- numeric(nrow(umap))
  
  for (i in seq_len(nrow(umap))) {
    x <- umap[i, 1]
    y <- umap[i, 2]
    
    lik <- numeric(length(celltypes))
    for (k in seq_along(celltypes)) {
      ct <- celltypes[k]
      kd <- kde_list[[ct]]
      xi <- nearest_idx(kd$x, x)
      yi <- nearest_idx(kd$y, y)
      lik[k] <- kd$likelihood[xi, yi]
    }
    
    post_unnorm <- priors * lik
    denom <- sum(post_unnorm)
    
    if (denom <= 0 || is.na(denom)) {
      pred_label[i] <- unknown_label
      odds[i] <- NA_real_
      top_p[i] <- NA_real_
      second_p[i] <- NA_real_
      next
    }
    
    post <- post_unnorm / denom
    ord <- order(post, decreasing = TRUE)
    best <- ord[1]
    second <- ord[2]
    
    top_p[i] <- post[best]
    second_p[i] <- post[second]
    odds[i] <- ifelse(second_p[i] > 0, top_p[i] / second_p[i], Inf)
    
    if (!is.finite(odds[i]) || is.na(odds[i])) {
      pred_label[i] <- unknown_label
    } else if (odds[i] >= odds_cutoff) {
      pred_label[i] <- celltypes[best]
    } else {
      pred_label[i] <- unknown_label
    }
  }
  
  data.frame(
    cellid = cell_ids,
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    pred = pred_label,
    odds = odds,
    top_p = top_p,
    second_p = second_p,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Metrics + plotting helpers
# ============================================================

one_vs_rest_metrics <- function(true_labels, pred_labels, positive_label) {
  true_pos <- (true_labels == positive_label)
  pred_pos <- (pred_labels == positive_label)
  
  TP <- sum(true_pos & pred_pos, na.rm = TRUE)
  FN <- sum(true_pos & !pred_pos, na.rm = TRUE)
  TN <- sum(!true_pos & !pred_pos, na.rm = TRUE)
  FP <- sum(!true_pos & pred_pos, na.rm = TRUE)
  
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  
  c(TP = TP, FN = FN, TN = TN, FP = FP, Sensitivity = sens, Specificity = spec)
}

per_celltype_sens_spec <- function(true_labels, pred_labels) {
  cts <- levels(droplevels(as.factor(true_labels)))
  out <- data.frame(celltype = cts, Sensitivity = NA_real_, Specificity = NA_real_)
  for (i in seq_along(cts)) {
    m <- one_vs_rest_metrics(true_labels, pred_labels, cts[i])
    out$Sensitivity[i] <- m["Sensitivity"]
    out$Specificity[i] <- m["Specificity"]
  }
  out
}

confusion_ratio <- function(true_labels, pred_labels) {
  true_levels <- levels(droplevels(as.factor(true_labels)))
  pred_levels <- levels(droplevels(as.factor(pred_labels)))
  
  mat <- matrix(0, nrow = length(pred_levels), ncol = length(true_levels))
  rownames(mat) <- pred_levels
  colnames(mat) <- true_levels
  
  for (t in true_levels) {
    idx <- which(true_labels == t)
    if (length(idx) == 0) next
    denom <- length(idx)
    for (p in pred_levels) {
      mat[p, t] <- sum(pred_labels[idx] == p, na.rm = TRUE) / denom
    }
  }
  mat
}

plot_integrated_umap <- function(combined, label_col, out_pdf, title = NULL) {
  n <- length(unique(combined@meta.data[[label_col]]))
  cols <- c("gray", get_palette(max(1, n)))
  p <- DimPlot(combined, reduction = "umap.result", group.by = label_col, cols = cols, label = FALSE) +
    ggtitle(title %||% "")
  ggsave(out_pdf, p, width = 6, height = 4)
}

plot_heatmap <- function(mat, out_pdf, out_csv) {
  write.csv(mat, out_csv)
  col_fun <- colorRamp2(breaks = c(0, 1), colors = c("gray", "red"))
  
  ht <- Heatmap(
    mat,
    name = "Ratio",
    row_title = "Predicted",
    column_title = "Truth",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    column_names_gp = grid::gpar(fontsize = 6),
    row_names_gp = grid::gpar(fontsize = 6)
  )
  
  pdf(out_pdf, height = 0.12 * ncol(mat) + 2, width = 0.35 * nrow(mat) + 2)
  draw(ht)
  dev.off()
}

# ============================================================
# Module 1: Subsampling Evaluation
# ============================================================

steve_module_subsampling <- function(obj,
                                     truth_col,
                                     ratios = seq(0.1, 0.9, by = 0.1),
                                     n_repeats = 10,
                                     outputdir = "./output_subsampling",
                                     integration = "harmony",
                                     dims = 1:20,
                                     cluster_resolution = 0.8,
                                     grid_n = 200,
                                     odds_cutoff = 2,
                                     min_cells = 10,
                                     seed = 123) {
  
  ensure_dir(outputdir)
  obj <- filter_min_cells(obj, truth_col, min_cells = min_cells)
  
  all_summary <- list()
  
  for (r in ratios) {
    for (iter in seq_len(n_repeats)) {
      tag <- sprintf("ratio_%02d_iter_%02d", round(r * 100), iter)
      out <- file.path(outputdir, tag)
      ensure_dir(out)
      set.seed(seed + iter + round(r * 1000))
      
      n <- ncol(obj)
      ref_n <- max(1, floor(n * r))
      ref_idx <- sample(seq_len(n), ref_n, replace = FALSE)
      query_idx <- setdiff(seq_len(n), ref_idx)
      
      ref_obj <- obj[, ref_idx]
      query_obj <- obj[, query_idx]
      
      combined <- steve_integrate(
        ref_obj, query_obj,
        integration = integration,
        dims = dims,
        cluster_resolution = cluster_resolution,
        seed = seed + iter
      )
      
      combined$STEVE_plot <- "query"
      ref_cells <- colnames(subset(combined, subset = orig.ident == "reference"))
      combined@meta.data[ref_cells, "STEVE_plot"] <- combined@meta.data[ref_cells, truth_col]
      combined$STEVE_plot <- factor(combined$STEVE_plot)
      
      plot_integrated_umap(
        combined, "STEVE_plot",
        file.path(out, "integrated_umap.pdf"),
        title = sprintf("Subsampling: ref=%.0f%%", r * 100)
      )
      
      model <- steve_build_kde_model(combined, ref_ident = "reference", ref_label_col = truth_col, grid_n = grid_n)
      pred <- steve_predict(combined, model, odds_cutoff = odds_cutoff, unknown_label = "unknown")
      
      query_cells <- colnames(subset(combined, subset = orig.ident == "query"))
      pred_q <- pred[pred$cellid %in% query_cells, , drop = FALSE]
      
      truth_q <- combined@meta.data[pred_q$cellid, truth_col]
      metrics <- per_celltype_sens_spec(truth_q, pred_q$pred)
      write.csv(metrics, file.path(out, "per_celltype_metrics.csv"), row.names = FALSE)
      
      macro_sens <- mean(metrics$Sensitivity, na.rm = TRUE)
      macro_spec <- mean(metrics$Specificity, na.rm = TRUE)
      
      mat <- confusion_ratio(truth_q, pred_q$pred)
      plot_heatmap(
        mat,
        out_pdf = file.path(out, "confusion_heatmap.pdf"),
        out_csv = file.path(out, "confusion_heatmap.csv")
      )
      
      pred_q$truth <- truth_q
      write.csv(pred_q, file.path(out, "predictions_query.csv"), row.names = FALSE)
      
      all_summary[[length(all_summary) + 1]] <- data.frame(
        ratio = r,
        iteration = iter,
        macro_sensitivity = macro_sens,
        macro_specificity = macro_spec,
        stringsAsFactors = FALSE
      )
      
      msg("Subsampling done: %s (macro sens=%.3f spec=%.3f)", tag, macro_sens, macro_spec)
    }
  }
  
  summary_df <- do.call(rbind, all_summary)
  write.csv(summary_df, file.path(outputdir, "summary_macro_metrics.csv"), row.names = FALSE)
  
  p <- ggplot(summary_df, aes(x = interaction(ratio, iteration, drop = TRUE))) +
    geom_point(aes(y = macro_sensitivity, color = "Sensitivity")) +
    geom_point(aes(y = macro_specificity, color = "Specificity")) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "run", y = "macro metric", title = "Subsampling: Macro Sensitivity/Specificity") +
    scale_y_continuous(limits = c(0, 1))
  ggsave(file.path(outputdir, "macro_metrics_plot.pdf"), p, width = 8, height = 3)
  
  summary_df
}

# ============================================================
# Module 2: Novel Cell Evaluation
# ============================================================

steve_module_novel_cell <- function(obj,
                                    truth_col,
                                    outputdir = "./output_novel",
                                    n_repeats = 5,
                                    integration = "harmony",
                                    dims = 1:20,
                                    cluster_resolution = 0.8,
                                    grid_n = 200,
                                    odds_cutoff = 2,
                                    min_cells = 10,
                                    seed = 123) {
  
  ensure_dir(outputdir)
  obj <- filter_min_cells(obj, truth_col, min_cells = min_cells)
  
  celltypes <- levels(droplevels(as.factor(obj@meta.data[[truth_col]])))
  all_runs <- list()
  
  for (iter in seq_len(n_repeats)) {
    set.seed(seed + iter)
    n <- ncol(obj)
    half <- floor(n / 2)
    ref_idx <- sample(seq_len(n), half, replace = FALSE)
    query_idx <- setdiff(seq_len(n), ref_idx)
    base_ref <- obj[, ref_idx]
    base_query <- obj[, query_idx]
    
    for (omit in celltypes) {
      out <- file.path(outputdir, sprintf("iter_%02d_omit_%s", iter, gsub("[^A-Za-z0-9_\\-]", "_", omit)))
      ensure_dir(out)
      
      ref_obj <- subset_by_meta(base_ref, truth_col, omit, invert = TRUE)
      query_obj <- base_query
      
      combined <- steve_integrate(
        ref_obj, query_obj,
        integration = integration,
        dims = dims,
        cluster_resolution = cluster_resolution,
        seed = seed + iter
      )
      
      model <- steve_build_kde_model(combined, ref_ident = "reference", ref_label_col = truth_col, grid_n = grid_n)
      pred <- steve_predict(combined, model, odds_cutoff = odds_cutoff, unknown_label = "unknown")
      
      query_cells <- colnames(subset(combined, subset = orig.ident == "query"))
      pred_q <- pred[pred$cellid %in% query_cells, , drop = FALSE]
      truth_q <- combined@meta.data[pred_q$cellid, truth_col]
      
      is_novel_true <- (truth_q == omit)
      is_novel_pred <- (pred_q$pred == "unknown")
      
      TP <- sum(is_novel_true & is_novel_pred, na.rm = TRUE)
      FN <- sum(is_novel_true & !is_novel_pred, na.rm = TRUE)
      TN <- sum(!is_novel_true & !is_novel_pred, na.rm = TRUE)
      FP <- sum(!is_novel_true & is_novel_pred, na.rm = TRUE)
      
      novel_sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
      novel_spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
      
      mis <- sort(table(pred_q$pred[is_novel_true]), decreasing = TRUE)
      mis_df <- data.frame(pred = names(mis), n = as.integer(mis), stringsAsFactors = FALSE)
      write.csv(mis_df, file.path(out, "novel_misclassification_counts.csv"), row.names = FALSE)
      
      pred_q$truth <- truth_q
      write.csv(pred_q, file.path(out, "predictions_query.csv"), row.names = FALSE)
      
      all_runs[[length(all_runs) + 1]] <- data.frame(
        iteration = iter,
        omitted_type = omit,
        novel_sensitivity = novel_sens,
        novel_specificity = novel_spec,
        stringsAsFactors = FALSE
      )
      
      msg("Novel eval iter=%d omit=%s (sens=%.3f spec=%.3f)", iter, omit, novel_sens, novel_spec)
    }
  }
  
  out_df <- do.call(rbind, all_runs)
  write.csv(out_df, file.path(outputdir, "novel_summary.csv"), row.names = FALSE)
  
  p <- ggplot(out_df, aes(x = omitted_type, y = novel_sensitivity)) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Novel Cell Evaluation: sensitivity (detect omitted type as unknown)",
         x = "omitted type", y = "sensitivity")
  ggsave(file.path(outputdir, "novel_sensitivity_plot.pdf"), p, width = 8, height = 3.5)
  
  out_df
}

# ============================================================
# Module 3: Annotation Benchmarking
# ============================================================

steve_module_benchmark <- function(obj,
                                   truth_col,
                                   candidate_col,
                                   outputdir = "./output_benchmark",
                                   n_repeats = 5,
                                   integration = "harmony",
                                   dims = 1:20,
                                   cluster_resolution = 0.8,
                                   grid_n = 200,
                                   odds_cutoff = 2,
                                   min_cells = 10,
                                   seed = 123) {
  
  ensure_dir(outputdir)
  if (!candidate_col %in% colnames(obj@meta.data)) {
    stop("Candidate annotation column not found in meta.data: ", candidate_col, "\nAvailable columns: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  }
  
  obj <- filter_min_cells(obj, truth_col, min_cells = min_cells)
  all_runs <- list()
  
  for (iter in seq_len(n_repeats)) {
    set.seed(seed + iter)
    n <- ncol(obj)
    half <- floor(n / 2)
    ref_idx <- sample(seq_len(n), half, replace = FALSE)
    query_idx <- setdiff(seq_len(n), ref_idx)
    
    ref_obj <- obj[, ref_idx]
    query_obj <- obj[, query_idx]
    
    combined <- steve_integrate(
      ref_obj, query_obj,
      integration = integration,
      dims = dims,
      cluster_resolution = cluster_resolution,
      seed = seed + iter
    )
    
    query_cells <- colnames(subset(combined, subset = orig.ident == "query"))
    cand_q <- combined@meta.data[query_cells, candidate_col]
    truth_q <- combined@meta.data[query_cells, truth_col]
    
    tool_metrics <- per_celltype_sens_spec(truth_q, cand_q)
    write.csv(tool_metrics, file.path(outputdir, sprintf("iter_%02d_candidate_vs_truth_per_celltype.csv", iter)),
              row.names = FALSE)
    
    cand_macro_sens <- mean(tool_metrics$Sensitivity, na.rm = TRUE)
    cand_macro_spec <- mean(tool_metrics$Specificity, na.rm = TRUE)
    
    model <- steve_build_kde_model(combined, ref_ident = "reference", ref_label_col = truth_col, grid_n = grid_n)
    pred <- steve_predict(combined, model, odds_cutoff = odds_cutoff, unknown_label = "unknown")
    pred_q <- pred[pred$cellid %in% query_cells, , drop = FALSE]
    
    steve_metrics <- per_celltype_sens_spec(truth_q, pred_q$pred)
    write.csv(steve_metrics, file.path(outputdir, sprintf("iter_%02d_steve_vs_truth_per_celltype.csv", iter)),
              row.names = FALSE)
    
    steve_macro_sens <- mean(steve_metrics$Sensitivity, na.rm = TRUE)
    steve_macro_spec <- mean(steve_metrics$Specificity, na.rm = TRUE)
    
    mat_tool <- confusion_ratio(truth_q, cand_q)
    plot_heatmap(
      mat_tool,
      out_pdf = file.path(outputdir, sprintf("iter_%02d_candidate_confusion_heatmap.pdf", iter)),
      out_csv = file.path(outputdir, sprintf("iter_%02d_candidate_confusion_heatmap.csv", iter))
    )
    
    all_runs[[length(all_runs) + 1]] <- data.frame(
      iteration = iter,
      candidate_macro_sensitivity = cand_macro_sens,
      candidate_macro_specificity = cand_macro_spec,
      steve_macro_sensitivity = steve_macro_sens,
      steve_macro_specificity = steve_macro_spec,
      stringsAsFactors = FALSE
    )
    
    msg("Benchmark iter=%d candidate(sens=%.3f spec=%.3f) steve(sens=%.3f spec=%.3f)",
        iter, cand_macro_sens, cand_macro_spec, steve_macro_sens, steve_macro_spec)
  }
  
  out_df <- do.call(rbind, all_runs)
  write.csv(out_df, file.path(outputdir, "benchmark_summary.csv"), row.names = FALSE)
  
  out_df
}

# ============================================================
# Module 4: Reference Transfer Annotation
# ============================================================

steve_module_transfer <- function(user_obj,
                                  ref_obj,
                                  ref_truth_col,
                                  outputdir = "./output_transfer",
                                  integration = "harmony",
                                  dims = 1:20,
                                  cluster_resolution = 0.8,
                                  grid_n = 200,
                                  odds_cutoff = 2,
                                  min_cells = 10,
                                  seed = 123,
                                  user_truth_col = NULL) {
  
  ensure_dir(outputdir)
  ref_obj <- filter_min_cells(ref_obj, ref_truth_col, min_cells = min_cells)
  
  combined <- steve_integrate(
    ref_obj, user_obj,
    integration = integration,
    dims = dims,
    cluster_resolution = cluster_resolution,
    seed = seed
  )
  
  combined$STEVE_plot <- "user"
  ref_cells <- colnames(subset(combined, subset = orig.ident == "reference"))
  combined@meta.data[ref_cells, "STEVE_plot"] <- combined@meta.data[ref_cells, ref_truth_col]
  combined$STEVE_plot <- factor(combined$STEVE_plot)
  
  plot_integrated_umap(
    combined, "STEVE_plot",
    file.path(outputdir, "integrated_umap.pdf"),
    title = "Reference Transfer Annotation"
  )
  
  model <- steve_build_kde_model(combined, ref_ident = "reference", ref_label_col = ref_truth_col, grid_n = grid_n)
  pred <- steve_predict(combined, model, odds_cutoff = odds_cutoff, unknown_label = "unknown")
  
  user_cells <- colnames(subset(combined, subset = orig.ident == "query"))
  pred_u <- pred[pred$cellid %in% user_cells, , drop = FALSE]
  write.csv(pred_u, file.path(outputdir, "transfer_predictions_user.csv"), row.names = FALSE)
  
  if (!is.null(user_truth_col)) {
    if (!user_truth_col %in% colnames(combined@meta.data)) {
      stop("user_truth_col not in user meta.data: ", user_truth_col, "\nAvailable columns: ",
           paste(colnames(combined@meta.data), collapse = ", "))
    }
    truth_u <- combined@meta.data[pred_u$cellid, user_truth_col]
    metrics <- per_celltype_sens_spec(truth_u, pred_u$pred)
    write.csv(metrics, file.path(outputdir, "transfer_metrics_user_vs_truth.csv"), row.names = FALSE)
    
    mat <- confusion_ratio(truth_u, pred_u$pred)
    plot_heatmap(
      mat,
      out_pdf = file.path(outputdir, "transfer_confusion_heatmap.pdf"),
      out_csv = file.path(outputdir, "transfer_confusion_heatmap.csv")
    )
  }
  
  pred_u
}

# ============================================================
# CLI
# ============================================================

option_list <- list(
  make_option(c("-m", "--module"), type = "character", default = "subsampling",
              help = "Module: subsampling | novel | benchmark | transfer [default %default]"),
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Working directory [default %default]"),
  make_option(c("-u", "--userdata"), type = "character", default = NULL,
              help = "User dataset (.seurat.rds). Required for all modules."),
  make_option(c("-g", "--groundtruthdata"), type = "character", default = NULL,
              help = "Reference dataset (.seurat.rds). Required for transfer module."),
  make_option(c("-o", "--outputdir"), type = "character", default = "./output",
              help = "Output directory [default %default]"),
  
  make_option(c("-f", "--reference_annotation"), type = "character", default = "groundtruth",
              help = "Ground-truth label column in meta.data [default %default]"),
  make_option(c("-a", "--user_annotation"), type = "character", default = "user_annotation",
              help = "Candidate annotation column in meta.data (used by benchmark) [default %default]"),
  make_option(c("--user_truth_annotation"), type = "character", default = NULL,
              help = "Optional: user truth label column (for transfer evaluation) [default %default]"),
  
  make_option(c("-i", "--integration"), type = "character", default = "harmony",
              help = "Integration: harmony | none [default %default]"),
  make_option(c("--dims"), type = "character", default = "1:20",
              help = "Dimensions for PCA/Harmony/UMAP, like '1:20' [default %default]"),
  make_option(c("-r", "--cluster_resolution"), type = "double", default = 0.8,
              help = "FindClusters resolution [default %default]"),
  make_option(c("--grid_n"), type = "integer", default = 200,
              help = "KDE grid resolution (kde2d n=) [default %default]"),
  make_option(c("--odds_cutoff"), type = "double", default = 2,
              help = "Posterior odds cutoff; else 'unknown' [default %default]"),
  make_option(c("-x", "--min_cells"), type = "integer", default = 10,
              help = "Filter cell types with <= min_cells [default %default]"),
  make_option(c("--seed"), type = "integer", default = 123,
              help = "Random seed [default %default]"),
  
  make_option(c("--ratios"), type = "character",
              default = "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
              help = "Subsampling reference fractions (comma-separated) [default %default]"),
  make_option(c("--n_repeats"), type = "integer", default = 10,
              help = "Number of repeats (subsampling/novel/benchmark) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

options(future.globals.maxSize = 8000 * 1024^2)

setwd(opt$workdir)
ensure_dir(opt$outputdir)

dims <- eval(parse(text = opt$dims))
ratios <- as.numeric(strsplit(opt$ratios, ",")[[1]])

if (is.null(opt$userdata)) stop("--userdata is required.")

msg("STEVE master starting. module=%s output=%s", opt$module, opt$outputdir)

if (opt$module %in% c("subsampling", "novel", "benchmark")) {
  obj <- read_seurat_rds(opt$userdata, assay = "RNA")
}

if (opt$module == "subsampling") {
  
  steve_module_subsampling(
    obj = obj,
    truth_col = opt$reference_annotation,
    ratios = ratios,
    n_repeats = opt$n_repeats,
    outputdir = opt$outputdir,
    integration = opt$integration,
    dims = dims,
    cluster_resolution = opt$cluster_resolution,
    grid_n = opt$grid_n,
    odds_cutoff = opt$odds_cutoff,
    min_cells = opt$min_cells,
    seed = opt$seed
  )
  
} else if (opt$module == "novel") {
  
  steve_module_novel_cell(
    obj = obj,
    truth_col = opt$reference_annotation,
    outputdir = opt$outputdir,
    n_repeats = opt$n_repeats,
    integration = opt$integration,
    dims = dims,
    cluster_resolution = opt$cluster_resolution,
    grid_n = opt$grid_n,
    odds_cutoff = opt$odds_cutoff,
    min_cells = opt$min_cells,
    seed = opt$seed
  )
  
} else if (opt$module == "benchmark") {
  
  steve_module_benchmark(
    obj = obj,
    truth_col = opt$reference_annotation,
    candidate_col = opt$user_annotation,
    outputdir = opt$outputdir,
    n_repeats = opt$n_repeats,
    integration = opt$integration,
    dims = dims,
    cluster_resolution = opt$cluster_resolution,
    grid_n = opt$grid_n,
    odds_cutoff = opt$odds_cutoff,
    min_cells = opt$min_cells,
    seed = opt$seed
  )
  
} else if (opt$module == "transfer") {
  
  if (is.null(opt$groundtruthdata)) stop("--groundtruthdata is required for module=transfer")
  
  user_obj <- read_seurat_rds(opt$userdata, assay = "RNA")
  ref_obj  <- read_seurat_rds(opt$groundtruthdata, assay = "RNA")
  
  steve_module_transfer(
    user_obj = user_obj,
    ref_obj = ref_obj,
    ref_truth_col = opt$reference_annotation,
    outputdir = opt$outputdir,
    integration = opt$integration,
    dims = dims,
    cluster_resolution = opt$cluster_resolution,
    grid_n = opt$grid_n,
    odds_cutoff = opt$odds_cutoff,
    min_cells = opt$min_cells,
    seed = opt$seed,
    user_truth_col = opt$user_truth_annotation
  )
  
} else {
  stop("Unknown --module: ", opt$module, " (use subsampling|novel|benchmark|transfer)")
}

msg("OK!")
