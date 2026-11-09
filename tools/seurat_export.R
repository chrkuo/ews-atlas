# Export a Seurat (v4 or v5) object into a folder that tools/folder_to_h5ad.py
# turns into an .h5ad for the atlas exporter.
#
# Usage (in R or RStudio):
#   source("tools/seurat_export.R")
#   obj <- readRDS("my_object.rds")
#   export_for_atlas(obj, out_dir = "export_tumor",
#                    reduction = "umap",                       # see Reductions(obj)
#                    meta_cols = c("cell_type", "timepoint"),  # columns of obj@meta.data
#                    colors = list(cell_type = my_palette))    # optional named vectors
#
# Writes: matrix.mtx (log-normalized "data" layer, genes x cells), genes.txt,
#         cells.txt, metadata.csv, embedding.csv, colors.json (if given)

export_for_atlas <- function(obj, out_dir, reduction = "umap", meta_cols,
                             assay = "RNA", colors = NULL) {
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  DefaultAssay(obj) <- assay
  # Seurat v5: split layers (data.1, data.2, ...) must be joined first
  if (inherits(obj[[assay]], "Assay5") && length(Layers(obj[[assay]], search = "data")) > 1) {
    obj[[assay]] <- JoinLayers(obj[[assay]])
  }
  if (inherits(obj[[assay]], "Assay5")) {
    if (!"data" %in% Layers(obj[[assay]])) stop("No 'data' layer: run NormalizeData() first.")
    X <- LayerData(obj, assay = assay, layer = "data")
  } else {
    X <- methods::slot(obj[[assay]], "data")   # Seurat v4-style Assay; works with any SeuratObject version
  }
  X <- as(X, "CsparseMatrix")

  if (!reduction %in% Reductions(obj))
    stop("Reduction '", reduction, "' not found. Available: ", paste(Reductions(obj), collapse = ", "))
  emb <- Embeddings(obj, reduction)[colnames(X), 1:2]

  missing <- setdiff(meta_cols, colnames(obj@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  meta <- obj@meta.data[colnames(X), meta_cols, drop = FALSE]
  # keep factor level order (it becomes legend order on the site)
  for (m in meta_cols) if (!is.factor(meta[[m]])) meta[[m]] <- factor(meta[[m]])

  message("Writing ", ncol(X), " cells x ", nrow(X), " genes to ", out_dir)
  writeMM(X, file.path(out_dir, "matrix.mtx"))
  writeLines(rownames(X), file.path(out_dir, "genes.txt"))
  writeLines(colnames(X), file.path(out_dir, "cells.txt"))
  write.csv(data.frame(cell = colnames(X), emb), file.path(out_dir, "embedding.csv"), row.names = FALSE)
  write.csv(data.frame(cell = colnames(X), meta, check.names = FALSE),
            file.path(out_dir, "metadata.csv"), row.names = FALSE)
  levs <- lapply(meta, levels)
  jsonlite::write_json(levs, file.path(out_dir, "levels.json"), auto_unbox = FALSE)
  if (!is.null(colors)) {
    # each entry: a named character vector, names = levels, values = hex colors
    jsonlite::write_json(lapply(colors, as.list), file.path(out_dir, "colors.json"), auto_unbox = TRUE)
  }
  message("Done.")
  invisible(out_dir)
}
