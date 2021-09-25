#' Find differentially expressed genes by species for each clusters in a dataset
#' @param object Seurat object
#' @param species.by A vector of features which is used to divid cells into two groups
#' @param clusters Clusters selected to do DGE
#' @param min.cells Minimum number of cells to do DGE
#' @return DGE table
#'
#' @author Wei Zhang
#' @export
scPathogenDGE <- function(object = NULL, species.by = NULL, 
                          clusters = NULL, min.cells = 20) {
  if (is.null(object)) {
    stop("Parameter 'object' must be specified!\n")
  }
  if (is.null(clusters)) {
    warning("All clusters will be evaluated!\n")
    clusters = levels(object)
  }
  clusters = as.character(clusters)
  if (!species.by %in% names(object@meta.data)) {
    stop("The feature ", species.by, " does not exist in MetaData slot!\n")
  }
  results <- list()
  seurat_clusters <- p_val_adj <- NULL
  for (feature in species.by) {
    message("### ", "Analysis feature ", feature, " ...\n")
    if (is.null(results[[feature]])) {
      results[[feature]] = list()
    }
    if (!feature %in% colnames(object@meta.data)) {
      message("all cells with 0 reads.\n")
      next
    }
    DE <- data.frame()
    for (cluster in clusters) {
      message("====== ", "Analysis cluster-", cluster, " ... ")
      Object <- subset(object, seurat_clusters == cluster)
      Object$group <- ifelse(Object@meta.data[, feature] > 0, "Pos", "Neg")
      Expr <- as.matrix(Seurat::GetAssayData(Object))
      if (length(Object$group[Object$group == "Pos"]) < min.cells) {
        message("too few cells to process.\n")
        next
      }
      if (length(Object$group[Object$group == "Neg"]) < min.cells) {
        message("too few cells to process.\n")
        next
      }
      tmpDE <- suppressMessages(suppressWarnings(Seurat::FindMarkers(object = Object,
                assay = "RNA", ident.1 = "Pos", ident.2 = "Neg", group.by = "group",
                only.pos = FALSE, min.pct = 0.25, test.use = "MAST", verbose = FALSE)))
      tmpDE$cluster <- cluster
      tmpDE$gene <- rownames(tmpDE)
      DE <- rbind(DE, tmpDE)
      message("done.\n")
    }
    results[[feature]] <- DE
  }
  # clean up results
  for (feature in species.by) {
    if (length(results[[feature]]) == 0) {
      results[[feature]] <- NULL
      next
    }
    for (cluster in clusters) {
      if (length(results[[feature]][[cluster]]) == 0) {
        next
      }
      results[[feature]][[cluster]] <- subset(results[[feature]][[cluster]],
        p_val_adj <= 1)
      if (nrow(results[[feature]][[cluster]]) <= 1) {
        results[[feature]][[cluster]] <- NULL
        next
      }
    }
    if (length(results[[feature]]) == 0) {
      results[[feature]] <- NULL
      next
    }
  }
  return(results)
}
