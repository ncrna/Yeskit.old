#' Differential gene expression analysis for each cluster
#' Perform differential analysis (DE) between two groups for each clusters
#' @param object Seurat object.
#' @param comparison Vector of nominal variables from group.by. Eg., comparison=c("Normal", "Tumor")
#' @param group.by Regroup cells before performing DGE, default group.by='group'
#' @param min.cells Minimum number of cells to perform DGE, default min.cells=20
#' @param logFC A positive value to set the cutoff of logFC, default logFC=0.25
#' @param clusters Vector of Clusters to perform DGE, default clusters=NULL,
#' which means all clusters will be performed
#' @return DGE Table
#'
#' @author rstatistics
#' @export
scDGE <- function(object = object, comparison = c("condA", "condB"), group.by=NULL, min.cells = 20, logFC = 0.25, clusters = NULL){
  results <- list()
  seurat_clusters <- group <- sample <- NULL
  if (is.null(clusters))
    clusters = levels(object)
  if (length(comparison) != 2){
    stop("Comparison must have 2 elements!")
  }
  condA <- comparison[1]
  condB <- comparison[2]
  if (is.null(group.by)){
    group.by <- "group"
  }
  clusters <- as.character(clusters)
  DE <- data.frame()
  for (cluster in clusters){
    cat(paste0("### ", "Comparing cluster-", cluster, " between ", condA, " and ", condB, " ...\n"))
    Object <- subset(object, seurat_clusters==eval(quote(cluster)))
    if (all(comparison %in% names(table(Object[[group.by]])))){
      Object <- Object[, Object@meta.data[[group.by]] %in% comparison]
    }else{
      warning("Cluster ", cluster, " has no cell in ", condA, " or ", condB, ". Ignored.\n")
      next
    }
    if (!all(as.data.frame(table(Object[[group.by]]))$Freq >= min.cells)){
      warning("Cell number in cluster ", cluster, " is less than ", min.cells, " cells. Ignored.\n")
      next
    }else if (!condA %in% names(table(Object$group))){
      warning("condA has no cell in cluster ", cluster, ". Ignored.\n")
      next
    }else if (!condB %in% names(table(Object$group))){
      warning("condB has no cell in cluster ", cluster, ". Ignored.\n")
      next
    }
    tmpDE <- suppressWarnings(
      Seurat::FindMarkers(object = Object, assay = 'RNA', ident.1 = condA, ident.2 = condB, group.by = group.by,
                  min.cells.group = min.cells, logfc.threshold = logFC, test.use = 'MAST', only.pos = FALSE))
    tmpDE$cluster <- cluster
    tmpDE$gene <- rownames(tmpDE)
    DE <- rbind(DE, tmpDE)
    cat(" done.\n")
  }
  return(DE)
}
