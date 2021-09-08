#' Average Heatmap for specific feature in all clusters
#' @title scHeatmap visualization
#' @param object Seurat object
#' @param key The key slot that stored the DGE table
#' @param cluster Cluster specified by the user
#' @param feature.by Feature which is used to divid cells into two groups
#' @param markers A vector of features to plot, default is top_n features
#' @param colors User defined colors
#' @param top_n Top n markers to plot
#' @param font.size Font size
#' @param assay Assay to pull from
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @return NULL
#'
#' @author rstatistics
#' @export
scHeatmap <- function(object = NULL, key = NULL, cluster = NULL, feature.by = NULL, markers = NULL,
                     colors = NULL, top_n = 5, font.size = 8, assay = NULL, slot = NULL){
  seurat_clusters <- avg_log2FC <- p_val_adj <- NULL
  if (is.null(colors)){
    colors <- c("white", "red")
  }
  if (is.null(assay)){
    assay = "RNA"
  }
  if (is.null(slot)){
    slot = "data"
  }
  Seurat::DefaultAssay(object = object) <- assay
  cluster = as.character(cluster)
  object <- subset(object, seurat_clusters %in% cluster)
  object$group <- ifelse(object@meta.data[, feature.by] > 0, "Pos", "Neg")
  data <- object@misc[[key]][[feature.by]]
  if ("avg_logFC" %in% colnames(data)){
    colnames(data)[which(colnames(data)=="avg_logFC")] <- "avg_log2FC"
  }
  if (is.null(markers)){
    markers.up <- rownames(data[data$cluster==cluster & data$avg_log2FC > 0 & data$p_val_adj < 1e-2, ])[1:top_n]
    markers.dn <- rownames(data[data$cluster==cluster & data$avg_log2FC < 0 & data$p_val_adj < 1e-2, ])[1:top_n]
    markers <- unique(c(markers.up, markers.dn))
    if (length(markers) < 1){
      warning("Some markers does not exist!")
    }else if (length(markers)==0){
      stop("No markers left for heatmap visualization!")
    }
  }else{
    markers_len <- length(unique(markers))
    markers = unique(intersect(markers, rownames(object)))
    if (length(markers) < markers_len){
      warning("Some markers does not exist!")
    }else if (length(markers)==0){
      stop("No markers left for heatmap visualization!")
    }
  }
  #markers <- intersect(markers, VariableFeatures(object = object))
  if(length(markers)==0){
    stop("All markers do not exist in the variable features!\n")
  }
  return(suppressMessages(Seurat::DoHeatmap(object, features = markers, group.by = "group", assay = "RNA", slot = "data") +
                            ggplot2::scale_fill_gradientn(colors = colors) + ggplot2::ggtitle(feature.by)))
}
