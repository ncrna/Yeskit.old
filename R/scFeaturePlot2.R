#' scFeaturePlot
#' wrapper for Seurat::FeaturePlot
#' @title scFeaturePlot
#' @param object Seurat object
#' @param reduction Dimensional reduction to use
#' @param features Vector of features to plot
#' @param pt.size Point size to plot
#' @param font.size Font size
#' @param colors Vector of two colors to form the gradient over
#' @param ncol Number of columns to combine features
#' @param combine Combine plots into a single patchworked ggplot object
#' @param order Determine whether to plot cells in order of expression
#' @return NULL
#'
#' @author rstatistics
#' @export

scFeaturePlot2 <- function(object=object, reduction=NULL, features=features, pt.size=NULL, font.size=NULL,
                         colors=NULL, ncol=NULL, combine=FALSE, order=FALSE){
  if (is.null(reduction)){
    if ("umap" %in% names(object)){
      reduction <- "umap"
    }else if ("tsne" %in% names(object)){
      reduction <- "tsne"
    }else{
      stop(paste0("The reduction parameter does not support! Please use 'umap', 'tsne', or 'pca' instead.\n"))
    }
  }
  if(is.null(pt.size)){
    pt.size = 0.1
  }
  if(is.null(colors)){
    colors = c("lightgrey","#FF0000")
  }
  if(is.null(font.size)){
    font.size = 9
  }
  features <- intersect(unique(features), rownames(object))
  if(is.null(ncol)){
    ncol = ceiling(sqrt(length(features)))
  }
  pp = Seurat::FeaturePlot(object = object, reduction = reduction, features = features, pt.size = pt.size,
                           cols = colors, combine = combine, order = order)
  plots <- lapply(X = pp, FUN = function(p){
    p + ggplot2::theme(axis.title = ggplot2::element_text(size = font.size),
                       axis.text = ggplot2::element_text(size = font.size),
                       plot.title = ggplot2::element_text(family = 'sans',face='italic',size=10),
                       legend.text = ggplot2::element_text(size = 10),
                       legend.key.height = ggplot2::unit(0.9,"line"),
                       legend.key.width = ggplot2::unit(0.6,"line"))})
  return(patchwork::wrap_plots(plots, ncol = ncol))
}
