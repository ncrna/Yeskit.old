#' @title scDensityPlot function
#' @description Visualize density plot in low dimensional embedding
#' @param object Seurat object
#' @param reduction Dimensional reduction to use, default reduction="umap"
#' @param title Figure title
#' @param split.by Name of a metadata column to split plot by
#' @param ncol Number of columns for display the plots
#' @param colors Vector of colours to use for n-colour gradient
#'
#'
#' @return ggplot object
#'
#' @import Seurat

#' @author Wei Zhang
#' @export
#'
scDensityPlot2 <- function(object=NULL, reduction="umap", title=NULL, split.by=NULL, ncol=NULL, colors=NULL){
    z.min <- z.max <- ..density.. <- pp <- NULL
    if (is.null(reduction)){
      if ("umap" %in% names(object)){
        reduction <- "umap"
      }else if ("tsne" %in% names(object)){
        reduction <- "tsne"
      }else{
        stop(paste0("The reduction parameter does not support! Please use 'umap', 'tsne', or 'pca' instead.\n"))
      }
    }
    pp <- function(data, z.min=z.min, z.max=z.max, title=NULL){
      if (is.null(colors)){
        colors <- c("white", "gray99", "#FFEDA0", "red", "darkred")
      }
      if (is.null(title)){
        title <- "Density plot"
      }
      p <- ggplot2::ggplot(data=data, ggplot2::aes(x=data[,1], y=data[,2]))
      p <- p + ggplot2::stat_density2d(geom = 'raster', interpolate=TRUE, ggplot2::aes(fill = ..density..), contour = FALSE) +
        ggplot2::scale_fill_gradientn(colors=colors, limits = c(z.min, z.max), breaks=round(seq(z.min, z.max, length=5), digits=5)) +
        cowplot::theme_cowplot() + ggplot2::ggtitle(title) + ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(size=5))) +
        ggplot2::labs(x=colnames(data)[1], y=colnames(data)[2])
      return(p)
    }
    Data <- Seurat::Embeddings(object = object[[reduction]])[, 1:2]
    Data <- as.data.frame(Data)
    m <- MASS::kde2d(x = Data[,1], y = Data[,2], h = c(MASS::bandwidth.nrd(Data[,1]), MASS::bandwidth.nrd(Data[,2])), n = 100, lims = c(range(Data[,1]), range(Data[,2])))
    z.min <- min(m$z)
    z.max <- max(m$z)
    if (is.null(x = split.by)){
      return(pp(Data, z.min, z.max, title))
    }
    plots <- list()
    z <- vector()
    if (! is.null(x = split.by)){
      Data[, split.by] <- object[[split.by, drop=FALSE]]
      if(is.null(ncol)){
        ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
      }
      for (s in unique(Data[, split.by])){
        data <- Data[Data[,3] == s,]
        m <- MASS::kde2d(x = data[,1], y = data[,2], h = c(MASS::bandwidth.nrd(data[,1]), MASS::bandwidth.nrd(data[,2])), n = 100, lims = c(range(data[,1]), range(data[,2])))
        z <- c(z, range(m$z))
      }
      z <- sort(z)
      z.min <- z[1]
      z.max <- z[length(z)]
      for (s in unique(Data[, split.by])){
        data <- Data[Data[,3] == s,]
        plots[[s]] <- pp(data, z.min, z.max, title = s)
      }
    }
    return(patchwork::wrap_plots(plots, ncol = ncol))
}
