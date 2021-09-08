#' MSigdb scoring DimPlot for single-cell clusters
#' @title scScoreDimPlot
#' @param object Seurat object
#' @param features Name of gene symbols
#' @param reduction	Which dimensionality reduction to use
#' @param cols Colors to plot
#' @param pt.size Adjust point size to plot, default pt.size=0.5
#' @param split.by Name of a metadata column to split plot by
#' @param title Title of the plot
#' @param ncol Number of columns for display the plots
#' @param raster Convert points to raster format, default is TRUE
#' @param scale Determine whether to scale the data, default is TRUE
#' @param col.min Minimum scaled average score threshold (smaller values will be set to this)
#' @param col.max Maximum scaled average score threshold (larger values will be set to this)
#' @return ggplot2 object
#' @export

scFeaturePlot <- function(object=NULL, features=NULL, reduction=NULL, cols=NULL, pt.size=NULL, split.by=NULL, title=NULL, ncol=NULL, raster=TRUE, scale=TRUE, col.min=NA, col.max=NA){
  if (is.null(features)){
    stop("Parameter 'features' must be specified!\n")
  }else if (length(features) != 1){
    stop("Parameter 'features' must be one pathway!\n")
  }
  if (is.null(reduction)){
    if ("umap" %in% names(object)){
      reduction <- "umap"
    }else if ("tsne" %in% names(object)){
      reduction <- "tsne"
    }else{
      stop(paste0("The reduction parameter does not support! Please use 'umap', 'tsne', or 'pca' instead.\n"))
    }
  }
  if (is.null(cols)){
    cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))(32)
  }
  if (is.null(pt.size)){
    pt.size <- 0.5
  }
  xmin <- xmax <- ymin <- ymax <- NULL
  ps <- function(data, min.value=min.value, max.value=max.value, title=NULL, legend_title=NULL){
    data <- suppressWarnings(data.frame(apply(data, 2, as.numeric)))
    p <- ggplot2::ggplot(data=data, ggplot2::aes(x=data[,1], y=data[,2], z=data[,3]))
    if (isTRUE(raster)){
      p <- p + ggrastr::rasterise(dpi = 300, ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size, na.rm = TRUE))
    }else{
      p <- p + ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size, na.rm = TRUE)
    }
    p <- p + ggplot2::scale_colour_gradientn(colours = cols, na.value = "white",
                                      breaks=seq(min.value, max.value, length=5),
                                      limit = c(min.value,max.value),
                                      guide = ggplot2::guide_colorbar(title = legend_title, order = 1,
                                                                      title.theme = ggplot2::element_text(size=8),
                                                                      label.theme = ggplot2::element_text(size=8)))
    p <- p + cowplot::theme_cowplot() + ggplot2::labs(x=colnames(data)[1], y=colnames(data)[2], title=title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin/5)*5, ceiling(xmax/5)*5,by = 5)) +
         ggplot2::scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin/5)*5, ceiling(ymax/5)*5,by = 5))
    return(p)
  }
  pm <- function(data, min.value=min.value, max.value=max.value, title=NULL, legend_title=NULL){
    data <- suppressWarnings(data.frame(apply(data, 2, as.numeric)))
    p <- ggplot2::ggplot(data=data, ggplot2::aes(x=data[,1], y=data[,2], z=data[,3]))
    if (isTRUE(raster)){
      p <- p + ggrastr::rasterise(dpi = 300, ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size, na.rm = TRUE))
    }else{
      p <- p + ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size, na.rm = TRUE)
    }
    p <- p + ggplot2::scale_colour_gradientn(colours = cols, na.value = "white",
                                             breaks=seq(min.value, max.value, length=5),
                                             limit = c(min.value,max.value))
    p <- p + cowplot::theme_cowplot() + ggplot2::labs(x=colnames(data)[1], y=colnames(data)[2], title=title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin/5)*5, ceiling(xmax/5)*5,by = 5)) +
         ggplot2::scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin/5)*5, ceiling(ymax/5)*5,by = 5))
    p <- p + ggplot2::theme(legend.position = "none")
    return(p)
  }
  reduction_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
  if (is.null(x = split.by)){
    Data <- object@meta.data[, c(reduction_ids, features)]
    #data <- Seurat::FetchData(object = object, vars = c(dims, "ident", features), cells = cells, slot = slot)
    Data <- object@meta.data[, reduction_ids]
    # 获取某个基因的表达量
    features <- intersect(features, rownames(Seurat::GetAssayData(object = object, assay = "RNA", slot = "scale.data")))
    Expr <- Seurat::GetAssayData(object = object, assay = "RNA", slot = "scale.data")[features, ]
    Data <- rbind(Data, t(Expr))
    Data[, features][is.na(Data[, features])] <- 0
    if (is.na(col.min)){
      col.min <- round(min(Data[, features]), 1)
    }
    if (is.na(col.max)){
      col.max <- round(max(Data[, features]), 1)
    }
    Data[, features][Data[, features] > col.max] <- col.max
    Data[, features][Data[, features] < col.min] <- col.min
    xmin <- min(Data[,1])
    xmax <- max(Data[,1])
    ymin <- min(Data[,2])
    ymax <- max(Data[,2])
    return(ps(data=Data, min.value = col.min, max.value = col.max, title=title, legend_title = features))
  }
  plots <- list()
  if (! is.null(x = split.by)){
    if (! split.by %in% colnames(object@meta.data)){
      stop(paste0("The parameter 'split.by' ", split.by, " does not exist in MetaData slot!\n"))
    }
    Data <- object@meta.data[, c(reduction_ids, features, split.by)]
    if (scale) {
      Data[, features] <- scale(x = Data[, features])
    }
    Data[, features][is.na(Data[, features])] <- 0
    if (is.na(col.min)){
      col.min <- round(min(Data[, features]), 1)
    }
    if (is.na(col.max)){
      col.max <- round(max(Data[, features]), 1)
    }
    Data[, features][Data[, features] > col.max] <- col.max
    Data[, features][Data[, features] < col.min] <- col.min
    xmin <- min(Data[,1])
    xmax <- max(Data[,1])
    ymin <- min(Data[,2])
    ymax <- max(Data[,2])
    if(is.null(ncol)){
      ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    }
    legend <- ps(data=Data, min.value = col.min, max.value = col.max, title=NULL, legend_title=features)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, 'box', trim=F)  # use trim depending on need
    for (s in unique(Data[, split.by])){
      data <- Data[Data[, split.by] == s,]
      plots[[s]] <- pm(data, min.value = col.min, max.value = col.max, title = s, legend_title = NULL)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | legend) + patchwork::plot_layout(widths = c(3,1)))
}
