#' DimPlot with rasterized point for single-cell visualization
#' @title scDimPlot
#' @param object Seurat object
#' @param cols Colors to plot, use ggplot2's default colors by default. We include a pallete called 'sc' which consists of 36 colors
#' @param pt.size Adjust point size to plot, default pt.size=0.5
#' @param reduction	Which dimensionality reduction to use
#' @param split.by Name of a metadata column to split plot by
#' @param label Whether to label the clusters
#' @param title Title of the plot
#' @param ncol Number of columns for display the plots
#' @param raster Convert points to raster format, default is TRUE
#' @return ggplot2 object
#' @export

scDimPlot <- function(object=NULL, cols=NULL, pt.size=NULL, reduction=NULL, split.by=NULL, label=TRUE, title=NULL, ncol=NULL, raster=TRUE){
  if (is.null(cols)){
    cols <- scales::hue_pal()(length(levels(Seurat::Idents(object))))
  }else if (cols == "sc"){
    if (length(levels(Seurat::Idents(object))) <= 36){
      cols <- c("#1660A7","#FF6A00","#219418","#CD0C18","#814BB2","#794339","#DC59B6","#CC79A7","#FF0000","#11B3C6",
                "#AFB400","#00FFFF", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00",
                "#CC79A7", "#00AFBB", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#4477AA",
                "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    }else{
      warning(paste0("Not enough colours provided for ", length(levels(Seurat::Idents(object))), " clusters! Use ggplot2's default colors instead\n"))
      cols <- scales::hue_pal()(length(levels(Seurat::Idents(object))))
    }
  }else if (length(levels(Seurat::Idents(object))) > length(cols)){
      stop(paste0("Not enough colours provided for ", length(levels(Seurat::Idents(object))), " clusters!"))
  }
  if (is.null(pt.size)){
    pt.size <- 0.5
  }
  if (is.null(reduction)){
    if ("umap" %in% names(object)){
      reduction <- "umap"
    }else if ("tsne" %in% names(object)){
      reduction <- "tsne"
    }else if ("pca" %in% names(object)){
      reduction <- "pca"
    }else{
      stop(paste0("The reduction parameter does not support! Please use 'umap', 'tsne', or 'pca' instead.\n"))
    }
  }
  xmin <- xmax <- ymin <- ymax <- NULL
  ps <- function(data, title=NULL, legend_title=NULL){
    p <- ggplot2::ggplot(data=data, ggplot2::aes(x=data[,1], y=data[,2], z=data[,3]))
    if (isTRUE(raster)){
      p <- p + ggrastr::rasterise(dpi = 300, ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size))
    }else{
      p <- p + ggplot2::geom_point(ggplot2::aes(colour=data[,3]), size = pt.size)
    }
    if (isTRUE(label)){
      p <- p + ggplot2::geom_text(ggplot2::aes(label=label), na.rm = TRUE)
    }
    p <- p + ggplot2::scale_colour_manual(values = cols, na.value = "white")
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title=NULL, override.aes = list(size = 3))) +
             ggplot2::labs(x=colnames(data)[1], y=colnames(data)[2], title=legend_title)
    p <- p + cowplot::theme_cowplot() + ggplot2::ggtitle(title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin/5)*5, ceiling(xmax/5)*5,by = 5)) +
             ggplot2::scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin/5)*5, ceiling(ymax/5)*5,by = 5))
    return(p)
  }
  pm <- function(data, title=NULL, legend_title=NULL){
    p <- ggplot2::ggplot(data=data, ggplot2::aes(x=data[,1], y=data[,2], z=data[,"color"]))
    if (isTRUE(raster)){
      p <- p + ggrastr::rasterise(dpi = 300, ggplot2::geom_point(colour=data[,"color"], size = pt.size))
    }else{
      p <- p + ggplot2::geom_point(colour=data[,"color"], size = pt.size)
    }
    if (isTRUE(label)){
      p <- p + ggplot2::geom_text(ggplot2::aes(label=label), na.rm = TRUE)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title=NULL, override.aes = list(size = 3))) +
             ggplot2::labs(x=colnames(data)[1], y=colnames(data)[2], title=title)
    p <- p + cowplot::theme_cowplot() + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin/5)*5, ceiling(xmax/5)*5,by = 5)) +
             ggplot2::scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin/5)*5, ceiling(ymax/5)*5,by = 5))
    return(p)
  }
  GetXYCenter <- function(data){
    cluster <- NULL
    result <- data.frame(x=NULL, y=NULL, ident=NULL)
    for (c in as.vector(unique(data[,"cluster"]))){
      min_x <- max_x <- min_y <- max_y <- center_x <- center_y <- NULL
      sub_data <- subset(data, cluster==c)
      center_x <- median(sub_data[,1])
      center_y <- median(sub_data[,2])
      label <- c
      result <- rbind(result, data.frame(x=center_x, y=center_y, label=label))
    }
    return(result)
  }
  reduction_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
  coord <- Seurat::Embeddings(object = object, reduction = reduction)[, 1:2]
  object <- Seurat::AddMetaData(object = object, metadata = coord, col.name = reduction_ids)
  object <- Seurat::AddMetaData(object = object, metadata = as.vector(Seurat::Idents(object)), col.name = "cluster")
  if (is.null(x = split.by)){
    Data <- object@meta.data[, c(reduction_ids, "cluster")]
    Data[, "cluster"] <- factor(Data[, "cluster"], levels=levels(Seurat::Idents(object)))
    Data[, "color"] <- Data[, "cluster"]
    levels(Data[, "color"]) <- cols[1:length(levels(Seurat::Idents(object)))]
    nearest.point <- RANN::nn2(data=Data[,1:2], query=GetXYCenter(Data)[,1:2], k=1)$nn.idx
    Data[, "label"] <- NA
    Data[nearest.point, "label"] <- as.vector(GetXYCenter(Data)[, "label"])
    xmin <- min(Data[,1])
    xmax <- max(Data[,1])
    ymin <- min(Data[,2])
    ymax <- max(Data[,2])
    return(ps(data=Data, title=title, legend_title = NULL))
  }
  plots <- list()
  if (! is.null(x = split.by)){
    if (! split.by %in% colnames(object@meta.data)){
      stop(paste0("The parameter 'split.by' ", split.by, " does not exist in MetaData slot!\n"))
    }
    Data <- object@meta.data[, c(reduction_ids, "cluster", split.by)]
    Data[, "cluster"] <- factor(Data[, "cluster"], levels=levels(Seurat::Idents(object)))
    Data[, "color"] <- Data[, "cluster"]
    levels(Data[, "color"]) <- cols[1:length(levels(Seurat::Idents(object)))]
    xmin <- min(Data[,1])
    xmax <- max(Data[,1])
    ymin <- min(Data[,2])
    ymax <- max(Data[,2])
    if(is.null(ncol)){
      ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    }
    legend <- ps(data=Data, title=NULL, legend_title=NULL)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, 'box', trim=F)  # use trim depending on need
    for (s in unique(Data[, split.by])){
      data <- Data[Data[, split.by] == s,]
      nearest.point <- RANN::nn2(data=data[,1:2], query=GetXYCenter(data)[,1:2], k=1)$nn.idx
      data[, "label"] <- NA
      data[nearest.point, "label"] <- as.vector(GetXYCenter(data)[, "label"])
      plots[[s]] <- pm(data, title = s, legend_title = NULL)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | legend) + patchwork::plot_layout(widths = c(3,1)))
}
