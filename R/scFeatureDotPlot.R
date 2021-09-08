#' scFeatureDotPlot
#' wrapper for Seurat::DotPlot
#' @title scFeatureDotPlot
#' @param object Seurat object
#' @param features Vector of features to plot
#' @param cols Vector of two colors to form the gradient over
#' @param dot.scale Scale the size of the dot
#' @param font.size Font size
#' @return NULL
#' @author Wei Zhang
#' @export
scFeatureDotPlot <- function(object, features, cols=NULL, dot.scale=NULL, font.size=NULL){
  features <- rev(features)
  features <- intersect(unique(features), rownames(object))
  if (is.null(cols)){
    cols <- c('#FFFFFF','#16388E')
  }
  if (is.null(dot.scale)){
    dot.scale <- 3
  }
  if (is.null(font.size)){
    font.size <- 8
  }
  p <- Seurat::DotPlot(object=object, features = features, cols = cols, dot.scale = dot.scale) +
       Seurat::RotatedAxis() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = font.size),
                                            axis.text.y = ggplot2::element_text(size = font.size)) +
       ggplot2::guides(color = ggplot2::guide_colorbar(title = 'Scaled Expression'),
                    size = ggplot2::guide_legend(title = 'Percent Expressed')) +
       ggplot2::theme(axis.line = ggplot2::element_line(size = 0.6)) + ggplot2::labs(x='', y='')
  return(p)
}
