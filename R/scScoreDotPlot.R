#' MSigdb Scoring DotPlot for single-cell clusters
#' @title scScoreDotPlot
#' @param object Seurat object
#' @param geneSets Vector of gene sets to display
#' @param cols Colors to plot, default cols=c('blue', 'orange', 'red')
#' @param col.min Minimum scaled average score threshold
#' @param col.max Maximum scaled average score threshold
#' @param dot.min Lower limit of cell fraction to draw the dots (default is 0)
#' @param dot.scale Scale the size of the dots
#' @param group.by Group the cells by
#' @param split.by Split the cells by
#' @param scale Whether to scale the data, default scale = TRUE
#' @param scale.min Lower limit to scale the data, default scale.min = NA
#' @param scale.max Upper limit to scale the data, default scale.max = NA
#' @return NULL
#'
#' @author rstatistics
#' @export
#'
scScoreDotPlot <- function (object=NULL, geneSets=NULL, cols = c('blue', 'orange', 'red'),
            col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL,
            split.by = NULL, scale = TRUE, scale.min = NA, scale.max = NA)
  {
    if (is.null(geneSets)){
      stop("geneSets must be specified!\n")
    }
    data.pathways <- object[[geneSets]]
    if (is.null(group.by)) {
      data.pathways$id <- Seurat::Idents(object = object)
    } else {
      data.pathways$id <- object[[group.by, drop = TRUE]]
    }
    id.levels <- levels(factor(data.pathways$id))
    data.pathways[, "id"] <- as.vector(data.pathways[, "id"])
    if (!is.null(split.by)) {
      split_elements <- object[[split.by, drop = TRUE]]
      if (length(unique(split_elements)) > length(cols)) { stop("Not enought colors for the number of split elements") }
      cols <- cols[1:length(unique(split_elements))]
      names(cols) <- unique(split_elements)
      data.pathways$id <- paste(data.pathways$id, split_elements, sep = "_")
      unique.split_elements <- unique(split_elements)
      id.levels <- paste0(rep(id.levels, each = length(unique.split_elements)),
                          "_", rep(unique(split_elements), times = length(id.levels)))
    }
    data.plot <- lapply(X = unique(data.pathways$id), FUN = function(ident) {
      data.use <- data.pathways[data.pathways$id == ident, 1:(ncol(data.pathways) - 1), drop = FALSE]
      avg.score <- apply(data.use, 2, function(x) { return(mean(expm1(x))) })
      pct.score <- apply(data.use, 2, function(x){ return(length(x[x > 0]) / length(x)) })
      return(list(avg.score = avg.score, pct.score = pct.score))
    })
    names(data.plot) <- unique(data.pathways$id)
    data.plot <- lapply(X = names(data.plot), FUN = function(i) {
      data.use <- as.data.frame(data.plot[[i]])
      data.use$pathways.plot <- rownames(data.use)
      data.use$id <- i
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(id.levels)) {
      data.plot$id <- factor(data.plot$id, levels = id.levels)
    }
    if (length(levels(data.plot$id)) == 1) {
      scale <- FALSE
      warning("Only one identity present, the pathway scores will not be scaled.")
    }
    avg.score.scaled <- sapply(X = unique(data.plot$pathways.plot), FUN = function(x) {
                               data.use <- data.plot[data.plot$pathways.plot == x, "avg.score"]
                               if (scale) { data.use <- scale(data.use)
                                            data.use <- Seurat::MinMax(data = data.use, min = col.min, max = col.max)
                               } else { data.use <- log(data.use) }
                               return(data.use)
                             })
    avg.score.scaled <- as.vector(t(avg.score.scaled))
    if (!is.null(split.by)) { avg.score.scaled <- as.numeric(cut(avg.score.scaled, breaks = 20)) }
    data.plot$avg.score.scaled <- avg.score.scaled
    data.plot$pathways.plot <- factor(data.plot$pathways.plot, levels = rev(geneSets))
    data.plot$pct.score[data.plot$pct.score < dot.min] <- NA
    data.plot$pct.score <- data.plot$pct.score * 100
    if (!is.null(split.by)) {
      split_elements.use <- vapply(X = as.character(data.plot$id), FUN = gsub, FUN.VALUE = character(length = 1L),
                                   pattern = paste0("^((", paste(sort(levels(object), decreasing = TRUE),
                                            collapse = "|"), ")_)"), replacement = "", USE.NAMES = FALSE)
      data.plot$colors <- mapply(FUN = function(color, value) { return(colorRampPalette(colors = c("grey", color))(20)[value])
                                                              }, color = cols[split_elements.use], value = avg.score.scaled)
    }
    color.by <- ifelse(test = is.null(split.by), yes = "avg.score.scaled", no = "colors")
    if (!is.na(scale.min)) { data.plot[data.plot$pct.score < scale.min, "pct.score"] <- scale.min }
    if (!is.na(scale.max)) { data.plot[data.plot$pct.score > scale.max, "pct.score"] <- scale.max }
    p <- ggplot2::ggplot(data = data.plot, ggplot2::aes_string(x = "pathways.plot", y = "id")) +
              ggplot2::geom_point(ggplot2::aes_string(size = "pct.score", color = color.by)) +
              ggplot2::scale_radius(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
              ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()) +
              ggplot2::guides(size = ggplot2::guide_legend(title = "Percent Expressed")) +
              cowplot::theme_cowplot() + ggplot2::labs(x = "Gene Sets", y = "Identity") + ggplot2::coord_flip() +
              Seurat::RotatedAxis()
    if (!is.null(split.by)) {
      p <- p + ggplot2::scale_color_identity()
    } else if (length(cols) == 1) {
      p <- p + ggplot2::scale_color_distiller(palette = cols)
    } else {
      p <- p + ggplot2::scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(split.by)) {
      p <- p + ggplot2::guides(color = ggplot2::guide_colorbar(title = "Average Score"))
    }
    return(p)
}
