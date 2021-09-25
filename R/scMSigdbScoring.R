#' MSigdb signature scoring for single-cell clusters
#' @title scMsigdbScoring
#' @param object Seurat object
#' @param category String of specific MSigDB category to use, 
#'   default category='H' stands for 'hallmark gene sets'
#' @param geneSets Vector of gene sets, default geneSets='NULL' 
#'   stands for all gene sets in category will be used
#' @return Seurat object
#'
#' @author Wei Zhang
#' @export

scMsigdbScoring <- function(object = NULL, category = NULL, geneSets = NULL) {
  if (is.null(category)) category <- "H"
  MSigDB <- system.file("extdata", "MSIGDB.gmt", package = "Yeskit")
  MSIGDB <- read.table(MSigDB, row.names = 1, col.names = paste0("V", 
    seq_len(max(count.fields(MSigDB, sep = "\t")))), header = FALSE, 
    sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  if (!category %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) {
    stop("category must be one of these choices: ", 
    "\"H\", \"C1\", \"C2\", \"C3\", \"C4\", \"C5\", \"C6\", \"C7\", \"C8\"")
  }
  MSIGDB <- MSIGDB[MSIGDB[, 1] == category, ]
  MSIGDB[, 1] <- NULL
  MSIGDB <- as.list(data.frame(apply(MSIGDB, 1, as.vector), 
    stringsAsFactors = FALSE), drop = TRUE)
  if (length(names(MSIGDB)) == 0) {
    stop("MSigDB has no valid entries!")
  }
  if (is.null(geneSets)) {
    geneSets = names(MSIGDB)
  } else {
    geneSets <- intersect(geneSets, names(MSIGDB))
    if (length(geneSets) == 0) {
      stop("geneSets ", geneSets, "does not exist in ", MSigDB, "!")
    }
  }
  for (item in geneSets) {
    if (item %in% names(object[[]])) {
      object[[item]] <- NULL
    }
    features <- MSIGDB[[item]]
    features <- features[grep("^$", features, invert = TRUE)]
    features <- list(Score = features)
    object.hallmark <- Seurat::AddModuleScore(object = object, 
      features = features, name = item, ctrl = 
        min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))))
    hallmark.columns <- grep(pattern = item, x = colnames(x = object.hallmark[[]]), value = TRUE)
    hallmark.scores <- object.hallmark[[hallmark.columns]]
    rm(object.hallmark)
    colnames(x = hallmark.scores) <- c(item)
    object[[colnames(x = hallmark.scores)]] <- hallmark.scores
  }
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions) {
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, reduction = reduction)[, c(1, 2)]
    object <- Seurat::AddMetaData(object = object, metadata = coord, col.name = meta_ids)
  }
  return(object)
}

