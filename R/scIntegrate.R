#' Integrate seurat objects into one
#' @title scIntegrate
#' @param object.list A list of seurat objects
#' @param object.names An array of seurat object names
#' @param nVariable Number of features to select as top variable features
#' @param nPC Number of PCs to use
#' @param res Resolution parameter to set. Default res=0.5
#' @param batch.rm Remove batch effect with 'seurat' or 'harmony'. Default batch.rm="harmony"
#' @return Seurat object.
#' @importFrom rlang %>%
#' @importFrom rlang %||%
#' @importFrom Seurat DefaultAssay

#' @author rstatistics
#' @export

scIntegrate <- function(object.list=NULL, object.names=NULL, nVariable=2000, nPC=30, res=0.5, batch.rm="harmony"){
  # Merge each seurat object
  meta.list <- list()
  meta.table <- object.list[[1]]@meta.data
  for (i in 2:length(object.list)){
      meta.table <- dplyr::bind_rows(meta.table, object.list[[i]]@meta.data)
  }
  meta.table[is.na(meta.table)] <- 0

  object <- merge(object.list[[1]],  y=object.list[2:length(object.list)],
                          add.cell.ids = as.character(object.names), project = "seurat")
  rownames(meta.table) <- rownames(object@meta.data)
  object@meta.data <- meta.table

  object <- Seurat::NormalizeData(object, verbose = TRUE)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = nVariable)
  object <- Seurat::ScaleData(object, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = TRUE)
  object <- Seurat::RunPCA(object, pc.genes = Seurat::VariableFeatures(object), npcs = nPC, verbose = TRUE)

  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- Seurat::DimPlot(object = object, reduction = "pca", pt.size = .1, group.by = "sample")
  p2 <- Seurat::VlnPlot(object = object, features = "PC_1", pt.size = .1, group.by = "sample")
  plot(cowplot::plot_grid(p1,p2))
  object <- harmony::RunHarmony(object = object, group.by.vars = "sample", plot_convergence = TRUE)

  if (batch.rm=="harmony"){
    object <- Seurat::FindNeighbors(object, reduction = "harmony", dims = 1:nPC)
    object <- Seurat::FindClusters(object, resolution = res)
    object <- Seurat::RunUMAP(object, reduction = "harmony", dims = 1:nPC)
    object <- Seurat::RunTSNE(object, reduction = "harmony", dims = 1:nPC, do.fast=TRUE)
  }else if (batch.rm=="seurat"){
    object <- Seurat::FindNeighbors(object, reduction = "pca", dims = 1:nPC)
    object <- Seurat::FindClusters(object, resolution = res)
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = 1:nPC)
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = 1:nPC)
  }else{
    stop("batch.rm must be 'harmony' or 'seurat'!")
  }

  cols <- NA
  if (length(levels(Seurat::Idents(object))) <= 36){
    cols = c("#1660A7","#FF6A00","#219418","#CD0C18","#814BB2","#794339","#DC59B6","#CC79A7","#FF0000","#11B3C6",
             "#AFB400","#00FFFF", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00",
             "#CC79A7", "#00AFBB", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#4477AA",
             "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    object@misc$cols = cols
  }
  plot(Seurat::DimPlot(object, reduction = "umap", pt.size = 0.1, cols = cols))
  # Add reduction coordinates into Meta.data slot
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions){
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, reduction = reduction)[, 1:2]
    object <- Seurat::AddMetaData(object = object, metadata = coord, col.name = meta_ids)
  }
  return(object)
}
