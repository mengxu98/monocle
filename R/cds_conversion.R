#' Export a monocle CellDataSet object to a Seurat object
#'
#' @param monocle_cds the Monocle CellDataSet you would like to export into a Seurat object
#' @param export_all Whether or not to export all the slots in Monocle and keep in the Seurat object. Default is FALSE (or only keep
#' minimal dataset). If export_all is setted to be true, the original monocle CellDataSet object will be keeped in the Seurat object too.
#' @return a new Seurat object.
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' seurat_lung <- exportCDS(lung)
#' seurat_lung_all <- exportCDS(lung, export_all = T)
#' }
exportCDS <- function(monocle_cds, export_all = FALSE) {
  requireNamespace("Seurat")
  data <- exprs(monocle_cds)
  ident <- colnames(monocle_cds)

  if (export_all) {
    monocle_cds@auxClusteringData$seurat <- NULL
    monocle_cds@auxClusteringData$scran <- NULL
    mist_list <- monocle_cds
  } else {
    mist_list <- list()
  }
  if ("use_for_ordering" %in% colnames(fData(monocle_cds))) {
    var.gene <- row.names(subset(fData(monocle_cds), use_for_ordering == TRUE))
  }

  srt <- Seurat::CreateSeuratObject(
    raw.data = data,
    normalization.method = "LogNormalize",
    do.scale = TRUE,
    do.center = TRUE,
    is.expr = monocle_cds@lowerDetectionLimit,
    project = "exportCDS",
    meta.data = pData(monocle_cds)
  )

  srt@misc <- mist_list
  srt@meta.data <- pData(monocle_cds)

  return(srt)
}

#' Import a Seurat object and convert it to a monocle CellDataSet object
#'
#' @param srt the Seurat object you would like to convert into a monocle CellDataSet object
#' @param import_all Whether or not to import all the slots in Seurat.
#' Default is FALSE (or only keep minimal dataset).
#' @return a new monocle CellDataSet object converted from a Seurat object.
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' seurat_lung <- exportCDS(lung)
#' seurat_lung_all <- exportCDS(lung, export_all = T)
#'
#' importCDS(seurat_lung)
#' importCDS(seurat_lung, import_all = T)
#' importCDS(seurat_lung_all)
#' importCDS(seurat_lung_all, import_all = T)
#' }
importCDS <- function(srt, import_all = FALSE) {
  requireNamespace("Seurat")
  data <- srt@raw.data
  pd <- tryCatch(
    {
      pd <- new("AnnotatedDataFrame", data = srt@meta.data)
      pd
    },
    error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)

      log_message("This Seurat object doesn't provide any meta data")
      pd
    }
  )

  if (length(setdiff(colnames(data), rownames(pd))) > 0) {
    data <- data[, rownames(pd)]
  }

  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  lowerDetectionLimit <- srt@is.expr

  if (all(data == floor(data))) {
    expressionFamily <- negbinomial.size()
  } else if (any(data < 0)) {
    expressionFamily <- uninormal()
  } else {
    expressionFamily <- tobit()
  }

  valid_data <- data[, row.names(pd)]

  monocle_cds <- newCellDataSet(data,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = lowerDetectionLimit,
    expressionFamily = expressionFamily
  )

  if (import_all) {
    if ("Monocle" %in% names(srt@misc)) {
      srt@misc$Monocle@auxClusteringData$seurat <- NULL
      srt@misc$Monocle@auxClusteringData$scran <- NULL

      monocle_cds <- srt@misc$Monocle
      mist_list <- srt
    } else {
      mist_list <- srt
    }
  } else {
    mist_list <- list()
  }

  if ("var.genes" %in% slotNames(srt)) {
    var.genes <- setOrderingFilter(monocle_cds, srt@var.genes)
  }
  monocle_cds@auxClusteringData$seurat <- mist_list

  return(monocle_cds)
}
