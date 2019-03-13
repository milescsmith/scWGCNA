#' @title GatherData
#'
#' @param object scRNA-seq data object to pull data from
#' @param assay Assay to pull data from
#' @param slot_use Slot to pull data from.  Ignored for SingleCellExperiment objects
#' @param ... Additional arguments
#' @importFrom purrr %||%
GatherData <- function(object, ...) {
  UseMethod("GatherData")
}

#' @rdname GatherData
#' @method GatherData Seurat
#' @importFrom Seurat GetAssayData
#' @return matrix
#' @export
GatherData.Seurat <- function(object,
                              assay,
                              slot_use,
                              ...) {
  assay <- assay %||% "RNA"
  slot_use <- slot_use %||% "data"
  obj_data <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot_use
  ) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays
#' @return matrix
#' @export
GatherData.SingleCellExperiment <- function(object,
                                            assay,
                                            ...) {
  assay <- assay %||% "logcounts"
  if (!assay %in% names(assays(object))) {
    stop(glue("Sorry, but {assay} is not present in the object.
              Please run the assay or choose a different assay slot."))
  }
  obj_data <- assay(
    object = object,
    i = assay
  ) %>%
    as.matrix()

  return(obj_data)
}
