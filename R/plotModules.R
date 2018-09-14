#' @title plotModules
#'
#' @description Use a bubbleplot to display expression levels and percentage for
#' genes identified as functioning within a network module.
#'
#' @param seuratObj
#' @param modules Gene expression modules to plot, as determined by
#' seuratClusterWGCNA
#' @param ... Additional arguments to pass to bubbleplot
#'
#' @import ggplot2
#' @import magrittr
#' @import patchwork
#' @importFrom seuratBubblePlot bubbleplot
#' @importFrom gridExtra grid.arrange
#' @importFrom purrr reduce
#'
#' @return ggplot2 plots
#' @export
#'
#' @examples
plotModules <- function(seuratObj, modules, ...) {
  modplots <- lapply(
    (1:length(modules)),
    function(x) {
      bp <- bubbleplot(
        seuratObj = seuratObj,
        genes.plot = modules[[x]],
        do.return = TRUE,
        ...
      ) +
        labs(title = names(modules)[[x]]) +
        xlab(NULL) +
        ylab(NULL) +
        theme(legend.position = "none")
    }
  )

  reduce(modplots, `+`)
}
